from copy import deepcopy
from pathlib import Path
from typing import Any
import sys
import time

import numpy as np
import json
import torch
import typer
from tqdm import tqdm
from eliot import start_action

from data.pdb_utils import VOCAB
from data.dataset import PDBDataset
from trainers.abs_trainer import Trainer
from models.prediction_model import PredictionModel

app = typer.Typer(help="Compute interaction scores for protein structures", add_completion=False)


def mask_block(data: dict[str, Any], block_idx: int) -> dict[str, Any]:
    """Mask a specific block in the data by replacing it with MASK token.
    
    Args:
        data: The input data dictionary
        block_idx: Index of the block to mask
        
    Returns:
        Modified data dictionary with the specified block masked
    """
    data = deepcopy(data)
    for key in data:
        if isinstance(data[key], np.ndarray):
            data[key] = data[key].tolist()
    data["B"][block_idx] = VOCAB.symbol_to_idx(VOCAB.MASK)
    block_start = sum(data["block_lengths"][:block_idx])
    block_end = block_start + data["block_lengths"][block_idx]
    data["block_lengths"][block_idx] = 1
    data["X"] = (
        data["X"][:block_start]
        + [np.mean(data["X"][block_start:block_end], axis=0).tolist()]
        + data["X"][block_end:]
    )
    data["A"] = (
        data["A"][:block_start]
        + [VOCAB.get_atom_mask_idx()]
        + data["A"][block_end:]
    )
    data["atom_positions"] = (
        data["atom_positions"][:block_start]
        + [VOCAB.get_atom_pos_mask_idx()]
        + data["atom_positions"][block_end:]
    )
    return data


def get_residue_model_score(
    model: PredictionModel, data: dict[str, Any], block_idx: int
) -> float:
    """Calculate cosine similarity score for a masked block.
    
    Args:
        model: The prediction model
        data: The input data dictionary
        block_idx: Index of the block to evaluate
        
    Returns:
        Cosine similarity score between original and masked predictions
    """
    with torch.no_grad():
        model.eval()
        masked_data = mask_block(data, block_idx)
        batch = PDBDataset.collate_fn([data, masked_data])
        batch = Trainer.to_device(batch, "cuda")
        output = model(
            batch["X"],
            batch["B"],
            batch["A"],
            batch["block_lengths"],
            batch["lengths"],
            batch["segment_ids"],
        )
        cos_distance = torch.nn.functional.cosine_similarity(
            output.graph_repr[0], output.graph_repr[1], dim=-1
        ).item()
    return cos_distance


def get_residue_model_scores(
    model: PredictionModel, data: dict[str, Any]
) -> tuple[list[float], list[int]]:
    """Calculate cosine similarity scores for all blocks in the data.
    
    Args:
        model: The prediction model
        data: The input data dictionary
        
    Returns:
        Tuple of (cosine_distances, block_indices) for non-GLB blocks
    """
    cos_distances: list[float] = []
    block_idx: list[int] = []
    for i in range(len(data["B"])):
        if data["B"][i] == VOCAB.symbol_to_idx(VOCAB.GLB):
            continue
        cos_distances.append(get_residue_model_score(model, data, i))
        block_idx.append(i)
    return cos_distances, block_idx


@app.command()
def compute_interact_scores(
    data_path: Path = typer.Option(
        ..., help="Path to the data file", exists=True, file_okay=True
    ),
    output_path: Path = typer.Option(
        ..., help="Output JSON file for importance scores", file_okay=True
    ),
    model_ckpt: Path = typer.Option(
        ..., help="Path to the model checkpoint", exists=True, file_okay=True
    ),
    start_idx: int = typer.Option(
        0, help="Start index for processing (useful for large files)"
    ),
    num_lines: int | None = typer.Option(
        None, help="Maximum number of lines to process (None = all lines)"
    ),
    summary_path: Path | None = typer.Option(
        None, help="Optional path to save processing summary report (JSON)"
    ),
) -> None:
    """Compute interaction scores (importance of each block) for protein structures.
    
    This tool evaluates how much each block (residue or molecule) contributes to
    the model's predictions by measuring the cosine similarity between original
    and block-masked predictions.
    
    For large .jsonl.gz files (gigabytes), use --start-idx and --num-lines to process
    specific slices:
    
    Example: Process only lines 1000-2000 from a large file:
        interact-profiler compute-interact-scores \\
          --data-path large_file.jsonl.gz \\
          --output-path output.jsonl \\
          --model-ckpt model.pt \\
          --start-idx 1000 \\
          --num-lines 1000 \\
          --summary-path summary.json
    """
    with start_action(
        action_type="compute_interact_scores",
        data_path=str(data_path),
        output_path=str(output_path),
        model_ckpt=str(model_ckpt),
        start_idx=start_idx,
        num_lines=num_lines,
        summary_path=str(summary_path) if summary_path else None,
    ) as action:
        # Clear CUDA cache before starting
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        
        model = PredictionModel.load_from_pretrained(str(model_ckpt))
        model = model.to("cuda")

        dataset = PDBDataset(str(data_path), start_idx=start_idx, num_lines=num_lines)
        
        times: list[float] = []
        peak_memory_per_structure: list[float] = []
        total_start = time.perf_counter()
        
        for i in tqdm(range(len(dataset)), total=len(dataset), desc="Processing structures"):
            structure_start = time.perf_counter()
            torch.cuda.reset_peak_memory_stats()
            
            cos_distances, block_idx = get_residue_model_scores(model, dataset[i])
            
            structure_time = time.perf_counter() - structure_start
            times.append(structure_time)
            
            # Get peak memory usage for this structure (in MB)
            peak_memory_mb = torch.cuda.max_memory_allocated() / 1024 / 1024
            peak_memory_per_structure.append(peak_memory_mb)
            
            output = {
                "id": dataset.indexes[i],
                "cos_distances": cos_distances,
                "block_idx": block_idx,
                "time_seconds": round(structure_time, 4),
                "peak_memory_mb": round(peak_memory_mb, 2),
            }
            with open(output_path, "a") as f:
                f.write(json.dumps(output) + "\n")
            
            # Log metrics periodically
            if (i + 1) % 10 == 0:
                avg_time = np.mean(times[-10:])
                avg_memory = np.mean(peak_memory_per_structure[-10:])
                tqdm.write(
                    f"[{i+1:6d}] Avg time (last 10): {avg_time:.2f}s, "
                    f"Avg memory (last 10): {avg_memory:.1f} MB"
                )

        total_time = time.perf_counter() - total_start
        
        # Compute statistics
        avg_time = np.mean(times)
        std_time = np.std(times)
        min_time = np.min(times)
        max_time = np.max(times)
        
        avg_memory = np.mean(peak_memory_per_structure)
        max_memory = np.max(peak_memory_per_structure)
        
        # Create summary report
        summary_report = {
            "total_structures": len(dataset),
            "total_time_seconds": round(total_time, 2),
            "total_time_minutes": round(total_time / 60, 2),
            "time_per_structure": {
                "mean": round(avg_time, 4),
                "std": round(std_time, 4),
                "min": round(min_time, 4),
                "max": round(max_time, 4),
            },
            "gpu_memory_mb": {
                "mean": round(avg_memory, 2),
                "max": round(max_memory, 2),
            },
            "output_file": str(output_path),
            "data_file": str(data_path),
            "start_idx": start_idx,
            "num_lines": num_lines,
        }
        
        # Print summary to console
        print("\n" + "="*70)
        print("PROCESSING SUMMARY")
        print("="*70)
        print(f"Total structures processed: {summary_report['total_structures']}")
        print(f"Total time: {summary_report['total_time_seconds']:.2f}s "
              f"({summary_report['total_time_minutes']:.2f} minutes)")
        print(f"Time per structure: {summary_report['time_per_structure']['mean']:.2f}s ± "
              f"{summary_report['time_per_structure']['std']:.2f}s "
              f"(min: {summary_report['time_per_structure']['min']:.2f}s, "
              f"max: {summary_report['time_per_structure']['max']:.2f}s)")
        print(f"Estimated remaining time per structure: ~{summary_report['time_per_structure']['mean']:.2f}s")
        print(f"\nGPU Memory per structure: {summary_report['gpu_memory_mb']['mean']:.1f} MB (avg), "
              f"{summary_report['gpu_memory_mb']['max']:.1f} MB (peak)")
        print("="*70)
        print(f"Output saved to: {output_path}")
        
        # Save summary report if requested
        if summary_path:
            with open(summary_path, "w") as f:
                json.dump(summary_report, f, indent=2)
            print(f"Summary report saved to: {summary_path}")
        print()
        
        action.add_success_fields(
            num_structures=len(dataset),
            total_time_seconds=round(total_time, 2),
            avg_time_per_structure=round(avg_time, 4),
            avg_memory_mb=round(avg_memory, 2),
            max_memory_mb=round(max_memory, 2),
            summary_saved=summary_path is not None,
        )


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()
            
