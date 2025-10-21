from copy import deepcopy
from pathlib import Path
from typing import Any, Optional, List
import sys
import time
import tempfile

import numpy as np
import json
import torch
import typer
from tqdm import tqdm
from eliot import start_action
import orjson
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx

from data.pdb_utils import VOCAB
from data.dataset import PDBDataset
from trainers.abs_trainer import Trainer
from models.prediction_model import PredictionModel

# Import from parent directory
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))
from pdb_converter import pdb_to_jsonl_item

app = typer.Typer(help="Compute interaction scores for protein structures", add_completion=False)


def is_pdb_format(file_path: Path) -> bool:
    """Check if file is PDB/CIF format (vs JSONL/pickle)."""
    suffix = file_path.suffix.lower()
    if suffix in ['.pdb', '.cif', '.ent']:
        return True
    if suffix == '.gz':
        stem_suffix = Path(file_path.stem).suffix.lower()
        return stem_suffix in ['.pdb', '.cif', '.ent']
    return False


def convert_pdb_to_jsonl_temp(
    pdb_path: Path,
    chains: Optional[List[str]] = None
) -> Path:
    """Convert PDB/CIF file to temporary JSONL file."""
    with start_action(action_type="convert_pdb_to_jsonl", pdb_path=str(pdb_path), chains=chains):
        temp_file = Path(tempfile.mktemp(suffix='.jsonl'))
        item = pdb_to_jsonl_item(pdb_path, chains=chains)
        
        # Convert integer keys to strings for orjson compatibility
        if 'block_to_pdb_indexes' in item and isinstance(item['block_to_pdb_indexes'], dict):
            item['block_to_pdb_indexes'] = {str(k): v for k, v in item['block_to_pdb_indexes'].items()}
        
        with open(temp_file, 'w') as f:
            f.write(orjson.dumps(item).decode('utf-8') + '\n')
        return temp_file


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
    model: PredictionModel, data: dict[str, Any], block_idx: int, device: str = "cuda"
) -> float:
    """Calculate cosine similarity score for a masked block.
    
    Args:
        model: The prediction model
        data: The input data dictionary
        block_idx: Index of the block to evaluate
        device: Device to use for computation ('cuda' or 'cpu')
        
    Returns:
        Cosine similarity score between original and masked predictions
    """
    with torch.no_grad():
        model.eval()
        masked_data = mask_block(data, block_idx)
        batch = PDBDataset.collate_fn([data, masked_data])
        batch = Trainer.to_device(batch, device)
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
    model: PredictionModel, data: dict[str, Any], device: str = "cuda"
) -> tuple[list[float], list[int]]:
    """Calculate cosine similarity scores for all blocks in the data.
    
    Args:
        model: The prediction model
        data: The input data dictionary
        device: Device to use for computation ('cuda' or 'cpu')
        
    Returns:
        Tuple of (cosine_distances, block_indices) for non-GLB blocks
    """
    cos_distances: list[float] = []
    block_idx: list[int] = []
    for i in range(len(data["B"])):
        if data["B"][i] == VOCAB.symbol_to_idx(VOCAB.GLB):
            continue
        cos_distances.append(get_residue_model_score(model, data, i, device))
        block_idx.append(i)
    return cos_distances, block_idx


def get_residue_info_from_pdb(
    pdb_file: Path,
    block_indices: list[int],
    scores: list[float]
) -> list[dict[str, Any]]:
    """Extract residue information from PDB file.
    
    Args:
        pdb_file: Path to PDB/CIF file
        block_indices: List of block indices
        scores: List of ATOMICA scores
        
    Returns:
        List of dictionaries with residue information
    """
    # Load structure
    if pdb_file.suffix == '.cif':
        pdbx_file = pdbx.CIFFile.read(str(pdb_file))
        structure = pdbx.get_structure(pdbx_file, model=1)
    else:
        pdb_file_obj = pdb.PDBFile.read(str(pdb_file))
        structure = pdb.get_structure(pdb_file_obj, model=1)
    
    # Get CA atoms (one per residue)
    ca_atoms = structure[structure.atom_name == 'CA']
    
    residues = []
    for idx, score in zip(block_indices, scores):
        if idx < len(ca_atoms):
            atom = ca_atoms[idx]
            residues.append({
                'block_idx': idx,
                'chain_id': atom.chain_id,
                'res_id': atom.res_id,
                'res_name': atom.res_name,
                'atomica_score': float(score),
                'importance_delta': float((1 - score) * 100),
            })
    return residues


def generate_pymol_commands(
    residues: list[dict[str, Any]],
    structure_id: str,
    top_n: int = 10
) -> str:
    """Generate PyMOL commands for visualization.
    
    Args:
        residues: List of residue dictionaries sorted by importance
        structure_id: Structure identifier
        top_n: Number of top residues to highlight
        
    Returns:
        PyMOL command script as string
    """
    top_res = residues[:top_n]
    res_ids = [r['res_id'] for r in top_res]
    chain = top_res[0]['chain_id']
    res_str = '+'.join([str(r) for r in res_ids])
    
    commands = f"""# ============================================================
# PyMOL Commands for {structure_id}
# ============================================================
# Copy and paste these commands into PyMOL
# ============================================================

# Load structure (adjust path as needed)
load {structure_id}.pdb

# Basic setup
hide all
show cartoon
color grey80, all

# Select and highlight critical residues (Top {top_n})
select critical_top{top_n}, resi {res_str} and chain {chain}
show sticks, critical_top{top_n}
color red, critical_top{top_n}
set stick_radius, 0.3, critical_top{top_n}

# Label critical residues
label critical_top{top_n} and name CA, "%s%s" % (resn,resi)
set label_size, 14
set label_color, red

# Create gradient coloring by importance (Top 5 in different reds)
"""
    
    # Individual residue selections with gradient colors
    colors = ['red', 'tv_red', 'salmon', 'lightsalmon', 'warmpink']
    for i, res in enumerate(top_res[:5]):
        commands += f"select critical_rank{i+1}, resi {res['res_id']} and chain {res['chain_id']}\n"
        commands += f"color {colors[i]}, critical_rank{i+1}\n"
    
    commands += f"""
# Show surface around critical residues
show surface, byres (critical_top{top_n} around 5)
set surface_color, white
set transparency, 0.5

# Center view on critical residues
zoom critical_top{top_n}

# ============================================================
# Additional useful commands:
# ============================================================
# To save session: save {structure_id}_critical.pse
# To save image: png {structure_id}_critical.png, dpi=300
# To highlight most critical (rank 1): select most_critical, critical_rank1
# ============================================================
"""
    return commands


def generate_tsv_report(
    residues: list[dict[str, Any]],
    structure_id: str,
    structure_name: str,
    scores: list[float]
) -> str:
    """Generate TSV report with all critical residues.
    
    Args:
        residues: List of residue dictionaries sorted by importance
        structure_id: Structure identifier
        structure_name: Human-readable structure name
        scores: List of all ATOMICA scores
        
    Returns:
        TSV formatted string with headers
    """
    # Header with metadata as comments
    tsv = f"# ATOMICA Critical Residue Analysis: {structure_name}\n"
    tsv += f"# Structure ID: {structure_id}\n"
    tsv += f"# Total residues analyzed: {len(residues)}\n"
    tsv += f"# Mean ATOMICA_SCORE: {np.mean(scores):.6f}\n"
    tsv += f"# Std Dev: {np.std(scores):.6f}\n"
    tsv += f"# Method: ATOMICA_SCORE (cosine similarity with masked residues)\n"
    tsv += f"# Lower scores = More critical for intermolecular interactions\n"
    tsv += "#\n"
    
    # Column headers
    tsv += "Rank\tChain_ID\tResidue_ID\tResidue_Name\tATOMICA_Score\tImportance_Delta_Percent\tBlock_Index\n"
    
    # Data rows
    for i, res in enumerate(residues, 1):
        tsv += f"{i}\t{res['chain_id']}\t{res['res_id']}\t{res['res_name']}\t{res['atomica_score']:.6f}\t{res['importance_delta']:.4f}\t{res['block_idx']}\n"
    
    return tsv


@app.command(name="interact-score")
def interact_score(
    input: Path = typer.Option(
        ..., help="Path to input file (JSON/JSONL/pickle or PDB/CIF file)", exists=True, file_okay=True
    ),
    output: Optional[Path] = typer.Option(
        None, help="Output JSON file for importance scores (default: output/{input_stem}_interact_scores.json)", file_okay=True
    ),
    model_ckpt: Optional[Path] = typer.Option(
        None, help="Path to the model checkpoint (or use --model-config and --model-weights)", exists=True, file_okay=True
    ),
    model_config: Optional[Path] = typer.Option(
        None, help="Path to model config JSON (default: downloads/ATOMICA_checkpoints/pretrain/pretrain_model_config.json)", exists=True, file_okay=True
    ),
    model_weights: Optional[Path] = typer.Option(
        None, help="Path to model weights (default: downloads/ATOMICA_checkpoints/pretrain/pretrain_model_weights.pt)", exists=True, file_okay=True
    ),
    chains: Optional[str] = typer.Option(
        None, help="Comma-separated chain IDs (e.g., 'A,B') - only for PDB/CIF input"
    ),
    start_idx: int = typer.Option(
        0, help="Start index for processing (useful for large files) - ignored for PDB input"
    ),
    num_lines: int | None = typer.Option(
        None, help="Maximum number of lines to process (None = all lines) - ignored for PDB input"
    ),
    summary_path: Path | None = typer.Option(
        None, help="Path to save processing summary report (default: output/{input_stem}_summary.json)"
    ),
    cpu: bool = typer.Option(
        False, help="Use CPU instead of GPU for computation"
    ),
) -> None:
    """Compute interaction scores (importance of each block) for protein structures.
    
    Supports both JSONL/pickle format and direct PDB/CIF input.
    
    This tool evaluates how much each block (residue or molecule) contributes to
    the model's predictions by measuring the cosine similarity between original
    and block-masked predictions.
    
    For large .jsonl.gz files (gigabytes), use --start-idx and --num-lines to process
    specific slices:
    
    Example: Process only lines 1000-2000 from a large file:
        interact-score \\
          --input large_file.jsonl.gz \\
          --output output.jsonl \\
          --model-ckpt model.pt \\
          --start-idx 1000 \\
          --num-lines 1000 \\
          --summary-path summary.json
    """
    # Set default paths if not provided
    if output is None:
        output_dir = Path("output")
        output_dir.mkdir(exist_ok=True)
        output = output_dir / f"{input.stem}_interact_scores.json"
    
    if summary_path is None:
        output_dir = Path("output")
        output_dir.mkdir(exist_ok=True)
        summary_path = output_dir / f"{input.stem}_summary.json"
    
    # Set default model paths if not provided
    if model_ckpt is None and (model_config is None or model_weights is None):
        project_root = Path(__file__).parent.parent
        default_config = project_root / "downloads/ATOMICA_checkpoints/pretrain/pretrain_model_config.json"
        default_weights = project_root / "downloads/ATOMICA_checkpoints/pretrain/pretrain_model_weights.pt"
        
        if model_config is None:
            if default_config.exists():
                model_config = default_config
            else:
                raise ValueError(
                    f"Default model config not found at {default_config}. "
                    "Please specify --model-config or --model-ckpt."
                )
        
        if model_weights is None:
            if default_weights.exists():
                model_weights = default_weights
            else:
                raise ValueError(
                    f"Default model weights not found at {default_weights}. "
                    "Please specify --model-weights or --model-ckpt."
                )
    
    # Parse chains if provided
    chain_list: Optional[List[str]] = None
    if chains:
        chain_list = [c.strip() for c in chains.split(",")]
    
    # Detect if input is PDB format and convert if needed
    temp_jsonl_file: Optional[Path] = None
    actual_data_path = input
    
    if is_pdb_format(input):
        typer.echo(f"🔄 PDB/CIF input detected, converting to JSONL...")
        temp_jsonl_file = convert_pdb_to_jsonl_temp(input, chains=chain_list)
        actual_data_path = temp_jsonl_file
        # For PDB input, ignore start_idx and num_lines
        start_idx = 0
        num_lines = None
        typer.echo(f"✓ Converted to temporary file")
    
    try:
        # Set device
        device = "cpu" if cpu else "cuda"
        
        with start_action(
            action_type="interact_score",
            input=str(actual_data_path),
            output=str(output),
            model_ckpt=str(model_ckpt),
            start_idx=start_idx,
            num_lines=num_lines,
            summary_path=str(summary_path) if summary_path else None,
            device=device,
        ) as action:
            # Clear CUDA cache before starting (only if using GPU)
            if not cpu:
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
            
            # Load model
            if model_ckpt:
                model = PredictionModel.load_from_pretrained(str(model_ckpt))
            elif model_config and model_weights:
                model = PredictionModel.load_from_config_and_weights(str(model_config), str(model_weights))
            else:
                raise ValueError("Either model_ckpt or both model_config and model_weights must be provided")
            
            model = model.to(device)

            dataset = PDBDataset(str(actual_data_path), start_idx=start_idx, num_lines=num_lines)
        
        times: list[float] = []
        peak_memory_per_structure: list[float] = []
        total_start = time.perf_counter()
        
        for i in tqdm(range(len(dataset)), total=len(dataset), desc="Processing structures"):
            structure_start = time.perf_counter()
            if not cpu:
                torch.cuda.reset_peak_memory_stats()
            
            cos_distances, block_idx = get_residue_model_scores(model, dataset[i], device)
            
            structure_time = time.perf_counter() - structure_start
            times.append(structure_time)
            
            # Get peak memory usage for this structure (in MB) - only for GPU
            if not cpu:
                peak_memory_mb = torch.cuda.max_memory_allocated() / 1024 / 1024
            else:
                peak_memory_mb = 0.0  # Not tracking CPU memory
            peak_memory_per_structure.append(peak_memory_mb)
            
            result = {
                "id": dataset.indexes[i],
                "cos_distances": cos_distances,
                "block_idx": block_idx,
                "time_seconds": round(structure_time, 4),
                "peak_memory_mb": round(peak_memory_mb, 2),
            }
            with open(output, "a") as f:
                f.write(json.dumps(result) + "\n")
            
            # Log metrics periodically
            if (i + 1) % 10 == 0:
                avg_time = np.mean(times[-10:])
                if not cpu:
                    avg_memory = np.mean(peak_memory_per_structure[-10:])
                    tqdm.write(
                        f"[{i+1:6d}] Avg time (last 10): {avg_time:.2f}s, "
                        f"Avg memory (last 10): {avg_memory:.1f} MB"
                    )
                else:
                    tqdm.write(f"[{i+1:6d}] Avg time (last 10): {avg_time:.2f}s")

        total_time = time.perf_counter() - total_start
        
        # Compute statistics
        avg_time = np.mean(times)
        std_time = np.std(times)
        min_time = np.min(times)
        max_time = np.max(times)
        
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
            "output_file": str(output),
            "input_file": str(input),
            "start_idx": start_idx,
            "num_lines": num_lines,
            "device": device,
        }
        
        # Add GPU memory stats only if using GPU
        if not cpu:
            avg_memory = np.mean(peak_memory_per_structure)
            max_memory = np.max(peak_memory_per_structure)
            summary_report["gpu_memory_mb"] = {
                "mean": round(avg_memory, 2),
                "max": round(max_memory, 2),
            }
        
        # Print summary to console
        typer.echo("\n" + "="*70)
        typer.echo("PROCESSING SUMMARY")
        typer.echo("="*70)
        typer.echo(f"Device: {device.upper()}")
        typer.echo(f"Total structures processed: {summary_report['total_structures']}")
        typer.echo(f"Total time: {summary_report['total_time_seconds']:.2f}s "
              f"({summary_report['total_time_minutes']:.2f} minutes)")
        typer.echo(f"Time per structure: {summary_report['time_per_structure']['mean']:.2f}s ± "
              f"{summary_report['time_per_structure']['std']:.2f}s "
              f"(min: {summary_report['time_per_structure']['min']:.2f}s, "
              f"max: {summary_report['time_per_structure']['max']:.2f}s)")
        typer.echo(f"Estimated remaining time per structure: ~{summary_report['time_per_structure']['mean']:.2f}s")
        if not cpu:
            typer.echo(f"\nGPU Memory per structure: {summary_report['gpu_memory_mb']['mean']:.1f} MB (avg), "
                  f"{summary_report['gpu_memory_mb']['max']:.1f} MB (peak)")
        typer.echo("="*70)
        typer.echo(f"Output saved to: {output}")
        
        # Save summary report
        with open(summary_path, "w") as f:
            json.dump(summary_report, f, indent=2)
        typer.echo(f"Summary report saved to: {summary_path}")
        
        # Generate extended analysis with residue information (only for PDB inputs)
        if is_pdb_format(input) and len(dataset) == 1:
            typer.echo("\n" + "="*70)
            typer.echo("GENERATING EXTENDED ANALYSIS")
            typer.echo("="*70)
            
            # Load the interaction scores
            with open(output, "r") as f:
                interact_data = json.loads(f.readline())
            
            scores = interact_data['cos_distances']
            block_indices = interact_data['block_idx']
            structure_id = Path(input).stem
            
            # Extract residue information from PDB
            try:
                residues = get_residue_info_from_pdb(input, block_indices, scores)
                
                # Sort by importance (lowest score = most critical)
                residues_sorted = sorted(residues, key=lambda x: x['atomica_score'])
                
                # Display top 10 critical residues
                typer.echo(f"\n🔴 TOP 10 MOST CRITICAL RESIDUES:")
                typer.echo(f"\n{'Rank':<6} {'Residue':<10} {'Chain':<7} {'Position':<10} {'ATOMICA_SCORE':<15} {'Importance Delta':<15}")
                typer.echo('-'*70)
                for i, res in enumerate(residues_sorted[:10], 1):
                    res_label = f"{res['res_name']}{res['res_id']}"
                    typer.echo(f"{i:<6} {res_label:<10} {res['chain_id']:<7} {res['res_id']:<10} {res['atomica_score']:.6f}      {res['importance_delta']:.4f}%")
                
                # Generate TSV report with all residues
                output_dir = output.parent
                report_file = output_dir / f"{structure_id}_critical_residues.tsv"
                tsv_text = generate_tsv_report(
                    residues_sorted, 
                    structure_id, 
                    interact_data['id'], 
                    scores
                )
                with open(report_file, 'w', encoding='utf-8') as f:
                    f.write(tsv_text)
                typer.echo(f"\n✅ Saved TSV report: {report_file} ({len(residues_sorted)} residues)")
                
                # Generate PyMOL commands
                pymol_file = output_dir / f"{structure_id}_pymol_commands.pml"
                pymol_cmds = generate_pymol_commands(residues_sorted, structure_id, top_n=10)
                with open(pymol_file, 'w', encoding='utf-8') as f:
                    f.write(pymol_cmds)
                typer.echo(f"✅ Saved PyMOL commands: {pymol_file}")
                
                typer.echo(f"\n💡 To visualize in PyMOL:")
                typer.echo(f"   1. Open structure: pymol {input}")
                typer.echo(f"   2. Run script: @{pymol_file}")
                typer.echo(f"   Or copy-paste commands from: {pymol_file}")
                typer.echo(f"\n📊 To analyze residues: Open {report_file} in Excel/LibreOffice")
                
            except Exception as e:
                typer.echo(f"⚠️  Could not generate extended analysis: {e}")
        
        typer.echo()
        
        success_fields = {
            "num_structures": len(dataset),
            "total_time_seconds": round(total_time, 2),
            "avg_time_per_structure": round(avg_time, 4),
            "summary_saved": True,
            "summary_path": str(summary_path),
            "device": device,
        }
        
        if not cpu:
            success_fields["avg_memory_mb"] = round(avg_memory, 2)
            success_fields["max_memory_mb"] = round(max_memory, 2)
        
        action.add_success_fields(**success_fields)
    finally:
        # Clean up temporary JSONL file if created
        if temp_jsonl_file and temp_jsonl_file.exists():
            temp_jsonl_file.unlink()


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()
            
