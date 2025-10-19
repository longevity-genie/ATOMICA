from typing import Optional, List, Dict, Any, Union
from pathlib import Path
import json
import sys
import time
import warnings

# Suppress PyTorch deprecation warnings
warnings.filterwarnings('ignore', category=UserWarning, module='torch')

import typer
import torch
from tqdm import tqdm
from eliot import start_action
import polars as pl

from data.dataset import PDBDataset, ProtInterfaceDataset
from models.prediction_model import PredictionModel
from models.pretrain_model import DenoisePretrainModel
from models.prot_interface_model import ProteinInterfaceModel
from trainers.abs_trainer import Trainer


app = typer.Typer(help="Generate embeddings from ATOMICA models", add_completion=False)


def load_model(
    model_ckpt: Optional[Path] = None,
    model_config: Optional[Path] = None,
    model_weights: Optional[Path] = None
) -> Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel]:
    """Load model from checkpoint or config/weights.
    
    If no model is provided, attempts to load the default pretrain model from downloads/ATOMICA_checkpoints/.
    """
    with start_action(action_type="load_model") as action:
        if model_ckpt:
            with start_action(action_type="loading_checkpoint", path=str(model_ckpt)):
                model: Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel] = torch.load(model_ckpt)
        elif model_config and model_weights:
            with start_action(action_type="loading_from_config", config=str(model_config), weights=str(model_weights)):
                with open(model_config, "r") as f:
                    model_config_dict: Dict[str, Any] = json.load(f)
                
                model_type: str = model_config_dict['model_type']
                if model_type in ('PredictionModel', 'DenoisePretrainModel'):
                    model: Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel] = PredictionModel.load_from_config_and_weights(model_config, model_weights)
                elif model_type == 'ProteinInterfaceModel':
                    model: Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel] = ProteinInterfaceModel.load_from_config_and_weights(model_config, model_weights)
                else:
                    raise NotImplementedError(f"Model type {model_type} not implemented")
        else:
            # Try to load default pretrain model from downloads/ATOMICA_checkpoints
            default_pretrain_dir = Path("downloads/ATOMICA_checkpoints/pretrain")
            default_config = default_pretrain_dir / "pretrain_model_config.json"
            default_weights = default_pretrain_dir / "pretrain_model_weights.pt"
            
            if default_config.exists() and default_weights.exists():
                with start_action(action_type="loading_default_model", config=str(default_config), weights=str(default_weights)):
                    with open(default_config, "r") as f:
                        model_config_dict: Dict[str, Any] = json.load(f)
                    model_type: str = model_config_dict.get('model_type', 'DenoisePretrainModel')
                    model: Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel] = PredictionModel.load_from_config_and_weights(str(default_config), str(default_weights))
            else:
                raise ValueError(f"Either model_ckpt or both model_config and model_weights must be provided. Also tried default model at {default_pretrain_dir} but it doesn't exist.")
        
        return model


def process_batch(
    model: Union[PredictionModel, DenoisePretrainModel],
    items: List[Dict[str, Any]],
    dataset: Union[PDBDataset, ProtInterfaceDataset],
    device: str = "cuda"
) -> List[Dict[str, Any]]:
    """Process a single batch of items to generate embeddings."""
    outputs: List[Dict[str, Any]] = []
    for item in items:
        outputs.append({"id": item["id"]})
    
    if isinstance(dataset, ProtInterfaceDataset):
        batch_items: List[Dict[str, Any]] = [item["prot_data"] for item in items]
    else:
        batch_items: List[Dict[str, Any]] = [item["data"] for item in items]
    
    batch: Dict[str, Any] = PDBDataset.collate_fn(batch_items)
    batch = Trainer.to_device(batch, device)
    return_obj: Any = model.infer(batch)
    
    curr_block: int = 0
    curr_atom: int = 0
    for i, item in enumerate(items):
        num_blocks: int = len(item["data"]["B"])
        num_atoms: int = len(item["data"]["A"])

        outputs[i]["graph_embedding"] = return_obj.graph_repr[i].detach().cpu().numpy()
        outputs[i]["block_embedding"] = return_obj.block_repr[curr_block: curr_block + num_blocks].detach().cpu().numpy()
        outputs[i]["atom_embedding"] = return_obj.unit_repr[curr_atom: curr_atom + num_atoms].detach().cpu().numpy()
        outputs[i]["block_id"] = item["data"]["B"]
        outputs[i]["atom_id"] = item["data"]["A"]

        curr_block += num_blocks
        curr_atom += num_atoms
    
    return outputs


def process_single_item(
    model: Union[PredictionModel, DenoisePretrainModel],
    item: Dict[str, Any],
    dataset: Union[PDBDataset, ProtInterfaceDataset],
    device: str = "cuda"
) -> Dict[str, Any]:
    """Process a single item to generate embeddings (fallback for OOM)."""
    output: Dict[str, Any] = {"id": item["id"]}
    data_item: Dict[str, Any] = item["prot_data"] if isinstance(dataset, ProtInterfaceDataset) else item["data"]
    batch: Dict[str, Any] = PDBDataset.collate_fn([data_item])
    batch = Trainer.to_device(batch, device)
    return_obj: Any = model.infer(batch)
    
    output["graph_embedding"] = return_obj.graph_repr[0].detach().cpu().numpy()
    output["block_embedding"] = return_obj.block_repr.detach().cpu().numpy()
    output["atom_embedding"] = return_obj.unit_repr.detach().cpu().numpy()
    output["block_id"] = item["data"]["B"]
    output["atom_id"] = item["data"]["A"]
    
    return output


def embeddings_generator(
    model: Union[PredictionModel, DenoisePretrainModel],
    dataset: Union[PDBDataset, ProtInterfaceDataset],
    batch_size: int,
    device: str
):
    """Generator that yields embeddings one at a time (memory efficient).
    
    Args:
        model: ATOMICA model for generating embeddings
        dataset: Dataset to process
        batch_size: Batch size for model inference
        device: Device for inference
        
    Yields:
        Dictionary with embedding data (with lists instead of numpy arrays)
    """
    import numpy as np
    
    total_items: int = len(dataset)
    total_batches: int = (total_items + batch_size - 1) // batch_size
    
    for batch_idx, idx in enumerate(tqdm(range(0, total_items, batch_size), desc="Processing batches", total=total_batches)):
        items: List[Dict[str, Any]] = dataset.data[idx:min(idx + batch_size, total_items)]
        
        with start_action(action_type="process_batch", batch_idx=batch_idx, batch_size=len(items), global_idx=idx):
            # Clear CUDA cache before processing batch
            if device == "cuda":
                torch.cuda.empty_cache()
            
            # Process batch
            try:
                outputs: List[Dict[str, Any]] = process_batch(model, items, dataset, device)
            except RuntimeError as e:
                if "CUDA out of memory" in str(e):
                    torch.cuda.empty_cache()
                    with start_action(action_type="oom_fallback", info="CUDA out of memory, processing items one by one"):
                        outputs = []
                        for item in items:
                            try:
                                output: Dict[str, Any] = process_single_item(model, item, dataset, device)
                                outputs.append(output)
                                # Clear after each single item in OOM fallback
                                if device == "cuda":
                                    torch.cuda.empty_cache()
                            except Exception as item_error:
                                with start_action(action_type="item_error", item_id=item['id'], error=str(item_error)):
                                    if device == "cuda":
                                        torch.cuda.empty_cache()
                else:
                    raise
            
            # Clear CUDA cache after processing batch (important!)
            if device == "cuda":
                torch.cuda.empty_cache()
            
            # Yield embeddings one at a time with numpy arrays converted to lists
            for emb in outputs:
                converted = {}
                for key, value in emb.items():
                    if isinstance(value, np.ndarray):
                        converted[key] = value.tolist()
                    else:
                        converted[key] = value
                yield converted


def stream_embeddings_to_parquet(
    model: Union[PredictionModel, DenoisePretrainModel],
    dataset: Union[PDBDataset, ProtInterfaceDataset],
    output_path: Path,
    batch_size: int = 4,
    device: str = "cuda",
    compression: str = "snappy",
    chunk_size: int = 1000
) -> int:
    """Stream embeddings directly to parquet file using Polars lazy frames (memory efficient).
    
    Collects lazy frames without materializing them, then uses sink_parquet for a single
    streaming write operation. Perfect for huge datasets.
    
    Memory Management:
    - batch_size controls CUDA/GPU memory (how many items on GPU at once)
    - chunk_size controls lazy frame accumulation (doesn't affect memory much)
    - CUDA cache is cleared before AND after each batch to minimize GPU memory usage
    - Polars handles streaming internally - data never fully loaded into memory
    
    Args:
        model: ATOMICA model for generating embeddings
        dataset: Dataset to process
        output_path: Path to save parquet file
        batch_size: Batch size for model inference (affects CUDA memory)
        device: Device for inference
        compression: Compression codec ('snappy', 'zstd', 'lz4', 'gzip', 'brotli', 'uncompressed')
        chunk_size: Number of embeddings per lazy frame (affects write granularity)
    
    Returns:
        Total number of embeddings processed
    """
    with start_action(
        action_type="stream_embeddings_to_parquet",
        path=str(output_path),
        chunk_size=chunk_size,
        compression=compression
    ) as action:
        # Use generator to avoid loading all embeddings into memory
        gen = embeddings_generator(model, dataset, batch_size, device)
        
        # Collect lazy frames (these are not materialized yet!)
        lazy_frames: List[pl.LazyFrame] = []
        embeddings_chunk: List[Dict[str, Any]] = []
        processed_items: int = 0
        
        for emb in gen:
            embeddings_chunk.append(emb)
            processed_items += 1
            
            # Create lazy frame when chunk reaches chunk_size
            if len(embeddings_chunk) >= chunk_size:
                df_chunk = pl.from_dicts(embeddings_chunk)
                lazy_frames.append(df_chunk.lazy())
                
                with start_action(action_type="chunk_collected", embeddings_in_chunk=len(embeddings_chunk), total_processed=processed_items, total_lazy_frames=len(lazy_frames)):
                    pass
                
                embeddings_chunk = []  # Clear chunk from CPU memory
        
        # Add any remaining embeddings as a lazy frame
        if embeddings_chunk:
            df_chunk = pl.from_dicts(embeddings_chunk)
            lazy_frames.append(df_chunk.lazy())
            
            with start_action(action_type="final_chunk_collected", embeddings_in_chunk=len(embeddings_chunk), total_lazy_frames=len(lazy_frames)):
                pass
        
        # Now concat all lazy frames and write once with sink_parquet (streaming write!)
        if lazy_frames:
            with start_action(action_type="streaming_write", total_lazy_frames=len(lazy_frames), total_embeddings=processed_items):
                if len(lazy_frames) == 1:
                    # Single frame, just sink it
                    lazy_frames[0].sink_parquet(output_path, compression=compression)
                else:
                    # Multiple frames, concat lazily and sink
                    combined = pl.concat(lazy_frames, how="vertical")
                    combined.sink_parquet(output_path, compression=compression)
        
        return processed_items


@app.command()
def main(
    data_path: Path = typer.Option(..., help="Path to the data file either in json or pickle format"),
    output_path: Path = typer.Option(..., help="Path to save the output embeddings (parquet format)"),
    model_ckpt: Optional[Path] = typer.Option(None, help="Path of the model checkpoint to load"),
    model_config: Optional[Path] = typer.Option(None, help="Path of the model config to load"),
    model_weights: Optional[Path] = typer.Option(None, help="Path of the model weights to load"),
    batch_size: int = typer.Option(4, help="Batch size for GPU/CUDA processing (controls GPU memory usage)"),
    device: str = typer.Option("cuda", help="Device to use for inference"),
    start_line: int = typer.Option(0, help="Starting line index for slicing JSONL data (0-based)"),
    num_lines: Optional[int] = typer.Option(None, help="Maximum number of lines to process (None = all remaining)"),
    chunk_size: Optional[int] = typer.Option(None, help="Number of embeddings to accumulate before writing (default: same as batch_size for max memory efficiency)"),
    compression: str = typer.Option("snappy", help="Parquet compression codec: 'snappy', 'zstd', 'lz4', 'gzip', 'brotli', 'uncompressed'"),
) -> None:
    """Generate embeddings from ATOMICA models for input data.
    
    Memory Efficient Streaming:
    - batch_size controls GPU memory: only one batch on GPU at a time
    - chunk_size controls CPU memory: embeddings accumulated in RAM before writing to disk
    - Default: chunk_size = batch_size (write after each batch, since CUDA is the bottleneck)
    - CUDA cache is cleared before AND after each batch to minimize GPU memory
    - If OOM occurs, falls back to processing items one-by-one
    - Parquet output with compression (snappy by default) for optimal storage
    
    Example:
        # Process with batch size 4 and write after each batch (most memory efficient)
        python get_embeddings.py --data-path proteins.json --output-path embeddings.parquet --batch-size 4
        
        # Or accumulate 100 embeddings before writing (if disk I/O is slow)
        python get_embeddings.py --data-path proteins.json --output-path embeddings.parquet --batch-size 4 --chunk-size 100
    """
    start_time: float = time.time()
    
    # Default chunk_size to batch_size for maximum memory efficiency
    if chunk_size is None:
        chunk_size = batch_size
    
    with start_action(action_type="generate_embeddings", data_path=str(data_path), output_path=str(output_path), start_line=start_line, num_lines=num_lines, batch_size=batch_size, chunk_size=chunk_size) as action:
        # Load model
        with start_action(action_type="load_model"):
            model: Union[PredictionModel, ProteinInterfaceModel, DenoisePretrainModel] = load_model(model_ckpt, model_config, model_weights)
        
        # Prepare dataset
        with start_action(action_type="prepare_dataset", data_file=str(data_path), start_idx=start_line, num_lines=num_lines):
            if isinstance(model, ProteinInterfaceModel):
                with start_action(action_type="extracting_prot_model"):
                    model = model.prot_model
                dataset: Union[PDBDataset, ProtInterfaceDataset] = ProtInterfaceDataset(str(data_path), start_idx=start_line, num_lines=num_lines)
            else:
                dataset: Union[PDBDataset, ProtInterfaceDataset] = PDBDataset(str(data_path), start_idx=start_line, num_lines=num_lines)
            
            # Convert DenoisePretrainModel to PredictionModel if needed
            if isinstance(model, DenoisePretrainModel) and not isinstance(model, PredictionModel):
                model = PredictionModel.load_from_pretrained(model_ckpt)
            
            model = model.to(device)
        
        # Stream embeddings to parquet (memory efficient)
        with start_action(action_type="streaming_mode", chunk_size=chunk_size, compression=compression):
            processed_items = stream_embeddings_to_parquet(
                model=model,
                dataset=dataset,
                output_path=output_path,
                batch_size=batch_size,
                device=device,
                compression=compression,
                chunk_size=chunk_size
            )
        
        # Calculate and log final benchmarks
        total_elapsed: float = time.time() - start_time
        total_seconds_per_embedding: float = total_elapsed / processed_items if processed_items > 0 else 0
        
        with start_action(
            action_type="benchmark_summary",
            total_embeddings=processed_items,
            total_time_seconds=round(total_elapsed, 2),
            seconds_per_embedding=round(total_seconds_per_embedding, 4),
            throughput_embeddings_per_second=round(processed_items / total_elapsed, 2) if total_elapsed > 0 else 0
        ):
            pass


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()