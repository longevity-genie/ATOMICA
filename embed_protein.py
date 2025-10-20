from typing import Optional, List, Dict, Any
from pathlib import Path
import tempfile
import json
import sys
import urllib.request
import urllib.error
import gzip
import shutil
import warnings

# Suppress PyTorch deprecation warnings
warnings.filterwarnings('ignore', category=UserWarning, module='torch')

import typer
import torch
import polars as pl
import numpy as np
from eliot import start_action, to_file
import biotite.structure.io.pdb as pdb
import biotite.structure.io.pdbx as pdbx
import biotite.database.rcsb as rcsb

from data.dataset import PDBDataset, ProtInterfaceDataset, blocks_to_data
from data.converter.pdb_to_list_blocks import pdb_to_list_blocks
from models.prediction_model import PredictionModel
from models.prot_interface_model import ProteinInterfaceModel
from trainers.abs_trainer import Trainer
from pdb_converter import pdb_to_jsonl_item


app = typer.Typer(
    help="Download PDB structure, process, and generate ATOMICA embeddings",
    add_completion=False
)


def download_alphafold_structure(
    uniprot_id: str,
    output_dir: Path,
    version: int = 4,
    use_pdb: bool = False
) -> Path:
    """Download AlphaFold predicted structure from AlphaFold DB (EBI).
    
    Note: Biotite does not have native AlphaFold DB support, so we use direct HTTP download.
    AlphaFold DB URL format: https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v{version}.{format}
    
    Args:
        uniprot_id: UniProt accession ID (e.g., 'P12345', 'Q9Y6K9')
        output_dir: Directory to save the downloaded file
        version: AlphaFold DB version (2, 3, or 4; default: 4)
        use_pdb: If True, download PDB format; if False, download mmCIF (recommended)
    
    Returns:
        Path to the downloaded structure file
    
    Raises:
        ValueError: If structure not found in AlphaFold DB
        
    Examples:
        >>> download_alphafold_structure("P12345", Path("/tmp"))
        Path('/tmp/AF-P12345-F1-model_v4.cif')
    """
    with start_action(action_type="download_alphafold", uniprot_id=uniprot_id, version=version):
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # AlphaFold DB URL format (EBI)
        file_format = "pdb" if use_pdb else "cif"
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.{file_format}"
        file_path = output_dir / f"AF-{uniprot_id}-F1-model_v{version}.{file_format}"
        
        try:
            # Download with better error handling
            with urllib.request.urlopen(url) as response:
                with open(file_path, 'wb') as out_file:
                    shutil.copyfileobj(response, out_file)
            return file_path
        except urllib.error.HTTPError as e:
            if e.code == 404:
                raise ValueError(
                    f"AlphaFold structure not found for UniProt ID '{uniprot_id}' (version {version}). "
                    f"Please check: 1) UniProt ID is valid, 2) Structure exists at https://alphafold.ebi.ac.uk/entry/{uniprot_id}"
                )
            elif e.code == 503:
                raise ValueError(f"AlphaFold DB service temporarily unavailable (503). Try again later.")
            else:
                raise ValueError(f"Failed to download from AlphaFold DB: HTTP {e.code} - {e.reason}")
        except urllib.error.URLError as e:
            raise ValueError(f"Network error downloading from AlphaFold DB: {e.reason}")


def download_pdb_structure(
    pdb_id: str,
    output_dir: Path,
    file_format: str = "cif",
    alphafold_fallback: bool = False
) -> tuple[Path, str]:
    """Download PDB structure from RCSB using biotite, with optional AlphaFold fallback.
    
    Args:
        pdb_id: 4-letter PDB ID or UniProt accession (if alphafold_fallback=True)
        output_dir: Directory to save the downloaded file
        file_format: 'pdb' or 'cif' (mmCIF format recommended)
        alphafold_fallback: If True and PDB download fails, try AlphaFold DB
    
    Returns:
        Tuple of (Path to the downloaded structure file, source name: 'rcsb' or 'alphafold')
    """
    with start_action(action_type="download_structure", pdb_id=pdb_id, format=file_format, alphafold_fallback=alphafold_fallback):
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # First, try RCSB PDB
        try:
            if file_format == "cif":
                file_path = output_dir / f"{pdb_id}.cif"
                cif_file = rcsb.fetch(pdb_id, "cif", target_path=str(file_path))
                return Path(cif_file), "rcsb"
            else:
                file_path = output_dir / f"{pdb_id}.pdb"
                pdb_file = rcsb.fetch(pdb_id, "pdb", target_path=str(file_path))
                return Path(pdb_file), "rcsb"
        except Exception as e:
            if not alphafold_fallback:
                raise
            
            # Try AlphaFold DB as fallback
            with start_action(action_type="alphafold_fallback", error=str(e)):
                try:
                    # Assume pdb_id is actually a UniProt ID if RCSB failed
                    af_file = download_alphafold_structure(pdb_id, output_dir)
                    return af_file, "alphafold"
                except Exception as af_error:
                    # Both failed, raise original error
                    raise ValueError(
                        f"Failed to download from RCSB PDB (error: {e}) "
                        f"and AlphaFold DB (error: {af_error}). "
                        f"For AlphaFold, provide a valid UniProt ID (e.g., P12345)."
                    )


def process_protein_structure(
    pdb_file: Path,
    pdb_id: str,
    chains: Optional[List[str]] = None,
    interface_dist_th: float = 8.0,
    fragmentation_method: Optional[str] = None
) -> Dict[str, Any]:
    """Process protein structure into ATOMICA format.
    
    Uses shared pdb_to_jsonl_item converter function.
    
    Args:
        pdb_file: Path to PDB/CIF file
        pdb_id: Identifier for this structure
        chains: List of chain IDs to process (None = all chains)
        interface_dist_th: Distance threshold for interface detection (unused, kept for compatibility)
        fragmentation_method: Fragmentation method for small molecules (unused, kept for compatibility)
    
    Returns:
        Dictionary with processed structure data
    """
    with start_action(
        action_type="process_structure",
        pdb_file=str(pdb_file),
        chains=chains
    ) as action:
        # Use shared converter function
        item = pdb_to_jsonl_item(pdb_file, pdb_id=pdb_id, chains=chains)
        action.log(message_type="processed_structure", num_blocks=len(item["data"]["B"]))
        return item


def generate_embeddings(
    processed_data: List[Dict[str, Any]],
    model_config: Path,
    model_weights: Path,
    device: str = "cuda",
    batch_size: int = 1
) -> List[Dict[str, Any]]:
    """Generate ATOMICA embeddings for processed structures.
    
    Args:
        processed_data: List of processed structure dictionaries
        model_config: Path to model config JSON
        model_weights: Path to model weights
        device: Device for inference
        batch_size: Batch size for processing
    
    Returns:
        List of embedding dictionaries
    """
    with start_action(
        action_type="generate_embeddings",
        num_structures=len(processed_data),
        device=device
    ) as action:
        # Load model
        with start_action(action_type="load_model"):
            with open(model_config, "r") as f:
                model_config_dict: Dict[str, Any] = json.load(f)
            
            model_type: str = model_config_dict.get("model_type", "PredictionModel")
            
            if model_type == "ProteinInterfaceModel":
                model: ProteinInterfaceModel = ProteinInterfaceModel.load_from_config_and_weights(
                    str(model_config), str(model_weights)
                )
                # Extract the protein-only model
                prot_model: PredictionModel = model.prot_model
            else:
                prot_model: PredictionModel = PredictionModel.load_from_config_and_weights(
                    str(model_config), str(model_weights)
                )
            
            prot_model = prot_model.to(device)
            prot_model.eval()
        
        embeddings: List[Dict[str, Any]] = []
        
        # Process in batches
        for idx in range(0, len(processed_data), batch_size):
            batch_items = processed_data[idx:min(idx + batch_size, len(processed_data))]
            
            with start_action(
                action_type="process_batch",
                batch_idx=idx // batch_size,
                batch_size=len(batch_items)
            ):
                # Prepare batch
                batch_data = [item["data"] for item in batch_items]
                batch = PDBDataset.collate_fn(batch_data)
                batch = Trainer.to_device(batch, device)
                
                # Generate embeddings
                with torch.no_grad():
                    return_obj = prot_model.infer(batch)
                
                # Extract embeddings for each item
                curr_block: int = 0
                curr_atom: int = 0
                
                for i, item in enumerate(batch_items):
                    num_blocks: int = len(item["data"]["B"])
                    num_atoms: int = len(item["data"]["A"])
                    
                    embedding = {
                        "id": item["id"],
                        "graph_embedding": return_obj.graph_repr[i].detach().cpu().numpy(),
                        "block_embedding": return_obj.block_repr[curr_block:curr_block + num_blocks].detach().cpu().numpy(),
                        "atom_embedding": return_obj.unit_repr[curr_atom:curr_atom + num_atoms].detach().cpu().numpy(),
                        "block_id": item["data"]["B"],
                        "atom_id": item["data"]["A"],
                    }
                    
                    embeddings.append(embedding)
                    
                    curr_block += num_blocks
                    curr_atom += num_atoms
        
        return embeddings


def save_embeddings_parquet(embeddings: List[Dict[str, Any]], output_path: Path, compression: str = "snappy") -> None:
    """Save embeddings to Parquet format with compression.
    
    Args:
        embeddings: List of embedding dictionaries
        output_path: Path to save parquet file
        compression: Compression codec ('snappy', 'zstd', 'lz4', 'gzip', 'brotli', 'uncompressed')
    """
    with start_action(action_type="save_embeddings", path=str(output_path), count=len(embeddings), compression=compression):
        # Convert numpy arrays to lists for parquet compatibility
        embeddings_converted = []
        for emb in embeddings:
            converted = {}
            for key, value in emb.items():
                if isinstance(value, np.ndarray):
                    converted[key] = value.tolist()
                else:
                    converted[key] = value
            embeddings_converted.append(converted)
        
        df = pl.from_records(embeddings_converted)
        df.write_parquet(output_path, compression=compression)


@app.command()
def from_pdb_id(
    pdb_id: str = typer.Argument(..., help="4-letter PDB ID or UniProt accession (with --alphafold-fallback)"),
    output_path: Path = typer.Option(..., "--output", "-o", help="Output path for embeddings parquet file"),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Comma-separated chain IDs (e.g., 'A,B'). If not specified, all chains are used."),
    model_config: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.json",
        "--model-config",
        help="Path to model config JSON"
    ),
    model_weights: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.pt",
        "--model-weights",
        help="Path to model weights"
    ),
    device: str = typer.Option("cuda", "--device", "-d", help="Device for inference (cuda/cpu)"),
    file_format: str = typer.Option("cif", "--format", "-f", help="File format to download (pdb/cif)"),
    alphafold_fallback: bool = typer.Option(False, "--alphafold-fallback", help="If PDB not found, try AlphaFold DB (provide UniProt ID)"),
    alphafold_version: int = typer.Option(4, "--alphafold-version", help="AlphaFold DB version (default: 4)"),
    temp_dir: Optional[Path] = typer.Option(None, "--temp-dir", help="Temporary directory for downloads (default: system temp)"),
    log_file: Optional[Path] = typer.Option(None, "--log-file", help="Path to save eliot logs"),
    compression: str = typer.Option("snappy", "--compression", help="Parquet compression codec: 'snappy', 'zstd', 'lz4', 'gzip', 'brotli', 'uncompressed'"),
) -> None:
    """Download a PDB structure, process it, and generate ATOMICA embeddings.
    
    Can fallback to AlphaFold predicted structures if PDB not available.
    
    Example:
        # Process PDB structure
        python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet
        
        # Try PDB first, fallback to AlphaFold if not found
        python embed_protein.py from-pdb-id P12345 -o embeddings.parquet --alphafold-fallback
        
        # Process specific chains
        python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet --chains A,B
        
        # Use AlphaFold structure directly (provide UniProt ID)
        python embed_protein.py from-pdb-id Q9Y6K9 -o embeddings.parquet --alphafold-fallback
    """
    if log_file:
        to_file(open(str(log_file), "w"))
    
    with start_action(action_type="embed_protein_from_id", identifier=pdb_id, alphafold_fallback=alphafold_fallback):
        # Parse chains
        chain_list: Optional[List[str]] = None
        if chains:
            chain_list = [c.strip() for c in chains.split(",")]
        
        # Create temp directory
        if temp_dir is None:
            temp_context = tempfile.TemporaryDirectory()
            temp_path = Path(temp_context.name)
        else:
            temp_path = temp_dir
            temp_path.mkdir(parents=True, exist_ok=True)
            temp_context = None
        
        try:
            # Step 1: Download structure
            typer.echo(f"📥 Downloading structure {pdb_id}...")
            pdb_file, source = download_pdb_structure(pdb_id, temp_path, file_format, alphafold_fallback)
            if source == "alphafold":
                typer.echo(f"✓ Downloaded AlphaFold structure to {pdb_file}")
            else:
                typer.echo(f"✓ Downloaded RCSB PDB structure to {pdb_file}")
            
            # Step 2: Process structure
            typer.echo(f"🔄 Processing structure...")
            processed = process_protein_structure(
                pdb_file,
                pdb_id,
                chains=chain_list
            )
            typer.echo(f"✓ Processed {len(processed['data']['B'])} blocks")
            
            # Step 3: Generate embeddings
            typer.echo(f"🧠 Generating embeddings with {device}...")
            embeddings = generate_embeddings(
                [processed],
                model_config,
                model_weights,
                device=device,
                batch_size=1
            )
            typer.echo(f"✓ Generated embeddings")
            
            # Step 4: Save to parquet
            typer.echo(f"💾 Saving to {output_path}...")
            save_embeddings_parquet(embeddings, output_path)
            typer.echo(f"✓ Saved embeddings to {output_path}")
            
            # Print summary
            emb = embeddings[0]
            typer.echo("\n" + "="*60)
            typer.echo("Summary:")
            typer.echo(f"  ID: {emb['id']}")
            typer.echo(f"  Graph embedding shape: {np.array(emb['graph_embedding']).shape}")
            typer.echo(f"  Block embeddings shape: {np.array(emb['block_embedding']).shape}")
            typer.echo(f"  Atom embeddings shape: {np.array(emb['atom_embedding']).shape}")
            typer.echo("="*60)
            
        finally:
            if temp_context:
                temp_context.cleanup()


@app.command()
def from_uniprot(
    uniprot_id: str = typer.Argument(..., help="UniProt accession ID (e.g., P12345, Q9Y6K9)"),
    output_path: Path = typer.Option(..., "--output", "-o", help="Output path for embeddings parquet file"),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Comma-separated chain IDs (e.g., 'A,B'). If not specified, all chains are used."),
    model_config: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.json",
        "--model-config",
        help="Path to model config JSON"
    ),
    model_weights: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.pt",
        "--model-weights",
        help="Path to model weights"
    ),
    device: str = typer.Option("cuda", "--device", "-d", help="Device for inference (cuda/cpu)"),
    alphafold_version: int = typer.Option(4, "--version", help="AlphaFold DB version (2, 3, or 4)"),
    temp_dir: Optional[Path] = typer.Option(None, "--temp-dir", help="Temporary directory for downloads (default: system temp)"),
    log_file: Optional[Path] = typer.Option(None, "--log-file", help="Path to save eliot logs"),
) -> None:
    """Download AlphaFold predicted structure and generate ATOMICA embeddings.
    
    Uses AlphaFold Database (https://alphafold.ebi.ac.uk/) for predicted structures.
    Perfect for proteins without experimental structures in PDB.
    
    Example:
        # Download AlphaFold structure for a UniProt ID
        python embed_protein.py from-uniprot P12345 -o embeddings.parquet
        
        # Use specific AlphaFold version
        python embed_protein.py from-uniprot Q9Y6K9 -o embeddings.parquet --version 3
        
        # Process specific chains
        python embed_protein.py from-uniprot P12345 -o embeddings.parquet --chains A
    """
    if log_file:
        to_file(open(str(log_file), "w"))
    
    with start_action(action_type="embed_protein_from_uniprot", uniprot_id=uniprot_id, alphafold_version=alphafold_version):
        # Parse chains
        chain_list: Optional[List[str]] = None
        if chains:
            chain_list = [c.strip() for c in chains.split(",")]
        
        # Create temp directory
        if temp_dir is None:
            temp_context = tempfile.TemporaryDirectory()
            temp_path = Path(temp_context.name)
        else:
            temp_path = temp_dir
            temp_path.mkdir(parents=True, exist_ok=True)
            temp_context = None
        
        try:
            # Step 1: Download AlphaFold structure
            typer.echo(f"📥 Downloading AlphaFold structure for {uniprot_id} (version {alphafold_version})...")
            af_file = download_alphafold_structure(uniprot_id, temp_path, version=alphafold_version)
            typer.echo(f"✓ Downloaded to {af_file}")
            
            # Step 2: Process structure
            typer.echo(f"🔄 Processing structure...")
            processed = process_protein_structure(
                af_file,
                uniprot_id,
                chains=chain_list
            )
            typer.echo(f"✓ Processed {len(processed['data']['B'])} blocks")
            
            # Step 3: Generate embeddings
            typer.echo(f"🧠 Generating embeddings with {device}...")
            embeddings = generate_embeddings(
                [processed],
                model_config,
                model_weights,
                device=device,
                batch_size=1
            )
            typer.echo(f"✓ Generated embeddings")
            
            # Step 4: Save to parquet
            typer.echo(f"💾 Saving to {output_path}...")
            save_embeddings_parquet(embeddings, output_path)
            typer.echo(f"✓ Saved embeddings to {output_path}")
            
            # Print summary
            emb = embeddings[0]
            typer.echo("\n" + "="*60)
            typer.echo("Summary:")
            typer.echo(f"  ID: {emb['id']}")
            typer.echo(f"  Source: AlphaFold DB v{alphafold_version}")
            typer.echo(f"  Graph embedding shape: {np.array(emb['graph_embedding']).shape}")
            typer.echo(f"  Block embeddings shape: {np.array(emb['block_embedding']).shape}")
            typer.echo(f"  Atom embeddings shape: {np.array(emb['atom_embedding']).shape}")
            typer.echo("="*60)
            
        finally:
            if temp_context:
                temp_context.cleanup()


@app.command()
def from_file(
    pdb_file: Path = typer.Argument(..., help="Path to PDB/CIF file"),
    output_path: Path = typer.Option(..., "--output", "-o", help="Output path for embeddings parquet file"),
    pdb_id: Optional[str] = typer.Option(None, "--id", help="Identifier for this structure (default: filename)"),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Comma-separated chain IDs (e.g., 'A,B'). If not specified, all chains are used."),
    model_config: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.json",
        "--model-config",
        help="Path to model config JSON"
    ),
    model_weights: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.pt",
        "--model-weights",
        help="Path to model weights"
    ),
    device: str = typer.Option("cuda", "--device", "-d", help="Device for inference (cuda/cpu)"),
    log_file: Optional[Path] = typer.Option(None, "--log-file", help="Path to save eliot logs"),
) -> None:
    """Process a local PDB/CIF file and generate ATOMICA embeddings.
    
    Example:
        # Process local PDB file
        python embed_protein.py from-file /path/to/protein.pdb -o embeddings.parquet
        
        # Process specific chains
        python embed_protein.py from-file /path/to/protein.cif -o embeddings.parquet --chains A
    """
    if log_file:
        to_file(open(str(log_file), "w"))
    
    # Use filename as ID if not provided
    if pdb_id is None:
        pdb_id = pdb_file.stem
    
    with start_action(action_type="embed_protein_from_file", pdb_file=str(pdb_file), pdb_id=pdb_id):
        # Parse chains
        chain_list: Optional[List[str]] = None
        if chains:
            chain_list = [c.strip() for c in chains.split(",")]
        
        # Step 1: Process structure
        typer.echo(f"🔄 Processing structure from {pdb_file}...")
        processed = process_protein_structure(
            pdb_file,
            pdb_id,
            chains=chain_list
        )
        typer.echo(f"✓ Processed {len(processed['data']['B'])} blocks")
        
        # Step 2: Generate embeddings
        typer.echo(f"🧠 Generating embeddings with {device}...")
        embeddings = generate_embeddings(
            [processed],
            model_config,
            model_weights,
            device=device,
            batch_size=1
        )
        typer.echo(f"✓ Generated embeddings")
        
        # Step 3: Save to parquet
        typer.echo(f"💾 Saving to {output_path}...")
        save_embeddings_parquet(embeddings, output_path)
        typer.echo(f"✓ Saved embeddings to {output_path}")
        
        # Print summary
        emb = embeddings[0]
        typer.echo("\n" + "="*60)
        typer.echo("Summary:")
        typer.echo(f"  ID: {emb['id']}")
        typer.echo(f"  Graph embedding shape: {np.array(emb['graph_embedding']).shape}")
        typer.echo(f"  Block embeddings shape: {np.array(emb['block_embedding']).shape}")
        typer.echo(f"  Atom embeddings shape: {np.array(emb['atom_embedding']).shape}")
        typer.echo("="*60)


@app.command()
def batch(
    input_csv: Path = typer.Argument(..., help="CSV file with columns: pdb_id, pdb_path (optional), chains (optional), source (optional)"),
    output_path: Path = typer.Option(..., "--output", "-o", help="Output path for embeddings parquet file"),
    model_config: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.json",
        "--model-config",
        help="Path to model config JSON"
    ),
    model_weights: Path = typer.Option(
        "downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.pt",
        "--model-weights",
        help="Path to model weights"
    ),
    device: str = typer.Option("cuda", "--device", "-d", help="Device for inference (cuda/cpu)"),
    batch_size: int = typer.Option(4, "--batch-size", "-b", help="Batch size for processing"),
    alphafold_fallback: bool = typer.Option(False, "--alphafold-fallback", help="If PDB not found, try AlphaFold DB"),
    temp_dir: Optional[Path] = typer.Option(None, "--temp-dir", help="Temporary directory for downloads"),
    log_file: Optional[Path] = typer.Option(None, "--log-file", help="Path to save eliot logs"),
    compression: str = typer.Option("snappy", "--compression", help="Parquet compression codec: 'snappy', 'zstd', 'lz4', 'gzip', 'brotli', 'uncompressed'"),
    chunk_size: Optional[int] = typer.Option(None, "--chunk-size", help="Number of proteins to process before writing (default: 1, write after each protein for max memory efficiency)"),
) -> None:
    """Process multiple proteins from a CSV file with memory-efficient streaming.
    
    CSV Format:
        pdb_id,pdb_path,chains,source
        1ABC,,A,rcsb
        P12345,,,alphafold
        2XYZ,/path/to/2xyz.pdb,A_B,
        3DEF,,,
    
    - pdb_id: PDB ID or UniProt accession (required)
    - pdb_path: Local file path (optional, if empty will download)
    - chains: Chain IDs separated by underscore (optional, if empty uses all)
    - source: 'rcsb', 'alphafold', or empty (optional, auto-detect or use --alphafold-fallback)
    
    Example:
        # Basic batch processing
        python embed_protein.py batch proteins.csv -o all_embeddings.parquet
        
        # With AlphaFold fallback for missing PDBs
        python embed_protein.py batch proteins.csv -o all_embeddings.parquet --alphafold-fallback
    """
    if log_file:
        to_file(open(str(log_file), "w"))
    
    # Default chunk_size to 1 for maximum memory efficiency (write after each protein)
    if chunk_size is None:
        chunk_size = 1
    
    with start_action(action_type="embed_proteins_batch", input_csv=str(input_csv), chunk_size=chunk_size):
        # Load model first (reuse for all proteins)
        with start_action(action_type="load_model"):
            with open(model_config, "r") as f:
                model_config_dict: Dict[str, Any] = json.load(f)
            
            model_type: str = model_config_dict.get("model_type", "PredictionModel")
            
            if model_type == "ProteinInterfaceModel":
                model: ProteinInterfaceModel = ProteinInterfaceModel.load_from_config_and_weights(
                    str(model_config), str(model_weights)
                )
                # Extract the protein-only model
                prot_model: PredictionModel = model.prot_model
            else:
                prot_model: PredictionModel = PredictionModel.load_from_config_and_weights(
                    str(model_config), str(model_weights)
                )
            
            prot_model = prot_model.to(device)
            prot_model.eval()
        
        # Read CSV lazily using scan_csv
        df_lazy = pl.scan_csv(input_csv)
        total_rows = df_lazy.select(pl.len()).collect().item()
        typer.echo(f"📋 Processing {total_rows} proteins from {input_csv}")
        
        # Create temp directory
        if temp_dir is None:
            temp_context = tempfile.TemporaryDirectory()
            temp_path = Path(temp_context.name)
        else:
            temp_path = temp_dir
            temp_path.mkdir(parents=True, exist_ok=True)
            temp_context = None
        
        try:
            # Collect lazy frames (not materialized yet!)
            lazy_frames: List[pl.LazyFrame] = []
            embeddings_chunk: List[Dict[str, Any]] = []
            
            # Process each protein (using lazy iteration)
            for row in df_lazy.collect().iter_rows(named=True):
                pdb_id = row["pdb_id"]
                pdb_path = row.get("pdb_path")
                chains_str = row.get("chains")
                source_hint = row.get("source", "").lower() if row.get("source") else None
                
                typer.echo(f"\n🔄 Processing {pdb_id}...")
                
                # Parse chains
                chain_list: Optional[List[str]] = None
                if chains_str and chains_str.strip():
                    chain_list = [c.strip() for c in chains_str.split("_")]
                
                # Get structure file
                if pdb_path and pdb_path.strip():
                    pdb_file = Path(pdb_path)
                    structure_source = "local"
                else:
                    # Download based on source hint or use fallback
                    if source_hint == "alphafold":
                        typer.echo(f"  📥 Downloading from AlphaFold DB...")
                        try:
                            pdb_file = download_alphafold_structure(pdb_id, temp_path)
                            structure_source = "alphafold"
                        except Exception as e:
                            typer.echo(f"  ✗ Error downloading from AlphaFold: {e}", err=True)
                            continue
                    else:
                        typer.echo(f"  📥 Downloading {pdb_id}...")
                        try:
                            pdb_file, structure_source = download_pdb_structure(
                                pdb_id, temp_path, "cif", alphafold_fallback
                            )
                            if structure_source == "alphafold":
                                typer.echo(f"  ℹ️  PDB not found, using AlphaFold structure")
                        except Exception as e:
                            typer.echo(f"  ✗ Error downloading {pdb_id}: {e}", err=True)
                            continue
                
                # Process structure
                try:
                    processed = process_protein_structure(
                        pdb_file,
                        pdb_id,
                        chains=chain_list
                    )
                    
                    # Generate embedding immediately and add to chunk (streaming mode)
                    batch_data = [processed["data"]]
                    batch = PDBDataset.collate_fn(batch_data)
                    batch = Trainer.to_device(batch, device)
                    
                    with torch.no_grad():
                        return_obj = prot_model.infer(batch)
                    
                    embedding = {
                        "id": processed["id"],
                        "graph_embedding": return_obj.graph_repr[0].detach().cpu().numpy().tolist(),
                        "block_embedding": return_obj.block_repr.detach().cpu().numpy().tolist(),
                        "atom_embedding": return_obj.unit_repr.detach().cpu().numpy().tolist(),
                        "block_id": processed["data"]["B"],
                        "atom_id": processed["data"]["A"],
                    }
                    
                    embeddings_chunk.append(embedding)
                    typer.echo(f"  ✓ Processed {len(processed['data']['B'])} blocks from {structure_source}")
                    
                    # Create lazy frame when chunk reaches chunk_size
                    if len(embeddings_chunk) >= chunk_size:
                        df_chunk = pl.from_dicts(embeddings_chunk)
                        lazy_frames.append(df_chunk.lazy())
                        typer.echo(f"  💾 Collected {len(embeddings_chunk)} embeddings (lazy frame {len(lazy_frames)})")
                        
                        embeddings_chunk = []  # Clear chunk from memory
                        
                        # Clear CUDA cache
                        if device == "cuda":
                            torch.cuda.empty_cache()
                        
                except Exception as e:
                    typer.echo(f"  ✗ Error processing {pdb_id}: {e}", err=True)
                    continue
            
            # Add any remaining embeddings as a lazy frame
            if embeddings_chunk:
                df_chunk = pl.from_dicts(embeddings_chunk)
                lazy_frames.append(df_chunk.lazy())
                typer.echo(f"  💾 Collected final {len(embeddings_chunk)} embeddings (lazy frame {len(lazy_frames)})")
            
            # Now write all lazy frames with a single sink_parquet (streaming write!)
            if lazy_frames:
                typer.echo(f"\n💾 Writing {len(lazy_frames)} lazy frames to {output_path}...")
                if len(lazy_frames) == 1:
                    # Single frame, just sink it
                    lazy_frames[0].sink_parquet(output_path, compression=compression)
                else:
                    # Multiple frames, concat lazily and sink
                    combined = pl.concat(lazy_frames, how="vertical")
                    combined.sink_parquet(output_path, compression=compression)
                
                # Read final count
                df_final = pl.read_parquet(output_path)
                total_embeddings = len(df_final)
                typer.echo(f"✓ Saved {total_embeddings} embeddings to {output_path}")
                
                # Print summary
                typer.echo("\n" + "="*60)
                typer.echo(f"Successfully processed {total_embeddings} proteins")
                typer.echo("="*60)
            else:
                typer.echo("❌ No structures were successfully processed!", err=True)
                raise typer.Exit(1)
            
        finally:
            if temp_context:
                temp_context.cleanup()


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()

