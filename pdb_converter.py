from typing import Optional, List, Dict, Any, Union
from pathlib import Path
import sys
import gzip
import orjson
import typer
from eliot import start_action
from tqdm import tqdm

from data.dataset import blocks_to_data
from data.converter.pdb_to_list_blocks import pdb_to_list_blocks


app = typer.Typer(help="Convert PDB/CIF files to JSONL format for ATOMICA", add_completion=False)


def pdb_to_jsonl_item(
    pdb_file: Union[str, Path],
    pdb_id: Optional[str] = None,
    chains: Optional[List[str]] = None,
    return_blocks: bool = False
) -> Dict[str, Any]:
    """Convert a single PDB/CIF file to JSONL item format.
    
    Args:
        pdb_file: Path to PDB or CIF file
        pdb_id: Identifier for this structure (default: filename)
        chains: List of chain IDs to process (None = all chains)
        return_blocks: If True, also return the blocks (for debugging)
    
    Returns:
        Dictionary with 'data', 'id', and optionally 'block_to_pdb_indexes'
    """
    with start_action(action_type="pdb_to_jsonl_item", pdb_file=str(pdb_file), chains=chains):
        pdb_file = Path(pdb_file)
        
        # Use filename as ID if not provided
        if pdb_id is None:
            pdb_id = pdb_file.stem
        
        # Extract blocks from PDB
        blocks, pdb_indexes = pdb_to_list_blocks(
            str(pdb_file),
            selected_chains=chains,
            return_indexes=True,
            use_model=0
        )
        
        # Flatten if single chain or treat all as one entity
        if isinstance(blocks[0], list) and len(blocks) == 1:
            blocks = blocks[0]
            pdb_indexes = pdb_indexes[0]
        elif isinstance(blocks[0], list):
            # Multiple chains - concatenate them
            blocks = sum(blocks, [])
            pdb_indexes = sum(pdb_indexes, [])
        
        if len(blocks) == 0:
            raise ValueError(f"No blocks extracted from {pdb_file}")
        
        # Convert blocks to data format
        data = blocks_to_data(blocks)
        
        # Create mapping from block index to PDB indexes
        pdb_indexes_map = dict(
            zip(range(1, len(blocks) + 1), pdb_indexes)  # +1 for global block
        )
        
        chain_str = "_".join(chains) if chains else "all"
        item = {
            "data": data,
            "block_to_pdb_indexes": pdb_indexes_map,
            "id": f"{pdb_id}_{chain_str}",
        }
        
        if return_blocks:
            item["_blocks"] = blocks  # For debugging
        
        return item


def convert_pdb_directory(
    input_dir: Path,
    output_file: Path,
    pattern: str = "*.pdb",
    chains: Optional[List[str]] = None,
    compress: bool = False
) -> int:
    """Convert all PDB files in a directory to a single JSONL file.
    
    Args:
        input_dir: Directory containing PDB files
        output_file: Output JSONL file path
        pattern: Glob pattern for PDB files (e.g., '*.pdb', '*.cif')
        chains: List of chain IDs to process (None = all chains)
        compress: If True, gzip compress the output
    
    Returns:
        Number of structures processed
    """
    with start_action(
        action_type="convert_pdb_directory",
        input_dir=str(input_dir),
        output_file=str(output_file),
        pattern=pattern
    ):
        pdb_files = list(input_dir.glob(pattern))
        
        if len(pdb_files) == 0:
            raise ValueError(f"No files matching '{pattern}' found in {input_dir}")
        
        # Open output file (optionally compressed)
        if compress or str(output_file).endswith('.gz'):
            f_out = gzip.open(output_file, 'wt')
        else:
            f_out = open(output_file, 'w')
        
        processed_count = 0
        
        try:
            for pdb_file in tqdm(pdb_files, desc="Converting PDBs"):
                try:
                    item = pdb_to_jsonl_item(pdb_file, chains=chains)
                    f_out.write(orjson.dumps(item).decode('utf-8') + '\n')
                    processed_count += 1
                except Exception as e:
                    with start_action(action_type="conversion_error", file=str(pdb_file), error=str(e)):
                        pass
        finally:
            f_out.close()
        
        return processed_count


@app.command()
def convert_file(
    pdb_file: Path = typer.Argument(..., help="PDB or CIF file to convert"),
    output: Path = typer.Option(..., "--output", "-o", help="Output JSONL file"),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Comma-separated chain IDs (e.g., 'A,B')"),
    pdb_id: Optional[str] = typer.Option(None, "--id", help="Custom structure ID (default: filename)"),
) -> None:
    """Convert a single PDB/CIF file to JSONL format.
    
    Example:
        python pdb_converter.py convert-file protein.pdb -o protein.jsonl
        python pdb_converter.py convert-file protein.cif -o protein.jsonl --chains A,B
    """
    # Parse chains
    chain_list: Optional[List[str]] = None
    if chains:
        chain_list = [c.strip() for c in chains.split(",")]
    
    typer.echo(f"Converting {pdb_file}...")
    
    # Convert
    item = pdb_to_jsonl_item(pdb_file, pdb_id=pdb_id, chains=chain_list)
    
    # Write output
    with open(output, 'w') as f:
        f.write(orjson.dumps(item).decode('utf-8') + '\n')
    
    typer.echo(f"✓ Converted to {output}")
    typer.echo(f"  ID: {item['id']}")
    typer.echo(f"  Blocks: {len(item['data']['B'])}")
    typer.echo(f"  Atoms: {len(item['data']['A'])}")


@app.command()
def convert_directory(
    input_dir: Path = typer.Argument(..., help="Directory containing PDB/CIF files"),
    output: Path = typer.Option(..., "--output", "-o", help="Output JSONL file"),
    pattern: str = typer.Option("*.pdb", "--pattern", "-p", help="File pattern (e.g., '*.pdb', '*.cif')"),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Comma-separated chain IDs (e.g., 'A,B')"),
    compress: bool = typer.Option(False, "--compress", help="Compress output with gzip"),
) -> None:
    """Convert all PDB/CIF files in a directory to a single JSONL file.
    
    Example:
        python pdb_converter.py convert-directory ./pdbs/ -o all_proteins.jsonl
        python pdb_converter.py convert-directory ./pdbs/ -o all_proteins.jsonl.gz --compress
        python pdb_converter.py convert-directory ./cifs/ -o all_proteins.jsonl --pattern "*.cif"
    """
    # Parse chains
    chain_list: Optional[List[str]] = None
    if chains:
        chain_list = [c.strip() for c in chains.split(",")]
    
    typer.echo(f"Converting PDB files from {input_dir}...")
    typer.echo(f"Pattern: {pattern}")
    
    # Convert
    count = convert_pdb_directory(input_dir, output, pattern=pattern, chains=chain_list, compress=compress)
    
    typer.echo(f"\n✓ Converted {count} structures to {output}")


@app.command()
def batch_convert(
    file_list: Path = typer.Argument(..., help="Text file with PDB file paths (one per line)"),
    output: Path = typer.Option(..., "--output", "-o", help="Output JSONL file"),
    chains: Optional[str] = typer.Option(None, "--chains", "-c", help="Comma-separated chain IDs (e.g., 'A,B')"),
    compress: bool = typer.Option(False, "--compress", help="Compress output with gzip"),
) -> None:
    """Convert PDB files listed in a text file to JSONL format.
    
    Example:
        # Create list: ls pdbs/*.pdb > file_list.txt
        python pdb_converter.py batch-convert file_list.txt -o all_proteins.jsonl
    """
    # Parse chains
    chain_list: Optional[List[str]] = None
    if chains:
        chain_list = [c.strip() for c in chains.split(",")]
    
    # Read file list
    with open(file_list) as f:
        pdb_files = [line.strip() for line in f if line.strip()]
    
    typer.echo(f"Converting {len(pdb_files)} PDB files...")
    
    # Open output file (optionally compressed)
    if compress or str(output).endswith('.gz'):
        f_out = gzip.open(output, 'wt')
    else:
        f_out = open(output, 'w')
    
    processed_count = 0
    
    try:
        for pdb_file in tqdm(pdb_files, desc="Converting"):
            try:
                item = pdb_to_jsonl_item(pdb_file, chains=chain_list)
                f_out.write(orjson.dumps(item).decode('utf-8') + '\n')
                processed_count += 1
            except Exception as e:
                typer.echo(f"✗ Error with {pdb_file}: {e}", err=True)
    finally:
        f_out.close()
    
    typer.echo(f"\n✓ Converted {processed_count}/{len(pdb_files)} structures to {output}")


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()

