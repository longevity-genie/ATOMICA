"""
CLI for downloading PDB structures with optional metadata resolution.

This module provides a command-line interface for downloading PDB files and their metadata
using biotite for PDB downloading and the metadata resolution functions from atomica-mcp.
"""
from typing import List, Optional
from pathlib import Path
import sys

import typer
from eliot import start_action, Logger
from pycomfort.logging import to_nice_file, to_nice_stdout

# Import biotite for PDB downloading
import biotite.database.rcsb as rcsb

def fetch_pdb_metadata(pdb_id: str, timeout: int = 10, retries: int = 3, use_tsv: bool = True, **kwargs):
    """
    Fetch PDB metadata including resolution, title, chains, and organism information.

    Uses biotite's RCSB API for reliable metadata fetching.

    Args:
        pdb_id: PDB identifier (e.g., '1tgr')
        timeout: Request timeout in seconds
        retries: Number of retry attempts
        use_tsv: If True, use local TSV files for data (not implemented in fallback)

    Returns:
        Dictionary with PDB metadata including:
        - pdb_id: PDB identifier
        - found: Whether PDB was found
        - source: "biotite"
        - title: Structure title
        - resolution: PDB resolution
        - experimental_method: Experimental method used
        - entities: List of polymer entities with chains, organism, UniProt IDs
        - error: Error message if not found
    """
    try:
        import json

        # Fetch entry information using biotite
        query = rcsb.BasicQuery(pdb_id)
        results = rcsb.search(query)

        if not results:
            return {"pdb_id": pdb_id, "found": False, "error": "PDB ID not found"}

        # Fetch detailed entry information
        entry_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
        import urllib.request
        with urllib.request.urlopen(entry_url, timeout=timeout) as response:
            entry_data = json.loads(response.read().decode())

        # Extract basic metadata
        metadata = {
            "pdb_id": pdb_id,
            "found": True,
            "source": "biotite",
            "title": entry_data.get("struct", {}).get("title", ""),
            "resolution": None,
            "experimental_method": None,
            "entities": []
        }

        # Get resolution and experimental method
        if "rcsb_entry_info" in entry_data:
            metadata["resolution"] = entry_data["rcsb_entry_info"].get("resolution_combined", [])
            metadata["experimental_method"] = entry_data["rcsb_entry_info"].get("experimental_method", "")

        # Get polymer entity information
        polymer_entity_ids = entry_data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
        for entity_id in polymer_entity_ids:
            try:
                entity_url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
                with urllib.request.urlopen(entity_url, timeout=timeout) as response:
                    entity = json.loads(response.read().decode())

                # Extract organism information
                organism_info = {"scientific_name": "Unknown", "taxonomy_id": None}
                if "entity_src_gen" in entity and entity["entity_src_gen"]:
                    src = entity["entity_src_gen"][0]
                    organism_info = {
                        "scientific_name": src.get("pdbx_gene_src_scientific_name", "Unknown"),
                        "taxonomy_id": src.get("pdbx_gene_src_ncbi_taxonomy_id", None),
                    }
                elif "entity_src_nat" in entity and entity["entity_src_nat"]:
                    src = entity["entity_src_nat"][0]
                    organism_info = {
                        "scientific_name": src.get("pdbx_organism_scientific", "Unknown"),
                        "taxonomy_id": src.get("pdbx_ncbi_taxonomy_id", None),
                    }
                elif "rcsb_entity_source_organism" in entity and entity["rcsb_entity_source_organism"]:
                    src = entity["rcsb_entity_source_organism"][0]
                    organism_info = {
                        "scientific_name": src.get("ncbi_scientific_name", "Unknown"),
                        "taxonomy_id": src.get("ncbi_taxonomy_id", None),
                    }

                entity_info = {
                    "entity_id": entity_id,
                    "description": entity.get("rcsb_polymer_entity", {}).get("pdbx_description", ""),
                    "chains": [c.strip() for c in entity.get("entity_poly", {}).get("pdbx_strand_id", "").split(",") if c.strip()],
                    "type": entity.get("entity_poly", {}).get("type", ""),
                    "organism": organism_info,
                    "uniprot_ids": entity.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids", [])
                }
                metadata["entities"].append(entity_info)
            except Exception:
                # Skip entity errors
                pass

        return metadata

    except Exception as e:
        return {"pdb_id": pdb_id, "found": False, "error": f"Metadata fetch failed: {str(e)}"}


app = typer.Typer(help="Download PDB structures with optional metadata resolution", add_completion=False)


@app.command()
def download(
    pdb_ids: List[str] = typer.Argument(..., help="PDB IDs to download (e.g., 1abc 2def)"),
    output_dir: Path = typer.Option(Path("downloads/pdbs"), help="Output directory for PDB files"),
    file_format: str = typer.Option("cif", help="File format: 'pdb' or 'cif' (mmCIF recommended)"),
    include_metadata: bool = typer.Option(True, help="Download and include metadata for each PDB"),
    metadata_dir: Optional[Path] = typer.Option(None, help="Directory to save metadata files (defaults to output_dir)"),
    log_to_file: bool = typer.Option(False, help="Log detailed output to files in ./logs directory"),
    log_dir: Path = typer.Option(Path("logs"), help="Directory for log files (only used if --log-to-file is set)"),
    log_file_name: str = typer.Option("pdb_downloader", help="Base name for log files without extension"),
    timeout: int = typer.Option(10, help="Request timeout in seconds for API calls"),
    retries: int = typer.Option(3, help="Number of retry attempts for failed API calls"),
    use_tsv: bool = typer.Option(True, help="Use local TSV files for metadata (if available)"),
) -> None:
    """
    Download PDB structures with optional metadata resolution.

    Downloads PDB files to the specified output directory and optionally fetches
    metadata including resolution, title, chains, organisms, and UniProt IDs.

    Default downloads go to downloads/pdbs/ which is created if it doesn't exist.
    Metadata files are saved as JSON files alongside PDB files.
    """
    # Configure logging
    if log_to_file:
        log_dir.mkdir(exist_ok=True)
        json_log = log_dir / f"{log_file_name}.json"
        rendered_log = log_dir / f"{log_file_name}.log"
        to_nice_file(json_log, rendered_log)
        # Keep only the last added destination (file logging) if flag is set
        if len(Logger._destinations._destinations) > 1:
            Logger._destinations._destinations = [Logger._destinations._destinations[-1]]
    else:
        to_nice_stdout()
        # Keep only the last added destination (stdout logging) if flag is set
        if len(Logger._destinations._destinations) > 1:
            Logger._destinations._destinations = [Logger._destinations._destinations[-1]]

    # Set default metadata directory
    if metadata_dir is None:
        metadata_dir = output_dir

    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir.mkdir(parents=True, exist_ok=True)

    # Counters for summary
    total_downloaded = 0
    total_metadata = 0
    failed_downloads = []

    with start_action(
        action_type="pdb_download_batch",
        pdb_ids=pdb_ids,
        output_dir=str(output_dir),
        format=file_format,
        include_metadata=include_metadata,
        metadata_dir=str(metadata_dir)
    ) as action:
        action.log(
            message_type="download_config",
            output_dir=str(output_dir),
            format=file_format,
            include_metadata=include_metadata,
            metadata_dir=str(metadata_dir)
        )

        for pdb_id in pdb_ids:
            pdb_id = pdb_id.lower().strip()

            with start_action(action_type="download_pdb", pdb_id=pdb_id) as pdb_action:
                try:
                    # Download PDB structure
                    if file_format == "cif":
                        file_path = output_dir / f"{pdb_id}.cif"
                    else:
                        file_path = output_dir / f"{pdb_id}.pdb"

                    # Use biotite to download
                    downloaded_file = rcsb.fetch(pdb_id, file_format, target_path=str(file_path))
                    actual_path = Path(downloaded_file)

                    total_downloaded += 1
                    pdb_action.log(
                        message_type="pdb_downloaded",
                        pdb_id=pdb_id,
                        file_path=str(actual_path),
                        format=file_format
                    )

                    # Optionally download metadata
                    if include_metadata:
                        try:
                            metadata = fetch_pdb_metadata(
                                pdb_id,
                                timeout=timeout,
                                retries=retries,
                                use_tsv=use_tsv
                            )

                            # Save metadata as JSON
                            metadata_file = metadata_dir / f"{pdb_id}_metadata.json"
                            import json
                            with open(metadata_file, 'w') as f:
                                json.dump(metadata, f, indent=2)

                            total_metadata += 1
                            pdb_action.log(
                                message_type="metadata_downloaded",
                                pdb_id=pdb_id,
                                metadata_file=str(metadata_file),
                                source=metadata.get("source", "unknown")
                            )

                        except Exception as e:
                            pdb_action.log(
                                message_type="metadata_download_failed",
                                pdb_id=pdb_id,
                                error=str(e)
                            )

                except Exception as e:
                    failed_downloads.append(pdb_id)
                    pdb_action.log(
                        message_type="pdb_download_failed",
                        pdb_id=pdb_id,
                        error=str(e)
                    )

    # Summary
    action.log(
        message_type="download_summary",
        total_requested=len(pdb_ids),
        total_downloaded=total_downloaded,
        total_metadata=total_metadata,
        failed_downloads=failed_downloads
    )

    # User feedback
    typer.echo(f"✅ Downloaded {total_downloaded}/{len(pdb_ids)} PDB files")
    if include_metadata:
        typer.echo(f"📋 Downloaded metadata for {total_metadata} PDB files")

    if failed_downloads:
        typer.echo(f"❌ Failed to download: {', '.join(failed_downloads)}", err=True)

    if total_downloaded == 0:
        typer.echo("❌ No PDB files were successfully downloaded", err=True)
        raise typer.Exit(code=1)


@app.command()
def download_from_file(
    input_file: Path = typer.Argument(..., help="Text file containing PDB IDs (one per line)"),
    output_dir: Path = typer.Option(Path("downloads/pdbs"), help="Output directory for PDB files"),
    file_format: str = typer.Option("cif", help="File format: 'pdb' or 'cif' (mmCIF recommended)"),
    include_metadata: bool = typer.Option(True, help="Download and include metadata for each PDB"),
    metadata_dir: Optional[Path] = typer.Option(None, help="Directory to save metadata files (defaults to output_dir)"),
    log_to_file: bool = typer.Option(False, help="Log detailed output to files in ./logs directory"),
    log_dir: Path = typer.Option(Path("logs"), help="Directory for log files (only used if --log-to-file is set)"),
    log_file_name: str = typer.Option("pdb_downloader", help="Base name for log files without extension"),
    timeout: int = typer.Option(10, help="Request timeout in seconds for API calls"),
    retries: int = typer.Option(3, help="Number of retry attempts for failed API calls"),
    use_tsv: bool = typer.Option(True, help="Use local TSV files for metadata (if available)"),
) -> None:
    """
    Download PDB structures from a file containing PDB IDs.

    Each line in the input file should contain one PDB ID.
    Empty lines and lines starting with '#' are ignored.
    """
    # Read PDB IDs from file
    pdb_ids = []
    try:
        with open(input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):
                    pdb_ids.append(line)
    except Exception as e:
        typer.echo(f"❌ Error reading input file: {e}", err=True)
        raise typer.Exit(code=1)

    if not pdb_ids:
        typer.echo("❌ No PDB IDs found in input file", err=True)
        raise typer.Exit(code=1)

    typer.echo(f"📄 Found {len(pdb_ids)} PDB IDs in {input_file}")

    # Call the main download function
    download(
        pdb_ids=pdb_ids,
        output_dir=output_dir,
        file_format=file_format,
        include_metadata=include_metadata,
        metadata_dir=metadata_dir,
        log_to_file=log_to_file,
        log_dir=log_dir,
        log_file_name=log_file_name,
        timeout=timeout,
        retries=retries,
        use_tsv=use_tsv,
    )


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()
