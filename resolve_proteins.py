from typing import Optional, List, Dict, Any, Set
from pathlib import Path
import gzip
import json
import sys
import csv
import time

import typer
from eliot import start_action, log_message
import biotite.database.rcsb as rcsb
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


app = typer.Typer(help="Resolve protein names from PDB IDs in JSONL.GZ files", add_completion=False)


def create_retry_session(
    retries: int = 3,
    backoff_factor: float = 0.5,
    timeout: int = 10
) -> requests.Session:
    """
    Create a requests session with retry logic and timeout.
    
    Args:
        retries: Number of retry attempts
        backoff_factor: Backoff factor for exponential retry delay
        timeout: Request timeout in seconds
    
    Returns:
        Configured requests Session
    """
    session = requests.Session()
    
    # Configure retry strategy
    retry_strategy = Retry(
        total=retries,
        backoff_factor=backoff_factor,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["HEAD", "GET", "OPTIONS"]
    )
    
    adapter = HTTPAdapter(max_retries=retry_strategy)
    session.mount("http://", adapter)
    session.mount("https://", adapter)
    
    return session


def parse_entry_id(entry_id: str) -> Dict[str, str]:
    """
    Parse entry ID to extract PDB ID and chain information.
    
    Format: PDB_ID_number_ChainA_ChainB
    Example: 2uxq_2_A_B -> {'pdb_id': '2uxq', 'chain1': 'A', 'chain2': 'B'}
    """
    parts = entry_id.split('_')
    if len(parts) >= 4:
        return {
            'pdb_id': parts[0].lower(),
            'chain1': parts[2],
            'chain2': parts[3]
        }
    elif len(parts) >= 3:
        return {
            'pdb_id': parts[0].lower(),
            'chain1': parts[2] if len(parts) > 2 else '',
            'chain2': ''
        }
    return {'pdb_id': parts[0].lower(), 'chain1': '', 'chain2': ''}


def fetch_pdb_metadata(pdb_id: str, timeout: int = 10, retries: int = 3) -> Dict[str, Any]:
    """
    Fetch PDB metadata including title, description, and chain information.
    
    Args:
        pdb_id: PDB identifier (e.g., '2uxq')
        timeout: Request timeout in seconds
        retries: Number of retry attempts
    
    Returns:
        Dictionary with PDB metadata
    """
    with start_action(action_type="fetch_pdb_metadata", pdb_id=pdb_id) as action:
        try:
            # Create retry session
            session = create_retry_session(retries=retries, timeout=timeout)
            
            # Use biotite to query RCSB PDB
            query = rcsb.BasicQuery(pdb_id)
            results = rcsb.search(query)
            
            if not results:
                action.log(message_type="pdb_not_found", pdb_id=pdb_id)
                return {"pdb_id": pdb_id, "found": False}
            
            # Fetch detailed information
            metadata = {
                "pdb_id": pdb_id,
                "found": True,
            }
            
            # Get structure title and details using fetch API
            try:
                # Fetch the summary from RCSB REST API with retry and timeout
                response = session.get(
                    f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}",
                    timeout=timeout
                )
                if response.status_code == 200:
                    data = response.json()
                    metadata["title"] = data.get("struct", {}).get("title", "")
                    metadata["description"] = data.get("struct", {}).get("pdbx_descriptor", "")
                    
                    # Get entity information (protein names)
                    if "rcsb_entry_info" in data:
                        metadata["experimental_method"] = data["rcsb_entry_info"].get("experimental_method", "")
                        metadata["resolution"] = data["rcsb_entry_info"].get("resolution_combined", [])
                
                # Get polymer entity IDs from entry data
                polymer_entity_ids = data.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
                
                # Fetch each polymer entity
                metadata["entities"] = []
                for entity_id in polymer_entity_ids:
                    try:
                        entity_response = session.get(
                            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}",
                            timeout=timeout
                        )
                        if entity_response.status_code == 200:
                            entity = entity_response.json()
                            chains_str = entity.get("entity_poly", {}).get("pdbx_strand_id", "")
                            chains = [c.strip() for c in chains_str.split(",")] if chains_str else []
                            
                            # Extract organism information
                            organism_info = {}
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
                            else:
                                organism_info = {
                                    "scientific_name": "Unknown",
                                    "taxonomy_id": None,
                                }
                            
                            # Extract UniProt IDs
                            uniprot_ids = entity.get("rcsb_polymer_entity_container_identifiers", {}).get("uniprot_ids", [])
                            
                            entity_info = {
                                "entity_id": entity_id,
                                "description": entity.get("rcsb_polymer_entity", {}).get("pdbx_description", ""),
                                "chains": chains,
                                "type": entity.get("entity_poly", {}).get("type", ""),
                                "organism": organism_info,
                                "uniprot_ids": uniprot_ids,
                            }
                            metadata["entities"].append(entity_info)
                            action.log(message_type="entity_parsed", entity_id=entity_id, 
                                      chains=chains, description=entity_info["description"],
                                      organism=organism_info.get("scientific_name"),
                                      uniprot_ids=uniprot_ids)
                    except Exception as entity_error:
                        action.log(message_type="entity_fetch_error", entity_id=entity_id, error=str(entity_error))
                        
            except requests.exceptions.Timeout as e:
                action.log(message_type="fetch_timeout", pdb_id=pdb_id, error=str(e))
                return {"pdb_id": pdb_id, "found": False, "error": f"Timeout after {timeout}s: {str(e)}"}
            except requests.exceptions.ConnectionError as e:
                action.log(message_type="connection_error", pdb_id=pdb_id, error=str(e))
                return {"pdb_id": pdb_id, "found": False, "error": f"Connection error: {str(e)}"}
            except requests.exceptions.RequestException as e:
                action.log(message_type="request_error", pdb_id=pdb_id, error=str(e))
                return {"pdb_id": pdb_id, "found": False, "error": f"Request error: {str(e)}"}
            except Exception as e:
                action.log(message_type="fetch_details_error", pdb_id=pdb_id, error=str(e))
            
            action.log(message_type="pdb_metadata_fetched", pdb_id=pdb_id)
            return metadata
            
        except requests.exceptions.Timeout as e:
            log_message(message_type="fetch_timeout", pdb_id=pdb_id, error=str(e))
            return {"pdb_id": pdb_id, "found": False, "error": f"Timeout after {timeout}s: {str(e)}"}
        except requests.exceptions.ConnectionError as e:
            log_message(message_type="connection_error", pdb_id=pdb_id, error=str(e))
            return {"pdb_id": pdb_id, "found": False, "error": f"Connection error: {str(e)}"}
        except Exception as e:
            log_message(message_type="fetch_error", pdb_id=pdb_id, error=str(e))
            return {"pdb_id": pdb_id, "found": False, "error": str(e)}


def classify_organism(taxonomy_id: Optional[int], scientific_name: str) -> str:
    """
    Classify organism type based on taxonomy ID or scientific name.
    
    Args:
        taxonomy_id: NCBI taxonomy ID
        scientific_name: Scientific name of organism
    
    Returns:
        Organism classification (e.g., 'Mammalian', 'Bacterial', 'Fungal', etc.)
    """
    scientific_name_lower = scientific_name.lower()
    
    # Common mammalian organisms
    mammals = ['homo sapiens', 'mus musculus', 'rattus norvegicus', 'bos taurus', 
               'sus scrofa', 'canis lupus', 'macaca', 'oryctolagus cuniculus', 'human', 'mouse', 'rat']
    
    # Common bacterial indicators
    bacteria_indicators = ['escherichia', 'staphylococcus', 'streptococcus', 'bacillus', 
                          'pseudomonas', 'salmonella', 'mycobacterium', 'clostridium',
                          'vibrio', 'yersinia', 'campylobacter', 'desulfotalea']
    
    # Common fungal organisms
    fungi = ['saccharomyces', 'candida', 'aspergillus', 'neurospora', 'yeast']
    
    # Viruses
    viruses = ['virus', 'phage', 'viral']
    
    # Plants
    plants = ['arabidopsis', 'nicotiana', 'solanum', 'oryza', 'zea mays']
    
    # Check mammals
    if any(mammal in scientific_name_lower for mammal in mammals):
        return "Mammalian"
    
    # Check bacteria
    if any(bacteria in scientific_name_lower for bacteria in bacteria_indicators):
        return "Bacterial"
    
    # Check fungi
    if any(fungus in scientific_name_lower for fungus in fungi):
        return "Fungal"
    
    # Check viruses
    if any(virus in scientific_name_lower for virus in viruses):
        return "Viral"
    
    # Check plants
    if any(plant in scientific_name_lower for plant in plants):
        return "Plant"
    
    # Use taxonomy ID ranges if available
    if taxonomy_id:
        try:
            tax_id_int = int(taxonomy_id)
            if 9604 <= tax_id_int <= 9606:  # Primates
                return "Mammalian"
            elif 40674 <= tax_id_int <= 40674:  # Mammalia
                return "Mammalian"
            elif tax_id_int == 2:  # Bacteria domain
                return "Bacterial"
            elif 2157 <= tax_id_int <= 2157:  # Archaea
                return "Archaeal"
            elif 4751 <= tax_id_int <= 4892:  # Fungi
                return "Fungal"
            elif 10239 <= tax_id_int <= 10239:  # Viruses
                return "Viral"
        except (ValueError, TypeError):
            pass  # If conversion fails, rely on name-based classification
    
    return "Other"


def get_chain_protein_name(metadata: Dict[str, Any], chain_id: str) -> str:
    """
    Get protein name for a specific chain from PDB metadata.
    
    Args:
        metadata: PDB metadata dictionary
        chain_id: Chain identifier (e.g., 'A', 'B')
    
    Returns:
        Protein name/description for the chain
    """
    if not metadata.get("found") or "entities" not in metadata:
        return "Unknown"
    
    for entity in metadata["entities"]:
        if chain_id in entity.get("chains", []):
            return entity.get("description", "Unknown")
    
    return "Unknown"


def get_chain_organism(metadata: Dict[str, Any], chain_id: str) -> Dict[str, Any]:
    """
    Get organism information for a specific chain from PDB metadata.
    
    Args:
        metadata: PDB metadata dictionary
        chain_id: Chain identifier (e.g., 'A', 'B')
    
    Returns:
        Dictionary with organism information including classification
    """
    if not metadata.get("found") or "entities" not in metadata:
        return {"scientific_name": "Unknown", "taxonomy_id": None, "classification": "Unknown"}
    
    for entity in metadata["entities"]:
        if chain_id in entity.get("chains", []):
            organism_info = entity.get("organism", {})
            scientific_name = organism_info.get("scientific_name", "Unknown")
            taxonomy_id = organism_info.get("taxonomy_id", None)
            
            return {
                "scientific_name": scientific_name,
                "taxonomy_id": taxonomy_id,
                "classification": classify_organism(taxonomy_id, scientific_name)
            }
    
    return {"scientific_name": "Unknown", "taxonomy_id": None, "classification": "Unknown"}


def get_chain_uniprot_ids(metadata: Dict[str, Any], chain_id: str) -> List[str]:
    """
    Get UniProt IDs for a specific chain from PDB metadata.
    
    Args:
        metadata: PDB metadata dictionary
        chain_id: Chain identifier (e.g., 'A', 'B')
    
    Returns:
        List of UniProt IDs for the chain
    """
    if not metadata.get("found") or "entities" not in metadata:
        return []
    
    for entity in metadata["entities"]:
        if chain_id in entity.get("chains", []):
            return entity.get("uniprot_ids", [])
    
    return []


def read_jsonl_gz_lines(file_path: Path, line_numbers: Set[int]) -> List[Dict[str, Any]]:
    """
    Read specific lines from a gzipped JSONL file.
    
    Args:
        file_path: Path to the .jsonl.gz file
        line_numbers: Set of line numbers to read (1-indexed)
    
    Returns:
        List of parsed JSON entries
    """
    entries = []
    with start_action(action_type="read_jsonl_gz", path=str(file_path)) as action:
        with gzip.open(file_path, 'rt', encoding='utf-8') as f:
            for line_num, line in enumerate(f, start=1):
                if line_num in line_numbers:
                    try:
                        entry = json.loads(line)
                        entries.append({"line_number": line_num, "entry": entry})
                        action.log(message_type="line_read", line_number=line_num, entry_id=entry.get("id"))
                    except json.JSONDecodeError as e:
                        action.log(message_type="json_decode_error", line_number=line_num, error=str(e))
                
                # Early exit if we've read all requested lines
                if len(entries) == len(line_numbers):
                    break
        
        action.log(message_type="read_complete", entries_read=len(entries))
    return entries


def matches_filter(result: Dict[str, Any], filter_organism: Optional[str], filter_classification: Optional[str]) -> bool:
    """
    Check if a result matches the organism or classification filter.
    
    Args:
        result: Result dictionary with chain_organisms information
        filter_organism: Organism name to filter (case-insensitive partial match)
        filter_classification: Classification to filter (exact match: Mammalian, Bacterial, etc.)
    
    Returns:
        True if result matches filter criteria
    """
    if not filter_organism and not filter_classification:
        return True
    
    chain_organisms = result.get("chain_organisms", {})
    
    for chain_key in ["chain1", "chain2"]:
        org_info = chain_organisms.get(chain_key, {})
        
        # Check organism filter (partial match, case-insensitive)
        if filter_organism:
            scientific_name = org_info.get("scientific_name", "").lower()
            if filter_organism.lower() in scientific_name:
                return True
        
        # Check classification filter (exact match)
        if filter_classification:
            classification = org_info.get("classification", "")
            if classification.lower() == filter_classification.lower():
                return True
    
    return False


def write_jsonl_output(results: List[Dict[str, Any]], original_entries: List[Dict[str, Any]], output_path: Path) -> None:
    """
    Write filtered results back to JSONL format with original entry data.
    
    Args:
        results: List of result dictionaries
        original_entries: List of original entries from input file
        output_path: Path to write JSONL output
    """
    with start_action(action_type="write_jsonl", output_path=str(output_path)) as action:
        # Create mapping from entry_id to original entry
        entry_map = {e["entry"]["id"]: e["entry"] for e in original_entries}
        
        with gzip.open(output_path, 'wt', encoding='utf-8') if output_path.suffix == '.gz' else open(output_path, 'w') as f:
            for result in results:
                entry_id = result["entry_id"]
                if entry_id in entry_map:
                    # Write original entry data
                    json.dump(entry_map[entry_id], f)
                    f.write('\n')
        
        action.log(message_type="jsonl_written", path=str(output_path), count=len(results))


def write_csv_summary(results: List[Dict[str, Any]], output_path: Path) -> None:
    """
    Write CSV summary with protein and species information.
    
    Args:
        results: List of result dictionaries
        output_path: Path to write CSV output
    """
    with start_action(action_type="write_csv", output_path=str(output_path)) as action:
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            writer.writerow([
                'entry_id', 'pdb_id', 'chain_id', 'protein_name', 
                'organism', 'taxonomy_id', 'classification', 'uniprot_ids'
            ])
            
            for result in results:
                entry_id = result["entry_id"]
                pdb_id = result["pdb_id"]
                
                # Write chain1
                chain1_protein = result.get("chain_proteins", {}).get("chain1", "Unknown")
                chain1_org = result.get("chain_organisms", {}).get("chain1", {})
                chain1_uniprot = result.get("chain_uniprot_ids", {}).get("chain1", [])
                writer.writerow([
                    entry_id,
                    pdb_id,
                    result["chains"]["chain1"],
                    chain1_protein,
                    chain1_org.get("scientific_name", "Unknown"),
                    chain1_org.get("taxonomy_id", ""),
                    chain1_org.get("classification", "Unknown"),
                    ";".join(chain1_uniprot) if chain1_uniprot else ""
                ])
                
                # Write chain2 if exists
                if result["chains"]["chain2"]:
                    chain2_protein = result.get("chain_proteins", {}).get("chain2", "Unknown")
                    chain2_org = result.get("chain_organisms", {}).get("chain2", {})
                    chain2_uniprot = result.get("chain_uniprot_ids", {}).get("chain2", [])
                    writer.writerow([
                        entry_id,
                        pdb_id,
                        result["chains"]["chain2"],
                        chain2_protein,
                        chain2_org.get("scientific_name", "Unknown"),
                        chain2_org.get("taxonomy_id", ""),
                        chain2_org.get("classification", "Unknown"),
                        ";".join(chain2_uniprot) if chain2_uniprot else ""
                    ])
        
        action.log(message_type="csv_written", path=str(output_path))


def parse_line_numbers(line_numbers_str: str) -> Set[int]:
    """
    Parse line numbers from string format.
    
    Supports:
    - Single numbers: "5"
    - Ranges: "1-10"
    - Lists: "1,5,10"
    - Mixed: "1-5,10,15-20"
    
    Args:
        line_numbers_str: String representation of line numbers
    
    Returns:
        Set of line numbers
    """
    line_nums = set()
    parts = line_numbers_str.split(',')
    
    for part in parts:
        part = part.strip()
        if '-' in part:
            start, end = part.split('-')
            line_nums.update(range(int(start), int(end) + 1))
        else:
            line_nums.add(int(part))
    
    return line_nums


@app.command()
def resolve(
    input_file: Path = typer.Option(..., help="Path to the .jsonl.gz file to process"),
    line_numbers: str = typer.Option(..., help="Line numbers to process (e.g., '1,5,10' or '1-10' or '1-5,10,15-20')"),
    output_file: Optional[Path] = typer.Option(None, help="Path to save the output JSON file (optional)"),
    output_jsonl: Optional[Path] = typer.Option(None, help="Path to save filtered entries as JSONL (supports .gz)"),
    output_csv: Optional[Path] = typer.Option(None, help="Path to save protein/species summary as CSV"),
    show_chains: bool = typer.Option(True, help="Show individual chain information"),
    filter_organism: Optional[str] = typer.Option(None, help="Filter by organism name (case-insensitive partial match)"),
    filter_classification: Optional[str] = typer.Option(None, help="Filter by classification (Mammalian, Bacterial, Fungal, Viral, Plant, Archaeal, Other)"),
    mammals_only: bool = typer.Option(False, help="Shorthand to filter only mammalian proteins (same as --filter-classification Mammalian)"),
    timeout: int = typer.Option(10, help="Request timeout in seconds for API calls"),
    retries: int = typer.Option(3, help="Number of retry attempts for failed API calls"),
) -> None:
    """
    Resolve protein names from PDB IDs in JSONL.GZ entries.
    
    Extracts PDB IDs from entries, queries RCSB PDB for metadata,
    resolves protein names for each chain, and optionally filters by organism or classification.
    """
    with start_action(action_type="resolve_proteins", input_file=str(input_file)) as action:
        # Handle mammals_only shorthand
        if mammals_only:
            filter_classification = "Mammalian"
            typer.echo("Filtering for mammalian proteins only...")
        
        # Parse line numbers
        line_nums = parse_line_numbers(line_numbers)
        action.log(message_type="processing_lines", count=len(line_nums), lines=sorted(list(line_nums)))
        
        # Log filter settings
        if filter_organism or filter_classification:
            action.log(message_type="filter_settings", 
                      organism=filter_organism, 
                      classification=filter_classification)
            if filter_organism:
                typer.echo(f"Filtering by organism: {filter_organism}")
            if filter_classification:
                typer.echo(f"Filtering by classification: {filter_classification}")
        
        # Read entries from file
        entries = read_jsonl_gz_lines(input_file, line_nums)
        
        if not entries:
            typer.echo("No entries found at the specified line numbers.", err=True)
            raise typer.Exit(1)
        
        # Process each entry
        results = []
        filtered_results = []
        processed_pdb_ids = {}  # Cache PDB metadata to avoid duplicate queries
        
        for entry_data in entries:
            line_num = entry_data["line_number"]
            entry = entry_data["entry"]
            entry_id = entry.get("id", "")
            
            with start_action(action_type="process_entry", line_number=line_num, entry_id=entry_id) as entry_action:
                # Parse entry ID to extract PDB and chain info
                parsed = parse_entry_id(entry_id)
                pdb_id = parsed["pdb_id"]
                
                # Fetch PDB metadata (use cache if available)
                if pdb_id not in processed_pdb_ids:
                    metadata = fetch_pdb_metadata(pdb_id, timeout=timeout, retries=retries)
                    processed_pdb_ids[pdb_id] = metadata
                else:
                    metadata = processed_pdb_ids[pdb_id]
                    entry_action.log(message_type="using_cached_metadata", pdb_id=pdb_id)
                
                # Build result
                result = {
                    "line_number": line_num,
                    "entry_id": entry_id,
                    "pdb_id": pdb_id,
                    "chains": {
                        "chain1": parsed["chain1"],
                        "chain2": parsed["chain2"],
                    },
                    "metadata": metadata,
                }
                
                # Resolve protein names, organism info, and UniProt IDs for chains
                if show_chains and metadata.get("found"):
                    result["chain_proteins"] = {
                        "chain1": get_chain_protein_name(metadata, parsed["chain1"]),
                        "chain2": get_chain_protein_name(metadata, parsed["chain2"]),
                    }
                    result["chain_organisms"] = {
                        "chain1": get_chain_organism(metadata, parsed["chain1"]),
                        "chain2": get_chain_organism(metadata, parsed["chain2"]),
                    }
                    result["chain_uniprot_ids"] = {
                        "chain1": get_chain_uniprot_ids(metadata, parsed["chain1"]),
                        "chain2": get_chain_uniprot_ids(metadata, parsed["chain2"]),
                    }
                
                # Always add to results for processing
                results.append(result)
                
                # Apply filters and add to filtered_results if matches
                passes_filter = matches_filter(result, filter_organism, filter_classification)
                if passes_filter:
                    filtered_results.append(result)
                
                # Display result only if it passes filter
                if passes_filter:
                    typer.echo(f"\n{'='*80}")
                    typer.echo(f"Line {line_num}: {entry_id}")
                    typer.echo(f"{'='*80}")
                    typer.echo(f"PDB ID: {pdb_id}")
                    
                    if metadata.get("found"):
                        if "title" in metadata:
                            typer.echo(f"Title: {metadata['title']}")
                    if show_chains and "chain_proteins" in result:
                        chain1_org = result.get("chain_organisms", {}).get("chain1", {})
                        chain2_org = result.get("chain_organisms", {}).get("chain2", {})
                        chain1_uniprot = result.get("chain_uniprot_ids", {}).get("chain1", [])
                        chain2_uniprot = result.get("chain_uniprot_ids", {}).get("chain2", [])
                        
                        typer.echo(f"\nChain {parsed['chain1']}: {result['chain_proteins']['chain1']}")
                        typer.echo(f"  Organism: {chain1_org.get('scientific_name', 'Unknown')} ({chain1_org.get('classification', 'Unknown')})")
                        if chain1_uniprot:
                            typer.echo(f"  UniProt: {', '.join(chain1_uniprot)}")
                        
                        if parsed['chain2']:
                            typer.echo(f"\nChain {parsed['chain2']}: {result['chain_proteins']['chain2']}")
                            typer.echo(f"  Organism: {chain2_org.get('scientific_name', 'Unknown')} ({chain2_org.get('classification', 'Unknown')})")
                            if chain2_uniprot:
                                typer.echo(f"  UniProt: {', '.join(chain2_uniprot)}")
                        if "entities" in metadata:
                            typer.echo(f"\nAll entities in structure:")
                            for entity in metadata["entities"]:
                                organism = entity.get('organism', {})
                                org_name = organism.get('scientific_name', 'Unknown')
                                classification = classify_organism(organism.get('taxonomy_id'), org_name)
                                uniprot_ids = entity.get('uniprot_ids', [])
                                typer.echo(f"  - {entity.get('description', 'N/A')}")
                                typer.echo(f"    Chains: {', '.join(entity.get('chains', []))}")
                                typer.echo(f"    Organism: {org_name} ({classification})")
                                if uniprot_ids:
                                    typer.echo(f"    UniProt: {', '.join(uniprot_ids)}")
                    else:
                        typer.echo("Status: PDB not found or error fetching metadata")
                        if "error" in metadata:
                            typer.echo(f"Error: {metadata['error']}")
        
        # Summary
        total_processed = len(results)
        total_passed_filter = len(filtered_results)
        typer.echo(f"\n\n{'='*80}")
        typer.echo(f"SUMMARY")
        typer.echo(f"{'='*80}")
        typer.echo(f"Total entries processed: {total_processed}")
        typer.echo(f"Entries passing filter: {total_passed_filter}")
        typer.echo(f"Unique PDB structures: {len(processed_pdb_ids)}")
        
        if filter_organism or filter_classification:
            typer.echo(f"Filtered out: {total_processed - total_passed_filter}")
        
        # Check if we have filtered results to save
        if not filtered_results and (output_jsonl or output_csv):
            typer.echo("\nWarning: No entries passed the filter. No output files will be written.", err=True)
        
        # Save JSON metadata if requested
        if output_file and filtered_results:
            with start_action(action_type="save_results", output_file=str(output_file)) as save_action:
                with open(output_file, 'w') as f:
                    json.dump(filtered_results, f, indent=2)
                save_action.log(message_type="results_saved", path=str(output_file), count=len(filtered_results))
                typer.echo(f"\nJSON metadata saved to: {output_file}")
        
        # Save JSONL output with original entry data if requested
        if output_jsonl and filtered_results:
            write_jsonl_output(filtered_results, entries, output_jsonl)
            typer.echo(f"Filtered JSONL data saved to: {output_jsonl}")
        
        # Save CSV summary if requested
        if output_csv and filtered_results:
            write_csv_summary(filtered_results, output_csv)
            typer.echo(f"CSV summary saved to: {output_csv}")


if __name__ == "__main__":
    # Show help if no arguments provided
    if len(sys.argv) == 1:
        app(["--help"])
    else:
        app()

