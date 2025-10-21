# PDB CLI

A unified command-line interface for PDB operations in ATOMICA.

## Features

- **Download**: PDB files in CIF or PDB format using biotite
- **Convert**: PDB/CIF files to JSONL format for ATOMICA models
- **Metadata**: Optional metadata resolution including organism information, UniProt IDs, resolution, etc.
- **Batch Processing**: Download and convert multiple files at once
- **Comprehensive logging**: Eliot-based structured logging
- **Default structure**: Downloads go to `downloads/pdbs/` directory
- **Clean output**: No nested directories - files go directly to output folder

## Installation

The PDB CLI is included with the ATOMICA project. Make sure all dependencies are installed:

```bash
uv sync
```

## Usage

The PDB CLI is accessed via the unified `pdb` command with subcommands:

```bash
uv run pdb --help
```

### Download PDB Files

Download a single PDB file:

```bash
uv run pdb download 1tgr
```

Download multiple PDB files:

```bash
uv run pdb download 1tgr 4xss 6pyh
```

### Download from File

Download PDB files from a text file (one PDB ID per line):

```bash
uv run pdb download-from-file pdb_ids.txt
```

Example `pdb_ids.txt`:
```
1tgr
4xss
# This is a comment and will be ignored
6pyh
```

### Convert PDB to JSONL

Convert PDB/CIF files to JSONL format for ATOMICA models:

```bash
# Convert all CIF files in downloads/pdbs
uv run pdb convert downloads/pdbs/*.cif --output-file structures.jsonl

# Convert specific files
uv run pdb convert 1tgr.cif 4xss.cif --output-file my_structures.jsonl

# Convert with specific chains
uv run pdb convert 1tgr.cif --output-file 1tgr.jsonl --chains A,B
```

### Download Options

#### Output Directory

Specify custom output directory:

```bash
uv run pdb download 1tgr --output-dir /path/to/custom/directory
```

Default: `downloads/pdbs/`

#### File Format

Choose between CIF (mmCIF) or PDB format:

```bash
uv run pdb download 1tgr --file-format pdb
```

Default: `cif` (recommended)

#### Metadata

Control metadata download:

```bash
# Disable metadata (download only PDB files)
uv run pdb download 1tgr --include-metadata false

# Specify custom metadata directory
uv run pdb download 1tgr --metadata-dir /path/to/metadata
```

Default: Downloads metadata alongside PDB files

#### Logging

Enable detailed logging to files:

```bash
uv run pdb download 1tgr --log-to-file --log-file-name my_download
```

This creates:
- `logs/my_download.json` (machine-readable logs)
- `logs/my_download.log` (human-readable logs)

#### Advanced Options

```bash
# Set timeout and retry options for API calls
uv run pdb download 1tgr --timeout 30 --retries 5

# Use local TSV files for metadata (if available)
uv run pdb download 1tgr --use-tsv true
```

## Output Structure

### Download Output

```
downloads/pdbs/
├── 1tgr.cif            # PDB structure file
├── 1tgr_metadata.json  # Metadata (organism, UniProt IDs, resolution, etc.)
├── 4xss.cif
└── 4xss_metadata.json
```

### Convert Output

The `convert` command creates JSONL.GZ files:

```bash
uv run pdb convert downloads/pdbs/*.cif --output-file structures.jsonl
```

This creates `structures.jsonl.gz` with one JSON object per line, each containing:
- `id`: PDB identifier with chain info
- `data`: Processed structure data for ATOMICA models
- `block_to_pdb_indexes`: Mapping of block indices to PDB atom indices

## Metadata Format

The metadata JSON file contains:

```json
{
  "pdb_id": "1tgr",
  "found": true,
  "source": "API",
  "title": "Structure Title",
  "resolution": 2.1,
  "experimental_method": "X-RAY DIFFRACTION",
  "entities": [
    {
      "entity_id": "1",
      "description": "Protein description",
      "chains": ["A", "B"],
      "type": "polymer",
      "organism": {
        "scientific_name": "Homo sapiens",
        "taxonomy_id": 9606
      },
      "uniprot_ids": ["P12345", "Q67890"]
    }
  ]
}
```

## Logging

The CLI uses Eliot for structured logging. When `--log-to-file` is enabled:

- **JSON logs**: Machine-readable structured logs for programmatic analysis
- **Human logs**: Formatted logs for easy reading

Example log output:
```
2025-10-21 15:55:56Z [INFO] Starting PDB download batch
2025-10-21 15:55:57Z [INFO] Downloaded PDB: 1tgr.cif
2025-10-21 15:55:57Z [INFO] Downloaded metadata: 1tgr_metadata.json
```

## Error Handling

The CLI handles common errors gracefully:

- **Network timeouts**: Automatic retries with configurable timeout
- **Invalid PDB IDs**: Logged and skipped, continues with other downloads
- **Permission errors**: Clear error messages for directory access issues

Failed downloads are reported in the summary:
```
❌ Failed to download: 2def, invalid_id
✅ Downloaded 1/3 PDB files
```

## Integration with ATOMICA

The PDB CLI integrates seamlessly with ATOMICA's workflow:

- **Download**: Uses biotite for reliable PDB downloading (same as embed_protein.py)
- **Convert**: Prepares structures in JSONL format for ATOMICA models
- **Metadata**: Self-contained resolution using biotite's RCSB API
- **Logging**: Follows ATOMICA's Eliot logging patterns

### Complete Workflow Example

```bash
# 1. Download PDB structures with metadata
uv run pdb download 1tgr 4xss 6pyh --include-metadata

# 2. Convert to JSONL for ATOMICA training
uv run pdb convert downloads/pdbs/*.cif --output-file training_data.jsonl

# 3. Use with ATOMICA models
# training_data.jsonl.gz can now be used with ATOMICA trainers
```

### Commands Overview

| Command | Description | Example |
|---------|-------------|---------|
| `pdb download` | Download PDB structures with metadata | `uv run pdb download 1tgr` |
| `pdb download-from-file` | Batch download from text file | `uv run pdb download-from-file ids.txt` |
| `pdb convert` | Convert PDB/CIF to JSONL for ATOMICA | `uv run pdb convert *.cif --output-file data.jsonl` |

### Metadata Implementation

The metadata resolution is implemented directly in the PDB CLI using biotite's RCSB API, providing:
- Structure titles and descriptions
- Resolution and experimental methods
- Chain information and organism details
- UniProt ID mappings
- No external dependencies beyond biotite
