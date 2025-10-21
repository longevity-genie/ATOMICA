# PDB Downloader CLI

A command-line interface for downloading PDB structures with optional metadata resolution.

## Features

- Download PDB files in CIF or PDB format using biotite
- Optional metadata resolution including organism information, UniProt IDs, resolution, etc.
- Batch download from command line arguments or text files
- Comprehensive logging with Eliot
- Default downloads go to `downloads/pdbs/` directory
- Metadata saved as JSON files alongside PDB files

## Installation

The PDB downloader is included with the ATOMICA project. Make sure all dependencies are installed:

```bash
uv sync
```

## Usage

### Basic Usage

Download a single PDB file:

```bash
uv run python -m pdb_downloader download 1tgr
```

Download multiple PDB files:

```bash
uv run python -m pdb_downloader download 1tgr 4xss 6pyh
```

### File Input

Download PDB files from a text file (one PDB ID per line):

```bash
uv run python -m pdb_downloader download-from-file pdb_ids.txt
```

Example `pdb_ids.txt`:
```
1tgr
4xss
# This is a comment and will be ignored
6pyh
```

### Options

#### Output Directory

Specify custom output directory:

```bash
uv run python -m pdb_downloader download 1tgr --output-dir /path/to/custom/directory
```

Default: `downloads/pdbs/`

#### File Format

Choose between CIF (mmCIF) or PDB format:

```bash
uv run python -m pdb_downloader download 1tgr --file-format pdb
```

Default: `cif` (recommended)

#### Metadata

Control metadata download:

```bash
# Disable metadata (download only PDB files)
uv run python -m pdb_downloader download 1tgr --include-metadata false

# Specify custom metadata directory
uv run python -m pdb_downloader download 1tgr --metadata-dir /path/to/metadata
```

Default: Downloads metadata alongside PDB files

#### Logging

Enable detailed logging to files:

```bash
uv run python -m pdb_downloader download 1tgr --log-to-file --log-file-name my_download
```

This creates:
- `logs/my_download.json` (machine-readable logs)
- `logs/my_download.log` (human-readable logs)

#### Advanced Options

```bash
# Set timeout and retry options for API calls
uv run python -m pdb_downloader download 1tgr --timeout 30 --retries 5

# Use local TSV files for metadata (if available)
uv run python -m pdb_downloader download 1tgr --use-tsv true
```

## Output Structure

```
downloads/pdbs/
├── 1tgr.cif/           # PDB structure file
│   └── 1tgr.cif
├── 1tgr_metadata.json  # Metadata (if enabled)
├── 4xss.cif/
│   └── 4xss.cif
└── 4xss_metadata.json
```

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

The PDB downloader integrates with ATOMICA's existing infrastructure:

- Uses biotite for reliable PDB downloading (same as embed_protein.py)
- Self-contained metadata resolution using biotite's RCSB API
- Follows ATOMICA's logging and CLI patterns
- Can be used as a preprocessing step for ATOMICA models

### Metadata Implementation

The metadata resolution is implemented directly in the PDB downloader using biotite's RCSB API, providing:
- Structure titles and descriptions
- Resolution and experimental methods
- Chain information and organism details
- UniProt ID mappings
- No external dependencies beyond biotite
