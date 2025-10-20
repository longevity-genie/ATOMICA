# PDB/CIF Direct Input Support

All three embedding scripts now support **direct PDB/CIF input** with automatic on-the-fly conversion to JSONL format.

## Quick Examples

### 1. Get Embeddings from PDB File
```bash
# Direct PDB input - automatically converted
python get_embeddings.py --data-path protein.pdb --output-path embeddings.parquet

# With specific chains
python get_embeddings.py --data-path protein.cif --output-path embeddings.parquet --chains A,B

# Still works with JSONL
python get_embeddings.py --data-path proteins.jsonl --output-path embeddings.parquet
```

### 2. Compute Interaction Scores from PDB
```bash
# Direct PDB input
python interaction_profiler/interact_score.py compute-interact-scores \
  --data-path protein.pdb \
  --output-path scores.jsonl \
  --model-ckpt model.pt

# With specific chains
python interaction_profiler/interact_score.py compute-interact-scores \
  --data-path protein.cif \
  --output-path scores.jsonl \
  --model-ckpt model.pt \
  --chains A
```

### 3. Embed Protein (unchanged, already supports PDB)
```bash
# from-file command accepts PDB/CIF
python embed_protein.py from-file protein.pdb -o embeddings.parquet
```

## Shared Converter CLI

A new `pdb-converter` CLI tool is available for batch PDB-to-JSONL conversion:

### Convert Single File
```bash
python pdb_converter.py convert-file protein.pdb -o protein.jsonl
python pdb_converter.py convert-file protein.cif -o protein.jsonl --chains A,B
```

### Convert Directory
```bash
# All PDB files in directory
python pdb_converter.py convert-directory ./pdbs/ -o all_proteins.jsonl

# All CIF files
python pdb_converter.py convert-directory ./cifs/ -o all_proteins.jsonl --pattern "*.cif"

# Compress output
python pdb_converter.py convert-directory ./pdbs/ -o all_proteins.jsonl.gz --compress

# Specific chains only
python pdb_converter.py convert-directory ./pdbs/ -o proteins_chainA.jsonl --chains A
```

### Batch Convert from File List
```bash
# Create file list
ls pdbs/*.pdb > file_list.txt

# Convert all
python pdb_converter.py batch-convert file_list.txt -o all_proteins.jsonl

# With compression
python pdb_converter.py batch-convert file_list.txt -o all_proteins.jsonl.gz --compress
```

## Installation

After pulling these changes:
```bash
uv sync
```

The following CLI commands will be available:
- `embed` - generate embeddings (now supports PDB)
- `interact-profiler` - interaction scores (now supports PDB)
- `pdb-converter` - standalone PDB-to-JSONL converter

## How It Works

1. **Automatic Detection**: Scripts detect `.pdb`, `.cif`, `.ent` extensions
2. **On-the-fly Conversion**: PDB files are converted to temporary JSONL
3. **Processing**: Temporary JSONL is processed normally
4. **Cleanup**: Temporary file is automatically deleted

## Features

- ✅ Supports PDB and mmCIF formats
- ✅ Chain selection support (`--chains A,B`)
- ✅ Automatic cleanup of temporary files
- ✅ Compatible with compressed files (`.pdb.gz`, `.cif.gz`)
- ✅ Backward compatible - JSONL/pickle still works
- ✅ Shared conversion code - consistent across all scripts

## Chain Selection

When using `--chains`:
- **Single chain**: `--chains A`
- **Multiple chains**: `--chains A,B,C`
- **No chains specified**: All chains are processed

## File Format Detection

Scripts automatically detect format based on extension:

**PDB Formats (auto-convert):**
- `.pdb`
- `.cif` (mmCIF)
- `.ent`
- `.pdb.gz` (compressed)
- `.cif.gz` (compressed)

**Data Formats (direct use):**
- `.json`
- `.jsonl`
- `.jsonl.gz`
- `.pkl` (pickle)

## JSONL Format

If you want to pre-convert PDB files to JSONL for reuse:

```bash
# Convert once
python pdb_converter.py convert-directory ./pdbs/ -o proteins.jsonl

# Use many times (faster, no conversion overhead)
python get_embeddings.py --data-path proteins.jsonl --output-path embeddings1.parquet
python get_embeddings.py --data-path proteins.jsonl --output-path embeddings2.parquet
```

Each line in JSONL contains:
```json
{
  "data": {
    "X": [...],       // atom coordinates
    "B": [...],       // block (residue) types
    "A": [...],       // atom types
    "atom_positions": [...],
    "block_lengths": [...],
    "segment_ids": [...]
  },
  "id": "1ABC_A",
  "block_to_pdb_indexes": {...}  // mapping to original PDB
}
```

## Performance Notes

- **Direct PDB input**: Slight conversion overhead (~1-2 seconds per file)
- **Pre-converted JSONL**: No conversion overhead, faster for repeated use
- **Batch conversion**: Best for processing many files multiple times

## Examples by Use Case

### One-off Analysis (Direct PDB)
```bash
# Fastest for single use
python get_embeddings.py --data-path protein.pdb --output-path embeddings.parquet
```

### Repeated Analysis (Pre-convert)
```bash
# Convert once
python pdb_converter.py convert-directory ./pdbs/ -o proteins.jsonl

# Use many times (faster)
python get_embeddings.py --data-path proteins.jsonl --output-path run1.parquet
python get_embeddings.py --data-path proteins.jsonl --output-path run2.parquet --batch-size 8
```

### Large Dataset Processing
```bash
# Convert with compression
python pdb_converter.py convert-directory ./pdbs/ -o proteins.jsonl.gz --compress

# Process compressed JSONL
python get_embeddings.py --data-path proteins.jsonl.gz --output-path embeddings.parquet
```

## Troubleshooting

### Error: "No blocks extracted from PDB"
- Check that the PDB file contains valid structure data
- Try specifying chains explicitly with `--chains`
- Verify the file isn't corrupted

### Conversion is Slow
- Pre-convert to JSONL for repeated use
- Use batch conversion with `pdb-converter batch-convert`
- Consider using compressed JSONL (`.jsonl.gz`)

---

**See also:**
- `EMBED_PROTEIN_USAGE.md` - embed_protein.py documentation
- `interaction_profiler/README.md` - interact_score.py documentation

