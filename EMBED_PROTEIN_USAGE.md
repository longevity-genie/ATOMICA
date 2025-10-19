# embed_protein.py - Quick Usage Guide

CLI tool to download PDB structures, process them, and generate ATOMICA embeddings.

## Features

- **Download from RCSB PDB** automatically using biotite
- **Download from AlphaFold DB** for predicted structures (no external libraries needed!)
- **AlphaFold fallback** - try PDB first, then AlphaFold if not found
- **Process local PDB/CIF files**
- **Batch processing** from CSV
- **Chain selection** support
- **GPU/CPU** inference
- **Parquet output** for efficient storage

## Installation

The tool uses existing ATOMICA dependencies. Make sure you have:
```bash
uv sync
```

## Note on AlphaFold Support

**Biotite does not have native AlphaFold DB support**, so we use direct HTTP download from the AlphaFold Database (https://alphafold.ebi.ac.uk/). This is a simple and robust solution that requires no additional dependencies.

## Commands

### 1. From PDB ID (Download from RCSB)

Download and embed a protein structure from RCSB PDB:

```bash
# Basic usage - entire structure
python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet

# Process specific chains only
python embed_protein.py from-pdb-id 6LLW --chains A -o 6llw_chainA.parquet

# Process multiple chains
python embed_protein.py from-pdb-id 2UXQ --chains A,B -o 2uxq_AB.parquet

# Use CPU instead of GPU
python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet --device cpu

# Use different model version
python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet \
    --model-config downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v2.json \
    --model-weights downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v2.pt

# Download as PDB format instead of CIF
python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet --format pdb

# Save logs
python embed_protein.py from-pdb-id 1ABC -o embeddings.parquet --log-file embedding.log
```

### 2. From UniProt ID (Download from AlphaFold DB)

Download and embed AlphaFold predicted structures:

```bash
# Download AlphaFold structure for a UniProt ID
python embed_protein.py from-uniprot P12345 -o embeddings.parquet

# Use specific AlphaFold version (2, 3, or 4)
python embed_protein.py from-uniprot Q9Y6K9 -o embeddings.parquet --version 3

# Process specific chains
python embed_protein.py from-uniprot P12345 -o embeddings.parquet --chains A

# Use CPU
python embed_protein.py from-uniprot P12345 -o embeddings.parquet --device cpu
```

**When to use this:**
- Protein has no experimental structure in PDB
- You want the AlphaFold predicted structure
- You have a UniProt accession ID

**Finding UniProt IDs:**
- Search at https://www.uniprot.org/
- Or use gene name: search "BRCA1 human" → P38398

### 3. From Local File

Process a local PDB or CIF file:

```bash
# Process local file
python embed_protein.py from-file /path/to/protein.pdb -o embeddings.parquet

# With custom ID
python embed_protein.py from-file /path/to/protein.cif -o embeddings.parquet --id my_protein

# Specific chains
python embed_protein.py from-file /path/to/protein.pdb --chains A,B -o embeddings.parquet
```

### 4. Batch Processing

Process multiple proteins from a CSV file:

```bash
python embed_protein.py batch proteins.csv -o all_embeddings.parquet --batch-size 8
```

**CSV Format:**
```csv
pdb_id,pdb_path,chains,source
1ABC,,A,rcsb
P12345,,,alphafold
2XYZ,/path/to/local.pdb,A_B,
Q9Y6K9,,,alphafold
3DEF,,,
```

- `pdb_id`: Required - PDB ID or UniProt accession
- `pdb_path`: Optional - if empty, downloads based on source
- `chains`: Optional - if empty, processes all chains (use underscore for multiple: A_B_C)
- `source`: Optional - 'rcsb', 'alphafold', or empty (auto-detect)

**Example batch.csv:**
```csv
pdb_id,pdb_path,chains,source
6LLW,,A,rcsb
2UXQ,,A_B,rcsb
P38398,,,alphafold
1ABC,/data/structures/1abc.pdb,A,
Q9Y6K9,,,alphafold
3I5X,,,rcsb
```

**With AlphaFold fallback:**
```bash
# If PDB not found, automatically try AlphaFold DB
python embed_protein.py batch proteins.csv -o all_embeddings.parquet --alphafold-fallback
```

## Output Format

The tool saves embeddings in Parquet format with the following columns:

- `id`: Structure identifier (e.g., "6LLW_A")
- `graph_embedding`: List - overall structure embedding (use this for similarity search!)
- `block_embedding`: List - per-residue embeddings
- `atom_embedding`: List - per-atom embeddings
- `block_id`: List - block type IDs
- `atom_id`: List - atom type IDs

### Reading the Output

```python
import polars as pl
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# Load embeddings
df = pl.read_parquet("embeddings.parquet")

# Extract graph embeddings for similarity search
graph_embeddings = np.array([np.array(emb) for emb in df['graph_embedding']])

# Compute pairwise similarities
similarities = cosine_similarity(graph_embeddings)

# Find most similar proteins to the first one
query_idx = 0
similar_indices = np.argsort(similarities[query_idx])[::-1][1:6]  # Top 5 (excluding self)
print(f"Most similar proteins to {df['id'][query_idx]}:")
for idx in similar_indices:
    print(f"  {df['id'][idx]}: similarity = {similarities[query_idx][idx]:.3f}")
```

## Examples

### Example 1: Mix PDB and AlphaFold structures

```bash
# Create CSV with both PDB and AlphaFold sources
cat > mixed_proteins.csv << EOF
pdb_id,pdb_path,chains,source
1HRC,,A,rcsb
P12345,,,alphafold
6LLW,,A,rcsb
Q9Y6K9,,,alphafold
EOF

# Process all
python embed_protein.py batch mixed_proteins.csv -o mixed_embeddings.parquet

# Analyze similarity
python << 'PYEND'
import polars as pl
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

df = pl.read_parquet("mixed_embeddings.parquet")
graph_embs = np.array([np.array(e) for e in df['graph_embedding']])

# Compare PDB vs AlphaFold embeddings
for i in range(len(df)):
    for j in range(i+1, len(df)):
        sim = cosine_similarity(graph_embs[i:i+1], graph_embs[j:j+1])[0][0]
        print(f"{df['id'][i]} vs {df['id'][j]}: {sim:.3f}")
PYEND
```

### Example 2: Analyze heme-binding proteins

```bash
# Create CSV with known heme-binding proteins
cat > heme_proteins.csv << EOF
pdb_id,pdb_path,chains
1A6N,,A
1HRC,,A
1MBN,,
2HHB,,A_B
EOF

# Generate embeddings
python embed_protein.py batch heme_proteins.csv -o heme_embeddings.parquet

# Analyze in Python
python << 'PYEND'
import polars as pl
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

df = pl.read_parquet("heme_embeddings.parquet")
graph_embs = np.array([np.array(e) for e in df['graph_embedding']])
sim_matrix = cosine_similarity(graph_embs)

print("Heme-binding protein similarity matrix:")
print(sim_matrix)
PYEND
```

### Example 3: AlphaFold-only analysis for novel proteins

```bash
# Create list of UniProt IDs for proteins without structures
cat > novel_proteins.csv << EOF
pdb_id,pdb_path,chains,source
P12345,,,alphafold
Q9Y6K9,,,alphafold
O15111,,,alphafold
P51608,,,alphafold
EOF

# Generate embeddings
python embed_protein.py batch novel_proteins.csv -o novel_embeddings.parquet

# Cluster them
python << 'PYEND'
import polars as pl
import numpy as np
from sklearn.cluster import KMeans

df = pl.read_parquet("novel_embeddings.parquet")
embeddings = np.array([np.array(e) for e in df['graph_embedding']])

# Cluster by interface similarity
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(embeddings)

for i, (protein_id, cluster) in enumerate(zip(df['id'], clusters)):
    print(f"{protein_id}: Cluster {cluster}")
PYEND
```

### Example 4: Find proteins similar to a query

```bash
# Embed your query protein
python embed_protein.py from-pdb-id 6LLW --chains A -o query.parquet

# Embed a database of proteins
python embed_protein.py batch protein_database.csv -o database.parquet

# Search for similar proteins in Python
python << 'PYEND'
import polars as pl
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

query_df = pl.read_parquet("query.parquet")
db_df = pl.read_parquet("database.parquet")

query_emb = np.array(query_df['graph_embedding'][0]).reshape(1, -1)
db_embs = np.array([np.array(e) for e in db_df['graph_embedding']])

similarities = cosine_similarity(query_emb, db_embs)[0]
top_10 = np.argsort(similarities)[::-1][:10]

print("Top 10 most similar proteins:")
for idx in top_10:
    print(f"{db_df['id'][idx]}: {similarities[idx]:.3f}")
PYEND
```

### Example 5: Use AlphaFold fallback for robustness

```bash
# Some IDs might be PDB, some might only be in AlphaFold
cat > mixed_ids.csv << EOF
pdb_id,pdb_path,chains
1ABC,,A
FAKE123,,A
P12345,,A
6LLW,,
EOF

# Try PDB first, fallback to AlphaFold (treating as UniProt ID)
python embed_protein.py batch mixed_ids.csv -o robust_embeddings.parquet --alphafold-fallback

# FAKE123 will fail (neither PDB nor UniProt)
# P12345 will use AlphaFold (not a PDB ID)
# 1ABC and 6LLW will use PDB
```

### Example 6: Cluster proteins by interface similarity

```bash
python embed_protein.py batch my_proteins.csv -o embeddings.parquet

python << 'PYEND'
import polars as pl
import numpy as np
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

df = pl.read_parquet("embeddings.parquet")
embeddings = np.array([np.array(e) for e in df['graph_embedding']])

# Cluster
kmeans = KMeans(n_clusters=5, random_state=42)
clusters = kmeans.fit_predict(embeddings)

# Visualize with t-SNE
tsne = TSNE(n_components=2, random_state=42)
coords_2d = tsne.fit_transform(embeddings)

plt.figure(figsize=(10, 8))
scatter = plt.scatter(coords_2d[:, 0], coords_2d[:, 1], c=clusters, cmap='viridis')
for i, protein_id in enumerate(df['id']):
    plt.annotate(protein_id, (coords_2d[i, 0], coords_2d[i, 1]), fontsize=8)
plt.colorbar(scatter)
plt.title("Protein Interface Clustering")
plt.savefig("protein_clusters.png", dpi=300, bbox_inches='tight')
PYEND
```

## Model Options

The tool uses ATOMICA-Interface by default. Three versions are available:

- **v1**: `atomica_interface_v1.json` / `atomica_interface_v1.pt`
- **v2**: `atomica_interface_v2.json` / `atomica_interface_v2.pt`
- **v3**: `atomica_interface_v3.json` / `atomica_interface_v3.pt` (default)

All versions are in `downloads/ATOMICA_checkpoints/prot_interface/`

## Performance Tips

1. **GPU vs CPU**: GPU is much faster for batch processing
2. **Batch size**: Increase `--batch-size` if you have enough GPU memory
3. **Chains**: Only process chains you need to save computation time
4. **Format**: CIF format is more reliable than legacy PDB format

## Troubleshooting

### CUDA Out of Memory
```bash
# Reduce batch size
python embed_protein.py batch proteins.csv -o output.parquet --batch-size 1

# Or use CPU
python embed_protein.py batch proteins.csv -o output.parquet --device cpu
```

### PDB Download Fails
```bash
# Try PDB format instead of CIF
python embed_protein.py from-pdb-id 1ABC -o output.parquet --format pdb

# Or download manually and use from-file
```

### Chain Not Found
```bash
# List available chains first (use biotite or PyMOL)
# Then specify the correct chain ID
python embed_protein.py from-pdb-id 1ABC --chains A -o output.parquet
```

## Advanced: Using Different Models

You can also use the pretrained model or ligand-specific models:

```bash
# Use pretrained model (works for complexes)
python embed_protein.py from-file complex.pdb -o embeddings.parquet \
    --model-config downloads/ATOMICA_checkpoints/pretrain/pretrain_model_config.json \
    --model-weights downloads/ATOMICA_checkpoints/pretrain/pretrain_model_weights.pt
```

## Integration with Your Pipeline

The tool is designed to be pipeline-friendly:

```bash
# Example: Process all PDBs in a directory
for pdb in /data/structures/*.pdb; do
    id=$(basename "$pdb" .pdb)
    python embed_protein.py from-file "$pdb" --id "$id" -o "embeddings/${id}.parquet"
done

# Combine all parquets
python << 'PYEND'
import polars as pl
from pathlib import Path

dfs = [pl.read_parquet(f) for f in Path("embeddings").glob("*.parquet")]
combined = pl.concat(dfs)
combined.write_parquet("all_embeddings.parquet")
PYEND
```

