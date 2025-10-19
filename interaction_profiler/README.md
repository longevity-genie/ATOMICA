# ATOMICA Interaction Profiler

Tools for analyzing protein-ligand interactions and computing interaction importance scores.

## Overview

The Interaction Profiler module provides two main components:

### 1. Interaction Scoring (`interact_score.py`)

Computes interaction importance scores for protein blocks using a trained ATOMICA model. This evaluates how much each block (residue or molecule) contributes to the model's predictions.

**How it works:**
- Takes a protein structure and a trained prediction model
- For each block in the structure, masks it with a special MASK token
- Measures cosine similarity between original and masked predictions
- Outputs importance scores (cosine distances) for each block
- Runs on GPU with CUDA acceleration

**Use case:**
Identify which parts of a protein structure are most important for accurate interaction predictions. Blocks with higher cosine similarity scores have greater impact on the model's predictions.

### 2. Molecular Preparation (`preparation.py`)

Utility functions for preparing and visualizing molecular structures.

**Features:**
- Convert atomic data to XYZ/SDF format
- Build PyBel molecule objects for interaction analysis
- Add hydrogens to molecules using OpenEye
- Visualize molecular segments from data

## Usage

### Output Files

The tool generates and saves the following outputs:

1. **Interaction Scores** (`--output-path`)
   - JSONL file with one result per line
   - Each line contains: id, cos_distances, block_idx, time_seconds, peak_memory_mb
   - Continuously written during processing
   
2. **Processing Summary** (optional, `--summary-path`)
   - JSON file with aggregated statistics
   - Total time, per-structure time (mean/std/min/max), GPU memory usage
   - File paths and processing parameters
   - Only saved if `--summary-path` is specified

3. **Console Output**
   - Real-time progress with tqdm progress bar
   - Periodic metrics (every 10 structures)
   - Final summary report printed to console
   - Eliot structured logs (if logging is configured)

### Command Line Interface

The interaction profiler is available as a CLI command:

```bash
# Compute interaction scores
interact-profiler compute-interact-scores \
  --data-path /path/to/data.jsonl \
  --output-path /path/to/output.jsonl \
  --model-ckpt /path/to/model.pt

# Save results AND processing summary
interact-profiler compute-interact-scores \
  --data-path /path/to/data.jsonl \
  --output-path /path/to/output.jsonl \
  --model-ckpt /path/to/model.pt \
  --summary-path /path/to/summary.json
```

**Example Summary File** (`summary.json`):
```json
{
  "total_structures": 100,
  "total_time_seconds": 128.45,
  "total_time_minutes": 2.14,
  "time_per_structure": {
    "mean": 1.2845,
    "std": 0.1523,
    "min": 0.8934,
    "max": 1.6723
  },
  "gpu_memory_mb": {
    "mean": 395.2,
    "max": 425.1
  },
  "output_file": "/path/to/output.jsonl",
  "data_file": "/path/to/data.jsonl",
  "start_idx": 0,
  "num_lines": null
}
```

**Options:**
- `--data-path`: Path to input protein structure data file (jsonl format)
- `--output-path`: Path to save interaction scores (jsonl format)
- `--model-ckpt`: Path to pretrained ATOMICA model checkpoint
- `--start-idx`: Start index for processing (default: 0, useful for large files)
- `--num-lines`: Maximum number of lines to process (default: None = all lines)
- `--summary-path`: Optional path to save processing summary report as JSON (default: None, only prints to console)

**Output format:**
Each line in the output file is a JSON object:
```json
{
  "id": "protein_id",
  "cos_distances": [0.95, 0.87, ...],
  "block_idx": [0, 1, ...],
  "time_seconds": 2.34,
  "peak_memory_mb": 412.5
}
```

**Output fields:**
- `id`: Structure identifier from the dataset
- `cos_distances`: List of cosine similarity scores (interaction importance per block)
- `block_idx`: List of block indices corresponding to the scores
- `time_seconds`: Time taken to process this structure (in seconds)
- `peak_memory_mb`: Peak GPU memory used for this structure (in MB)

### Processing Large Files (Slicing)

For gigabyte-sized `.jsonl.gz` files, you can process specific slices:

```bash
# Process only lines 1000-2000 from a large file
interact-profiler compute-interact-scores \
  --data-path large_file.jsonl.gz \
  --output-path output_slice.jsonl \
  --model-ckpt model.pt \
  --start-idx 1000 \
  --num-lines 1000

# Process first 100 lines
interact-profiler compute-interact-scores \
  --data-path huge_file.jsonl.gz \
  --output-path output_first_100.jsonl \
  --model-ckpt model.pt \
  --num-lines 100

# Skip first 5000 lines, process remaining
interact-profiler compute-interact-scores \
  --data-path huge_file.jsonl.gz \
  --output-path output_from_5000.jsonl \
  --model-ckpt model.pt \
  --start-idx 5000
```

**Use cases:**
- Test on a small subset before processing entire gigabyte files
- Parallelize processing across multiple machines (split the file)
- Resume interrupted processing
- Memory-constrained environments

### Performance Monitoring

The tool provides real-time monitoring of:
- **Time per structure**: Processing time for each protein complex
- **GPU memory usage**: Peak CUDA memory allocated per structure
- **Running statistics**: Average times and memory for the last 10 structures

**Example output during processing:**

```
Processing structures: 42%|████▎     | 42/100 [00:53<01:15, 0.77/s]
[    10] Avg time (last 10): 1.23s, Avg memory (last 10): 387.5 MB
[    20] Avg time (last 10): 1.18s, Avg memory (last 10): 392.1 MB
[    30] Avg time (last 10): 1.25s, Avg memory (last 10): 401.3 MB
...
======================================================================
PROCESSING SUMMARY
======================================================================
Total structures processed: 100
Total time: 128.45s (2.14 minutes)
Time per structure: 1.28s ± 0.15s (min: 0.89s, max: 1.67s)
Estimated remaining time per structure: ~1.28s

GPU Memory per structure: 395.2 MB (avg), 425.1 MB (peak)
======================================================================
```

### Analyzing Results

Parse timing and memory metrics from output:

```python
import json

# Analyze performance metrics
times = []
memories = []

with open("output/interaction_scores.jsonl") as f:
    for line in f:
        result = json.loads(line)
        times.append(result["time_seconds"])
        memories.append(result["peak_memory_mb"])

import statistics
print(f"Mean time: {statistics.mean(times):.2f}s")
print(f"Median memory: {statistics.median(memories):.1f} MB")
print(f"Max memory used: {max(memories):.1f} MB")

# Find slowest structures
slow_structures = [(r["id"], r["time_seconds"]) 
                   for line in open("output/interaction_scores.jsonl")
                   for r in [json.loads(line)]
                   if r["time_seconds"] > 2.0]
print(f"Slow structures (>2s): {slow_structures}")
```

### Python API

```python
from interaction_profiler.interact_score import get_residue_model_scores
from models.prediction_model import PredictionModel
from data.dataset import PDBDataset

# Load model and data
model = PredictionModel.load_from_pretrained("path/to/model.pt")
dataset = PDBDataset("path/to/data.jsonl")

# Get scores for a structure
data = dataset[0]
cos_distances, block_idx = get_residue_model_scores(model, data)
```

#### Preparation functions:

```python
from interaction_profiler.preparation import build_pybel_mols, visualize_data_sdf

# Build PyBel molecules from data
pybel_mols = build_pybel_mols(data, with_hydrogens=True)

# Visualize segments
ob_mols, sdf_strings = visualize_data_sdf(data, save_dir="./molecules")
```

## Technical Details

### Model Requirements
- ATOMICA prediction model (trained on protein-ligand structures)
- GPU memory: Depends on structure size (typically 4-8GB for medium proteins)

### Dependencies
- PyTorch with CUDA support
- OpenBabel for molecular format conversion
- OpenEye Toolkit for hydrogen addition
- Typer for CLI framework
- Eliot for structured logging

### Performance
- GPU-accelerated inference with CUDA
- Processes multiple structures with progress bar
- Logs structured actions for debugging

## Implementation Details

### Type Safety
All functions include complete type hints for better IDE support and error detection.

### Logging
Uses Eliot structured logging with action contexts:
```python
with start_action(action_type="compute_interact_scores", ...) as action:
    # computation
    action.add_success_fields(num_structures=N)
```

### CLI Framework
Uses Typer for intuitive command-line interface with:
- Automatic help generation
- Type-based validation
- Auto-completion support
- Path validation

## Examples

### Compute Scores for a Dataset
```bash
interact-profiler compute-interact-scores \
  --data-path downloads/dataverse_files/pretraining_data/PL.jsonl.gz \
  --output-path output/interaction_scores.jsonl \
  --model-ckpt downloads/ATOMICA_checkpoints/prot_interface/atomica_interface_v3.pt
```

### Process Results
```python
import json

# Read and analyze scores
with open("output/interaction_scores.jsonl") as f:
    for line in f:
        result = json.loads(line)
        protein_id = result["id"]
        importance_scores = result["cos_distances"]
        print(f"{protein_id}: average importance = {sum(importance_scores)/len(importance_scores):.3f}")
```

## Future Enhancements

- Batch processing with configurable batch sizes
- Parallel processing across multiple GPUs
- Export results to different formats (CSV, HDF5)
- Visualization utilities for importance heatmaps
- Integration with molecular dynamics
