# AlphaFold Database Integration

## Summary

**Biotite does NOT have native AlphaFold DB support**, so we implemented direct HTTP download from the AlphaFold Database (https://alphafold.ebi.ac.uk/). This is a simple, robust solution requiring no additional dependencies.

## Why This Approach?

After researching available Python libraries:
- `biotite.database` only supports RCSB PDB, not AlphaFold DB
- No mature Python packages exist for AlphaFold DB downloads
- Direct HTTP download is:
  - Simple and reliable
  - No extra dependencies needed
  - Officially documented by AlphaFold DB
  - Used by the research community

## AlphaFold Database API

The AlphaFold Database provides a simple REST API:

```
https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.{format}
```

- `uniprot_id`: UniProt accession (e.g., P12345, Q9Y6K9)
- `version`: AlphaFold DB version (2, 3, or 4)
- `format`: `cif` (mmCIF, recommended) or `pdb`

### Example URLs:
```
https://alphafold.ebi.ac.uk/files/AF-P12345-F1-model_v4.cif
https://alphafold.ebi.ac.uk/files/AF-Q9Y6K9-F1-model_v4.pdb
```

## Implementation

Our implementation in `embed_protein.py`:

```python
def download_alphafold_structure(
    uniprot_id: str,
    output_dir: Path,
    version: int = 4,
    use_pdb: bool = False
) -> Path:
    """Download AlphaFold predicted structure from AlphaFold DB."""
    # Direct HTTP download using urllib
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v{version}.{format}"
    
    with urllib.request.urlopen(url) as response:
        with open(file_path, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
    
    return file_path
```

**Features:**
- ✅ Robust error handling (404, 503, network errors)
- ✅ Helpful error messages
- ✅ Support for multiple AlphaFold versions
- ✅ Both CIF and PDB formats
- ✅ No external dependencies beyond Python standard library

## Usage Examples

### 1. Direct AlphaFold Download

```bash
# Download AlphaFold structure for BRCA1 (P38398)
python embed_protein.py from-uniprot P38398 -o brca1_embeddings.parquet

# Use specific version
python embed_protein.py from-uniprot P38398 -o brca1_embeddings.parquet --version 3
```

### 2. AlphaFold Fallback

```bash
# Try PDB first, fallback to AlphaFold if not found
python embed_protein.py from-pdb-id P38398 -o embeddings.parquet --alphafold-fallback
```

### 3. Batch Processing with Mixed Sources

```csv
pdb_id,pdb_path,chains,source
1ABC,,A,rcsb
P38398,,,alphafold
6LLW,,A,rcsb
Q9Y6K9,,,alphafold
```

```bash
python embed_protein.py batch proteins.csv -o all_embeddings.parquet
```

## Finding UniProt IDs

### Method 1: UniProt Website
1. Go to https://www.uniprot.org/
2. Search by gene name: "BRCA1 human"
3. Get UniProt ID: P38398

### Method 2: From PDB
If you know a PDB ID but want the AlphaFold structure:
1. Go to https://www.rcsb.org/structure/1ABC
2. Look for "UniProt" in the sequence tab
3. Use that UniProt ID with `from-uniprot`

### Method 3: Programmatically
```python
import biotite.database.uniprot as uniprot

# Search by gene name
query = uniprot.SimpleQuery("gene", "BRCA1") & uniprot.SimpleQuery("organism_id", "9606")
uniprot_ids = uniprot.search(query)
print(uniprot_ids[0])  # P38398
```

## AlphaFold Coverage

AlphaFold DB contains predicted structures for:
- **~214 million** protein structures (as of 2024)
- All proteins in:
  - Human proteome (~20,000 proteins)
  - Mouse proteome
  - Other model organisms
  - Swiss-Prot reviewed proteins
  - Many more species

Coverage is MUCH larger than PDB (~200,000 experimental structures).

## When to Use AlphaFold vs PDB

### Use PDB (RCSB) when:
- ✅ Experimental structure exists
- ✅ Need bound ligands/cofactors in structure
- ✅ Need specific experimental conditions
- ✅ Need high confidence in structure accuracy

### Use AlphaFold when:
- ✅ No experimental structure available
- ✅ Novel or understudied proteins
- ✅ Need structure for all isoforms
- ✅ Studying protein domains in isolation
- ✅ High-throughput studies

### Use Both:
- ✅ Compare predicted vs experimental
- ✅ Validation studies
- ✅ Use AlphaFold with `--alphafold-fallback` for robustness

## Confidence Scores

AlphaFold structures include per-residue confidence scores (pLDDT):
- **90-100**: Very high confidence (often as good as experimental)
- **70-90**: High confidence
- **50-70**: Low confidence
- **<50**: Very low confidence (disordered regions)

You can check confidence in the downloaded CIF files or at:
https://alphafold.ebi.ac.uk/entry/{uniprot_id}

## Limitations

1. **AlphaFold predicts apo structures** (without ligands)
   - For binding site analysis with ligands, use PDB when available
   - ATOMICA-Interface can still predict interface regions

2. **Flexible/disordered regions** have low confidence
   - Check pLDDT scores
   - These regions may not be reliable

3. **Multimeric states** 
   - AlphaFold predicts monomers by default
   - For protein complexes, check AlphaFold-Multimer predictions

4. **Not all proteins** are in AlphaFold DB
   - Coverage is excellent but not 100%
   - Newly discovered proteins may not be included

## Alternative Libraries Considered

We investigated these options:

1. **`biotite.database.alphafold`** - Does not exist
2. **`alphafoldmodel`** - Only for parsing, not downloading
3. **`af-analysis`** - For local AlphaFold runs, not DB downloads
4. **`alphapulldown`** - For AlphaFold-Multimer, not DB downloads
5. **Direct HTTP** - ✅ **Chosen**: Simple, reliable, no dependencies

## References

- AlphaFold Database: https://alphafold.ebi.ac.uk/
- API Documentation: https://alphafold.ebi.ac.uk/api-docs
- AlphaFold Paper: Jumper et al., Nature (2021)
- AlphaFold Database Paper: Varadi et al., Nucleic Acids Research (2022)

## Future Considerations

If AlphaFold DB access becomes more complex, we could:
1. Add retry logic with exponential backoff
2. Implement local caching
3. Add batch download optimization
4. Support AlphaFold-Multimer structures

For now, the simple HTTP approach is sufficient and maintainable.

