# IGF1 PDB Structures

This directory contains downloaded IGF1 (Insulin-like Growth Factor 1) PDB structures.

## Files

### 1TGR.pdb - Best IGF1 Model
- **Resolution:** 1.42 Ångströms (highest resolution available)
- **Chains:** A, B
- **Description:** Crystal structure of mini-IGF-1 (residues 1-52)
- **Size:** 115 KB
- **PDB Link:** https://www.rcsb.org/structure/1TGR

**Best for:**
- Structural visualization of IGF1
- Understanding the IGF1 fold and domain organization
- Backbone and secondary structure analysis

### 4XSS.pdb - IGF1 + Receptor Complex
- **Resolution:** 3.00 Ångströms
- **Chains:** B (IGF1), E (Insulin Receptor), F (IGF1R alpha-CT peptide)
- **Description:** IGF1 in complex with hybrid insulin/IGF1 receptor
- **Size:** 555 KB
- **PDB Link:** https://www.rcsb.org/structure/4XSS

**Best for:**
- Understanding IGF1-receptor binding interactions
- Studying receptor activation mechanism
- Functional analysis of IGF1 signaling

## Usage Examples

### Load with Biotite
```python
import biotite.structure.io as strucio

# Load IGF1 structure
igf1 = strucio.load_structure("1TGR.pdb")
print(igf1)

# Load complex structure
complex_struct = strucio.load_structure("4XSS.pdb")
print(complex_struct)
```

### Load with BioPython
```python
from Bio.PDB import PDBParser

parser = PDBParser(QUIET=True)
igf1_structure = parser.get_structure("1TGR", "1TGR.pdb")

# Access chains and atoms
for model in igf1_structure:
    for chain in model:
        print(f"Chain {chain.id}: {len(chain)} residues")
```

### Visualize with PyMOL
```bash
pymol 1TGR.pdb
pymol 4XSS.pdb
```

## Additional Resources

For other IGF1 structures, visit:
- https://www.rcsb.org/search?q=IGF1

Related structures:
- **7WRQ:** IGF1/IGFBP3/ALS ternary complex (3.60 Å)
- **5U8Q:** IGF1/IGF1R complex (3.27 Å)
- **1WQJ:** IGF1/Integrin complex (1.60 Å)
