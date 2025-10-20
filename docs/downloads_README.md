# ATOMICA: Learning Universal Representations of Intermolecular Interactions

This repo contains the trained model weights and configs for the ATOMICA models. 

ATOMICA is a geometric AI model that learns universal representations of molecular interactions at an atomic scale. The model is pretrained on 2,037,972 molecular interaction interfaces from the Protein Data Bank and Cambridge Structural Database, this includes protein-small molecule, protein-ion, small molecule-small molecule, protein-protein, protein-peptide, protein-RNA, protein-DNA, and nucleic acid-small molecule complexes. Embeddings of ATOMICA can be generated with the open source model weights and code to be used for various downstream tasks. In the paper, we demonstrate the utility of ATOMICA embeddings for studying the human interfaceome network with ATOMICANets and for annotating ions and small molecules to proteins in the dark proteome.

[Preprint](https://www.biorxiv.org/content/10.1101/2025.04.02.646906v1) | [Project Website](https://zitniklab.hms.harvard.edu/projects/ATOMICA) | [GitHub](https://github.com/mims-harvard/ATOMICA)

### Model Checkpoints
The following models are available:
* ATOMICA model
* Pretrained ATOMICA-Interface model for construction of ATOMICANets
* Finetuned ATOMICA-Ligand prediction models for the following ligands:
    * metal ions: Ca, Co, Cu, Fe, K, Mg, Mn, Na, Zn
    * small molecules: ADP, ATP, GTP, GDP, FAD, NAD, NAP, NDP, HEM, HEC, CIT, CLA

### Setup Instructions
1. Install the huggingface cli `pip install -U "huggingface_hub[cli]"`
2. Download the checkpoints with `hf download ada-f/ATOMICA`
3. Known issue: `ATOMICA_checkpoints/ligand/small_molecules/NAD/NAD_v2.pt` has a [HuggingFace server-side issue](https://github.com/mims-harvard/ATOMICA/issues/8) where the uploaded and downloaded file does not match. In the interim, please use the checkpoint provided on [Google Drive](https://drive.google.com/file/d/1Dwajwx7hgOCEZYN2qwl6H8vJsnwcZSov/view?usp=sharing).

---
license: cc-by-4.0
---
