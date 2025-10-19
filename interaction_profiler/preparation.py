from pathlib import Path
from typing import Any
import os
import numpy as np
import tempfile

from openbabel import openbabel
from openbabel import pybel
from openeye import oechem

from data.pdb_utils import VOCAB


def generated_to_xyz(data: tuple[int, list[str], list[list[float]]]) -> str:
    """Convert generated atomic data to XYZ format string.
    
    Args:
        data: Tuple of (num_atoms, atom_types, atom_coords)
        
    Returns:
        XYZ format string representation of the molecule
    """
    num_atoms, atom_type, atom_coords = data
    xyz = "%d\n\n" % (num_atoms,)
    for i in range(num_atoms):
        symb = atom_type[i]
        x, y, z = atom_coords[i]
        xyz += "%s %.8f %.8f %.8f\n" % (symb, x, y, z)
    return xyz


def generated_to_sdf(
    data: tuple[int, list[str], list[list[float]]]
) -> tuple[str, openbabel.OBMol]:
    """Convert generated atomic data to SDF format.
    
    Args:
        data: Tuple of (num_atoms, atom_types, atom_coords)
        
    Returns:
        Tuple of (sdf_string, openbabel_molecule)
    """
    xyz = generated_to_xyz(data)
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "sdf")

    mol = openbabel.OBMol()
    obConversion.ReadString(mol, xyz)
    sdf = obConversion.WriteString(mol)
    return sdf, mol


def visualize_mol_sdf(
    node_type: list[str],
    coords: list[list[float]],
    save_sdf_path: Path | str | None = None,
) -> tuple[openbabel.OBMol, str]:
    """Visualize a molecule from node types and coordinates.
    
    Args:
        node_type: List of atomic symbols
        coords: List of atomic coordinates
        save_sdf_path: Optional path to save the SDF file
        
    Returns:
        Tuple of (openbabel_molecule, sdf_string)
    """
    num_atoms = len(node_type)
    data = (num_atoms, node_type, coords)
    sdf_string, ob_mol = generated_to_sdf(data)
    if save_sdf_path is None:
        return ob_mol, sdf_string
    with open(save_sdf_path, "w") as f:
        f.write(sdf_string)
    return ob_mol, sdf_string


def visualize_data_sdf(
    data: dict[str, Any], save_dir: Path | str | None = None
) -> tuple[dict[int, openbabel.OBMol], dict[int, str]]:
    """Visualize all segments in the data as SDF files.
    
    Args:
        data: Data dictionary containing molecular information
        save_dir: Optional directory to save SDF files
        
    Returns:
        Tuple of (dict of molecules, dict of SDF strings) keyed by segment_id
    """
    ob_mols: dict[int, openbabel.OBMol] = {}
    sdf_strings: dict[int, str] = {}

    if save_dir is not None and not os.path.exists(save_dir):
        os.makedirs(save_dir)
    segment_start = 0
    block_segment_ids = np.array(data["segment_ids"])
    block_len = np.array(data["block_lengths"])
    for segment_id in sorted(set(data["segment_ids"])):
        segment_len = block_len[block_segment_ids == segment_id].sum()
        segment_end = segment_start + segment_len
        if VOCAB.idx_to_atom(data["A"][segment_start]) == VOCAB.atom_global:
            segment_start += 1
        atom_types = [VOCAB.idx_to_atom(x) for x in data["A"][segment_start:segment_end]]
        coords = data["X"][segment_start:segment_end]
        ob_mol, sdf_string = visualize_mol_sdf(
            atom_types,
            coords,
            f"{save_dir}/segment_{segment_id}.sdf" if save_dir is not None else None,
        )
        segment_start = segment_end
        ob_mols[segment_id] = ob_mol
        sdf_strings[segment_id] = sdf_string
    return ob_mols, sdf_strings


def addh_to_mol(input_sdf_path: Path | str, output_sdf_path: Path | str) -> None:
    """Add hydrogens to a molecule using OpenEye.
    
    OpenEye is used instead of OpenBabel as the latter doesn't work reliably
    for hydrogen addition in some cases.
    
    Args:
        input_sdf_path: Path to input SDF file
        output_sdf_path: Path to output SDF file with added hydrogens
    """
    # Create an input stream for reading the SDF file
    ifs = oechem.oemolistream()
    if not ifs.open(str(input_sdf_path)):
        print(f"Could not open {input_sdf_path} for reading.")
        exit()

    # Create an output stream for writing the modified SDF file
    ofs = oechem.oemolostream()
    if not ofs.open(str(output_sdf_path)):
        print(f"Could not open {output_sdf_path} for writing.")
        exit()

    # Create a molecule object to store each molecule
    mol = oechem.OEGraphMol()

    # Loop through the molecules in the SDF file
    while oechem.OEReadMolecule(ifs, mol):
        initial_formula = oechem.OEMolecularFormula(mol)
        # Add hydrogens to the molecule
        oechem.OEAssignImplicitHydrogens(mol)
        oechem.OEAddExplicitHydrogens(mol)
        final_formula = oechem.OEMolecularFormula(mol)
        # Write the molecule with added hydrogens to the output file
        oechem.OEWriteMolecule(ofs, mol)
        print(
            "Hydrogens have been added with OpenEye. Initial={}, Final={}".format(
                initial_formula, final_formula
            )
        )

    # Close the input and output streams
    ifs.close()
    ofs.close()


def build_pybel_mols(
    data: dict[str, Any], with_hydrogens: bool = True, tmpfile_dir: str = "./"
) -> dict[int, pybel.Molecule]:
    """Build PyBel molecule objects from data.
    
    Args:
        data: Data dictionary containing molecular information
        with_hydrogens: Whether to add hydrogens to molecules
        tmpfile_dir: Directory for temporary files
        
    Returns:
        Dictionary of PyBel molecules keyed by segment_id
    """
    ob_mols, sdf_strings = visualize_data_sdf(data)
    pybel_mols: dict[int, pybel.Molecule] = {}
    if not with_hydrogens:
        for k, v in ob_mols.items():
            pybel_mols[k] = pybel.Molecule(v)
    else:
        for segment_id in ob_mols:
            # Create a temporary file with the .sdf suffix
            fd, tmpfile_name = tempfile.mkstemp(suffix=".sdf", dir=tmpfile_dir)
            os.close(fd)
            with open(tmpfile_name, "w") as f:
                f.write(sdf_strings[segment_id])

            fd, tmpfile_with_hydrogens_name = tempfile.mkstemp(
                suffix=".sdf", dir=tmpfile_dir
            )
            os.close(fd)
            print(f"Created tmpfiles: {tmpfile_name}, {tmpfile_with_hydrogens_name}")

            addh_to_mol(tmpfile_name, tmpfile_with_hydrogens_name)
            for molecule in pybel.readfile("sdf", tmpfile_with_hydrogens_name):
                molecule.OBMol.AddNewHydrogens(0, True, 7.4)
                pybel_mols[segment_id] = molecule
            os.remove(tmpfile_name)
            os.remove(tmpfile_with_hydrogens_name)
            print(f"Removed tmpfiles: {tmpfile_name}, {tmpfile_with_hydrogens_name}")
    return pybel_mols

