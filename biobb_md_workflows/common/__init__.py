#!/usr/bin/env python3

"""Shared helpers for the biobb_md_workflows CLI pipelines."""

from typing import Any, Dict, List

import numpy as np
from Bio.PDB import PDBParser


# Constants
# Possible ion names in the input structure not recognized by GROMACS default "Ion" group
ions_library: List[str] = ["K+", "CL-", "MG", "Cl-", "Na+"]

# All other solvent names in the input structure not recognized by GROMACS default "SOL" group
solvent_library: List[str] = ["WAT", "SOL"]

# Group names used for T-coupling and trajectory post-processing. Shared by md_gromacs (which
# builds the index) and traj_postprocessing (whose YAML references these names) - they must match.
solvent_group: str = "Solvent_group"
solute_group: str = "Solute_group"
output_group: str = "Output_group"


def to_yaml(value: Any) -> str:
    """Render a Python value as a YAML scalar for injection into a config template.

    Ensures ``None`` becomes ``null`` (not the string ``"None"``), booleans become
    lowercase ``true``/``false``, and lists become YAML flow sequences. Everything
    else is rendered with ``str()``.
    """
    if value is None:
        return "null"
    # bool must be checked before int/float (bool is a subclass of int)
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (list, tuple)):
        return "[" + ", ".join(to_yaml(v) for v in value) + "]"
    return str(value)


# GROMACS index (.ndx) files
def read_groups(index_path: str) -> List[str]:
    """Read all groups in an index file and return a list of their names."""

    # Read index file
    with open(index_path, 'r') as f:
        lines = f.readlines()

    # Extract group names
    group_names = []
    for line in lines:
        if line.strip().startswith("[") and line.strip().endswith("]"):
            group_name = line.strip()[1:-1].strip()
            group_names.append(group_name)

    return group_names


def read_group_indices(index_path: str) -> Dict[str, List[int]]:
    """
    Parse a GROMACS index file into a dictionary mapping each group name to its atom indices.

    Used to compare the size of two groups (e.g. Output_group vs System) without re-running
    make_ndx. Repeated group names keep the last occurrence, as GROMACS does.
    """
    groups: Dict[str, List[int]] = {}
    current_group = None

    with open(index_path, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                current_group = stripped[1:-1].strip()
                groups[current_group] = []
            elif current_group is not None and stripped:
                groups[current_group].extend(int(index) for index in stripped.split())

    return groups


def rename_last_ndx_group(ndx_path: str, new_name: str) -> None:
    """
    Reads a GROMACS index (.ndx) file and renames the last group in the file.

    Args:
        ndx_path (str): Path to the index file.
        new_name (str): The new name for the group (e.g., 'Protein_Ligand').
    """
    with open(ndx_path, 'r') as f:
        lines = f.readlines()

    # Iterate backwards through the file to find the last group header
    for i in range(len(lines) - 1, -1, -1):
        line = lines[i].strip()
        # Group headers in GROMACS look like: [ Group_Name ]
        if line.startswith("[") and line.endswith("]"):
            lines[i] = f"[ {new_name} ]\n"
            break

    # Overwrite the file with the corrected name
    with open(ndx_path, 'w') as f:
        f.writelines(lines)


def add_group(atom_indices: List, group_name: str, old_ndx_path: str, new_ndx_path: str) -> None:
    """
    Add all the atom indices to a new group at the end of the old ndx file. Save the new ndx file onto
    a different path.

    Args:
        atom_indices (List): List of atom indices to add to the new group.
        group_name (str): Name of the new group.
        old_ndx_path (str): Path to the existing index file.
        new_ndx_path (str): Path to save the new index file.

    Returns:
        None
    """
    COLUMNS = 15

    # Format the group header
    group_block = f"\n[ {group_name} ]\n"

    # Format indices in rows of 15, right-justified to match GROMACS style
    for i, idx in enumerate(atom_indices):
        group_block += f"{idx:>4}"
        if (i + 1) % COLUMNS == 0:
            group_block += "\n"
        else:
            group_block += " " if (i + 1) < len(atom_indices) else ""

    # Ensure the block ends with a newline
    if len(atom_indices) % COLUMNS != 0:
        group_block += "\n"

    # Read old ndx file and append the new group
    with open(old_ndx_path, 'r') as f:
        old_content = f.read()

    with open(new_ndx_path, 'w') as f:
        f.write(old_content)
        f.write(group_block)


def build_solvent_selection(solvent_names: List[str], ion_names: List[str]) -> str:
    """
    Create a GROMACS selection string based on the Default solvent and ion groups in
    GROMACS + the specific solvent and ion names provided.

    These groups will be useful for Temperature coupling and post-processing of the traj
    """

    # Find the selection of default and additional solvent molecules
    if solvent_names:
        solvent_selection = f'"SOL" | {" | ".join(f"r {solvent}" for solvent in solvent_names)}'
    else:
        solvent_selection = '"SOL"'

    # Find the selection of default and additional ions
    if ion_names:
        ions_selection = f'"Ion" | {" | ".join(f"a {ion}" for ion in ion_names)}'
    else:
        ions_selection = '"Ion"'

    # Join both selections
    return f'{solvent_selection} | {ions_selection}'


# Biopython helpers
def get_residue_types(pdb_path: str, target_resnames: List[str]) -> List[str]:
    """
    Load the pdb path with BioPython and find if any residue names from the
    target_resnames list exist in the structure.

    Inputs
    ------
        pdb_path (str): Path to the pdb file.
        target_resnames (List[str]): List of residue names to search for (e.g., ions or solvents).

    Returns
    -------
        list
            list of unique residue names found in the structure that match the target list.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)

    found_residues = set()

    # We convert the target library to a set for O(1) lookup speed
    # This is important if your structure has 50,000+ water molecules
    target_set = set(target_resnames)

    for residue in structure.get_residues():
        res_name = residue.get_resname().strip()

        if res_name in target_set:
            found_residues.add(res_name)

    return list(found_residues)


def get_atom_types(pdb_path: str, target_atom_names: List[str]) -> List[str]:
    """
    Load the pdb path and find if any atom names from the target list exist
    in the structure.

    Inputs
    ------
        pdb_path (str): Path to the pdb file.
        target_atom_names (List[str]): List of atom names to search for.

    Returns
    -------
        list
            list of unique atom names found in the structure.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_path)

    found_atoms = set()
    # Convert list to set for faster lookup (O(1))
    target_set = set(target_atom_names)

    # .get_atoms() is a generator that recursively yields every Atom in the structure
    for atom in structure.get_atoms():
        # Atom names in PDB are 4 chars. " CA " becomes "CA" after strip()
        atom_name = atom.get_name().strip()

        if atom_name in target_set:
            found_atoms.add(atom_name)

    return list(found_atoms)


def get_central_atom_index(pdb_path: str) -> int:
    """
    Find the index of the atom closest to the geometric center of a PDB file.

    Args:
        pdb_path: Path to the PDB file

    Returns:
        The atom index (1-based) of the atom closest to the geometric center
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Collect all atoms and their coordinates
    atoms = list(structure.get_atoms())
    coords = np.array([atom.get_vector().get_array() for atom in atoms])

    # Calculate geometric center
    center = coords.mean(axis=0)

    # Find atom closest to center
    distances = np.linalg.norm(coords - center, axis=1)
    central_atom_index = int(np.argmin(distances)) + 1

    return central_atom_index
