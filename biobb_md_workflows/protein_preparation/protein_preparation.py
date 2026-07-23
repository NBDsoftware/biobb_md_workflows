#!/usr/bin/env python3

# Importing all the needed libraries
from typing import List, Dict, Union, Optional, Literal
from pathlib import Path
import propka.run
import argparse
import shutil
import time
import os
import yaml

from Bio.SeqIO.PdbIO import PdbSeqresIterator
from Bio.PDB import PDBParser, PDBIO
from Bio import SeqIO

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_md_workflows.common import to_yaml
from biobb_md_workflows import __version__
from biobb_io.api.canonical_fasta import canonical_fasta
from biobb_model.model.fix_backbone import fix_backbone
from biobb_model.model.fix_side_chain import fix_side_chain
from biobb_model.model.fix_ssbonds import fix_ssbonds
from biobb_model.model.fix_altlocs import fix_altlocs
from biobb_model.model.fix_amides import fix_amides
from biobb_model.model.fix_chirality import fix_chirality
from biobb_model.model.mutate import mutate
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_structure_utils.utils.renumber_structure import renumber_structure
from biobb_pdb_tools.pdb_tools.biobb_pdb_tofasta import biobb_pdb_tofasta
from biobb_amber.pdb4amber.pdb4amber_run import pdb4amber_run
from biobb_chemistry.ambertools.reduce_remove_hydrogens import reduce_remove_hydrogens


# Constants

# Titratable residues in GROMACS
gmx_titra_resnames = {
    'LYS': ['LYSN', 'LYS'],                 # deprotonated, protonated
    'ARG': ['ARGN', 'ARG'],                 # deprotonated, protonated
    'ASP': ['ASP', 'ASPH'],                 # deprotonated, protonated
    'GLU': ['GLU', 'GLUH'],                 # deprotonated, protonated
    'HIS': ['HISD', 'HISE', 'HISH', 'HIS1'] # delta, epsilon, protonated, bound to HEME
}

# Titratable residues in AMBER
amber_titra_resnames = {
    'LYS': ['LYN', 'LYS'],                  # deprotonated, protonated
    'ASP': ['ASP', 'ASH'],                  # deprotonated, protonated
    'GLU': ['GLU', 'GLH'],                  # deprotonated, protonated  
    'HIS': ['HID', 'HIE', 'HIP']            # delta, epsilon, protonated
}

# Mapping between residue names and their protonation states
titra_mapping = {
    "amber": amber_titra_resnames,
    "gromacs": gmx_titra_resnames
}

# List of GROMACS 4-letter residue names
gmx_4letter_resnames = [
    'LYSN',
    'ARGN',
    'ASPH',
    'GLUH',
    'HISD',
    'HISE',
    'HISH',
    'HIS1'
]

# Reverse index mapping any titratable variant resname to its (base residue, protonation index).
# Built from both format maps so amber (HID, ...) and gromacs (HISD, ...) names are recognized.
titra_variant_index = {}
for _fmt_resnames in (amber_titra_resnames, gmx_titra_resnames):
    for _base, _names in _fmt_resnames.items():
        for _idx, _name in enumerate(_names):
            titra_variant_index[_name] = (_base, _idx)

# Disulfide-bond cysteine markers and their per-format target name
ss_bond_resnames = {'CYX', 'CYS2'}
ss_target_resname = {'amber': 'CYX', 'gromacs': 'CYS2'}

# Non-standard titratable variant names (i.e. those that differ from the plain base name).
# These are the names captured before standardization and restored at the end.
nonstd_titra_resnames = {
    'LYN', 'ASH', 'GLH', 'HID', 'HIE', 'HIP',                       # amber variants
    'LYSN', 'ARGN', 'ASPH', 'GLUH', 'HISD', 'HISE', 'HISH', 'HIS1'  # gromacs variants
}

# Biopython helpers
def highest_occupancy_altlocs(pdb_file, global_log) -> List[str]:
    """
    Reads a PDB file and returns a list of the highest occupancy alternative locations
    for each residue that has multiple conformations (alternative locations). 

    The output is a list where each element is a string in the format: 
    "<chain_id><residue_number>:<altLoc>", representing the chain, residue number, 
    and the alternative location identifier with the highest occupancy.

    Args:
        pdb_file (str): Path to the PDB file to be parsed.

    Returns:
        List[str]: A list of strings where each string represents the residue's chain ID, 
                   residue number, and the highest occupancy alternative location identifier.
                   The format of each string is "<chain_id><residue_number>:<altLoc>".

    Example:
        If a residue with ID 339 in chain 'A' has two alternative locations 'A' and 'B',
        and 'A' has a higher occupancy, the output for this residue would be "A339:A".

        >>> highest_occupancy_altlocs('example.pdb')
        ["A339:A", "A171:B", "A768:A"]
    """
    
    parser = PDBParser(QUIET=True)
    
    # Check if the file exists
    if not os.path.exists(pdb_file):
        global_log.error(f"File {pdb_file} not found")
        return []
    
    structure = parser.get_structure('structure', pdb_file)

    altloc_residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                
                altloc_dict = {}
                for atom in residue:
                    altloc = atom.get_altloc()
                    
                    if altloc and altloc != " ":  # Check if the atom has an alternative location
                        # Keep track of highest occupancy for each altloc
                        if altloc not in altloc_dict or atom.get_occupancy() > altloc_dict[altloc]['occupancy']:
                            altloc_dict[altloc] = {'occupancy': atom.get_occupancy(), 'residue': residue}
                            
                # If there is any alternative location
                if altloc_dict:
                    # Find the altloc with the highest occupancy
                    best_altloc = max(altloc_dict, key=lambda x: altloc_dict[x]['occupancy'])
                    res_id = residue.get_id()[1]
                    chain_id = chain.get_id()
                    altloc_residues.append(f"{chain_id}{res_id}:{best_altloc}")
    
    # Log result
    if altloc_residues:
        global_log.info(f"Found residues with alternative locations: {altloc_residues}")
    else:
        global_log.info("No residues with alternative locations found")
        
    return altloc_residues
      
def fasta_from_pdb(input_pdb_path: str, output_fasta_path: str, global_log) -> bool:
    """
    Try to obtain the FASTA sequence using the SEQRES records in the PDB file with Biopython. If the SEQRES records are available, 
    write the FASTA sequence to the output file and return True. If the SEQRES records are not available, return False.
    
    Inputs
    ------
    
        input_pdb_path (str): Path to the input PDB file.
        output_fasta_path (str): Path to the output FASTA file.
    
    Returns
    -------
    
        bool: Whether the FASTA sequence was obtained from SEQRES records or not.
    """
    
    # Open the PDB file and use the PdbSeqresIterator to extract sequences
    with open(input_pdb_path, "r") as handle:
        sequences = list(PdbSeqresIterator(handle))  # Extract all sequences into a list
        
        if not sequences:
            global_log.warning(f"PDB doesn't contain SEQRES records")
            return False
        
        global_log.info(f"PDB does contain SEQRES records {sequences}")
        
        # Find parent folder of output_fasta_path
        parent_folder = os.path.dirname(output_fasta_path)
        
        # Create parent folder if it does not exist
        if parent_folder and not os.path.exists(parent_folder):
            os.makedirs(parent_folder)
        
        # Write sequences to a FASTA file
        with open(output_fasta_path, "w") as fasta_out:
            SeqIO.write(sequences, fasta_out, "fasta")
    
    return True

def get_pdb_code(pdb_file: str) -> Union[str, None]:
    """
    Retrieve the PDB code from a PDB file, if available.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file.

    Returns
    -------
    str or None:
        The PDB code if available, otherwise None.
    """
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    # Access the header information
    header = structure.header

    # Retrieve the PDB code (idcode)
    pdb_code = header.get('idcode', None)
    
    # If empty, set to None
    if pdb_code is not None:
        pdb_code = pdb_code.strip()
        if len(pdb_code) == 0:
            pdb_code = None

    return pdb_code

# Other helpers
def rename_ss_bonds(pdb_file: str, format: Literal['standard', 'gromacs', 'amber']) -> str:
    """ 
    Read a PDB file and convert the nomenclature of cysteine residues involved 
    in disulfide bonds according to the specified target format.
    
    Format Conversions:
        - amber:    Converts 'CYS2' to 'CYX'
        - gromacs:  Converts 'CYX' to 'CYS2'
        - standard: Converts 'CYS2' or 'CYX' to 'CYS'
    
    Note:
        CYX residues are recognized as cysteine residues with a disulfide bond 
        only from GROMACS version 2024. Older versions require CYS2 or CYS 
        as the cysteine name to form bonds in the topology. 
        See https://gitlab.com/gromacs/gromacs/-/work_items/4929
    
    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file.
    format : Literal['standard', 'gromacs', 'amber']
        Target formatting convention for the disulfide-bonded CYS residues.
    
    Returns
    -------
    str
        Path to the newly generated PDB file with the renamed residues.
    """
    
    # Check requested format
    if format not in ['standard', 'gromacs', 'amber']:
        raise ValueError("Format must be 'standard', 'gromacs' or 'amber'")
    
    parent_path = Path(pdb_file).parent
    pdb_name = Path(pdb_file).stem
    new_pdb_file = os.path.join(parent_path, f"{pdb_name}_CYX.pdb")
    
    # Parse the PDB manually
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    # Replace residues based on format
    with open(new_pdb_file, 'w') as f:
        for line in lines:
            # Check if line is an ATOM/HETATM record to avoid modifying HEADER/REMARK lines
            if len(line) > 20 and (line.startswith("ATOM") or line.startswith("HETATM")):
                # Read columns 18-21 (4 characters)
                resname = line[17:21]
                
                if format == 'amber':
                    if resname == 'CYS2':
                        line = line[:17] + 'CYX ' + line[21:]
                elif format == 'gromacs':
                    if resname.strip() == 'CYX':
                        line = line[:17] + 'CYS2' + line[21:]
                elif format == 'standard':
                    if resname in ['CYS2', 'CYX ']:
                        line = line[:17] + 'CYS ' + line[21:]
                    
            f.write(line)
    
    return new_pdb_file

def rename_his(pdb_file: str, format: Literal['standard', 'gromacs', 'amber']) -> str:
    """ 
    Read a pdb file and rename HIS residues to the corresponding format.
    
    Formats:
    
        - gromacs:
            HID -> HISD
            HIE -> HISE
            HIP -> HISH
        - amber:
            HISD -> HID  
            HISE -> HIE 
            HISH -> HIP 
            HIS1 -> HIS
        - standard:
            HID -> HIS
            HIE -> HIS
            HIP -> HIS
            HISD -> HIS
            HISE -> HIS
            HISH -> HIS
            
        HIS -> HIS
    
    Parameters
    ----------
    
    pdb_file : str
        Path to the PDB file.
    format : str
        Format to rename the HIS residues. Can be 'standard', 'gromacs' or 'amber'
    
    Outputs
    -------
    
    new_pdb_file : str
        Path to the new PDB file with the renamed HIS residues.
    """
    
    # Pad Amber names to 4 characters to maintain strict slicing alignment
    gmx_names   = ['HISD', 'HISE', 'HISH', 'HIS1']
    amber_names = ['HID ', 'HIE ', 'HIP ', 'HIS ']
    
    if format not in ['standard', 'gromacs', 'amber']:
        raise ValueError("Format must be 'standard', 'gromacs' or 'amber'")
    
    parent_path = Path(pdb_file).parent
    pdb_name = Path(pdb_file).stem
    new_pdb_file = os.path.join(parent_path, f"{pdb_name}_HIS.pdb")
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    with open(new_pdb_file, 'w') as f:
        for line in lines:
            if len(line) > 20 and (line.startswith("ATOM") or line.startswith("HETATM")):
                pdb_resname = line[17:21]
                stripped_resname = pdb_resname.strip()
                
                if format == 'amber':
                    if pdb_resname in gmx_names:
                        new_name = amber_names[gmx_names.index(pdb_resname)]
                        line = line[:17] + new_name + line[21:]
                elif format == 'gromacs':
                    # Allow matching of both "HID " and "HID"
                    if pdb_resname in amber_names:
                        new_name = gmx_names[amber_names.index(pdb_resname)]
                        line = line[:17] + new_name + line[21:]
                    elif stripped_resname in [n.strip() for n in amber_names]:
                        idx = [n.strip() for n in amber_names].index(stripped_resname)
                        line = line[:17] + gmx_names[idx] + line[21:]
                elif format == 'standard':
                    if (pdb_resname in gmx_names) or (stripped_resname in [n.strip() for n in amber_names]):
                        line = line[:17] + 'HIS ' + line[21:]
                        
            f.write(line)
            
    return new_pdb_file

def rename_ter(pdb_file: str, format: Literal['standard', 'gromacs', 'amber']) -> str:
    """ 
    Rename atoms from terminal residues (ACE/NME) to the corresponding format.
    
    Formats: 
        - amber: 
            ACE: CH3, C, O
            NME: N, CH3
        - standard: 
            ACE: CA, C, O
            NME: N, CA 
            
    Parameters
    ----------
    
    pdb_file : str
        Path to the PDB file.
    format : str
        Format to rename the terminal residues. Can be 'amber' or 'standard'.
    
    Outputs
    -------
    
    new_pdb_file : str
        Path to the new PDB file with the renamed terminal residues.
    """
    
    if format not in ['standard', 'gromacs', 'amber']:
        raise ValueError("Format must be 'standard', 'gromacs' or 'amber'")
        
    parent_path = Path(pdb_file).parent
    pdb_name = Path(pdb_file).stem
    new_pdb_file = os.path.join(parent_path, f"{pdb_name}_{format}.pdb")
    
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    
    with open(new_pdb_file, 'w') as f:
        for line in lines:
            if len(line) > 20 and (line.startswith("ATOM") or line.startswith("HETATM")):
                # Residue name is strictly 3 chars for ACE/NME
                pdb_resname = line[17:20].strip()
                # Atom name sits in columns 13-16
                pdb_atomname = line[12:16].strip()
                
                if pdb_resname in ["ACE", "NME"]:
                    if format in ["amber", "gromacs"]:
                        if pdb_atomname == "CA":
                            # Replace strictly within indices 12:16. Pad to length 4.
                            line = line[:12] + " CH3" + line[16:]
                    elif format == "standard":
                        if pdb_atomname == "CH3":
                            line = line[:12] + " CA " + line[16:]
                            
            f.write(line)

    return new_pdb_file

# Preservation of user-provided special residue names (SS bonds / protonation states)
def capture_special_resnames(pdb_file: str) -> Dict:
    """
    Record the residues in a PDB file that carry a user-provided special residue name:
    disulfide cysteines (CYX/CYS2) and non-standard titratable variants (HID/HIE/HIP,
    HISD/HISE/HISH/HIS1, ASH/ASPH, GLH/GLUH, LYN/LYSN, ARGN).

    These names are captured before the workflow standardizes them for the Modeller /
    structure-checking steps, so they can be restored at the end.

    Parameters
    ----------
    pdb_file : str
        Path to the PDB file to parse.

    Returns
    -------
    Dict
        Mapping {(chain_id, resnum): resname} keyed by residue identity, using the residue
        numbering present in the input file.
    """
    captured = {}
    with open(pdb_file, 'r') as f:
        for line in f:
            if len(line) > 26 and (line.startswith("ATOM") or line.startswith("HETATM")):
                resname = line[17:21].strip()
                if resname in ss_bond_resnames or resname in nonstd_titra_resnames:
                    chain_id = line[21]
                    resnum = line[22:26].strip()
                    if resnum:
                        captured[(chain_id, int(resnum))] = resname
    return captured

def translate_special_resnames(captured: Dict, mapping_json_path: str, global_log) -> Dict:
    """
    Translate the residue keys captured by `capture_special_resnames` from the original
    numbering to the renumbered one produced by `renumber_structure` (step8).

    Parameters
    ----------
    captured : Dict
        Mapping {(chain_id, original_resnum): resname}.
    mapping_json_path : str
        Path to the mapping.json written by renumber_structure. Its "residues" entry maps
        chain -> {original_resnum(str): new_resnum(str)}.
    global_log :
        Logger object for logging messages.

    Returns
    -------
    Dict
        Mapping {(chain_id, new_resnum): resname}. Residues absent from the mapping are dropped.
    """
    import json
    with open(mapping_json_path, 'r') as f:
        residue_mapping = json.load(f).get("residues", {})

    translated = {}
    for (chain_id, resnum), resname in captured.items():
        new_resnum = residue_mapping.get(chain_id, {}).get(str(resnum))
        if new_resnum is None:
            global_log.warning(f"Residue {chain_id}{resnum} ({resname}) not found in renumbering map. Its name will not be preserved.")
            continue
        translated[(chain_id, int(new_resnum))] = resname
    return translated

def standardize_titra_variants(pdb_file: str) -> str:
    """
    Convert the ASP/GLU/LYS/ARG protonation variants to their standard 3-letter base name,
    so the Modeller / structure-checking steps (which do not recognize the gromacs 4-letter
    names) process these residues. HIS and CYS variants are handled separately by
    `rename_his` and `rename_ss_bonds`.

    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file.

    Returns
    -------
    str
        Path to the new PDB file with the standardized residue names.
    """
    variants = {'ASH', 'ASPH', 'GLH', 'GLUH', 'LYN', 'LYSN', 'ARGN'}

    parent_path = Path(pdb_file).parent
    pdb_name = Path(pdb_file).stem
    new_pdb_file = os.path.join(parent_path, f"{pdb_name}_stdtitra.pdb")

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(new_pdb_file, 'w') as f:
        for line in lines:
            if len(line) > 20 and (line.startswith("ATOM") or line.startswith("HETATM")):
                resname = line[17:21].strip()
                if resname in variants:
                    base = titra_variant_index[resname][0]
                    # Write the 3-letter base name padded to 4 chars (cols 18-21)
                    line = line[:17] + f"{base:<4}" + line[21:]
            f.write(line)

    return new_pdb_file

def restore_special_resnames(pdb_file: str, resmap: Dict, output_format: Literal['amber', 'gromacs'],
                             restore_ss: bool, restore_titra: bool, global_log) -> str:
    """
    Restore the user-provided special residue names captured earlier, expressed in the
    requested output format. Only the categories enabled by the flags are restored.

    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file (with residue numbering matching `resmap`).
    resmap : Dict
        Mapping {(chain_id, resnum): original_resname} of residues to restore.
    output_format : Literal['amber', 'gromacs']
        Target naming convention for the restored names.
    restore_ss : bool
        Whether to restore disulfide-bond cysteine names (CYX/CYS2).
    restore_titra : bool
        Whether to restore titratable protonation variant names.
    global_log :
        Logger object for logging messages.

    Returns
    -------
    str
        Path to the new PDB file with the restored residue names.
    """
    parent_path = Path(pdb_file).parent
    pdb_name = Path(pdb_file).stem
    new_pdb_file = os.path.join(parent_path, f"{pdb_name}_restored.pdb")

    with open(pdb_file, 'r') as f:
        lines = f.readlines()

    with open(new_pdb_file, 'w') as f:
        for line in lines:
            if len(line) > 26 and (line.startswith("ATOM") or line.startswith("HETATM")):
                chain_id = line[21]
                resnum = line[22:26].strip()
                key = (chain_id, int(resnum)) if resnum else None
                original = resmap.get(key)
                if original is not None:
                    target = None
                    if original in ss_bond_resnames:
                        if restore_ss:
                            target = ss_target_resname[output_format]
                    elif original in titra_variant_index:
                        if restore_titra:
                            base, idx = titra_variant_index[original]
                            target_names = titra_mapping[output_format].get(base)
                            if target_names is not None and idx < len(target_names):
                                target = target_names[idx]
                            else:
                                global_log.warning(f"Cannot express {original} ({chain_id}{resnum}) in {output_format} format. Leaving it as is.")
                    if target is not None:
                        # Write the target name padded to 4 chars (cols 18-21), preserving the rest
                        line = line[:17] + f"{target:<4}" + line[21:]
            f.write(line)

    return new_pdb_file

# Propka
def biobb_propka(input_structure_path: str, output_summary_path: str, properties: dict, global_log) -> None:
    """ 
    A mock function to simulate the behaviour of a biobb that uses propka
    
    Inputs:
    -------
        input_structure_path (str): Path to the input PDB file with the protein structure to predict the pKa of its aminoacids
        output_summary_path (str): Path to the output summary file with the pKa predictions
        protperties (dict): Dictionary with the properties of the function
        global_log: Logger object for logging messages
    
    NOTE: This will be changed by the true Biobb function when available
    """
    
    # Create output directory if it does not exist
    if not os.path.exists(properties["path"]):
        os.makedirs(properties["path"])
        
    # Check input structure path
    if not os.path.exists(input_structure_path):
        raise FileNotFoundError(f"Input structure file {input_structure_path} not found")
    
    # Run propka on the input structure
    mol = propka.run.single(input_structure_path)
    
    # Find the file name of the input structure
    input_structure_name = Path(input_structure_path).stem
    
    # Find the path to the propka summary file
    propka_summary_path = f"{input_structure_name}.pka"
    
    # Check the file exists
    if os.path.exists(propka_summary_path):
        # Move the propka summary file to the output path
        shutil.move(propka_summary_path, output_summary_path)
    else:
        global_log.error(f"Propka summary file {propka_summary_path} not found")
        
    return

def propka_summary(filepath) -> Dict:  
    """
    Parse the pKa summary section from a PROPKA output file.
    
    Inputs
    ------
    
        filepath (str): Path to the PROPKA output file.
    
    Returns
    -------
    
        Dict: A dictionary containing the parsed pKa data.
        
            {
                "resnum:chain": {
                    "pKa": float,
                    "model_pKa": float
                },
                ...
            }
    """
    
    results = {}
    in_summary = False

    with open(filepath) as f:
        for line in f:
            line = line.strip()

            # Detect start of summary
            if line.startswith("SUMMARY OF THIS PREDICTION"):
                in_summary = True
                continue
            
            # Skip header line
            if line.startswith("Group"):
                continue

            # Detect end of summary (first dashed line after the summary)
            if in_summary and line.startswith("----"):
                break

            if in_summary and line:
                # Format: RESNAME  RESNUM CHAIN  pKa  model-pKa
                parts = line.split()
                if len(parts) >= 5:
                    resname = parts[0]
                    resnum = int(parts[1])
                    chain = parts[2]
                    pKa = float(parts[3])
                    model_pKa = float(parts[4])
                    results.update({f"{resnum}:{chain}": { "pKa": pKa, "model_pKa": model_pKa }})
                
    return results

def add_reduce_his_resnames(pKa_data: Dict, pdb_file: str, global_log) -> Dict:
    """ 
    Read the histidine resnames in the PDB file (HID, HIE, HIP) and add it to the corresponding residue in the pKa data.
    
    The pdb_file is assumed to be the output of the reduce program in AmberTools. Which optimizes the H-bonds 
    of HIS to place the protons.
    
    Input
    -----
    
        pKa_data : 
            Dictionary with the parsed pKa data.
        pdb_file : 
            Path to the PDB file with the HIS resnames.
        global_log :
            Logger object for logging messages
            
    Returns
    -------
    
        pKa_data : 
            Updated list of dictionaries with the pKa data.
    """
    # Possible Histidine names coming from pdb4amber
    his_names = ['HIE', 'HID', 'HIP']
    
    # Parse the PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)
    
    any_missing_his = False
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Find the residue name
                res_name = residue.get_resname()
                # Find the residue number
                res_num = residue.get_id()[1]
                # Check if the residue is a histidine
                if res_name in his_names: 
                    res_key = f"{res_num}:{chain.get_id()}"
                    if res_key in pKa_data: 
                        pKa_data[res_key]['reduce_resname'] = res_name
                    else:
                        any_missing_his = True
                        
    if any_missing_his:
        global_log.warning("Some HIS residues found by reduce (AmberTools) were not found by propka. Check the residue names.")
            
    return pKa_data
 
def biobb_titrate(input_structure_path: str, output_structure_path: str, properties: dict, global_log) -> None:
    """ 
    A mock function to simulate the behaviour of a biobb that uses pKa estimations and the pH to determine the 
    protonation state of titratable residues. The names of such residues are changed according to their predicted
    protonation state.
    
    NOTE: Any residues with a four character name will keep the first three characters only - fix this if needed.
    
    Parameters 
    ----------
    
        input_structure_path  (str): Path to the original PDB file
        output_structure_path (str): Path to the new PDB file with the residue names changed
        properties           (dict): Dictionary with the properties of the function
        global_log                 : Logger object for logging messages
    """

    # Create output directory if it does not exist
    if not os.path.exists(properties["path"]):
        os.makedirs(properties["path"])
    
    # Check input structure path
    if not os.path.exists(input_structure_path):
        raise FileNotFoundError(f"Input structure file {input_structure_path} not found")
    
    # Properties
    pH = properties.get("pH", 7.0)
    manual_his = properties.get("his", None)
    pdb_format = properties.get("pdb_format", 'amber')
    pKa_data = properties.get("pKa_results", None)
    titra_resnames = titra_mapping[pdb_format]
    amber_his_resnames = titra_mapping['amber']['HIS']
    
    # Read the input PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_structure_path)
    
    his_count = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Identify residue
                res_name = residue.get_resname()
                res_num = residue.get_id()[1]
                chain_id = chain.get_id()
                res_key = f"{res_num}:{chain_id}"
                
                # If residue has a pKa value ...
                if res_key in pKa_data:
                    pKa = pKa_data[res_key]['pKa']
                    
                    # ... and is titratable for this format
                    if res_name in titra_resnames:
                        
                        if res_name == 'HIS':    
                            if manual_his is not None:
                                # Safely handle manual selection of protonation states for HIS
                                if his_count < len(manual_his):
                                    new_res_name = titra_resnames['HIS'][int(manual_his[his_count])]
                                else:
                                    global_log.warning(f"Not enough manual HIS states provided. Falling back to default for HIS {res_num}.")
                                    new_res_name = titra_resnames['HIS'][1] # Fallback to epsilon
                                his_count += 1
                            else:
                                protonated_his_resname = titra_resnames['HIS'][2]
                                # Automatic based on pKa
                                if pH < pKa:
                                    # Protonated
                                    new_res_name = protonated_his_resname
                                if pH >= pKa:
                                    # Deprotonated
                                    reduce_resname = pKa_data[res_key].get('reduce_resname')
                                    if (reduce_resname is not None) and (reduce_resname != protonated_his_resname):
                                        # Given by reduce: epsilon or delta
                                        new_res_name = titra_resnames['HIS'][amber_his_resnames.index(reduce_resname)]
                                    else:
                                        # Default to epsilon protonated
                                        new_res_name = titra_resnames['HIS'][1]
                        else:
                            if pH < pKa:
                                # Protonated
                                new_res_name = titra_resnames[res_name][1]
                            else:
                                # Deprotonated
                                new_res_name = titra_resnames[res_name][0]
                        
                        # Update the residue name
                        residue.resname = new_res_name
    
    # Write the modified structure to the output file
    io = PDBIO()
    io.set_structure(structure)
    tmp_structure_path = os.path.join(properties["path"], "tmp.pdb")
    io.save(tmp_structure_path, preserve_atom_numbering=True)          
    
    # Check format issues with 4 letter residue names
    with open(tmp_structure_path, 'r') as f:
        lines = f.readlines()
        
    with open(output_structure_path, 'w') as f:
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Biopython expands the %3s formatter to 4 chars for 4-letter resnames, 
                # which pushes the chain ID and everything after it 1 column to the right.
                # Columns 18-21 (indices 17:21) contain the resname.
                if line[17:21].strip() in gmx_4letter_resnames:
                    # Remove the shifted space at index 21 to pull the rest of the line back into strict PDB alignment
                    line = line[:21] + line[22:]
            f.write(line)
    
    # Remove the temporary file
    os.remove(tmp_structure_path)
    return
           
# YML construction
def config_contents(
    input_pdb_path : str,
    pdb_chains: Optional[List],
    mutation_list: Optional[List],
    modeller_key: Optional[str],
    cap_ter: bool = False,
    pH: float = 7.0,
    his: Optional[List] = None,
    keep_hs: bool = False,
    output_format: Literal['amber', 'gromacs'] = 'amber',
    restart: bool = False,
    debug: bool = False
    ) -> str:
    """
    Returns the contents of the YAML configuration file as a string.

    The YAML file contains the configuration for the protein preparation workflow.

    Returns
    -------
    str
        The contents of the YAML configuration file.
    """

    molecule_type = "all"
    step1_chains_property = ""

    # If chains are given as argument
    if pdb_chains is not None:
        molecule_type = "chains"
        step1_chains_property = f"chains : {to_yaml(pdb_chains)}"

    # Only emitted when mutations are requested (the mutate tool has no null default)
    step3_mutations_property = ""
    if mutation_list is not None:
        step3_mutations_property = f"mutation_list : {','.join(mutation_list)}"

    # Modeller key is injected into every step that supports it (null when not given);
    # use_modeller is enabled only where a key is available.
    use_modeller = to_yaml(modeller_key is not None)
    
    # Convert into absolute path
    input_pdb_path = os.path.abspath(input_pdb_path)

    return f"""
# Global properties (common for all steps)
global_properties:
  working_dir_path: output                                          # Workflow default output directory
  can_write_console_log: false                                      # Verbose writing of log information
  restart: {to_yaml(restart)}                                       # Skip steps already performed
  remove_tmp: {to_yaml(not debug)}                                  # Remove temporal files

# Step 1: extract atoms from input PDB file
step1_extractAtoms:
  tool: extract_molecule
  paths:
    input_structure_path: {input_pdb_path}
    output_molecule_path: main_structure.pdb
  properties:
    molecule_type: {molecule_type}
    {step1_chains_property}

# Step 2: fix alternative locations of residues if any with biobb_structure_checking and the Modeller suite (if key property is given)
step2_fixaltlocs:
  tool: fix_altlocs
  paths:
    input_pdb_path: dependency/step1_extractAtoms/output_molecule_path
    output_pdb_path: fixaltlocs.pdb
  properties:
    altlocs: null                                                    # Format: ["Chain Residue_number:Altloc"], e.g. # ["A339:A", "A171:B", "A768:A"]
    modeller_key : {to_yaml(modeller_key)}                           # MODELLER license key

# Step 3: Mutate residues in the structure if needed
step3_mutations:
  tool: mutate
  paths:
    input_pdb_path: dependency/step2_fixaltlocs/output_pdb_path
    output_pdb_path: mutated.pdb
  properties:
    modeller_key : {to_yaml(modeller_key)}
    use_modeller : {use_modeller}
    {step3_mutations_property}

# Step 4: Download a FASTA file with the canonical sequence of the protein
# It requires internet connection and a PDB code
step4_canonical_fasta:
  tool: canonical_fasta
  paths:
    output_fasta_path: canonicalFasta.fasta
  properties:
    pdb_code: null                                                    # Resolved at runtime (user arg or PDB header)

# Step 2 C: Extract the residue sequence from the PDB file to FASTA format
step5_pdb_tofasta:
  tool: biobb_pdb_tofasta
  paths:
    input_file_path: dependency/step1_extractAtoms/input_structure_path
    output_file_path: pdbFasta.fasta
  properties:
    multi: true

# Step 2 D: Model missing backbone atoms with biobb_structure_checking and the Modeller suite
# It requires a MODELLER license 
step6_fixbackbone:
  tool: fix_backbone
  paths:
    input_pdb_path: dependency/step3_mutations/output_pdb_path
    input_fasta_canonical_sequence_path: dependency/step4_canonical_fasta/output_fasta_path
    output_pdb_path: fixbackbone.pdb
  properties:
    add_caps: {to_yaml(cap_ter)}
    modeller_key : {to_yaml(modeller_key)}

# Step 2 E: Model missing side chain atoms with biobb_structure_checking and the Modeller suite (if key property is given)
step7_fixsidechain:
  tool: fix_side_chain
  paths:
    input_pdb_path: dependency/step6_fixbackbone/output_pdb_path
    output_pdb_path: fixsidechain.pdb
  properties:
    modeller_key : {to_yaml(modeller_key)}
    use_modeller : {use_modeller}

step8_renumberstructure:
  tool: renumber_structure
  paths:
    input_structure_path: dependency/step7_fixsidechain/output_pdb_path
    output_structure_path: renumbered.pdb
    output_mapping_json_path: mapping.json
  properties:
    renumber_residues: true
    renumber_residues_per_chain: false

# Step 2 G: Flip clashing amides with biobb_structure_checking and the Modeller suite
# Optional step (activate from command line with --fix_amides)
step9_fixamides:
  tool: fix_amides
  paths:
    input_pdb_path: dependency/step8_renumberstructure/output_structure_path
    output_pdb_path: fixamides.pdb
  properties:
    modeller_key : {to_yaml(modeller_key)}

step10_fixchirality:
  tool: fix_chirality
  paths:
    input_pdb_path: dependency/step9_fixamides/output_pdb_path
    output_pdb_path: fixchirality.pdb
  properties:
    modeller_key : {to_yaml(modeller_key)}

# Step 2 F: Fix disulfide bonds with biobb_structure_checking (CYS -> CYX for cysteines involved in disulfide bonds)
# Optional step (activate from command line with --fix_ss)
step11_fixssbonds:
  tool: fix_ssbonds
  paths:
    input_pdb_path: dependency/step10_fixchirality/output_pdb_path
    output_pdb_path: fixssbonds.pdb
  properties:
    modeller_key : {to_yaml(modeller_key)}  # MODELLER license key

step12_remove_hs:
  tool: remove_hydrogens
  paths:
    input_path: dependency/step11_fixssbonds/output_pdb_path
    output_path: remove_hs.pdb

# Use propka to predict the pKa of the titratable residues
step13_propka:
  paths: 
    input_structure_path: dependency/step12_remove_hs/output_path
    output_summary_path: summary.pka

# Use reduce to optimize H-bonds of HIS residues
step14_his_hbonds:
  tool: pdb4amber
  paths:
    input_pdb_path: dependency/step12_remove_hs/output_path
    output_pdb_path: his_hbonds.pdb
  properties:
    reduce: true

# Rename titratable residues according to protonation state
step15_titrate:
  paths:
    input_structure_path: dependency/step12_remove_hs/output_path
    output_structure_path: titrate.pdb
  properties:
    pH: {pH}
    his: {to_yaml(his)}
    pdb_format: {output_format}

# Put back hydrogens
step16_pdb4amber:
  tool: pdb4amber_run
  paths:
    input_pdb_path: dependency/step15_titrate/output_structure_path
    output_pdb_path: pdb4amber.pdb
  properties:
    reduce: {to_yaml(not keep_hs)}
"""

def create_config_file(output_path: str,
                       **config_args) -> str:
    """
    Create a YAML configuration file for the workflow in the output path.
    Return the path to the configuration file.
    
    Parameters
    ----------
    output_path : str
        Path to the output folder
    config_args : dict
        Arguments to be used in the configuration file.
    
    Returns
    -------
    
    str
        Path to configuration file
    """
    
    config_path = os.path.join(output_path, 'config.yml')
    
    # Write the contents to the file
    with open(config_path, 'w') as f:
        f.write(config_contents(**config_args))

    return config_path

# Main workflow
def protein_preparation(
    input_pdb_path: str,
    pdb_code: Optional[str] = None, 
    pdb_chains: Optional[List] = None, 
    mutation_list: Optional[List] = None, 
    skip_bc_fix: bool = False, 
    modeller_key: Optional[str] = None,
    cap_ter: bool = False,
    skip_sc_fix: bool = False, 
    skip_ss_bonds: bool = False, 
    skip_amides_flip: bool = False, 
    pH: float = 7.0,
    his: Optional[List] = None, 
    keep_hs: bool = False, 
    output_format: Literal['amber', 'gromacs'] = 'amber',
    restart: bool = False,
    debug: bool = False,
    output_path: Optional[str] = None
    ):
    '''
    Protein preparation workflow. Can be used to fix some PDB defects of the structure, add specific mutations 
    and protonate the titratable residues at a given pH.
    
    Inputs
    ------

        input_pdb_path: 
            path to input PDB file
        pdb_code: 
            PDB code to be used to get the canonical FASTA sequence
        pdb_chains: 
            list of chains to be extracted from the PDB file and fixed
        mutation_list: 
            list of mutations to be introduced in the structure
        skip_bc_fix: 
            whether to skip the fix of the backbone atoms. Default: False.
        modeller_key:
            Modeller key to be used to model the missing atoms in the structure. 
            If None, the Modeller suite will not be used. Default: None.
        cap_ter: 
            whether to add caps to the terminal residues. Default: False.
        skip_sc_fix:
            Skip the side chain modeling of missing atoms. Otherwise the 
            missing atoms in the side chains of the PDB structure will be modeled using 
            'biobb_structure_checking' and the 'Modeller suite' (if the Modeller key is given). 
            Default: False
        skip_ss_bonds:
            Skip the automatic distance-based detection of disulfide bonds. When set, disulfide
            cysteines already marked in the input PDB (CYX for amber, CYS2 for gromacs) are preserved
            in the output; if none are marked, no disulfide bonds are formed. Default: False.
        skip_amides_flip: 
            whether to flip clashing amides to relieve the clashes
        pH: 
            pH of the system. Used together with a pKa estimation (propka) to determine the 
            protonation state of titratable residues. Default: 7.0
        his: 
            Manual selection of histidine protonation states (HID: 0, HIE: 1, HIP:2). 
            If given, the pKa estimation and the pH won't be used to protonate histidine residues. 
            Default: None. Example: '0 1 1'
        keep_hs:
            Keep the hydrogen atoms in the input PDB file (otherwise they are removed). When set, the
            protonation state determination (propka/reduce/titrate) is skipped and the titratable
            protonation residue names marked in the input are preserved in the output. Default: False
        output_format:
            Output format of the PDB file. Can be 'amber' or 'gromacs'. Default: 'amber'.
        debug:
            whether to keep temporary files. Default: False
        output_path:
            path to output folder

    Outputs
    -------

        /output folder
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''
    
    start_time = time.time()
    
    # Determine final output path
    output_path = fu.get_working_dir_path(output_path, restart=restart)
    
    # Initialize a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)
    global_log.info(f"biobb_md_workflows version {__version__}")
    
    # Create the configuration file
    config_args = {
        'input_pdb_path' : input_pdb_path,
        'pdb_chains' : pdb_chains,
        'mutation_list' : mutation_list,
        'modeller_key' : modeller_key,
        'cap_ter': cap_ter,
        'pH': pH,
        'his': his,
        'keep_hs': keep_hs,
        'output_format' : output_format,
        'restart': restart,
        'debug': debug
    }
    configuration_path = create_config_file(output_path, **config_args)

    # Receiving the input configuration file (YAML)
    conf = settings.ConfReader(config=configuration_path)
    conf.working_dir_path = output_path

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    # STEP 1: extract main structure of interest while removing water and ligands (heteroatoms)
    global_log.info("step1_extractAtoms: extract chain of interest")
    extract_molecule(**global_paths["step1_extractAtoms"], properties=global_prop["step1_extractAtoms"])
    
    # STEP 2: Fix alternative locations
    global_log.info("step2_fixaltlocs: Fix alternative locations")
    global_prop["step2_fixaltlocs"]["altlocs"] = highest_occupancy_altlocs(global_paths["step1_extractAtoms"]["input_structure_path"], global_log)
    fix_altlocs(**global_paths["step2_fixaltlocs"], properties=global_prop["step2_fixaltlocs"])

    # Capture user-provided special residue names (SS bonds / protonation variants) before
    # standardizing, so they can be restored at the end when --skip_ss_bonds / --keep_hs are set.
    special_resnames = {}
    if skip_ss_bonds or keep_hs:
        special_resnames = capture_special_resnames(global_paths["step2_fixaltlocs"]["output_pdb_path"])

    # Standardize residue names for the Modeller / structure-checking steps, which only
    # recognize standard 3-letter names (gromacs 4-letter variants would be skipped).
    last_pdb_path = rename_his(global_paths["step2_fixaltlocs"]["output_pdb_path"], 'standard')
    last_pdb_path = rename_ss_bonds(last_pdb_path, 'standard')
    last_pdb_path = standardize_titra_variants(last_pdb_path)
    
    # STEP 3: Add mutations if requested
    global_paths["step3_mutations"]["input_pdb_path"] = last_pdb_path
    mutate(**global_paths["step3_mutations"], properties=global_prop["step3_mutations"])
    last_pdb_path = global_paths["step3_mutations"]["output_pdb_path"]
    
    # Model the backbone atoms
    if not skip_bc_fix:
        
        # STEP 4: Try to get the FASTA sequence to model the backbone
        fasta_available = False
        
        # Find the PDB code of the input PDB file
        file_pdb_code = get_pdb_code(global_paths["step1_extractAtoms"]["input_structure_path"])
        
        # Check for conflicts but strictly enforce user precedence
        if pdb_code is not None and file_pdb_code is not None:
            if pdb_code.upper() != file_pdb_code.upper():
                global_log.warning(f"step4_canonical_fasta: Provided PDB code ({pdb_code}) differs from the one in the input PDB file header ({file_pdb_code}).")
                global_log.warning(f"step4_canonical_fasta: User-provided PDB code ({pdb_code}) will take precedence.")
                # pdb_code remains exactly what the user provided
                
        # If user didn't provide one, fallback to the file's metadata
        elif pdb_code is None and file_pdb_code is not None:
            global_log.info(f"step4_canonical_fasta: Using PDB code ({file_pdb_code}) extracted from the input PDB file.")
            pdb_code = file_pdb_code
            
        # Optional: Log if the user provided one and the file had none
        elif pdb_code is not None and file_pdb_code is None:
            global_log.info(f"step4_canonical_fasta: Using user-provided PDB code ({pdb_code}). No PDB code found in file header.")
            
        # If we have a PDB code, get the FASTA sequence from an http request to the PDB
        if pdb_code is not None:
            try:
                global_log.info("step4_canonical_fasta: Get canonical FASTA")
                global_prop["step4_canonical_fasta"]["pdb_code"] = pdb_code
                canonical_fasta(**global_paths["step4_canonical_fasta"], properties=global_prop["step4_canonical_fasta"])
                fasta_available = True
            except:
                global_log.warning("step4_canonical_fasta: Could not get canonical FASTA. Check the internet connection on the machine running the workflow. Trying to get the canonical FASTA from the PDB file...")
                fasta_available = False
        else:
            global_log.warning("step4_canonical_fasta: No PDB code found. Trying to get the FASTA from the PDB file...")
            fasta_available = False
        
        # ... from SEQRES records in the PDB file
        if not fasta_available:
            global_log.info("step5_pdb_tofasta: Get FASTA from SEQRES of PDB file")
            fasta_available = fasta_from_pdb(global_paths["step1_extractAtoms"]["input_structure_path"], global_paths["step5_pdb_tofasta"]["output_file_path"], global_log)
            global_paths['step6_fixbackbone']['input_fasta_canonical_sequence_path'] = global_paths['step5_pdb_tofasta']['output_file_path']
            
        # ... from the residues in the PDB file (not canonical)
        if not fasta_available:
            global_log.info("step5_pdb_tofasta: Get FASTA from PDB file")
            # Update the input file path
            global_paths['step5_pdb_tofasta']['input_file_path'] = global_paths["step1_extractAtoms"]["input_structure_path"]
            # Only existing residues in the PDB file are included
            biobb_pdb_tofasta(**global_paths["step5_pdb_tofasta"], properties=global_prop["step5_pdb_tofasta"])
            global_paths['step6_fixbackbone']['input_fasta_canonical_sequence_path'] = global_paths['step5_pdb_tofasta']['output_file_path']
            fasta_available = True
            
        # STEP 6: Model missing heavy atoms of backbone
        if fasta_available:
            global_log.info("step6_fixbackbone: Modeling the missing heavy atoms in the structure side chains")
            fix_backbone(**global_paths["step6_fixbackbone"], properties=global_prop["step6_fixbackbone"])
            last_pdb_path = global_paths["step6_fixbackbone"]["output_pdb_path"]
        else:
            global_log.warning("step6_fixbackbone: Could not get FASTA sequence. Skipping modeling of the missing heavy atoms in the backbone.")
    else:
        global_log.info("step6_fixbackbone: Skipping modeling of the missing heavy atoms in the backbone")

    # STEP 7: model missing heavy atoms of side chains NOTE: does this step erase unknown atoms?
    if not skip_sc_fix:
        global_log.info("step7_fixsidechain: Modeling the missing heavy atoms in the structure side chains")
        global_paths['step7_fixsidechain']['input_pdb_path'] = last_pdb_path
        fix_side_chain(**global_paths["step7_fixsidechain"], properties=global_prop["step7_fixsidechain"])
        last_pdb_path = global_paths["step7_fixsidechain"]["output_pdb_path"]
    else:
        global_log.info("step7_fixsidechain: Skipping modeling of the missing heavy atoms in the side chains")

    # STEP 8: Renumber structure atoms and residues
    global_log.info("step8_renumberstructure: renumber structure")
    global_paths['step8_renumberstructure']['input_pdb_path'] = last_pdb_path
    renumber_structure(**global_paths["step8_renumberstructure"], properties=global_prop["step8_renumberstructure"])
    last_pdb_path = global_paths["step8_renumberstructure"]["output_structure_path"]

    # Translate captured special-residue keys to the renumbered identities produced by step8
    if special_resnames:
        special_resnames = translate_special_resnames(special_resnames, global_paths["step8_renumberstructure"]["output_mapping_json_path"], global_log)

    # STEP 9: Flip amide groups to relieve clashes
    if not skip_amides_flip:
        global_log.info("step9_fixamides: fix clashing amides")
        global_paths['step9_fixamides']['input_pdb_path'] = last_pdb_path
        fix_amides(**global_paths["step9_fixamides"], properties=global_prop["step9_fixamides"])
        last_pdb_path = global_paths["step9_fixamides"]["output_pdb_path"]
    else:
        global_log.info("step9_fixamides: Skipping fixing clashing amides")

    # STEP 10: Fix chirality
    global_paths['step10_fixchirality']['input_pdb_path'] = last_pdb_path
    global_log.info("step10_fixchirality: fix chirality of residues")
    fix_chirality(**global_paths["step10_fixchirality"], properties=global_prop["step10_fixchirality"])
    
    # STEP 11: model SS bonds (CYS -> CYX)
    if not skip_ss_bonds:
        global_log.info("step11_fixssbonds: Fix SS bonds")
        global_paths['step11_fixssbonds']['input_pdb_path'] = last_pdb_path
        fix_ssbonds(**global_paths["step11_fixssbonds"], properties=global_prop["step11_fixssbonds"])
        last_pdb_path = global_paths["step11_fixssbonds"]["output_pdb_path"]
    else:
        global_log.info("step11_fixssbonds: Skipping modeling of the SS bonds")
        
    # STEP 12: Remove hydrogens
    if not keep_hs:
        global_log.info("step12_remove_hs: Remove hydrogens")
        global_paths['step12_remove_hs']['input_path'] = last_pdb_path
        reduce_remove_hydrogens(**global_paths["step12_remove_hs"], properties=global_prop["step12_remove_hs"])
        last_pdb_path = global_paths["step12_remove_hs"]["output_path"]
    else:
        global_log.info("step12_remove_hs: Skipping removal of hydrogens")
        
    # Rename atoms from terminal residues (ACE/NME) to amber format
    last_pdb_path = rename_ter(last_pdb_path, format='amber')

    # STEPS 13-15: Determine protonation states of titratable residues.
    # Skipped when keeping input hydrogens, so the input protonation is left unchanged.
    if not keep_hs:
        # STEP 13: Estimate pKa of titratable residues with propka
        # We need standardized residue names for propka! Otherwise it will skip them
        global_log.info("step13_propka: Estimate protonation state of titratable residues from empirical pKa calculation with propka")
        global_paths["step13_propka"]["input_structure_path"] = last_pdb_path
        biobb_propka(**global_paths["step13_propka"], properties=global_prop["step13_propka"], global_log=global_log)
        pKa_results = propka_summary(global_paths["step13_propka"]["output_summary_path"])

        # STEP 14: Run reduce from Amber Tools to estimate optimal proton placement in histidines
        global_log.info("step14_his_hbonds: Estimate optimal proton placement in Histidines")
        global_paths["step14_his_hbonds"]["input_pdb_path"] = last_pdb_path
        pdb4amber_run(**global_paths["step14_his_hbonds"], properties=global_prop["step14_his_hbonds"])
        pKa_results = add_reduce_his_resnames(pKa_results, global_paths["step14_his_hbonds"]["output_pdb_path"], global_log)

        # STEP 15: Choose protonation states for titratable residues
        global_log.info("step15_titrate: Choose protonation states for titratable residues")
        global_paths["step15_titrate"]["input_structure_path"] = last_pdb_path
        global_prop["step15_titrate"]["pKa_results"] = pKa_results
        biobb_titrate(**global_paths["step15_titrate"], properties=global_prop["step15_titrate"], global_log=global_log)
        last_pdb_path = global_paths["step15_titrate"]["output_structure_path"]

        # NOTE: Manually add Glutamine protonation states - assume model pKa value
    else:
        global_log.info("step13_propka/step14_his_hbonds/step15_titrate: Keeping input hydrogens, skipping protonation state determination")

    # Restore user-provided special residue names: disulfides (CYX/CYS2) when --skip_ss_bonds,
    # and titratable protonation variants when --keep_hs. Done here so the numbering still
    # matches step8; restored names survive the downstream pdb4amber / gromacs renaming.
    if special_resnames and (skip_ss_bonds or keep_hs):
        global_log.info("Restoring user-provided residue names (SS bonds / protonation states)")
        last_pdb_path = restore_special_resnames(last_pdb_path, special_resnames, output_format,
                                                 restore_ss=skip_ss_bonds, restore_titra=keep_hs,
                                                 global_log=global_log)

    if output_format == 'amber':
        # STEP 16: run pdb4amber to generate final PDB file (rebuilds hydrogens unless keep_hs)
        global_log.info("step16_pdb4amber: Generate final PDB file with hydrogens")
        # Point pdb4amber at the last produced structure: the config dependency references
        # step15_titrate, which is skipped under --keep_hs and bypassed by the restore step.
        global_paths["step16_pdb4amber"]["input_pdb_path"] = last_pdb_path
        pdb4amber_run(**global_paths["step16_pdb4amber"], properties=global_prop["step16_pdb4amber"])
        last_pdb_path = global_paths["step16_pdb4amber"]["output_pdb_path"]
    elif output_format == 'gromacs':
        # Rename atoms to gromacs format
        last_pdb_path = rename_ss_bonds(last_pdb_path, format='gromacs')
    
    # NOTE: We should make sure that PDB complies with the PDB format (while considering 4 letter resnames)
    
    # Copy the final PDB file to the output path
    final_pdb_path = os.path.join(output_path, f"{Path(input_pdb_path).stem}.pdb")
    shutil.copy(last_pdb_path, final_pdb_path)

    # Write a stable output manifest for external consumers (see manifest.yaml in output_path)
    manifest_outputs = {
        "final_structure": {
            "format": output_format,
            "pdb": os.path.relpath(final_pdb_path, output_path),
        },
    }
    mapping_json_path = global_paths["step8_renumberstructure"]["output_mapping_json_path"]
    if os.path.exists(mapping_json_path):
        manifest_outputs["renumbering_map"] = {
            "json": os.path.relpath(mapping_json_path, output_path),
        }
    with open(os.path.join(output_path, "manifest.yaml"), "w") as manifest_file:
        yaml.safe_dump({"schema_version": 1, "outputs": manifest_outputs}, manifest_file, sort_keys=False)

    # Print timing information to log file
    elapsed_time = time.time() - start_time
    global_log.info('')
    global_log.info('')
    global_log.info('Execution successful: ')
    global_log.info('  Workflow_path: %s' % output_path)
    global_log.info('  Config File: %s' % configuration_path)
    global_log.info('')
    global_log.info('Elapsed time: %.1f minutes' % (elapsed_time/60))
    global_log.info('')

    return global_paths, global_prop


def main():

    parser = argparse.ArgumentParser("Protein preparation with biobb_structure_checking and Modeller")

    parser.add_argument('--input_pdb', dest='input_pdb_path', type=str,
                        help="Input PDB file.",
                        required=True)

    parser.add_argument('--pdb_code', dest='pdb_code', type=str,
                        help="""PDB code to get the canonical FASTA sequence of the input PDB file. 
                        If not given the workflow will look for it in the HEADER of the PDB. Default: None""",
                        required=False)

    parser.add_argument('--pdb_chains', nargs='+', dest='pdb_chains',
                        help="Protein PDB chains to be extracted from PDB file and fixed. Default: all chains in PDB. Example: A B C",
                        required=False)

    parser.add_argument('--mutation_list', nargs='+', dest='mutation_list',
                        help="List of mutations to be introduced in the protein. Default: None. Example: A:Arg220Ala B:Ser221Gly",
                        required=False)

    parser.add_argument('--skip_bc_fix', action='store_true', dest='skip_bc_fix', 
                        help="""Skip the modeling of missing backbone atoms. Otherwise the missing atoms in the backbone 
                        of the PDB structure will be modeled using 'biobb_structure_checking' and the 'Modeller suite' 
                        (if the Modeller key is given). Note that modeling of missing loops is only possible if the Modeller 
                        key is provided. To obtain one register at: https://salilab.org/modeller/registration.html. Default: False""",
                        required=False, default=False)
    
    parser.add_argument('--modeller_key', dest='modeller_key', type=str,
                        help="""Modeller key to be used for the backbone modeling of missing atoms. Note that modeling of missing 
                        loops is only possible if the Modeller key is provided.
                        To obtain one register at: https://salilab.org/modeller/registration.html""",
                        required=False, default=None)
    
    parser.add_argument('--cap_ter', action='store_true', dest='cap_ter',
                        help="Cap the termini of the PDB file with ACE and NME residues, preserving existing atoms. Default: False",
                        required=False, default=False)
    
    parser.add_argument('--skip_sc_fix', action='store_true', dest='skip_sc_fix',
                        help="""Skip the side chain modeling of missing atoms. Otherwise the missing atoms in the side chains 
                        of the PDB structure will be modeled using 'biobb_structure_checking' and the 'Modeller suite' 
                        (if the Modeller key is given). Default: False""",
                        required=False, default=False)
    
    parser.add_argument('--skip_ss_bonds', action='store_true',
                        help="""Skip the automatic detection of disulfide bonds based on a distance criteria.
                        When set, any disulfide cysteines the user already marked in the input PDB
                        (CYX for AMBER, CYS2 for GROMACS) are preserved in the output; if none are marked,
                        no disulfide bonds are formed. Default: False""",
                        required=False, default=False)

    parser.add_argument('--skip_amides_flip', action='store_true', dest='skip_amides_flip',
                        help="""Skip the fliping of clashing amide groups in ASP or GLU residues.
                        Otherwise the amide orientations will be changed if needed to relieve clashes using
                        'biobb_structure_checking'. Note that amide group orientations coming from PDB structures is 
                        not reliable in general due to symmetries in the electron density. Default: False""",
                        required=False, default=False)

    parser.add_argument('--ph', dest='ph', type=float,
                        help="""pH of the system. Used together with a pKa estimation (with propka) to determine the 
                        protonation state of titratable residues. Default: 7.0""",
                        required=False, default=7.0)
    
    parser.add_argument('--his', nargs='+', dest='his',
                        help="""Manual selection of histidine protonation states (delta: 0, epsilon: 1, fully protonated: 2, 
                        bound to heme: 3). If given, the pKa estimation and the pH won't be used to protonate histidine residues. 
                        Default: None. Example: 0 1 1""",
                        required=False)
    
    parser.add_argument('--keep_hs', action='store_true',
                        help="""Keep the hydrogen atoms in the input PDB file (otherwise they are discarded).
                        When set, the protonation state determination (propka/reduce/titrate) is skipped and the
                        titratable protonation residue names marked in the input (e.g. HID/HIE/HIP, ASH, GLH, LYN
                        for AMBER; HISD/HISE/HISH, ASPH, GLUH, LYSN, ARGN for GROMACS) are preserved in the output.
                        Default: False""",
                        required=False)

    parser.add_argument('--output_format', dest='output_format', type=str,
                        help="""PDB format to be used. Options: amber, gromacs. Default: 'amber'""",
                        required=False, default='amber')

    parser.add_argument('--restart', action='store_true',
                        help="Restart the workflow from the last completed step. Default: False",
                        required=False, default=False)

    parser.add_argument('--debug', action='store_true',
                        help="Keep temporary files. Default: False",
                        required=False, default=False)

    parser.add_argument('--output', dest='output_path', type=str,
                        help="Output path. Default: 'output' in the current working directory",
                        required=False, default='output')

    args = parser.parse_args()
    
    protein_preparation(
        input_pdb_path=args.input_pdb_path, 
        pdb_code=args.pdb_code, 
        pdb_chains=args.pdb_chains, 
        mutation_list=args.mutation_list, 
        skip_bc_fix=args.skip_bc_fix, 
        modeller_key=args.modeller_key,
        cap_ter=args.cap_ter,
        skip_sc_fix=args.skip_sc_fix, 
        skip_ss_bonds=args.skip_ss_bonds, 
        skip_amides_flip=args.skip_amides_flip, 
        pH=args.ph,
        his=args.his,
        keep_hs=args.keep_hs, 
        output_format=args.output_format,
        restart=args.restart,
        debug=args.debug,
        output_path=args.output_path)


if __name__ == '__main__':
    main()
