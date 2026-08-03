#!/usr/bin/env python3

from typing import List, Optional
from pathlib import Path
import argparse
import time
import os
import yaml

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_md_workflows.common import (
    add_group,
    build_solvent_selection,
    get_atom_types,
    get_central_atom_index,
    get_residue_types,
    ions_library,
    output_group,
    read_group_indices,
    read_groups,
    rename_last_ndx_group,
    solute_group,
    solvent_group,
    solvent_library,
    to_yaml,
)
from biobb_md_workflows import __version__

from biobb_gromacs.gromacs.make_ndx import make_ndx
from biobb_gromacs.gromacs.editconf import editconf
from biobb_gromacs.gromacs.convert_tpr import convert_tpr
from biobb_analysis.gromacs.gmx_image import gmx_image
from biobb_analysis.gromacs.gmx_trjconv_trj import gmx_trjconv_trj
from biobb_analysis.gromacs.gmx_trjconv_str import gmx_trjconv_str


def get_input_pdb(input_structure: str, gmx_bin: str, output_path: str) -> str:
    """     
    Extract a PDB structure from the input structure file (GRO or TPR) using editconf, and return its path.
    """
    
    # Create a directory for the pre-processing step
    step_dir = os.path.join(output_path, '0_input_pdb')
    os.makedirs(step_dir, exist_ok=True)
    
    # Extract PDB from the input structure (GRO or TPR)
    output_pdb = os.path.join(step_dir, 'input_structure.pdb')
    editconf(
        input_gro_path=input_structure,
        output_gro_path=output_pdb,
        properties={'binary_path': gmx_bin}
    )

    return output_pdb

def build_dry_tpr(full_tpr: str, full_ndx: str, group_name: str, gmx_bin: str, output_path: str) -> str:
    """
    Trim the full-system tpr down to the atoms of `group_name` (e.g. Output_group), producing a
    "dry" tpr whose atom set and numbering match the stripped trajectory. Returns its path.
    """
    step_dir = os.path.join(output_path, 'dry_topology')
    os.makedirs(step_dir, exist_ok=True)
    dry_tpr = os.path.join(step_dir, 'dry.tpr')
    convert_tpr(
        input_tpr_path=full_tpr,
        output_tpr_path=dry_tpr,
        input_ndx_path=full_ndx,
        properties={'binary_path': gmx_bin, 'output_group': group_name}
    )
    return dry_tpr

def build_dry_index(dry_tpr: str, solvent_selection: str, gmx_bin: str, output_path: str,
                    central_index: Optional[int]) -> str:
    """
    Build an index file in the dry (stripped) atom numbering, exposing the group names the
    post-processing steps expect: `Solute_group`, `Output_group` and (optionally) `Center`.

    make_ndx is run on the dry tpr so the numbering matches the Output_group dry trajectory.
    `Output_group` is always the whole dry system ("System"). `Solute_group` (used for fit/cluster)
    is the solute only: when the dry system still contains solvent/ions (e.g. residues retained via
    --keep_residues) it is the complement of the detected Solvent_group; otherwise (default case,
    no solvent left after stripping) it is the whole system. If `central_index` is None (dry
    structure extraction failed) the Center group is omitted; callers fall back to Solute_group.
    """
    step_dir = os.path.join(output_path, 'dry_index')
    os.makedirs(step_dir, exist_ok=True)

    # Default groups on the dry system
    base_ndx = os.path.join(step_dir, 'base.ndx')
    make_ndx(input_structure_path=dry_tpr, output_ndx_path=base_ndx,
             properties={'binary_path': gmx_bin})
    default_groups = read_groups(base_ndx)

    # Try to detect leftover solvent/ions (present only when residues were retained)
    solvent_ndx = os.path.join(step_dir, 'solvent.ndx')
    make_ndx(input_structure_path=dry_tpr, input_ndx_path=base_ndx, output_ndx_path=solvent_ndx,
             properties={'binary_path': gmx_bin, 'selection': solvent_selection})
    if len(read_groups(solvent_ndx)) > len(default_groups):
        rename_last_ndx_group(solvent_ndx, solvent_group)
        solute_selection = f'! "{solvent_group}"'
    else:
        solute_selection = '"System"'

    # Solute group (fit/cluster selection): solute only
    solute_ndx = os.path.join(step_dir, 'solute.ndx')
    make_ndx(input_structure_path=dry_tpr, input_ndx_path=solvent_ndx, output_ndx_path=solute_ndx,
             properties={'binary_path': gmx_bin, 'selection': solute_selection})
    rename_last_ndx_group(solute_ndx, solute_group)

    # Output group: the whole dry system (the dry trajectory is the full Output_group set)
    output_ndx = os.path.join(step_dir, 'output.ndx')
    make_ndx(input_structure_path=dry_tpr, input_ndx_path=solute_ndx, output_ndx_path=output_ndx,
             properties={'binary_path': gmx_bin, 'selection': '"System"'})
    rename_last_ndx_group(output_ndx, output_group)

    if central_index is not None:
        final_ndx = os.path.join(step_dir, 'index.ndx')
        add_group([central_index], 'Center', output_ndx, final_ndx)
        return final_ndx
    return output_ndx

def index_config_contents(
    gmx_bin: str = 'gmx',
    solvent_selection: str = '"SOL" | "Ion"',
    output_selection: str = f'"{solute_group}"',
    structure_path: str = ''
) -> str:
    """
    Steps building the index file (Solvent_group, Solute_group, Output_group) from the input
    structure. Only included in the config when no usable index file is provided by the caller.
    """
    return f"""
###################################################
# Section 1: Index group creation from structure  #
###################################################

# Create base index file
step0_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}

# Add solvent and ions to index file
step1_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    input_ndx_path: dependency/step0_make_ndx/output_ndx_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '{solvent_selection}'

# Add solute group (complement of solvent) to index file
step2_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    input_ndx_path: dependency/step1_make_ndx/output_ndx_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '! "{solvent_group}"'

# Add output group to index file
step3_make_ndx:
  tool: make_ndx
  paths:
    input_structure_path: {structure_path}
    input_ndx_path: dependency/step2_make_ndx/output_ndx_path
    output_ndx_path: index.ndx
  properties:
    binary_path: {gmx_bin}
    selection: '{output_selection}'
"""

def common_config_contents(
    gmx_bin: str = 'gmx',
    debug: bool = False,
    solvent_selection: str = '"SOL" | "Ion"',
    output_selection: str = f'"{solute_group}"',
    structure_path: str = '',
    input_topology_path: str = '',
    input_traj_path: str = '',
    restart: bool = False,
    index_path: str = 'dependency/step3_make_ndx/output_ndx_path'

) -> str:
    global_properties = f"""
# Global properties (common for all steps)
global_properties:
  can_write_console_log: false
  restart: {to_yaml(restart)}
  remove_tmp: {to_yaml(not debug)}
"""

    # Build the index in this workflow unless the caller already provided one
    if index_path == 'dependency/step3_make_ndx/output_ndx_path':
        index_steps = index_config_contents(gmx_bin, solvent_selection, output_selection, structure_path)
    else:
        index_steps = f"""
###################################################
# Section 1: Index groups provided by the caller  #
###################################################
# Index file: {index_path}
"""

    return global_properties + index_steps + f"""
###################################################
# Section 2: Structure and trajectory processing  #
###################################################

# Extract the output structure (solute + any kept residues), matching the dry trajectory
step4_dry_str:
  tool: gmx_trjconv_str
  paths:
    input_structure_path: {structure_path}
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_str_path: dry_structure.pdb
  properties:
    binary_path: {gmx_bin}
    selection: "{output_group}"
    center: false
    pbc: none

# Step 5 (Center group creation) is done in Python after extracting the
# dry structure in step 4. The Center group is needed for centering
# the trajectory in step 7.

# Extract dry (or full) trajectory
step6_dry_traj:
  tool: gmx_trjconv_trj
  paths:
    input_traj_path: {input_traj_path}
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: dry_traj.xtc
  properties:
    binary_path: {gmx_bin}
    selection: "{output_group}"
"""

def complete_postprocessing(
    gmx_bin: str = 'gmx',
    input_topology_path: str = '',
    index_path: str = 'dependency/step3_make_ndx/output_ndx_path'
) -> str:
    return f"""
# Make the molecules whole 
step7_whole:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step6_dry_traj/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: whole_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: false
    pbc: whole

# Cluster the molecules
step8_cluster:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step7_whole/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: cluster_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    cluster_selection: "{solute_group}"
    center: false
    pbc: cluster

# Extract the first frame to use as ref
step9_extract_ref:
  tool: gmx_trjconv_trj
  paths:
    input_traj_path: dependency/step8_cluster/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: first_frame.gro
  properties:
    binary_path: {gmx_bin}
    selection: "{output_group}"
    dump: 0

# Use as ref with pbc nojump 
step10_nojump:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step8_cluster/output_traj_path
    input_top_path: dependency/step9_extract_ref/output_traj_path
    input_index_path: {index_path}
    output_traj_path: nojump_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: false
    pbc: nojump

# Center the system - whole solute group / central atom
step11_center:
  tool: gmx_image
  paths: 
    input_traj_path: dependency/step10_nojump/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: centered_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center_selection: "Center"
    center: true
    ur: compact
    pbc: none

# Image the trajectory to put all molecules back in the box
step12_image:
  tool: gmx_image
  paths: 
    input_traj_path: dependency/step11_center/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: centered_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: false
    ur: compact
    pbc: mol

# Fit the trajectory by rotation and translation
step13_fit:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step12_image/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: fitted_traj.xtc
  properties:
    binary_path: {gmx_bin}
    fit_selection: "{solute_group}"
    output_selection: "{output_group}"
    center: false
    fit: rot+trans
"""

def fast_postprocessing(
    gmx_bin: str = 'gmx',
    input_topology_path: str = '',
    index_path: str = 'dependency/step3_make_ndx/output_ndx_path'
) -> str:
    return f"""
# Center the trajectory on central atom
step7_center:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step6_dry_traj/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: center_traj.xtc
  properties:
    binary_path: {gmx_bin}
    center_selection: "Center"
    output_selection: "{output_group}"
    center: true
    ur: compact
    pbc: none

# Apply periodic boundary conditions imaging
step8_image:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step7_center/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: imaged_traj.xtc
  properties:
    binary_path: {gmx_bin}
    output_selection: "{output_group}"
    center: false
    ur: compact
    pbc: mol

# Fit the trajectory by rotation and translation
step9_fit:
  tool: gmx_image
  paths:
    input_traj_path: dependency/step8_image/output_traj_path
    input_top_path: {input_topology_path}
    input_index_path: {index_path}
    output_traj_path: fitted_traj.xtc
  properties:
    binary_path: {gmx_bin}
    fit_selection: "{solute_group}"
    output_selection: "{output_group}"
    center: false
    fit: rot+trans
"""

def config_contents(
    gmx_bin: str = 'gmx',
    debug: bool = False,
    solvent_selection: str = '"SOL" | "Ion"',
    output_selection: str = f'"{solute_group}"',
    structure_path: str = '',
    input_topology_path: str = '',
    input_traj_path: str = '',
    fast: bool = False,
    restart: bool = False,
    index_path: str = 'dependency/step3_make_ndx/output_ndx_path'
) -> str:
    common_contents = common_config_contents(gmx_bin,
                                            debug,
                                            solvent_selection,
                                            output_selection,
                                            structure_path,
                                            input_topology_path,
                                            input_traj_path,
                                            restart,
                                            index_path)

    if fast:
        return common_contents + fast_postprocessing(gmx_bin, input_topology_path, index_path)
    else:
        return common_contents + complete_postprocessing(gmx_bin, input_topology_path, index_path)

def create_config_file(output_path: str, **config_args) -> str:
    """Write YAML config to output_path/config.yml and return its path."""
    config_path = os.path.join(output_path, 'config.yml')
    with open(config_path, 'w') as f:
        f.write(config_contents(**config_args))
    return config_path


def traj_postprocessing(
    input_traj_path: str,
    input_topology_path: str,
    input_structure_path: str,
    input_index_path: Optional[str] = None,
    gmx_bin: str = 'gmx',
    keep_solvent: bool = False,
    residues_to_keep: Optional[List[int]] = None,
    extra_ions: List[str] = [],
    extra_solvents: List[str] = [],
    solvent_selection: Optional[str] = None,
    fast: bool = False,
    debug: bool = False,
    restart: bool = False,
    output_path: str = 'output',
    output_traj_path='trajectory.xtc',
    output_str_path='structure.pdb'
):
    '''
    Post-process a GROMACS MD trajectory: strip solvent, center, image, and fit.

    Inputs
    ------
        input_traj_path:
            path to the input trajectory file (.xtc). Required.
        input_topology_path:
            path to the binary run input file (.tpr). Required.
        input_structure_path:
            path to an input structure file (.gro or .pdb).
            Used to define solvent/output index groups and to find the center group for centering.
            Make sure the structure is not broken due to PBC.
        input_index_path:
            path to an index file (.ndx) already exposing the Solute_group and Output_group groups
            (e.g. the one built by md_gromacs). If given and usable, the index-building steps
            (step0-step3) are skipped and this file is used instead.
            Default: None (build the index from the input structure)
        gmx_bin:
            GROMACS binary path. Default: gmx
        keep_solvent:
            include solvent and ions in the output structure and trajectory.
            Default: False (dry output)
        residues_to_keep:
            residue indices to retain in the output besides the solute.
            Default: None (only solute)
        extra_ions:
            additional ion atom names to include in the solvent group (e.g. --ions NA+ CA2+). Default: []
        extra_solvents:
            additional solvent residue names to include in the solvent group (e.g. --solvents TIP3 TIP4). Default: []
        solvent_selection:
            GROMACS selection string for the solvent and ions. If given, the solvent/ion detection
            on the input structure is skipped (extra_ions/extra_solvents are then ignored).
            Default: None (detect solvent and ions in the input structure)
        debug:
            keep intermediate files. Default: False
        output_path:
            directory for workflow output. Default: output

    Outputs
    -------
        output/
            dry_structure.pdb   — processed structure
            fitted_traj.xtc     — processed trajectory (dry, centered, imaged, fitted)
        global_paths (dict), global_prop (dict)
    '''

    start_time = time.time()

    # Determine final output path
    output_path = fu.get_working_dir_path(output_path, restart=restart)

    # Initialize a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)
    global_log.info(f"biobb_md_workflows version {__version__}")

    ###########################################
    # Build GROMACS selection strings for ndx #
    ###########################################

    # Check whether the provided index file (if any) can be reused
    provided_groups = {}
    build_index = True
    if input_index_path and os.path.exists(input_index_path):
        provided_groups = read_group_indices(input_index_path)
        missing_groups = [group for group in (solute_group, output_group) if group not in provided_groups]
        if missing_groups:
            global_log.warning(f"Provided index file lacks the {missing_groups} group(s), building a new index file")
        else:
            global_log.info(f"Reusing the provided index file: {input_index_path}")
            build_index = False
    elif input_index_path:
        global_log.warning(f"Provided index file {input_index_path} not found, building a new index file")

    # Construct solvent selection, solute selection will be the the rest of the system
    if solvent_selection is None:

        # Convert provided GRO to PDB for solvent/ion detection if needed
        if Path(input_structure_path).suffix.lower() != '.pdb':
            pdb_structure_path = get_input_pdb(input_structure_path, gmx_bin, output_path)
        else:
            pdb_structure_path = input_structure_path

        solvent_names = get_residue_types(pdb_structure_path, solvent_library + list(extra_solvents))
        ion_names = get_atom_types(pdb_structure_path, ions_library + list(extra_ions))
        solvent_selection = build_solvent_selection(solvent_names, ion_names)

    # Determine output selection based on user options
    if keep_solvent:
        output_selection = '"System"'
    elif residues_to_keep:
        residues_sel = f"ri {' '.join(str(r) for r in residues_to_keep)}"
        output_selection = f'"{solute_group}" | {residues_sel}'
    else:
        output_selection = f'"{solute_group}"'

    #################################
    # Create workflow configuration #
    #################################

    config_path = create_config_file(
        output_path,
        gmx_bin=gmx_bin,
        debug=debug,
        solvent_selection=solvent_selection,
        output_selection=output_selection,
        structure_path=os.path.abspath(input_structure_path),
        input_topology_path=os.path.abspath(input_topology_path),
        input_traj_path=os.path.abspath(input_traj_path),
        fast=fast,
        restart=restart,
        index_path=('dependency/step3_make_ndx/output_ndx_path' if build_index
                    else os.path.abspath(input_index_path))
    )
    global_log.info(f"Configuration file: {config_path}")

    conf = settings.ConfReader(config=config_path, system=None)
    conf.working_dir_path = output_path
    prop = conf.get_prop_dic()
    paths = conf.get_paths_dic()

    ######################
    # Create index files #
    ######################

    if build_index:

        # Find default groups in the input structure
        global_log.info("step0_make_ndx: Create base index file")
        make_ndx(**paths['step0_make_ndx'], properties=prop['step0_make_ndx'])
        default_groups = read_groups(paths['step0_make_ndx']['output_ndx_path'])

        # Try to add a solvent group
        global_log.info("step1_make_ndx: Create solvent group in index file")
        make_ndx(**paths['step1_make_ndx'], properties=prop['step1_make_ndx'])
        all_groups = read_groups(paths['step1_make_ndx']['output_ndx_path'])

        # Check if a Solvent group was added correctly
        solvent_created = len(all_groups) > len(default_groups)
        if solvent_created:
            # If there is a solvent group, rename it and leave solute/output selections as they are
            global_log.info(f"Renaming last created group to {solvent_group}")
            rename_last_ndx_group(paths['step1_make_ndx']['output_ndx_path'], solvent_group)
        else:
            # If there is no solvent group, change the solute/output selections to System
            global_log.info(f"No Solvent group was created.")
            global_log.info(f"If your input structure contains solvent molecules or ions not recognized by GROMACS, ")
            global_log.info(f"make sure you include them with the --ions or --solvent arguments")
            prop['step2_make_ndx']['selection'] = '"System"'
            prop['step3_make_ndx']['selection'] = '"System"'

        global_log.info("step2_make_ndx: Create solute group in index file")
        make_ndx(**paths['step2_make_ndx'], properties=prop['step2_make_ndx'])
        rename_last_ndx_group(paths['step2_make_ndx']['output_ndx_path'], solute_group)

        global_log.info("step3_make_ndx: Create output group in index file")
        make_ndx(**paths['step3_make_ndx'], properties=prop['step3_make_ndx'])
        rename_last_ndx_group(paths['step3_make_ndx']['output_ndx_path'], output_group)

        input_ndx_path = paths['step3_make_ndx']['output_ndx_path']

        # The trajectory is only reduced (renumbered) when a solvent group was stripped. When nothing
        # is stripped (--keep_solvent, or no solvent group detected -> Output_group == System) the full
        # tpr + full index already match the trajectory, so the dry-tpr rewiring below is skipped.
        needs_trim = solvent_created and not keep_solvent

    else:
        # Reuse the index provided by the caller. The trajectory is reduced (renumbered) whenever the
        # Output_group is a strict subset of the system, which is read directly from the index file.
        input_ndx_path = os.path.abspath(input_index_path)
        if 'System' in provided_groups:
            needs_trim = len(provided_groups[output_group]) < len(provided_groups['System'])
        else:
            global_log.warning("Provided index file has no System group, assuming the trajectory is stripped")
            needs_trim = not keep_solvent

    ########################################
    # Extract solute and find central atom #
    ########################################

    centering_step = 'step7_center' if fast else 'step11_center'

    # Extract the solute from the input structure
    try:
        global_log.info("step4_dry_str: Center dry structure")
        gmx_trjconv_str(**paths['step4_dry_str'], properties=prop['step4_dry_str'])
        dry_str_ok = True
    except SystemExit as e:
        global_log.error(f"Structure post-processing failed (SystemExit, code={e.code})")
        dry_str_ok = False
    except Exception:
        global_log.exception("Structure post-processing failed with unexpected exception")
        dry_str_ok = False

    # Extract central atom index from the dry structure (dry numbering) to center the trajectory
    if dry_str_ok:
        central_index = get_central_atom_index(paths['step4_dry_str']['output_str_path'])
        global_log.info(f"Central atom index for centering: {central_index}")
    else:
        central_index = None
        global_log.warning("Structure post-processing failed falling back to Solute_group for centering")

    if needs_trim:
        # Build a dry tpr (Output_group subset) and a dry index so the post-strip steps operate in the stripped trajectory
        global_log.info("Building dry topology (convert-tpr) and dry index for post-processing")
        dry_tpr = build_dry_tpr(os.path.abspath(input_topology_path), input_ndx_path,
                                output_group, gmx_bin, output_path)
        dry_ndx = build_dry_index(dry_tpr, solvent_selection, gmx_bin, output_path, central_index)

        # Rewire post-strip steps to the dry tpr (-s) and dry index (-n)
        if fast:
            post_steps = ['step7_center', 'step8_image', 'step9_fit']
        else:
            post_steps = ['step7_whole', 'step8_cluster', 'step9_extract_ref', 'step10_nojump',
                          'step11_center', 'step12_image', 'step13_fit']
        # step10_nojump already uses a dry .gro reference as -s step9_extract_ref/output_traj_path
        for step in post_steps:
            paths[step]['input_index_path'] = dry_ndx
            if step != 'step10_nojump':
                paths[step]['input_top_path'] = dry_tpr

        if central_index is None:
            # Center group could not be built; fall back to centering on the whole solute group
            prop[centering_step]['center_selection'] = solute_group
    else:
        # No stripping: keep the original behavior (full tpr + full index, Center on the full index)
        if dry_str_ok:
            os.makedirs(os.path.join(output_path, 'step5_center_group'), exist_ok=True)
            center_ndx_path = os.path.join(output_path, 'step5_center_group', 'center.ndx')
            add_group([central_index], 'Center', input_ndx_path, center_ndx_path)
            paths[centering_step]['input_index_path'] = center_ndx_path
        else:
            prop[centering_step]['center_selection'] = solute_group

    #########################################################################
    # Process trajectory: strip, whole, nojump, cluster, center, image, fit #
    #########################################################################
        
    final_step = 'step9_fit' if fast else 'step13_fit'

    try:
        global_log.info("step6_dry_traj: Extract dry trajectory")
        gmx_trjconv_trj(**paths['step6_dry_traj'], properties=prop['step6_dry_traj'])
        
        if fast:
            global_log.info("Running in fast mode: skipping whole, nojump, and cluster steps")
            
            global_log.info("step7_center: Center the trajectory using the Center group")
            gmx_image(**paths['step7_center'], properties=prop['step7_center'])
            if not debug:
                os.remove(paths['step6_dry_traj']['output_traj_path'])
            
            global_log.info("step8_image: Image the trajectory to put all molecules back in the box")
            gmx_image(**paths['step8_image'], properties=prop['step8_image'])
            if not debug:
                os.remove(paths['step7_center']['output_traj_path'])
                
            global_log.info("step9_fit: Fit the trajectory by rotation and translation")
            gmx_image(**paths['step9_fit'], properties=prop['step9_fit'])
        else:
            global_log.info("Running in complete mode: performing whole, nojump, and cluster steps for better results")
            
            global_log.info("step7_whole: Make the molecules whole in the trajectory")
            gmx_image(**paths['step7_whole'], properties=prop['step7_whole'])
            if not debug:
                os.remove(paths['step6_dry_traj']['output_traj_path'])

            global_log.info("step8_cluster: Cluster the molecules in the trajectory")
            gmx_image(**paths['step8_cluster'], properties=prop['step8_cluster'])
            if not debug:
                os.remove(paths['step7_whole']['output_traj_path'])

            global_log.info("step9_extract_ref: Extract the first frame to use as reference")
            gmx_trjconv_trj(**paths['step9_extract_ref'], properties=prop['step9_extract_ref'])

            global_log.info("step10_nojump: Make the trajectory whole with nojump PBC")
            gmx_image(**paths['step10_nojump'], properties=prop['step10_nojump'])
            if not debug:
                os.remove(paths['step8_cluster']['output_traj_path'])

            global_log.info("step11_center: Center the trajectory using the Center group")
            gmx_image(**paths['step11_center'], properties=prop['step11_center'])
            if not debug: 
                os.remove(paths['step10_nojump']['output_traj_path'])
            
            global_log.info("step12_image: Image the trajectory to put all molecules back in the box")
            gmx_image(**paths['step12_image'], properties=prop['step12_image'])
            if not debug:
                os.remove(paths['step11_center']['output_traj_path'])
            
            global_log.info("step13_fit: Fit the trajectory by rotation and translation")
            gmx_image(**paths['step13_fit'], properties=prop['step13_fit'])

    except SystemExit as e:
        global_log.error(f"Trajectory post-processing failed (SystemExit, code={e.code})")
    except Exception:
        global_log.exception("Trajectory post-processing failed with unexpected exception")
        
    # Move final outputs to user-specified paths. Both steps are inside try/except blocks above, so
    # the files may be missing when post-processing failed - do not raise here, the manifest and the
    # log below still have to be written (this function is also called from md_gromacs).
    final_traj_path = paths[final_step]['output_traj_path']
    final_str_path = paths['step4_dry_str']['output_str_path']
    for source_path, destination_path in ((final_traj_path, output_traj_path),
                                          (final_str_path, output_str_path)):
        if os.path.exists(source_path):
            os.rename(source_path, destination_path)
        else:
            global_log.error(f"Expected output {source_path} not found, skipping move to {destination_path}")

    # Write a stable output manifest for external consumers (see manifest.yaml in output_path)
    manifest_outputs = {}
    if os.path.exists(output_str_path):
        manifest_outputs["structure"] = {"pdb": os.path.relpath(output_str_path, output_path)}
    if os.path.exists(output_traj_path):
        manifest_outputs["trajectory"] = {"xtc": os.path.relpath(output_traj_path, output_path)}
    with open(os.path.join(output_path, "manifest.yaml"), "w") as manifest_file:
        yaml.safe_dump({"schema_version": 1, "outputs": manifest_outputs}, manifest_file, sort_keys=False)

    elapsed = time.time() - start_time
    global_log.info('')
    global_log.info('Execution successful:')
    global_log.info(f'  Workflow path: {output_path}')
    global_log.info(f'  Config file:   {config_path}')
    global_log.info(f'  Elapsed time:  {elapsed/60:.1f} minutes')

    return paths, prop


def main():

    parser = argparse.ArgumentParser(
        description="Post-process a GROMACS MD trajectory: strip solvent, center, image, fit."
    )

    ###############
    # Input files #
    ###############

    parser.add_argument(
        '--input_traj', dest='input_traj_path', type=str, required=True,
        help="Input trajectory file (.xtc). Required."
    )
    parser.add_argument(
        '--input_top', dest='input_topology_path', type=str, required=True,
        help="Input binary run input file (.tpr). Required."
    )
    parser.add_argument(
        '--input_structure', dest='input_structure_path', type=str, required=True,
        help=("Input structure file (.gro or .pdb). Used to define solvent/output "
              "index groups and to find the center group for centering. "
              "Make sure the structure is not broken due to PBC.")
    )
    parser.add_argument(
        '--input_index', dest='input_index_path', type=str, required=False,
        help=(f"Input index file (.ndx) already containing the {solute_group} and {output_group} "
              "groups. If given and usable, the index creation steps are skipped. "
              "Default: None (build the index from the input structure)")
    )

    #########################
    # Configuration options #
    #########################

    parser.add_argument(
        '--gmx_bin', type=str, required=False, default='gmx',
        help="GROMACS binary path. Default: gmx"
    )
    parser.add_argument(
        '--keep_solvent', action='store_true', required=False, default=False,
        help="Keep solvent and ions in the output structure and trajectory. Default: False"
    )
    parser.add_argument(
        '--keep_residues', dest='residues_to_keep', type=int, nargs='+',
        required=False,
        help=("Residue indices to keep in the output besides the solute "
              "(e.g. --keep_residues 15 23 105). Default: None")
    )
    parser.add_argument(
        '--ions', dest='extra_ions', type=str, nargs='+',
        required=False, default=[],
        help=("Additional ion atom names to include in the solvent group (e.g. --ions NA+ CA2+). Default: []" )
    )
    parser.add_argument(
        '--solvents', dest='extra_solvents', type=str, nargs='+',
        required=False, default=[],
        help=("Additional solvent residue names to include in the solvent group (e.g. --solvents TIP3 TIP4). Default: []" )
    )
    parser.add_argument(
        '--fast', action='store_true', required=False, default=False,
        help="Skip making solute whole, removing jumps and clustering. Default: False"
    )
    parser.add_argument(
        '--debug', action='store_true', required=False, default=False,
        help="Keep intermediate files. Default: False"
    )
    parser.add_argument(
        '--restart', action='store_true', required=False, default=False,
        help="Restart the workflow from the last completed step. Default: False"
    )
    parser.add_argument(
        '--output', dest='output_path', type=str, required=False, default='output',
        help="Output directory path for the workflow where the steps will be written. Default: output"
    )
    parser.add_argument(
        '--output_traj', dest='output_traj_path', type=str, required=False, default='trajectory.xtc',
        help="Output trajectory file name (e.g. processed_traj.xtc). Default: trajectory.xtc"
    )
    parser.add_argument(
        '--output_str', dest='output_str_path', type=str, required=False, default='structure.pdb',
        help="Output structure file name (e.g. processed_structure.pdb). Default: structure.pdb"
    )

    args = parser.parse_args()

    traj_postprocessing(
        input_traj_path=args.input_traj_path,
        input_topology_path=args.input_topology_path,
        input_structure_path=args.input_structure_path,
        input_index_path=args.input_index_path,
        gmx_bin=args.gmx_bin,
        keep_solvent=args.keep_solvent,
        residues_to_keep=args.residues_to_keep,
        extra_ions=args.extra_ions,
        extra_solvents=args.extra_solvents,
        fast = args.fast,
        debug=args.debug,
        restart=args.restart,
        output_path=args.output_path,
        output_traj_path=args.output_traj_path,
        output_str_path=args.output_str_path
    )


if __name__ == '__main__':
    main()
