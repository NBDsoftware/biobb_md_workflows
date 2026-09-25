#!/usr/bin/env python3

# Importing all the needed libraries
from typing import Dict, Literal, Optional, Tuple
import argparse
import logging
import os
import time
import yaml

from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu

from biobb_md_workflows.common import to_yaml
from biobb_md_workflows import __version__

from biobb_pdb_tools.pdb_tools.biobb_pdb_fetch import biobb_pdb_fetch
from biobb_structure_utils.utils.extract_molecule import extract_molecule
from biobb_godmd.godmd.godmd_prep import godmd_prep
from biobb_godmd.godmd.godmd_run import godmd_run
from biobb_analysis.ambertools.cpptraj_convert import cpptraj_convert


def check_inputs(
    global_log: logging.Logger,
    origin_pdb_path: Optional[str],
    origin_pdb_code: Optional[str],
    target_pdb_path: Optional[str],
    target_pdb_code: Optional[str]
) -> Tuple[Literal['input_pdb', 'pdb_code'], Literal['input_pdb', 'pdb_code']]:
    """
    Check the inputs of the workflow.

    Origin and target are resolved independently - one can be a local file while the
    other is fetched by PDB code.

    Inputs
    ------
        global_log: Logger object for logging messages.
        origin_pdb_path (str): Path to the origin PDB file.
        origin_pdb_code (str): PDB code of the origin structure.
        target_pdb_path (str): Path to the target PDB file.
        target_pdb_code (str): PDB code of the target structure.

    Returns
    -------
        (origin_mode, target_mode): each 'input_pdb' or 'pdb_code'.
    """

    def check_one(label: str, pdb_path: Optional[str], pdb_code: Optional[str]) -> Literal['input_pdb', 'pdb_code']:
        if pdb_path:
            if not os.path.exists(pdb_path):
                global_log.error(f"{label} PDB file {pdb_path} not found.")
                raise FileNotFoundError(f"{label} PDB file {pdb_path} not found.")
            return 'input_pdb'
        elif pdb_code:
            return 'pdb_code'
        else:
            global_log.error(f"Either --{label.lower()}_pdb or --{label.lower()}_pdb_code must be given.")
            raise ValueError(f"Either --{label.lower()}_pdb or --{label.lower()}_pdb_code must be given.")

    origin_mode = check_one('origin', origin_pdb_path, origin_pdb_code)
    target_mode = check_one('target', target_pdb_path, target_pdb_code)

    return origin_mode, target_mode


# YML construction
def config_contents(
    origin_mode: Literal['input_pdb', 'pdb_code'],
    target_mode: Literal['input_pdb', 'pdb_code'],
    origin_pdb_path: Optional[str],
    origin_pdb_code: Optional[str],
    origin_chain: Optional[str],
    target_pdb_path: Optional[str],
    target_pdb_code: Optional[str],
    target_chain: Optional[str],
    temperature: float = 300,
    random_seed: int = 2525,
    num_blocks: int = 10000,
    snapshot_interval: int = 500,
    rmsd_tolerance: float = 1.5,
    gapopen: float = 12.0,
    gapextend: float = 2.0,
    restart: bool = False,
    debug: bool = False,
    ) -> str:
    """
    Returns the contents of the YAML configuration file as a string.

    The YAML file contains the configuration for the GOdMD workflow.

    Returns
    -------
    str
        The contents of the YAML configuration file.
    """

    # Resolve the input path each extract_molecule step should read from: a dependency on the
    # fetch step when the corresponding structure is being downloaded, or the (absolute) user path.
    if origin_mode == 'pdb_code':
        origin_extract_input = 'dependency/step1_origin_fetch/output_file_path'
    else:
        origin_extract_input = os.path.abspath(origin_pdb_path)

    if target_mode == 'pdb_code':
        target_extract_input = 'dependency/step3_target_fetch/output_file_path'
    else:
        target_extract_input = os.path.abspath(target_pdb_path)

    # Chain selection: restrict to the given chain when provided, otherwise keep all
    # (extract_molecule still strips waters/ligands in "all" mode).
    origin_molecule_type = "all"
    origin_chains_property = ""
    if origin_chain is not None:
        origin_molecule_type = "chains"
        origin_chains_property = f"chains : {to_yaml([origin_chain])}"

    target_molecule_type = "all"
    target_chains_property = ""
    if target_chain is not None:
        target_molecule_type = "chains"
        target_chains_property = f"chains : {to_yaml([target_chain])}"

    return f"""
# Global properties (common for all steps)
global_properties:
  can_write_console_log: false                                      # Verbose writing of log information
  restart: {to_yaml(restart)}                                       # Skip steps already performed
  remove_tmp: {to_yaml(not debug)}                                  # Remove temporal files

# Step 1 (Optional): fetch origin PDB structure from RCSB website
step1_origin_fetch:
  tool: biobb_pdb_fetch
  paths:
    output_file_path: origin_raw.pdb
  properties:
    biounit: false
    pdbid: {origin_pdb_code}

# Step 2: extract chain(s) of interest from the origin structure, stripping waters/ligands
step2_origin_extract:
  tool: extract_molecule
  paths:
    input_structure_path: {origin_extract_input}
    output_molecule_path: origin.pdb
  properties:
    molecule_type: {origin_molecule_type}
    {origin_chains_property}

# Step 3 (Optional): fetch target PDB structure from RCSB website
step3_target_fetch:
  tool: biobb_pdb_fetch
  paths:
    output_file_path: target_raw.pdb
  properties:
    biounit: false
    pdbid: {target_pdb_code}

# Step 4: extract chain(s) of interest from the target structure, stripping waters/ligands
step4_target_extract:
  tool: extract_molecule
  paths:
    input_structure_path: {target_extract_input}
    output_molecule_path: target.pdb
  properties:
    molecule_type: {target_molecule_type}
    {target_chains_property}

# Step 5: build a one-to-one residue mapping (origin <-> target) from a pairwise sequence
# alignment (EMBOSS water), required by GOdMD to move atoms from origin to target positions
step5_godmd_prep:
  tool: godmd_prep
  paths:
    input_pdb_orig_path: dependency/step2_origin_extract/output_molecule_path
    input_pdb_target_path: dependency/step4_target_extract/output_molecule_path
    output_aln_orig_path: origin.aln
    output_aln_target_path: target.aln
  properties:
    gapopen: {gapopen}
    gapextend: {gapextend}

# Step 6: run GOdMD, computing the conformational transition trajectory from origin to target
step6_godmd_run:
  tool: godmd_run
  paths:
    input_pdb_orig_path: dependency/step2_origin_extract/output_molecule_path
    input_pdb_target_path: dependency/step4_target_extract/output_molecule_path
    input_aln_orig_path: dependency/step5_godmd_prep/output_aln_orig_path
    input_aln_target_path: dependency/step5_godmd_prep/output_aln_target_path
    output_log_path: godmd.log
    output_ene_path: godmd_energy.out
    output_trj_path: transition.mdcrd
    output_pdb_path: transition.pdb
  properties:
    godmdin:
      TEMP: {temperature}
      SEED: {random_seed}
      NBLOC: {num_blocks}
      TSNAP: {snapshot_interval}
      ERRORACCEPTABLE: {rmsd_tolerance}

# Step 7: convert the GOdMD trajectory (mdcrd) to xtc for downstream analysis/visualization
step7_convert:
  tool: cpptraj_convert
  paths:
    input_top_path: dependency/step6_godmd_run/output_pdb_path
    input_traj_path: dependency/step6_godmd_run/output_trj_path
    output_cpptraj_path: transition.xtc
  properties:
    format: xtc
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
def godmd_simulation(
    origin_pdb_path: Optional[str] = None,
    origin_pdb_code: Optional[str] = None,
    origin_chain: Optional[str] = None,
    target_pdb_path: Optional[str] = None,
    target_pdb_code: Optional[str] = None,
    target_chain: Optional[str] = None,
    temperature: float = 300,
    random_seed: int = 2525,
    num_blocks: int = 10000,
    snapshot_interval: int = 500,
    rmsd_tolerance: float = 1.5,
    gapopen: float = 12.0,
    gapextend: float = 2.0,
    debug: bool = False,
    restart: bool = False,
    output_path: str = 'output'
    ) -> Tuple[Dict, Dict]:
    '''
    Compute an approximate conformational transition between an origin and a target
    structure of the same protein using GOdMD (discrete MD with GO-like biasing potentials).

    Inputs
    ------

        origin_pdb_path:
            path to a local PDB file for the origin structure. Mutually exclusive with origin_pdb_code.
        origin_pdb_code:
            PDB code to fetch the origin structure from RCSB. Mutually exclusive with origin_pdb_path.
        origin_chain:
            chain of interest in the origin structure. If None, all chains are kept.
        target_pdb_path:
            path to a local PDB file for the target structure. Mutually exclusive with target_pdb_code.
        target_pdb_code:
            PDB code to fetch the target structure from RCSB. Mutually exclusive with target_pdb_path.
        target_chain:
            chain of interest in the target structure. If None, all chains are kept.
        temperature:
            GOdMD simulation temperature (K).
        random_seed:
            random number generator seed.
        num_blocks:
            number of simulation blocks (total simulated time = num_blocks * snapshot_interval).
        snapshot_interval:
            output structure interval, in fs.
        rmsd_tolerance:
            minimum acceptable RMSD (Angstrom) to consider the transition converged.
        gapopen:
            EMBOSS water pairwise alignment gap-open penalty.
        gapextend:
            EMBOSS water pairwise alignment gap-extend penalty.
        debug:
            whether to run the workflow in debug mode (keep temporary files).
        restart:
            whether to restart the workflow from the last completed step.
        output_path:
            path to output folder

    Outputs
    -------

        /output folder
            transition.xtc, transition.pdb, godmd.log, godmd_energy.out, manifest.yaml
        global_paths    (dict): dictionary with all workflow paths
        global_prop     (dict): dictionary with all workflow properties

    '''

    ###########################
    # Workflow initialization #
    ###########################

    start_time = time.time()

    # Determine final output path
    output_path = fu.get_working_dir_path(output_path, restart=restart)

    # Initialize a global log file
    global_log, _ = fu.get_logs(path=output_path, light_format=True)
    global_log.info(f"biobb_md_workflows version {__version__}")

    # Check inputs
    origin_mode, target_mode = check_inputs(
        global_log,
        origin_pdb_path,
        origin_pdb_code,
        target_pdb_path,
        target_pdb_code
    )

    # Create the configuration file
    config_args = {
        'origin_mode': origin_mode,
        'target_mode': target_mode,
        'origin_pdb_path': origin_pdb_path,
        'origin_pdb_code': origin_pdb_code,
        'origin_chain': origin_chain,
        'target_pdb_path': target_pdb_path,
        'target_pdb_code': target_pdb_code,
        'target_chain': target_chain,
        'temperature': temperature,
        'random_seed': random_seed,
        'num_blocks': num_blocks,
        'snapshot_interval': snapshot_interval,
        'rmsd_tolerance': rmsd_tolerance,
        'gapopen': gapopen,
        'gapextend': gapextend,
        'restart': restart,
        'debug': debug
    }
    configuration_path = create_config_file(output_path, **config_args)
    conf = settings.ConfReader(config=configuration_path)
    conf.working_dir_path = output_path

    # Parsing the input configuration file (YAML);
    # Dividing it in global paths and global properties
    global_prop = conf.get_prop_dic(global_log=global_log)
    global_paths = conf.get_paths_dic()

    ##################
    # Fetch/prepare  #
    ##################

    if origin_mode == 'pdb_code':
        global_log.info("step1_origin_fetch: Fetch origin PDB structure")
        biobb_pdb_fetch(**global_paths["step1_origin_fetch"], properties=global_prop["step1_origin_fetch"])

    global_log.info("step2_origin_extract: Extract origin chain(s) of interest")
    extract_molecule(**global_paths["step2_origin_extract"], properties=global_prop["step2_origin_extract"])

    if target_mode == 'pdb_code':
        global_log.info("step3_target_fetch: Fetch target PDB structure")
        biobb_pdb_fetch(**global_paths["step3_target_fetch"], properties=global_prop["step3_target_fetch"])

    global_log.info("step4_target_extract: Extract target chain(s) of interest")
    extract_molecule(**global_paths["step4_target_extract"], properties=global_prop["step4_target_extract"])

    global_log.info("step5_godmd_prep: Compute origin-target residue mapping")
    godmd_prep(**global_paths["step5_godmd_prep"], properties=global_prop["step5_godmd_prep"])

    ##############
    # GOdMD run  #
    ##############

    try:
        global_log.info("step6_godmd_run: Compute the conformational transition")
        godmd_run(**global_paths["step6_godmd_run"], properties=global_prop["step6_godmd_run"])
    except SystemExit as e:
        global_log.error(f"GOdMD run failed (SystemExit, code={e.code})")
    except Exception:
        global_log.exception("GOdMD run failed with unexpected exception")

    ##########################
    # Trajectory conversion  #
    ##########################

    # mdcrd -> xtc only makes sense if the run actually produced a trajectory
    if os.path.exists(global_paths["step6_godmd_run"]["output_trj_path"]) and \
       os.path.exists(global_paths["step6_godmd_run"]["output_pdb_path"]):
        global_log.info("step7_convert: Convert transition trajectory to xtc")
        cpptraj_convert(**global_paths["step7_convert"], properties=global_prop["step7_convert"])
    else:
        global_log.warning("step7_convert: Skipped, no trajectory produced by step6_godmd_run")

    #####################
    # Final output copy #
    #####################

    output_traj_path = os.path.join(output_path, 'transition.xtc')
    output_str_path = os.path.join(output_path, 'transition.pdb')
    output_log_path = os.path.join(output_path, 'godmd.log')
    output_ene_path = os.path.join(output_path, 'godmd_energy.out')

    if os.path.exists(global_paths["step7_convert"]["output_cpptraj_path"]):
        os.rename(global_paths["step7_convert"]["output_cpptraj_path"], output_traj_path)
    if os.path.exists(global_paths["step6_godmd_run"]["output_pdb_path"]):
        os.rename(global_paths["step6_godmd_run"]["output_pdb_path"], output_str_path)
    if os.path.exists(global_paths["step6_godmd_run"]["output_log_path"]):
        os.rename(global_paths["step6_godmd_run"]["output_log_path"], output_log_path)
    if os.path.exists(global_paths["step6_godmd_run"]["output_ene_path"]):
        os.rename(global_paths["step6_godmd_run"]["output_ene_path"], output_ene_path)

    # Write a stable output manifest for external consumers (see manifest.yaml in output_path)
    manifest_outputs = {}
    if os.path.exists(output_traj_path):
        manifest_outputs["trajectory"] = {"xtc": os.path.relpath(output_traj_path, output_path)}
    if os.path.exists(output_str_path):
        manifest_outputs["structure"] = {"pdb": os.path.relpath(output_str_path, output_path)}
    if os.path.exists(output_log_path) or os.path.exists(output_ene_path):
        manifest_outputs["energy"] = {}
        if os.path.exists(output_log_path):
            manifest_outputs["energy"]["log"] = os.path.relpath(output_log_path, output_path)
        if os.path.exists(output_ene_path):
            manifest_outputs["energy"]["ene"] = os.path.relpath(output_ene_path, output_path)
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

    parser = argparse.ArgumentParser(
        description="Compute an approximate conformational transition between an origin and a "
                     "target structure of a protein using GOdMD (discrete MD)."
    )

    ###############
    # Input files #
    ###############

    parser.add_argument(
        '--origin_pdb', dest='origin_pdb_path', type=str, required=False,
        help="Local PDB file for the origin structure. Mutually exclusive with --origin_pdb_code."
    )
    parser.add_argument(
        '--origin_pdb_code', type=str, required=False,
        help="PDB code to fetch the origin structure from RCSB. Mutually exclusive with --origin_pdb."
    )
    parser.add_argument(
        '--origin_chain', type=str, required=False,
        help="Chain of interest in the origin structure. Default: all chains kept."
    )
    parser.add_argument(
        '--target_pdb', dest='target_pdb_path', type=str, required=False,
        help="Local PDB file for the target structure. Mutually exclusive with --target_pdb_code."
    )
    parser.add_argument(
        '--target_pdb_code', type=str, required=False,
        help="PDB code to fetch the target structure from RCSB. Mutually exclusive with --target_pdb."
    )
    parser.add_argument(
        '--target_chain', type=str, required=False,
        help="Chain of interest in the target structure. Default: all chains kept."
    )

    #########################
    # Configuration options #
    #########################

    parser.add_argument(
        '--temp', dest='temperature', type=float, required=False, default=300,
        help="GOdMD simulation temperature in K. Default: 300"
    )
    parser.add_argument(
        '--seed', dest='random_seed', type=int, required=False, default=2525,
        help="Random seed for the simulation. Default: 2525"
    )
    parser.add_argument(
        '--nbloc', dest='num_blocks', type=int, required=False, default=10000,
        help="Number of simulation blocks. Total simulated time = nbloc * tsnap. Default: 10000"
    )
    parser.add_argument(
        '--tsnap', dest='snapshot_interval', type=int, required=False, default=500,
        help="Output structure interval, in fs. Default: 500"
    )
    parser.add_argument(
        '--rmsd_tol', dest='rmsd_tolerance', type=float, required=False, default=1.5,
        help="Minimum acceptable RMSD (Angstrom) to consider the transition converged. Default: 1.5"
    )
    parser.add_argument(
        '--gapopen', type=float, required=False, default=12.0,
        help="EMBOSS water pairwise alignment gap-open penalty. Default: 12.0"
    )
    parser.add_argument(
        '--gapextend', type=float, required=False, default=2.0,
        help="EMBOSS water pairwise alignment gap-extend penalty. Default: 2.0"
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

    args = parser.parse_args()

    godmd_simulation(
        origin_pdb_path=args.origin_pdb_path,
        origin_pdb_code=args.origin_pdb_code,
        origin_chain=args.origin_chain,
        target_pdb_path=args.target_pdb_path,
        target_pdb_code=args.target_pdb_code,
        target_chain=args.target_chain,
        temperature=args.temperature,
        random_seed=args.random_seed,
        num_blocks=args.num_blocks,
        snapshot_interval=args.snapshot_interval,
        rmsd_tolerance=args.rmsd_tolerance,
        gapopen=args.gapopen,
        gapextend=args.gapextend,
        debug=args.debug,
        restart=args.restart,
        output_path=args.output_path
    )


if __name__ == '__main__':
    main()
