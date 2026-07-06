# Protein preparation

This workflow can be used to fix PDB defects, choose protonation states for tritatable residues and prepare the system for simulation. 

## Quick installation and run

---

Install the repo's conda environment:

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f ../../environment.yml
conda activate biobb_md
```


### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.

```bash
conda activate biobb_md
protein_preparation --help
```

```
usage: MD Simulation with GROMACS [-h] [--input_pdb INPUT_PDB_PATH] [--ligands_folder LIGANDS_TOP_FOLDER] [--input_gro INPUT_GRO_PATH] [--input_top INPUT_TOP_PATH]
                                  [--input_tpr INPUT_TPR_PATH] [--input_cpt INPUT_CPT_PATH] [--input_plumed_path INPUT_PLUMED_PATH] [--input_plumed_folder INPUT_PLUMED_FOLDER]
                                  [--config CONFIG_PATH] [--gmx_bin GMX_BIN] [--mpi_bin MPI_BIN] [--mpi_np MPI_NP] [--num_threads_mpi NUM_THREADS_MPI] [--num_threads_omp NUM_THREADS_OMP]
                                  [--use_gpu] [--restart] [--forcefield FORCEFIELD] [--ions_concentration IONS_CONCENTRATION] [--temp TEMPERATURE] [--seed RANDOM_SEED] [--setup_only]
                                  [--skip_restraints] [--skip_solvation] [--dt DT] [--equil_time EQUIL_TIME] [--equil_frames EQUIL_FRAMES] [--equil_only] [--prod_time PROD_TIME]
                                  [--prod_frames PROD_FRAMES] [--remove_raw_traj] [--keep_solvent] [--keep_residues RESIDUES_TO_KEEP [RESIDUES_TO_KEEP ...]] [--debug]
                                  [--output OUTPUT_PATH]

options:
  -h, --help            show this help message and exit
  --input_pdb INPUT_PDB_PATH
                        Input PDB file. The workflow assumes the protonation state specified by the residue names is the correct one. Default: None
  --ligands_folder LIGANDS_TOP_FOLDER
                        Path to folder with .itp and .gro files for the ligands that should be included in the simulation. Compatible with '--input_pdb' and '--input_gro'/'--input_top'.
                        Default: None
  --input_gro INPUT_GRO_PATH
                        Input structure file (.gro). Use together with '--input_top'. Restraints and ligands can be added; use '--skip_solvation' if the system is already solvated.
                        Default: None
  --input_top INPUT_TOP_PATH
                        Input compressed topology file (.zip). Use together with '--input_gro'. Default: None
  --input_tpr INPUT_TPR_PATH
                        Input portable binary run input file (.tpr) to restart a simulation. Use together with '--input_cpt'. Default: None
  --input_cpt INPUT_CPT_PATH
                        Input checkpoint file (.cpt) to restart a simulation. Use together with '--input_tpr'. Default: None
  --input_plumed_path INPUT_PLUMED_PATH
                        Path to the main PLUMED input file (plumed.dat). If provided, PLUMED will be used during the production run. Default: None
  --input_plumed_folder INPUT_PLUMED_FOLDER
                        Path to the folder with all files needed by the main PLUMED input file, see input_plumed_path. Default: None
  --config CONFIG_PATH  Configuration file (YAML)
  --gmx_bin GMX_BIN     Path to GROMACS binary (gmx for single node and gmx_mpi for multi-node). Default: gmx
  --mpi_bin MPI_BIN     Path to MPI binary. Default: null
  --mpi_np MPI_NP       Number of MPI processes given to the mpi_bin. Default: None
  --num_threads_mpi NUM_THREADS_MPI
                        Number of MPI threads. Default: 0 (Let GROMACS guess)
  --num_threads_omp NUM_THREADS_OMP
                        Number of OpenMP threads. Default: 0 (Let GROMACS guess)
  --use_gpu             Calculate non-bonding interactions and particle-mesh ewald in GPU by adding '-nb gpu -pme gpu' to mdrun call. If not used, gmx will still use a GPU for these
                        calculations if available. Default: False
  --restart             Restart the workflow from the last completed step. Default: False
  --forcefield FORCEFIELD
                        Forcefield to use. Default: amber99sb-ildn
  --ions_concentration IONS_CONCENTRATION
                        Concentration of ions in the system in mol/L (M). Default: 0.15 M
  --temp TEMPERATURE    Temperature of the system in K. Default: 300
  --seed RANDOM_SEED    Random seed for the simulations. If given, new velocities will be generated with this seed. Default: -1
  --setup_only          Only setup the system. Default: False
  --skip_restraints     Skip adding chain backbone position restraints to the topology. Only used for input_pdb and input_gro_top modes. Ligand restraints are always added if a ligands
                        folder is provided. Default: False
  --skip_solvation      Skip adding simulation box, solvent and ions (editconf, solvate, grompp, genion). Use when the input structure is already solvated. Only used for input_pdb and
                        input_gro_top modes. Default: False
  --dt DT               Time step in fs. Default: 2 fs
  --equil_time EQUIL_TIME
                        Time of each equilibration step in ns. Default: 1.0 ns
  --equil_frames EQUIL_FRAMES
                        Number of frames to save during the equilibration steps. Default: 500 frames
  --equil_only          Only run the equilibration steps. Default: False
  --prod_time PROD_TIME
                        Total time of the production simulation in ns. Default: 100.0 ns
  --prod_frames PROD_FRAMES
                        Number of frames to save during the production steps. Default: 2000 frames
  --remove_raw_traj     Delete the heavy, raw production trajectories after post-processing is complete to save disk space. Default: False
  --keep_solvent        Keep solvent and ions in the final post-processed trajectory. Default: False
  --keep_residues RESIDUES_TO_KEEP [RESIDUES_TO_KEEP ...]
                        List of specific residue indices to keep in the final post-processed trajectory (e.g., --keep_residues 15 23 105). Default: None
  --debug               Activate debug mode with more verbose logging. Default: False
  --output OUTPUT_PATH  Output path. Default: 'output' in the current working directory
```

## Description

1. **Extraction of structure from PDB**

2. **Fix PDB defects (A-I)**
    Steps to fix different possible defects in the input pdb structure. See below.

    1. **Fix alternative locations** 
    Provide a list with the choices of alternative locations to keep in the final structure. If no list is given (_null_ value) it will select the alternative location with the highest occupancy (the workflow will use Biopython to do so). 

    2. **Mutate initial pdb structure** 
    Mutations can be requested through the mutation_list command line argument. Where each mutation is defined by string with the following format: "Chain:Wild_type_residue_name Residue_number Mutated_type_residue_name". The residue name should be a 3 letter code starting with capital letters, e.g. "A:Arg220Ala".

    3. **Obtain the Sequence in FASTA format** 
    The sequence is then used to model missing backbone atoms in the next step. The workflow first tries to download the canonical FASTA (including all residues for that protein) from the Protein Data Bank. If there is no internet connection, it will try to obtain the sequence from the _SEQRES_ records in the PDB. If there are no _SEQRES_, then only the residues that contain at least one atom in the structure will be included. This step can be skipped including the ```--skip_bc_fix``` option.  

    4. **Model missing backbone atoms**
    Add missing backbone heavy atoms using _biobb_structure_checking_ and Modeller suite. A modeller license key and the previous FASTA file are required for this step. This step can be skipped including the ```--skip_bc_fix``` option.  

    5. **Model missing side chain atoms**
    Add missing side chain atoms using _biobb_structure_checking_ (and Modeller suite if a license key is provided).

    6. **Renumber atomic indices**
    So they start at 1.

    7. **Relieve clashes flipping amide groups**
    It flips the clashing amide groups to relieve clashes.

    8. **Fix chirality of residues**
    Creates a new PDB file fixing stereochemical errors in residue side-chains changing it's chirality when needed.

    9. **Add missing disulfide bonds**
    It changes CYS for CYX to mark cysteines residues pertaining to a [di-sulfide bond](https://en.wikipedia.org/wiki/Disulfide). It uses a distance criteria to determine if nearby cysteines are part of a di-sulfide bridge (_check_structure getss_). Use carefully, this step can be skipped using ```--skip_ss_bonds```

    10. **Remove all hydrogens**

    11. **Estimate the pKa of titratable residues with propka**
    Use an empirical method to estimate the pKa of residues considering the local environment. 

    12. **Optimize the hydrogen placement in Histidines with reduce**
    Use reduce from AmberTools to optimize the hydrogen bonds of histidines.

    13. **Select protonation state for titratable residues** 
    Taking into account previous two steps and pH

    14. **Add hydrogens back**
    If output format is amber.