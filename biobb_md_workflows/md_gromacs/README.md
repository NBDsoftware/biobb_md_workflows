# MD with GROMACS

This workflow uses BioBBs to fix PDB defects, prepare the MD simulations, equilibrate and execute production runs, do basic analysis of the trajectories and post-process the trajectories.

![alt text](../../img/MD_setup.png?raw=true)

## Quick installation and run

---

Install the repo's conda environment (running in Nostrum's cluster use the already installed environments located in */shared/work/BiobbWorkflows/envs*)

```bash
export KEY_MODELLER="HERE YOUR MODELLER KEY"
conda env create -f ../../environment.yml
conda activate biobb_md
```

See [biobb documentation](https://mmb.irbbarcelona.org/biobb/documentation/source) for additional properties not included in the YAML configuration file.

To run a single call to the workflow in an HPC environment use:

```bash
sbatch run_HPC.sl
```

To run a long MD simulation while respecting the time limit of the HPC jobs use:

```bash
./launch_long_MD.sh
```

## Inputs

---

### Configuration file

Take a look at the YAML configuration file to see the different properties that can be set.

```bash
vi input.yml
```

Specially important are: the binary path of GROMACS and the MODELLER key. Make sure the binary path specified and the module loaded in the run file (HPC only) agree between them.

### Command line arguments

The command line arguments can be used to provide some inputs and settings that will be prioritized over those in the YAML configuration file.


```bash
conda activate biobb_md
md_gromacs --help
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

This workflow has several steps. The input for the workflow can be (1) a prepared PDB file or (2) a gromacs structure file and .zip topology files ready to be minimized.

3. **Preparation of topology and coordinates for MD**

    A. **pdb2gmx** 
    Uses pdb2gmx to obtain a gromacs structure file (.gro) and topology file from the fixed PDB. Hydrogen atoms will be added in this step, one can choose to ignore the hydrogens in the original structure or not (```ignh``` property). The protonation state of histidines can be provided (```his``` property) in the form of a list of numbers see below. A force field and water model are chosen here.
    
    For the ```his``` property include a string with the protonation states '0 0 1 1 0 0 0', where:

        - 0 : H on ND1 only (HID)
        - 1 : H on NE2 only (HIE)
        - 2 : H on ND1 and NE2 (HIP)
        - 3 : Coupled to Heme (HIS1)
    
    NOTE: default behavior is to add charged termini - if one wants ACE and NME it should be provided already with the correct atom names - look at the force field being used: /eb/x86_64/software/GROMACS/2023.3-foss-2022a-CUDA-11.7.0-PLUMED-2.9.0/share/gromacs/top/amber99sb-ildn.ff/aminoacids.rtp
    
    NOTE: this will be automatized within the workflow in the future, right now pdb4amber can be used externally to obtain a guess

    B. **Generate cofactors topology**: if there are any cofactors with parameters in the `--ligand_parameters` folder, use _tleap_ to build the corresponding AMBER topology and coordinates file.

    C. **Convert AMBER topology and coordinates to GMX**: convert AMBER topology and coordinate files to GROMACS format.

    D-F: **Merge cofactors and structure**: if any parameterized cofactors are present, convert the structures to PDB format and concatenate them.

    G-H. **Generate restraints for cofactor heavy atoms**.

    I. **Append cofactor topology to structure topology**.

    J. **Generate simulation box** 
    Generate a simulation box around the structure. Some box types are more efficient than others (octahedron for globular proteins)

    K. **Solvation** 
    Generate a box of solvent around the structure.

    L-M. **Add ions** 
    Randomly replace solvent molecules with monoatomic ions. 

    N. **Convert final topology to PDB**
    This will be used in subsequent post-processing steps.

To prepare the system externally, use the ```--input_gro``` and ```--input_top``` command line arguments.

To make sure the system has been correctly prepared before minimizing or running MD, launch the workflow adding the ```--setup_only``` command line option. This will stop the workflow before the energy minimization. 

4. **Minimize and equilibrate the initial configuration**

    A-D. **Energy minimization** (including position restraints on protein/DNA heavy atoms)

    E-G. **NVT equilibration** (including position restraints on protein/DNA heavy atoms)

    H-J. **NPT equilibration** (including position restraints on protein/DNA heavy atoms)

5. **Production run**

    Launch several production trajectories from equilibrated structure (see `--num_parts` or `--num_replica` command line argument).

6. **Basic analysis**

    Computation of RMSD (with fitting) with respect to experimental structure and with respect to equilibrated structure (protein backbone atoms). Computation of Radius of gyration (protein backbone atoms) and RMSF (protein heavy atoms).

7. **Post-processing**

    Image, dry and fitting.

Note that re-launching the workflow will skip the previously successful steps if restart is True and the output folder is the same. 

