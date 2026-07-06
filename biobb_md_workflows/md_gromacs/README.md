# MD with GROMACS

Run a full GROMACS MD pipeline — system setup, equilibration, production, and trajectory post-processing — starting from a prepared PDB or from pre-built structure/topology files.

## Description

**Input modes** (mutually exclusive, resolved at runtime):

- `--input_pdb` — build the system from a prepared PDB (runs `pdb2gmx`).
- `--input_gro` + `--input_top` — start from a pre-built structure (`.gro`) and topology (`.zip`), skipping `pdb2gmx`.
- `--input_tpr` + `--input_cpt` — resume a simulation from a checkpoint.

Ligand topologies produced by `ligand_parameterization` can be added with `--ligands_folder`.

**Pipeline sections:**

1. **Setup** — `pdb2gmx`, ligand insertion + restraints, simulation box, solvation, and ion generation. Skip solvation for an already-solvated input with `--skip_solvation`; skip backbone restraints with `--skip_restraints`. Stop here with `--setup_only`.
2. **Equilibration** — energy minimization, then NVT and NPT equilibration with position restraints on solute heavy atoms. Stop after this with `--equil_only`.
3. **Production** — production `mdrun` from the equilibrated structure. If `--input_plumed_path` (and optionally `--input_plumed_folder`) is given, the run uses PLUMED.
4. **Analysis & post-processing** — RMSD, radius of gyration, and RMSF, followed by trajectory cleanup (dry/center/image/fit), shared with `traj_postprocessing`.

**Execution:** single-node (`gmx` with thread-MPI/OpenMP) or multi-node (`gmx_mpi` with `--mpi_bin`/`--mpi_np`, e.g. `srun`/`mpirun`). GPU offload is enabled with `--use_gpu`.

## Usage

```bash
conda activate biobb_md
# From a prepared PDB, 100 ns
md_gromacs --input_pdb structure.pdb --prod_time 100 --temp 300 --output output

# From pre-built GROMACS files, using the GPU
md_gromacs --input_gro system.gro --input_top system.zip --use_gpu --output output
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step. Run `md_gromacs --help` for the full option list.

## Options

### Inputs

| Flag | Default | Description |
|------|---------|-------------|
| `--input_pdb` | `None` | Prepared PDB; protonation is taken from the residue names. |
| `--input_gro` | `None` | Input structure (`.gro`); use with `--input_top`. |
| `--input_top` | `None` | Input compressed topology (`.zip`); use with `--input_gro`. |
| `--input_tpr` | `None` | `.tpr` to restart a simulation; use with `--input_cpt`. |
| `--input_cpt` | `None` | `.cpt` checkpoint; use with `--input_tpr`. |
| `--ligands_folder` | `None` | Folder of ligand `.itp` + `.gro` files to include. |
| `--input_plumed_path` | `None` | Main PLUMED input file; enables PLUMED in production. |
| `--input_plumed_folder` | `None` | Folder of files referenced by the PLUMED input. |

### Execution

| Flag | Default | Description |
|------|---------|-------------|
| `--gmx_bin` | `gmx` | GROMACS binary (`gmx` single-node, `gmx_mpi` multi-node). |
| `--mpi_bin` | `null` | MPI binary path (e.g. `srun`, `mpirun`). |
| `--mpi_np` | `None` | Number of MPI processes for `mpi_bin`. |
| `--num_threads_mpi` | `0` | thread-MPI ranks (`0` = let GROMACS guess). |
| `--num_threads_omp` | `0` | OpenMP threads (`0` = let GROMACS guess). |
| `--use_gpu` | `False` | Add `-nb gpu -pme gpu` to `mdrun`. |
| `--debug` | `False` | Verbose logging and keep temporary files. |

### Simulation

| Flag | Default | Description |
|------|---------|-------------|
| `--forcefield` | `amber99sb-ildn` | pdb2gmx force field. Available force fields depend on the GROMACS version |
| `--ions_concentration` | `0.15` | Salt concentration (mol/L). |
| `--temp` | `300` | Temperature (K). |
| `--seed` | `-1` | Velocity-generation random seed (`-1` = random). |
| `--dt` | `2` | Time step (fs; 1–4). |
| `--equil_time` | `1.0` | Time per equilibration step (ns). |
| `--equil_frames` | `500` | Frames saved during equilibration. |
| `--prod_time` | `100.0` | Total production time (ns). |
| `--prod_frames` | `2000` | Frames saved during production. |
| `--setup_only` | `False` | Only set up the system (stop before minimization). |
| `--equil_only` | `False` | Only run setup + equilibration. |
| `--skip_restraints` | `False` | Skip backbone position restraints (input_pdb / input_gro+top modes). |
| `--skip_solvation` | `False` | Skip box/solvent/ions (input already solvated). |
| `--remove_raw_traj` | `False` | Delete raw production trajectories after post-processing. |
| `--keep_solvent` | `False` | Keep solvent and ions in the post-processed trajectory. |
| `--keep_residues` | `None` | Extra residue indices to keep in the post-processed trajectory (e.g. `15 23 105`). |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `output` | Output directory. |

## Outputs

Written into `--output`, organized by section:

- `1_setup/`, `2_equil/`, `3_prod/`, `4_analysis/` — per-section working directories.
- Post-processed trajectory and structure under `4_analysis/` (e.g. `fitted_traj.xtc`, `dry_structure.pdb`), plus the RMSD/Rgyr/RMSF analysis outputs.
- `config.yml` and `log.out` for inspection.
