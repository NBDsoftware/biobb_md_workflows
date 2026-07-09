# MD with GROMACS

Run a full GROMACS MD pipeline — system setup, equilibration, production, and trajectory post-processing — starting from a prepared PDB or from pre-built structure/topology files.

## Description

**Input modes** (mutually exclusive, resolved at runtime):

- `--input_pdb` — build the system from a prepared PDB (runs `pdb2gmx`).
- `--input_gro` + `--input_top` — start from a pre-built structure (`.gro`) and topology (`.zip`), skipping `pdb2gmx`.
- `--input_tpr` + `--input_cpt` — resume a simulation from a checkpoint.

Ligand topologies produced by `ligand_parameterization` can be added with `--ligands_folder`.

**Pipeline sections:**

1. **Setup** — generation of protonated GROMACS topology and coordinates with `pdb2gmx`, addition of restraints to all chains (protein and NA), insertion of ligand to the topology, generation of restraints for the ligand, creation of simulation box, solvation of the system, and generation of ions. Skip solvation for an already-solvated input with `--skip_solvation`; skip backbone restraints with `--skip_restraints` if the `.top` already includes the restraints you need. Stop here with `--setup_only`.
2. **Equilibration** — energy minimization, then NVT and NPT equilibration with position restraints on solute heavy atoms. Stop after this with `--equil_only`.
3. **Production** — production `mdrun` from the equilibrated structure. If `--input_plumed_path` (and optionally `--input_plumed_folder`) is given, the run uses PLUMED.
4. **Analysis & post-processing** — RMSD, radius of gyration, and RMSF, followed by trajectory cleanup (dry/center/image/fit), shared with `traj_postprocessing`.

**Execution:** single-node (`gmx` with thread-MPI/OpenMP) or multi-node (`gmx_mpi` with `--mpi_bin`/`--mpi_np`, e.g. `srun`/`mpirun`). GPU offload of the non-bonded and PME work is enabled with `--use_gpu` (`-nb gpu -pme gpu`).

## Usage

```bash
conda activate biobb_md
# From a prepared PDB, 100 ns
md_gromacs --input_pdb structure.pdb --prod_time 100 --temp 300 --output output

# From pre-built GROMACS files, using the GPU
md_gromacs --input_gro system.gro --input_top system.zip --use_gpu --output output
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `md_gromacs --help` for the full option list.

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


## Recommendations

### Improving performance

GROMACS auto-detects the hardware and, left alone, makes near-optimal use of it: the
defaults (`--num_threads_mpi 0`, `--num_threads_omp 0`) let `mdrun` choose the rank/thread
split. Tune only when you need to fill a specific allocation, and only scale out as far as
the system justifies — small systems saturate quickly, and adding ranks past the scaling
limit wastes the allocation.

- **Single node (default).** Uses the built-in thread-MPI. `--num_threads_mpi` sets the number of thread-MPI ranks (`-ntmpi`), `--num_threads_omp` the OpenMP threads per rank (`-ntomp`); their product should equal the cores you were given.
- **GPU (`--use_gpu`).** Offloads the non-bonded and PME work to the GPU (`-nb gpu -pme gpu`); minimization always runs on CPU. For a single GPU one MPI rank is usually fastest (`--num_threads_mpi 1`) but you can add OpenMP threads; with N GPUs the number of MPI ranks must equal N. 
- **Multi-node.** Requires a GROMACS built with external MPI: set `--gmx_bin gmx_mpi` and launch it through `--mpi_bin` (`srun`/`mpirun`). One rank per core scales well down to ~200 particles/core; beyond that, add OpenMP threads per rank. Match `--num_threads_mpi` and `--num_threads_omp` to your SLURM `--ntasks` and `--cpus-per-task` respectively.

### Simulation parameters

The production protocol uses the leap-frog algorithm for integrating Newton’s equations of 
motion. The neighbor search is done with the Verlet cut-off scheme. Electrostatics are computed 
with the Fast smooth Particle-Mesh Ewald (SPME) algorithm. Van der Waals and Coulomb have a 1.0 nm, 
cut-off, long-range dispersion correction for energy and pressure are applied togehter with LINCS
constraints on bonds to hydrogen. The V-rescale thermostat and Parrinello-Rahman barostat are used. 
See a copy of the input .mdp files below. 

- **Force field (`--forcefield`, default `amber99sb-ildn`).** Must be one your GROMACS build
  ships (run `pdb2gmx` to list them). The **water model is hard-coded to TIP3P** and the
  non-bonded settings are fixed, thus a force field that expects a different water model 
  or non-bonded treatment is inconsistent here (e.g. CHARMM* or ff19SB**).
- **Time step (`--dt`, default 2 fs).** Only bonds to hydrogen are constrained, which makes
  2 fs safe. The accepted range is 1–4 fs, but there is no hydrogen-mass repartitioning, so
  **4 fs is not stable** here — stay at 2 fs.
- **Temperature (`--temp`, default 300 K).** Sets both the thermostat reference and the
  initial velocity generation. Solute and solvent use different temperature-coupling groups
  and thus different V-rescale thermostats; pressure is held at 1 bar (Parrinello-Rahman)
  during NPT and production.
- **Salt (`--ions_concentration`, default 0.15 mol/L).** Added after neutralizing the
  system's net charge.
- **Equilibration length (`--equil_time`, default 1 ns per phase).** Short by default (1 ns
  NVT + 1 ns NPT, heavy-atom position restraints on the solute); increase it for large,
  membrane, or slowly-relaxing systems before trusting the production run.

*: CHARMM employs a specific Lennard-Jones potential that uses a force-switching function, the .mdp must be carefully configured to replicate these physics correctly - not supported yet.

**: needs OPC water model and matching ions

## Output

Files:

Written into `--output`, organized by section: `1_setup/`, `2_equil/`, `3_prod/`, `4_analysis/`.

- Post-processed trajectory and structure under `4_analysis/` (e.g. `step10_fit_traj/fitted_traj.xtc`, `step6_dry_str/dry_structure.pdb`), plus the RMSD/Rgyr/RMSF analysis outputs.
- `config.yml` and `log.out` for inspection.

## Limitations

- **Fixed protocol.** Thermostat (V-rescale), barostat (Parrinello-Rahman, isotropic 1 bar),
  cut-offs (1.0 nm), Verlet/PME scheme, dispersion correction, and the box (truncated
  octahedron, 1.0 nm padding) cannot be changed for now.
- **Water model fixed to TIP3P**, and the non-bonded settings are AMBER-style; other water
  models / force-field families will be supported.
- **No hydrogen-mass repartitioning**, so the time step is effectively capped at 2 fs.
- **Restraints are all-or-nothing.** Solute heavy atoms are restrained at full strength
  through both equilibration phases and then fully released for production; there is no
  selective or progressively-released restraint schedule.
- **PLUMED runs only in production**, not during equilibration.

## Reference

### MDP files used

For the minimization

```
;Neighbour searching
cutoff-scheme = Verlet
ns-type = grid
rcoulomb = 1.0
vdwtype = cut-off
rvdw = 1.0
nstlist = 10
rlist = 1

;Eletrostatics
coulombtype = PME

;Periodic boundary conditions
pbc = xyz
ld-seed = 1
```

For the NVT equilibration:

```
;Bond parameters
constraint-algorithm = lincs
constraints = h-bonds
lincs-iter = 1
lincs-order = 4
continuation = no

;Neighbour searching
cutoff-scheme = Verlet
ns-type = grid
rcoulomb = 1.0
vdwtype = cut-off
rvdw = 1.0
nstlist = 10
rlist = 1

;Eletrostatics
coulombtype = PME
pme-order = 4
fourierspacing = 0.12
fourier-nx = 0
fourier-ny = 0
fourier-nz = 0
ewald-rtol = 1e-5

;Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1	  0.1
ref-t = 300 	  300

;Pressure coupling
pcoupl = no

;Dispersion correction
DispCorr = EnerPres

;Velocity generation
gen-vel = yes
gen-temp = 300
gen-seed = -1

;Periodic boundary conditions
pbc = xyz
ld-seed = 1
``` 

For the NPT equilibration:

```
;Bond parameters
constraint-algorithm = lincs
constraints = h-bonds
lincs-iter = 1
lincs-order = 4
continuation = yes

;Neighbour searching
cutoff-scheme = Verlet
ns-type = grid
rcoulomb = 1.0
vdwtype = cut-off
rvdw = 1.0
nstlist = 10
rlist = 1

;Eletrostatics
coulombtype = PME
pme-order = 4
fourierspacing = 0.12
fourier-nx = 0
fourier-ny = 0
fourier-nz = 0
ewald-rtol = 1e-5

;Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1	  0.1
ref-t = 300 	  300

;Pressure coupling
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 1.0
ref-p = 1.0
compressibility = 4.5e-5
refcoord-scaling = com

;Dispersion correction
DispCorr = EnerPres

;Velocity generation
gen-vel = no

;Periodic boundary conditions
pbc = xyz
ld-seed = 1
```

For the production run

```
;Bond parameters
constraint-algorithm = lincs
constraints = h-bonds
lincs-iter = 1
lincs-order = 4
continuation = yes

;Neighbour searching
cutoff-scheme = Verlet
ns-type = grid
rcoulomb = 1.0
vdwtype = cut-off
rvdw = 1.0
nstlist = 10
rlist = 1

;Eletrostatics
coulombtype = PME
pme-order = 4
fourierspacing = 0.12
fourier-nx = 0
fourier-ny = 0
fourier-nz = 0
ewald-rtol = 1e-5

;Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1	  0.1
ref-t = 300 	  300

;Pressure coupling
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 1.0
ref-p = 1.0
compressibility = 4.5e-5
refcoord-scaling = com

;Dispersion correction
DispCorr = EnerPres

;Velocity generation
gen-vel = no

;Periodic boundary conditions
pbc = xyz
ld-seed = 1
```



