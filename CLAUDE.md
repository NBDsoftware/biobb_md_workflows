# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Style

Be extremely concise in interactions and commit messages — sacrifice grammar for concision. At the end of any plan, list unresolved questions (if any), extremely concisely.

## Project overview

BioBB Workflows for MD simulations with GROMACS: four CLI pipelines built on top of BioExcel Building Blocks (BioBB), developed for the European BioExcel project. Each pipeline lives in its own package under `biobb_md_workflows/` and is exposed as a standalone CLI command:

| Command | Module | Purpose |
|---------|--------|---------|
| `protein_preparation` | `biobb_md_workflows/protein_preparation/protein_preparation.py` | Clean/protonate/fix a raw PDB structure for MD |
| `ligand_parameterization` | `biobb_md_workflows/ligand_parameterization/ligand_parameterization.py` | Generate GROMACS/AMBER topology+coords for ligands/cofactors (GAFF via acpype, or custom AMBER params via LEaP) |
| `md_gromacs` | `biobb_md_workflows/md_gromacs/md_gromacs.py` | Full GROMACS MD: setup → equilibration → production → trajectory post-processing |
| `traj_postprocessing` | `biobb_md_workflows/traj_postprocessing/traj_postprocessing.py` | Standalone post-processing (strip solvent, center, image, fit) of an existing raw trajectory |

## Setup / install

```bash
conda env create -f environment.yml   # creates the `biobb_md` env, installs this repo with `pip install -e .`
conda activate biobb_md
```

`pyproject.toml` requires exactly Python 3.11. Two dependencies (`biobb_gromacs`, `biobb_analysis`) are installed from NBDsoftware forks on GitHub, not PyPI — required by `md_gromacs` and `traj_postprocessing`.

After editing `pyproject.toml` entry points or dependencies, reinstall with `pip install .`.

Important note regarding installation. The `modeller` package is only intended for academic use and requires the export of a key before the installation of the environment. See the README file.

## Running / testing

There is no automated test suite (no pytest). `tests/<workflow>/` contains SLURM job scripts (`run.sl`), example configs (`input.yml`), and recorded output trees from real runs — used as manual smoke tests, not CI. To exercise a workflow, adapt and run its `run.sl` or invoke the CLI directly, e.g.:

```bash
protein_preparation --input_pdb data/1r9o.pdb --output output --ph 7 --cap_ter --output_format gromacs
protein_preparation --help   # every workflow supports --help
```

The `config.yml` is always auto-generated from the CLI arguments (there is no `--config` override). CLI options are injected into the generated YAML template wherever possible; the generated file is written into the output directory for inspection/reproducibility.

## Architecture

All four workflows share the same shape, built around `biobb_common`:

1. **YAML config-driven steps.** Each workflow defines a `config_contents()` returning a YAML string (an f-string that interpolates CLI-derived values, using `common.to_yaml()` to render nullable/bool/list scalars safely), always written to `<output>/config.yml` by `create_config_file(output_path, **config_args)`. The YAML has a `global_properties` block (`working_dir_path`, `restart`, `remove_tmp`, `can_write_console_log`) and one block per step named `stepN_name` (or prefixed `sectionprefix_stepN_name`, see below), each with a `tool`, a `paths` dict, and a `properties` dict.
2. **Dependency chaining.** Step `paths` reference earlier steps' outputs as `dependency/stepX_name/output_key` strings; `settings.ConfReader` resolves these into real paths at runtime via `conf.get_paths_dic()` / `conf.get_prop_dic()`.
3. **Workflow function body.** The main function (`protein_preparation()`, `ligand_parameterization()`, `md_gromacs()`, `traj_postprocessing()`) reads the config, then calls each BioBB building block in sequence as `tool_fn(**global_paths["stepN"], properties=global_prop["stepN"])`. Prefer injecting CLI options directly into the generated YAML in `config_contents()`; only mutate `global_paths`/`global_prop` after reading the config for values that cannot be known at config-generation time — runtime-derived data (protonation states parsed from the structure, propka pKa results, altlocs) or per-item iteration (per-ligand, per-chain, per-residue). Steps are logged individually via `global_log`.
4. **Restart support.** `restart: True` in `global_properties` + `biobb_common`'s own step-skipping means re-running a workflow against the same `output_path` resumes from the last completed step rather than redoing everything.
5. **Prefixed sub-configs.** `ligand_parameterization` calls `conf.get_prop_dic(prefix=ligand_name)` per ligand so multiple ligands share one YAML template but get independent step namespaces. `md_gromacs` does the same per pipeline section (`1_setup`, `2_equil`, `3_prod`, `4_analysis`).
6. **CLI = thin argparse wrapper.** Each module's `main()` just parses args and calls the workflow function; the workflow function is also importable directly (returns `(global_paths, global_prop)`).

### protein_preparation specifics

Fixes structure defects (altlocs, backbone/side-chain gaps via Modeller — needs `--modeller_key`, chirality, amides, SS bonds), applies mutations, then determines protonation states of titratable residues (LYS/ARG/ASP/GLU/HIS) via `propka` pKa estimation + AmberTools `reduce` for HIS H-bond optimization, and renames residues accordingly. `biobb_propka` and `biobb_titrate` (top of the file) are **mock functions standing in for not-yet-published real BioBBs** — replace them once those land upstream. Residue-naming format conversion between `standard`/`gromacs`/`amber` (CYS/CYX, HIS variants, ACE/NME terminal atoms) is done with hand-rolled fixed-column PDB line rewriting (`rename_his`, `rename_ss_bonds`, `rename_ter`) — column offsets are load-bearing, be careful editing them.

### ligand_parameterization specifics

Per ligand found in the input PDB (heteroatom records), branches on whether a custom AMBER parameter set (`.frcmod`+`.prep`, matched by ligand name) is supplied:
- **Custom params → path A**: `leap_gen_top` (AMBER) then optionally `acpype_convert_amber_to_gmx`.
- **No custom params → path B**: protonate (`ambertools`/`obabel`/`none`) → minimize (`babel_minimize`) → `acpype_params_gmx` or `acpype_params_ac` (GAFF).

Output `.top` files get converted to `.itp` via `gmx_top2itp` (strips `[ defaults ]`/`[ molecules ]` sections) since `.itp` files are meant to be `#include`d into a master topology, not standalone.

### md_gromacs specifics

Three mutually-exclusive input modes, resolved by `check_inputs()`: `input_pdb` (build from scratch), `input_gro_top` (pre-built structure+topology, skips `pdb2gmx`), `restart_simulation` (`.tpr`+`.cpt`, resumes an in-flight run). Pipeline sections: `1_setup` (pdb2gmx, ligand insertion + restraints, solvation, ion generation) → `2_equil` (minimization, NVT, NPT) → `3_prod` (production `mdrun`, optionally via `mdrun_plumed` if `--input_plumed_path`/`--input_plumed_folder` given) → `4_analysis` (RMSD/Rgyr/RMSF, then dry/center/image/fit trajectory post-processing, shared with `traj_postprocessing`). `--setup_only` / `--equil_only` stop the pipeline early. Supports single-node (`gmx` + thread-MPI/OpenMP) and multi-node (`gmx_mpi` + `mpi_bin`/`mpi_np`, e.g. `srun`/`mpirun`) execution, and GPU offload via `--use_gpu`.

### traj_postprocessing specifics

Standalone version of `md_gromacs`'s stage 4 trajectory cleanup, for post-processing a trajectory produced outside this repo. Builds GROMACS index-group selections by auto-detecting solvent/ion residue names from the structure (extendable via `--extra_ions`/`--extra_solvents`), then runs dry → center → image → fit. `--fast` picks a cheaper post-processing path (`fast_postprocessing`) vs the full one (`complete_postprocessing`).

## Gotchas

- AMBER and GROMACS use different residue-naming conventions for protonation/disulfide state (e.g. `HID`/`HIE`/`HIP` vs `HISD`/`HISE`/`HISH`, `CYX` vs `CYS2`); workflows convert between them at specific points — check `output_format`/`format` args before assuming a PDB's naming convention.
- PDB line rewriting in `protein_preparation.py` (`rename_his`, `rename_ss_bonds`, `rename_ter`, and the 4-letter-resname fixup in `biobb_titrate`) relies on fixed-column slicing (PDB columns 18-21, 13-16, etc.) — a malformed/nonstandard input PDB can silently produce wrong output here.
- `tests/*/output/` directories are recorded real-run outputs used for reference/debugging, not fixtures regenerated by a test runner — don't assume they're currently reproducible from the exact commands in the accompanying `run.sl` without checking data/versions.
