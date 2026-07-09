# Trajectory post-processing

Dry, center, image, and fit a GROMACS MD trajectory. It follows the trajectory post-processing recommended by the GROMACS user manual and forums.

## Description

The workflow auto-detects solvent and ion residue names from the input structure to build the GROMACS index-group selections, then cleans the trajectory:

- **Complete** (default): make the solute whole → cluster → extract reference → remove jumps → center → image → fit.
- **Fast** (`--fast`): center → image → fit only (skips the whole/nojump/cluster steps).

Solvent and ions are stripped unless `--keep_solvent` is given; extra residues can be retained with `--keep_residues`. The default solvent/ion molecules recognized by GROMACS as Solvent/Ion groups can be extended with `--ions` and `--solvents`.

Because centering and imaging results are system-dependent, always inspect the result — this workflow will probably not be adequate to all systems.

## Usage

```bash
conda activate biobb_md
traj_postprocessing --input_traj traj.xtc --input_top run.tpr --input_structure structure.gro --output output
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `traj_postprocessing --help` for the full option list.

## Options

### Inputs

| Flag | Default | Description |
|------|---------|-------------|
| `--input_traj` | *required* | Input trajectory (`.xtc`). |
| `--input_top` | *required* | Binary run input file (`.tpr`). |
| `--input_structure` | *required* | Structure (`.gro`/`.pdb`) defining the solvent/output index groups and center group. Must not be broken across PBC. |

### Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--gmx_bin` | `gmx` | GROMACS binary path. |
| `--fast` | `False` | Skip making solute whole, removing jumps, and clustering. |
| `--keep_solvent` | `False` | Keep solvent and ions in the output. |
| `--keep_residues` | `None` | Extra residue indices to keep besides the solute (e.g. `15 23 105`). |
| `--ions` | `[]` | Extra ion atom names for the solvent group (e.g. `NA+ CA2+`). |
| `--solvents` | `[]` | Extra solvent residue names for the solvent group (e.g. `TIP3 TIP4`). |
| `--debug` | `False` | Keep intermediate files. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `output` | Output directory. |
| `--output_traj` | `trajectory.xtc` | Output trajectory file name. |
| `--output_str` | `structure.pdb` | Output structure file name. |


## Options recommendations

### Complete vs. fast processing

:::{admonition} 🚧 To be written
:class: caution
When to use the default `complete` path (whole → cluster → nojump → center →
image → fit) versus `--fast` (center → image → fit only), and **why centering and
imaging are system-dependent** — plus how to visually validate the processed
trajectory (e.g. inspect in VMD/PyMOL, check the solute stays whole and centered).
:::

### Selecting what to keep

:::{admonition} 🚧 To be written
:class: caution
Guidance on `--keep_solvent`, `--keep_residues`, and extending the auto-detected
groups with `--ions` / `--solvents` when solvent or ion residue names are not
recognized automatically.
:::

## Output

Checklist:

Files:

Written into `--output`:

- The post-processed trajectory (`--output_traj`, default `trajectory.xtc`) and structure (`--output_str`, default `structure.pdb`).
- `config.yml`, `log.out`, and per-step directories for inspection.

## Limitations

- Centering and imaging are **system-dependent**. The recommended steps are not
  appropriate for every system — always visually inspect the result.
