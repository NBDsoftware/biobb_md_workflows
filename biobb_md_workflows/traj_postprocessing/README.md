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


## Recommendations

### Complete vs. fast processing

Each operation is a separate `gmx trjconv` call because PBC treatment and fitting cannot be combined in
one pass — order matters, and the default **complete** path is the robust general case:

1. **whole** (`-pbc whole`) — reconnect molecules broken across the box edge.
2. **cluster** (`-pbc cluster`) — wrap all atoms around the solute center of mass (iteratively-updated).
3. **extract reference** then **nojump** (`-pbc nojump`) — undo box jumps for a continuous trajectory.
   The reference is the *clustered first frame* (dumped from step 2), not the raw run input: `nojump`
   removes jumps relative to its `-s` reference, so a non-whole/non-clustered reference would undo
   steps 1–2. This is the key ordering subtlety in the GROMACS suggested workflow.
4. **center** (`-center`, `-ur compact`) on the atom nearest the solute's geometric center — this shifts
   the system, so no `nojump` may follow.
5. **image** (`-pbc mol`, `-ur compact`) — put every molecule's center of mass back in the box.
6. **fit** (`-fit rot+trans` on the solute) — remove global rotation/translation.

`--fast` keeps only center → image → fit. Use it when the solute is already whole and does not jump
across the box (small single-chain solute, or a trajectory already made whole upstream).

**Always inspect the result** in a viewer (VMD/PyMOL): `-pbc cluster` is only meaningful if you genuinely have a cluster and the molecules are not already broken; `-pbc nojump` keeps the trajectory continuous but lets molecules diffuse out of the box; and `-ur compact` is intended for the truncated-octahedron/rhombic-dodecahedron boxes `md_gromacs` produces. Confirm the solute stays whole and centered and that imaging did not split it.

### Selecting what to keep

By default only the solute is written and solvent/ions are stripped — the main disk-space saving of this
step, and the recommended dry `.xtc` for analysis. The Solvent group is built by auto-detecting solvent
and ion names from `--input_structure`; GROMACS's default `SOL`/`Ion` groups plus a few common ions
(`K+`, `CL-`, `MG`) are recognized. If your solvent or ions carry **unusual names**, add them with
`--solvents` (residue names, e.g. `TIP4 SPC`) / `--ions` (atom names, e.g. `NA+ CA2+`), otherwise they
are counted as part of the solute and never stripped. To keep them, use `--keep_solvent` (writes the
whole `System`); to retain a few specific residues (a ligand, key waters) alongside the solute, list
their residue indices with `--keep_residues`.

## Output

Checklist:

Files:

Written into `--output`:

- The post-processed trajectory (`--output_traj`, default `trajectory.xtc`) and structure (`--output_str`, default `structure.pdb`).
- `config.yml`, `log.out`, and per-step directories for inspection.

## Limitations

- **System-dependent PBC handling.** The recommended whole/cluster/nojump/center/image sequence is not
  appropriate for every system (multiple solutes, membranes, dissociating complexes) — always visually
  inspect the result.
- **Whole reference structure.** `--input_structure` is the PBC reference for the whole/nojump
  steps and the source of the centering atom; if it is itself broken across periodic boundaries, the
  output will be wrong.
- **Solvent/ion detection.** Solvent and ions are recognized by residue/atom name; anything
  not in the GROMACS defaults or added via `--solvents`/`--ions` is treated as solute and not stripped.
- **Fitting.** The `rot+trans` fit onto the solute discards overall
  rotation/translation. A single reference frame is used (not progressive fitting), so very large conformational transitions may not fit
  cleanly.
