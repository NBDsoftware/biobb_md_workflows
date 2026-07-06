# Ligand parameterization

Generate force-field parameters (topology + coordinates) for the ligands and cofactors in a PDB, for use in GROMACS or AMBER MD simulations.

## Description

For each ligand found in the input PDB, the workflow branches on whether a custom AMBER parameter set is supplied:

- **Custom parameters** (`--custom_parameters` folder with matching `.frcmod` + `.prep` files, e.g. from the [Amber Parameter Database (Manchester)](http://amber.manchester.ac.uk/)): LEaP builds the AMBER topology, optionally converted to GROMACS. The ligand keeps the charge and protonation state of its template.
- **No custom parameters:** the ligand is protonated (`--protonation_tool`), minimized (unless `--skip_min`), and parameterized with GAFF via antechamber/acpype.

The output is a folder of topology + coordinate files: `.gro` + `.itp` for GROMACS, or `.frcmod` + `.prep`/`.lib` for AMBER.

Note the per-format caveat: for **GROMACS**, the `.gro` and `.itp` must agree, so the workflow must be re-run for every new PDB (the coordinates change per system) before merging with the protein topology after `pdb2gmx`. For **AMBER**, `tleap` can reconstruct missing atoms from the `.prep`/`.lib` files, so existing `.frcmod` + `.prep`/`.lib` sets can be reused across PDBs without re-running.

## Usage

```bash
conda activate biobb_md
ligand_parameterization --input_pdb data/complex.pdb --ligands JZ4 --chains A --format gromacs --output output
```

The `config.yml` is auto-generated from the CLI arguments into the output folder. `--restart` resumes from the last completed step. Run `ligand_parameterization --help` for the full option list.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input_pdb` | *required* | Input PDB file containing the ligands to parameterize. |
| `--ligands` | all ligands | Ligand names (from the PDB) to parameterize. |
| `--chains` | `A` | Chain IDs to extract ligands from. |
| `--model` | `0` | Model number to extract ligands from. |
| `--format` | `gromacs` | Output topology format: `gromacs` or `amber`. |
| `--charges` | guessed | Per-ligand charges as `name:charge` (e.g. `JZ4:-2 FLP:1`). GAFF/acpype path only; by default acpype guesses from the protonation state. |
| `--protonation_tool` | `ambertools` | Ligand protonation tool: `ambertools`, `obabel`, or `none`. |
| `--skip_min` | `False` | Skip the minimization step. |
| `--custom_parameters` | `None` | Folder with `.frcmod` + `.prep` files (named by ligand); triggers the LEaP path. |
| `--forcefields` | `protein.ff14SB DNA.bsc1 gaff` | Force fields for LEaP (only used with custom parameters). See `$AMBERHOME/dat/leap/cmd/`. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output_top_path` | 'topologies' inside `output` | Output folder for the ligand topologies/coordinates. |
| `--output` | 'output' | Output directory. |

## Outputs

- A `topologies/` folder with one topology + coordinate pair per ligand: `<LIG>.gro` + `<LIG>.itp` (GROMACS) or `<LIG>.frcmod` + `<LIG>.prep`/`.lib` (AMBER).
- Per-ligand working directories and `log.out` for inspection.
