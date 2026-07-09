# Ligand parameterization

Generate force-field parameters (topology + coordinates) for the ligands and cofactors in a PDB, for use in GROMACS or AMBER MD simulations.

## Description

**Input modes** (resolved at runtime. Can be combined if there are different ligands. E.g. a cofactor with a parameter set and a ligand):

- `--input_pdb` + `--custom_parameters` — Build the topology using the parameter set ( `.frcmod` + `.prep` files). The ligands present in custom_parameters, keep the charge and protonation state of their template. The name of the files should be the name of the ligand they refer to. 

- `--input_pdb` — when no custom parameter set is available, the ligand is protonated with (`--protonation_tool`), minimized (unless `--skip_min`), and parameterized with GAFF via antechamber and acpype.

Note: for **GROMACS**, the coordinates `.gro` and topology `.itp` must agree, so the workflow must be re-run for every new PDB (the coordinates change for every system even if the ligand is the same). For **AMBER**, `tleap` can reconstruct missing atoms from the `.prep`/`.lib` files, so `.frcmod` + `.prep`/`.lib` sets can be reused across PDBs without re-running.

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


## Recommendations

### Parameters

Prefer a curated custom AMBER parameter set when one exists: pass it with `--custom_parameters`
(`<LIG>.frcmod` + `<LIG>.prep`, named after the ligand). These are literature-validated and keep the
template's charge and protonation, so they are more reliable than automatically generated parameters.
Common ligands and cofactors are available in the [Amber Parameter Database, Manchester](http://amber.manchester.ac.uk/).

Otherwise the GAFF path parameterizes the ligand automatically (antechamber atom types + AM1-BCC
charges) — adequate for most drug-like organic molecules. Both paths can run in the same job (e.g. a
cofactor with custom parameters alongside a GAFF-parameterized ligand).

### Charges and protonation

Get protonation right first: it fixes the ligand's net charge, and GAFF/acpype derives the AM1-BCC
partial charges from the protonated 3D structure. Choose the tool with `--protonation_tool`:

- `ambertools` (reduce, default) — completes hydrogens using the wwPDB HET connectivity dictionary and
  optimizes OH/SH/NH₃ and Asn/Gln/His orientations. Models a **neutral environment: it ignores pH and
  does not ionize acids/bases** (e.g. carboxylates stay protonated). Reliable for standard PDB ligands;
  may skip or misplace hydrogens for novel molecules absent from the dictionary or with nonstandard atom names.
- `obabel` — perceives bonds from the 3D coordinates and protonates for **pH 7.4** using tabulated
  per-group pKa rules. More robust for arbitrary/novel small molecules and assigns a physiological charge state.
- `none` — keep the input protonation unchanged; use when the ligand is already correctly protonated
  (e.g. from an upstream preparation tool).

acpype guesses the net charge from the protonation state; check it is reasonable and override with
`--charges LIG:<q>` (GAFF path only) if it is wrong.

## Output

Checklist:

For GROMACS the
`.gro`/`.itp` pair is system-specific and must be regenerated per PDB, whereas
AMBER `.frcmod`/`.prep` sets are reusable. Note how these files are consumed by
`md_gromacs` via `--ligands_folder`.
Check the guessed charge of the ligand is reasonable when using GAFF/acpype to parameterize

Files:

- A `topologies/` folder with one topology + coordinate pair per ligand: `<LIG>.gro` + `<LIG>.itp` (GROMACS) or `<LIG>.frcmod` + `<LIG>.prep`/`.lib` (AMBER).
- Per-ligand working directories and `log.out` for inspection.


## Limitations

- **Approximate protonation.** No rigorous pKa calculation is performed: `ambertools`/reduce assumes a
  neutral environment (won't ionize groups) and `obabel` applies fixed tabulated pKa rules at pH 7.4.
  Inspect ionizable groups and tautomers manually.
- **GAFF / AM1-BCC quality.** The automatic path (antechamber + acpype) uses general atom types and
  semi-empirical charges. It covers most drug-like organic molecules but is unreliable for metals and
  metal coordination, uncommon elements, and unusual chemistries — prefer curated custom parameters there.
- **Input atom names.** Nonstandard or duplicate atom names in the input PDB can break reduce and later steps.
- **GROMACS outputs are system-specific.** `.gro`/`.itp` pairs must be regenerated per PDB; AMBER
  `.frcmod`/`.prep` sets are reusable.


