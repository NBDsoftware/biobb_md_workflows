# Protein preparation

Fix PDB defects, choose protonation states for titratable residues, and prepare a protein structure for MD simulation.

## Description

Starting from a raw PDB, the workflow:

1. **Extracts** the requested chains from the input structure.
2. **Fixes PDB defects:**
   - **Alternative locations** — keep the highest-occupancy alternative location.
   - **Mutations** — apply the mutations given in `--mutation_list` (e.g. `A:Arg220Ala`).
   - **Sequence (FASTA)** — obtain the canonical sequence for backbone modeling (through the `--pdb_code` or `SEQRES` records).
   - **Missing backbone atoms** — model missing atoms with `biobb_structure_checking` + Modeller for missing loops.
   - **Missing side-chain atoms** — model missing atoms with `biobb_structure_checking`. 
   - **Amide flips** — relieve clashing ASP/GLU amide groups.
   - **Chirality** — fix stereochemical errors in side chains.
   - **Disulfide bonds** — mark CYS→CYX/CYS2 with a distance criterion.
3. **Assigns protonation states:** assign protonation states to protein residues using propka heurisitc method and a pH (unless `--keep_hs` is selected)
   - Removes existing hydrogens
   - Estimates pKa of titratable residues with `propka`
   - Optimizes histidine protonation state with AmberTools `reduce`
   - Renames residues according to their protonation state and `--output_format` (`amber` or `gromacs`)
   - Specific histidine protonation states can be provided through `--his`.

## Usage

```bash
conda activate biobb_md
protein_preparation --input_pdb data/1r9o.pdb --ph 7 --cap_ter --output_format gromacs --output output
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `protein_preparation --help` for the full option list.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input_pdb` | *required* | Input PDB file. |
| `--pdb_code` | `None` | PDB code used to fetch the canonical FASTA sequence. If omitted, the PDB HEADER is used. |
| `--pdb_chains` | all chains | Chains to extract and fix (e.g. `A B C`). |
| `--mutation_list` | `None` | Mutations to introduce (e.g. `A:Arg220Ala B:Ser221Gly`). |
| `--skip_bc_fix` | `False` | Skip modeling of missing backbone atoms. |
| `--modeller_key` | `None` | Modeller license key ([register here](https://salilab.org/modeller/registration.html)); required for loop modeling. Only for academic use. |
| `--cap_ter` | `False` | Cap termini with ACE/NME residues. |
| `--skip_sc_fix` | `False` | Skip modeling of missing side-chain atoms. |
| `--skip_ss_bonds` | `False` | Skip addition of disulfide bonds. |
| `--skip_amides_flip` | `False` | Skip flipping of clashing ASP/GLU amide groups. |
| `--ph` | `7.0` | pH used with the propka pKa estimation to set protonation states. |
| `--his` | `None` | Force histidine states (delta `0`, epsilon `1`, protonated `2`, heme-bound `3`), e.g. `0 1 1`. Overrides the pKa/pH estimation for HIS. |
| `--keep_hs` | `False` | Keep input hydrogen atoms instead of discarding them. |
| `--output_format` | `amber` | Residue-naming convention: `amber` or `gromacs`. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `output` | Output directory. |

## Recommendations

### General

:::{admonition} 🚧 To be written
:class: caution
Make sure you select the chains you want with `--pdb_chains`; how to specify `--mutation_list`
when to use `--cap_ter`; when and how to force
histidine states with `--his`; Explain how accurate can propka be expected to be. Avoid providing a modeller key if the use is not academic, in this case modelling of missing loop is not yet supported. Beware that when working with PDB files, especially those from NMR structures, it’s often more convenient to disregard hydrogen (H) atoms. This is because, if hydrogen atoms are included, they need to be named exactly as specified by the force fields in GROMACS. Since naming conventions for hydrogen atoms vary, this can sometimes cause complications. If you want to keep the original coordinates of hydrogen atoms you'll probably need to rename them.
:::

### Preserving input protonation states with `--keep_hs`

:::{note}
When `--keep_hs` is set, the original hydrogen atoms are kept and the automatic protonation
state determination (propka + reduce + titration) is **skipped**. The titratable protonation
residue names marked in the input PDB are therefore respected and written unchanged (converted
only to the requested `--output_format`): `HID`/`HIE`/`HIP`, `ASH`, `GLH`, `LYN` for AMBER, and
`HISD`/`HISE`/`HISH`, `ASPH`, `GLUH`, `LYSN`, `ARGN` for GROMACS. Residues left with a standard
name (e.g. plain `HIS`) are not reassigned. Downstream, `md_gromacs`'s `pdb2gmx` reads these
residue names to set the protonation state.
:::

### Disulfide bonds

Same thing here with `--skip_ss_bonds`. To specify them manually instead, mark the
bonded cysteines in the input PDB (`CYS2` for GROMACS — safer for old versions — or `CYX` for AMBER)
and run with `--skip_ss_bonds`: automatic detection is turned off and the marked cysteines are
preserved in the output. With `--skip_ss_bonds` and no marked cysteines, no disulfide bonds are formed.


## Output

Checklist:

How to verify the result: check the chosen protonation states in the output PDB,
confirm modelled loops/side chains look reasonable. Confirm the SS bonds are created correctly (CYS2 for GROMACS and CYX for AMBER)

Files:

Written into `--output`, organized by steps:

- The prepared PDB structure (named after the input), ready for `md_gromacs`: `<output_folder>/<pdb_name>.pdb`
- `config.yml` — the generated configuration used for the run.
- `log.out` and per-step directories (`step*_*`) for inspection.

## Limitations

- **Limited protonation accuracy.** The pKa estimation uses propka, an heuristic algorithm that determines protonation states based on neighboring residues.
- **Modeling of gaps in proteins requires Modeller**. Modeller is **academic-use only** and needs a
  license key. See [Installation](installation.md).
- **No support for custom force fields**.
- **Non-standard residues** are not supported.
- **Fixed-column PDB rewriting.** Residue renaming for protonation/disulfide
  states relies on fixed-column PDB slicing. A malformed or non-standard input
  PDB can silently produce wrong output here.
- **Disulfide bonds.** By default the pipeline detects disulfide bonds automatically with
  `biobb_structure_checking` using a distance criterion.
