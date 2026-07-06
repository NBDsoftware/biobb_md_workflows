# Protein preparation

Fix PDB defects, choose protonation states for titratable residues, and prepare a protein structure for MD simulation.

## Description

Starting from a raw PDB, the workflow:

1. **Extracts** the requested chains from the input structure.
2. **Fixes PDB defects:**
   - **Alternative locations** — keep a chosen altloc, or the highest-occupancy one by default.
   - **Mutations** — apply the mutations given in `--mutation_list` (e.g. `A:Arg220Ala`).
   - **Sequence (FASTA)** — obtain the canonical sequence (from the PDB, `--pdb_code`, or `SEQRES` records) to guide backbone modeling.
   - **Missing backbone atoms** — model with `biobb_structure_checking` + Modeller (loop modeling needs `--modeller_key`). Skip with `--skip_bc_fix`.
   - **Missing side-chain atoms** — model with `biobb_structure_checking`. Skip with `--skip_sc_fix`.
   - **Amide flips** — relieve clashing ASP/GLU amide groups. Skip with `--skip_amides_flip`.
   - **Chirality** — fix stereochemical errors in side chains.
   - **Disulfide bonds** — mark CYS→CYX by distance criterion. Skip with `--skip_ss_bonds`.
3. **Assigns protonation states:** removes existing hydrogens, estimates pKa of titratable residues with `propka`, optimizes histidine hydrogen placement with AmberTools `reduce`, then selects and renames residues for the requested `--output_format` (`amber` or `gromacs`). Histidine states can be forced with `--his`.

## Usage

```bash
conda activate biobb_md
protein_preparation --input_pdb data/1r9o.pdb --ph 7 --cap_ter --output_format gromacs --output output
```

The `config.yml` is auto-generated from the CLI arguments into `--output` for reproducibility. `--restart` resumes from the last completed step when re-run against the same output folder. Run `protein_preparation --help` for the full option list.

## Options

| Flag | Default | Description |
|------|---------|-------------|
| `--input_pdb` | *required* | Input PDB file. |
| `--pdb_code` | `None` | PDB code used to fetch the canonical FASTA sequence. If omitted, the PDB HEADER is used. |
| `--pdb_chains` | all chains | Chains to extract and fix (e.g. `A B C`). |
| `--mutation_list` | `None` | Mutations to introduce (e.g. `A:Arg220Ala B:Ser221Gly`). |
| `--skip_bc_fix` | `False` | Skip modeling of missing backbone atoms. |
| `--modeller_key` | `None` | Modeller license key ([register here](https://salilab.org/modeller/registration.html)); required for loop modeling. Only for academic use. |
| `--cap_ter` | `False` | Cap termini with ACE/NME residues, preserving existing atoms. |
| `--skip_sc_fix` | `False` | Skip modeling of missing side-chain atoms. |
| `--skip_ss_bonds` | `False` | Skip addition of disulfide bonds. |
| `--skip_amides_flip` | `False` | Skip flipping of clashing ASP/GLU amide groups. |
| `--ph` | `7.0` | pH used with the propka pKa estimation to set protonation states. |
| `--his` | `None` | Force histidine states (delta `0`, epsilon `1`, protonated `2`, heme-bound `3`), e.g. `0 1 1`. Overrides the pKa/pH estimation for HIS. |
| `--keep_hs` | `False` | Keep input hydrogen atoms instead of discarding them. |
| `--output_format` | `amber` | Residue-naming convention: `amber` or `gromacs`. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `output` | Output directory. |

## Outputs

Written into `--output`:

- The prepared PDB structure (named after the input, e.g. `1r9o.pdb`), ready for `md_gromacs`.
- `config.yml` — the generated configuration used for the run.
- `log.out` and per-step directories (`step*_*`) for inspection.
