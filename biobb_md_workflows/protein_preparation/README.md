# Protein preparation

Fix PDB defects, choose protonation states for titratable residues, and prepare a protein structure to run MD simulations.

## Description

Starting from a raw PDB, the workflow:

1. **Extracts** the requested chains from the input structure.
2. **Fixes PDB defects:**
   - **Alternative locations** — keep the highest-occupancy alternative location.
   - **Mutations** — apply the mutations given in `--mutation_list`.
   - **Sequence (FASTA)** — obtain the canonical sequence for backbone modeling through the `--pdb_code` or `SEQRES` records.
   - **Missing backbone atoms** — model missing atoms with `biobb_structure_checking` + Modeller for missing loops.
   - **Missing side-chain atoms** — model missing atoms with `biobb_structure_checking`. 
   - **Amide flips** — relieve clashing ASP/GLU amide groups due to a wrong fit of the electronic density.
   - **Chirality** — fix stereochemical errors in side chains.
   - **Disulfide bonds** — mark CYS→CYX/CYS2 with a distance criterion.
3. **Assigns protonation states:** assign protonation states to protein residues using propka heurisitc method and a pH, unless `--keep_hs` is selected.
   - Removes existing hydrogens
   - Estimates pKa of titratable residues with `propka`
   - Optimizes histidine protonation state with AmberTools `reduce`
   - Renames residues according to their protonation state and `--output_format`
   - Specific histidine protonation states can be provided through `--his`.

## Usage

```bash
conda activate biobb_md
protein_preparation --input_pdb data/1r9o.pdb --ph 7 --cap_ter --output_format gromacs --output output
```

Run `protein_preparation --help` for the full option list.

## Options

| Flag&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Default | Description |
|------|---------|-------------|
| `--input_pdb` | *required* | Input PDB file. |
| `--pdb_code` | `None` | PDB code used to fetch the canonical FASTA sequence. If omitted, the PDB HEADER is used. |
| `--pdb_chains` | all chains | Chains to extract and fix (e.g. `A B C`). |
| `--mutation_list` | `None` | Mutations to introduce in the  `Chain:WtResnumMut` format (e.g. `A:Arg220Ala B:Ser221Gly`). |
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

### Structure fixing

Fix only the chains you intend to simulate: select them with `--pdb_chains`;
everything else — extra copies, crystallographic partners, waters, ligands — is discarded at extraction.
Fixing runs *before* protonation, a complete structure makes the downstream pKa estimate more meaningful.

- **Missing atoms vs loops.** Missing side-chain and backbone *atoms* are modelled without a Modeller
  license; modelling missing *loops* (gaps in the chain) needs a `--modeller_key`. The
  key is **academic-use only**. Skip modelling entirely with `--skip_bc_fix` / `--skip_sc_fix`.
- **Mutations (`--mutation_list`).** Space-separated mutations,
  e.g. `A:Arg220Ala B:Ser221Gly`.
- **Capping termini (`--cap_ter`).** Add ACE/NME caps when the input is a fragment or truncated construct,
  so the exposed backbone ends are not treated as charged N-/C-termini. Skip it when the termini are the
  genuine chain ends.
- **Output format (`--output_format`).** Pick `amber` or `gromacs` to match the downstream tool — it fixes
  the residue-naming convention for protonation/disulfide states (`HID`/`HIE`/`HIP` vs `HISD`/`HISE`/`HISH`, `CYX` vs `CYS2`).

### Protonation states

Titratable residues (ASP, GLU, LYS, ARG, HIS) are protonated for the chosen `--ph` (default 7.0) from a
`propka` pKa estimate — a residue is protonated when the pH is below its predicted pKa. propka is an
**empirical** predictor: it estimates every pKa in a medium-sized protein in seconds and, for the vast
majority of residues, is as reliable as far more expensive Poisson–Boltzmann or FEP methods. Two things to keep in mind:

- **It reads a single static structure** — no conformational sampling — the estimate is only as good as
  the structure. Be wary of buried or strongly-shifted residues in flexible regions, where the empirical model is least accurate.
- **Histidine is the hardest case**. The workflow combines the propka pKa
  with an AmberTools `reduce` H-bond-network optimization to choose the delta/epsilon tautomer, but when
  you know the correct state — from the literature, a metal/heme contact, or inspection — set it explicitly
  with `--his` (delta `0`, epsilon `1`, protonated `2`, heme-bound `3`; one value per HIS in order, e.g.
  `0 1 1`). This overrides the pKa/pH estimate for histidine.

### Keeping input hydrogens (`--keep_hs`)

With `--keep_hs` the input hydrogens are kept and the whole protonation determination (propka + reduce +
titration) is **skipped**: the titratable residue names already marked in the input (`HID`/`HIE`/`HIP`,
`ASH`, `GLH`, `LYN` for AMBER; `HISD`/`HISE`/`HISH`, `ASPH`, `GLUH`, `LYSN`, `ARGN` for GROMACS) are
preserved and only converted to `--output_format`; residues left with a standard name (plain `HIS`) are
not reassigned. Use it when the protonation is already correct (e.g. from an upstream tool). Otherwise,
for X-ray/NMR inputs it is usually easier to **discard** the input hydrogens (the default) and let the
workflow re-place them: GROMACS requires hydrogen names to match the force field exactly and kept hydrogens often have to be renamed by hand to avoid `pdb2gmx` errors downstream.

### Disulfide bonds

Disulfide bonds are detected automatically from a distance criterion and the paired cysteines renamed
(`CYS2` for GROMACS — safer for GROMACS < 2024 — or `CYX` for AMBER). To set them manually instead, mark
the bonded cysteines in the input PDB with `CYS2` or `CYX` and run with `--skip_ss_bonds`: automatic
detection is turned off and the marked cysteines are preserved. With `--skip_ss_bonds` and no marked
cysteines, no disulfide bonds are formed.


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

- **Approximate protonation.** propka is an empirical, single-structure pKa predictor: no conformational
  sampling, typical errors below ~1 pKa unit but larger for buried, strongly-coupled, and histidine
  residues.Verify critical residues manually (see `--his` for histidines).
- **Loop modelling requires Modeller.** Modeller is **academic-use only** and needs a license key
  (see [Installation](https://nbdsoftware.github.io/biobb_md_workflows/installation.html)); without it, missing loops are left unmodelled.
- **Standard protein residues only.** Waters, ions, ligands, cofactors, and modified/non-standard
  residues are removed at extraction — parameterize ligands and cofactors separately with the
  `ligand_parameterization` workflow.
- **Output naming is AMBER or GROMACS only.** Other force-field residue-naming conventions are not
  supported.
- **Fixed-column PDB rewriting.** Residue renaming for protonation/disulfide states relies on fixed-column
  PDB slicing; a malformed or non-standard input PDB can silently produce wrong output here.
