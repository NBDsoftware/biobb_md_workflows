# GOdMD simulation workflow

Run a GO-like potential discrete MD simulation, a coarse grain method to trace approximately large conformational transitions.

## Description

Computes an approximate **conformational transition** between two known structural states (origin and target) of the same protein, using **GOdMD**: discrete MD with GO-like biasing potentials combining essential dynamics and Maxwell-Demon sampling (Sfriso et al., *Bioinformatics* 2013). Each structure can be a local PDB file or fetched by PDB code from RCSB; chain selection strips other chains, and waters/ligands are always removed since GOdMD works on the protein backbone/side-chain representation only.

Pipeline:

1. Fetch (if given by PDB code) and extract the chain of interest from the origin structure.
2. Fetch and extract the chain of interest from the target structure.
3. Build a one-to-one origin↔target residue mapping from a pairwise sequence alignment (EMBOSS `water`).
4. Run GOdMD (`discrete` binary) to compute the transition trajectory.
5. Convert the resulting `.mdcrd` trajectory to `.xtc` for downstream analysis/visualization.

## Usage

```bash
conda activate biobb_md

# From PDB codes and chains
godmd_simulation --origin_pdb_code 1ake --origin_chain A --target_pdb_code 4ake --target_chain A --output output

# From PDB files
godmd_simulation --origin_pdb 1ake.pdb --origin_chain A --target_pdb target.pdb --target_chain A --output output
```

The `config.yml` is auto-generated from the CLI arguments into `--output`. `--restart` resumes from the last completed step when re-run against the same output folder. Run `godmd_simulation --help` for the full option list.

## Options

### Inputs

| Flag&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Default | Description |
|------|---------|-------------|
| `--origin_pdb` | — | Local PDB file for the origin structure. Mutually exclusive with `--origin_pdb_code`. |
| `--origin_pdb_code` | — | PDB code to fetch the origin structure from RCSB. Mutually exclusive with `--origin_pdb`. |
| `--origin_chain` | `None` | Chain of interest in the origin structure. If not given, all chains are kept. |
| `--target_pdb` | — | Local PDB file for the target structure. Mutually exclusive with `--target_pdb_code`. |
| `--target_pdb_code` | — | PDB code to fetch the target structure from RCSB. Mutually exclusive with `--target_pdb`. |
| `--target_chain` | `None` | Chain of interest in the target structure. If not given, all chains are kept. |

One of `--origin_pdb`/`--origin_pdb_code` and one of `--target_pdb`/`--target_pdb_code` is required. Origin and target can independently be a local file or a fetched code.

### Parameters

| Flag&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; | Default | Description |
|------|---------|-------------|
| `--temp` | `300` | GOdMD simulation temperature (K). |
| `--seed` | `2525` | Random number generator seed. |
| `--nbloc` | `10000` | Number of simulation blocks. Total simulated time = `nbloc * tsnap`. |
| `--tsnap` | `500` | Output structure interval, in fs. |
| `--rmsd_tol` | `1.5` | Minimum acceptable RMSD (Å) to consider the transition converged. |
| `--gapopen` | `12.0` | EMBOSS `water` pairwise alignment gap-open penalty. |
| `--gapextend` | `2.0` | EMBOSS `water` pairwise alignment gap-extend penalty. |
| `--debug` | `False` | Keep intermediate files. |
| `--restart` | `False` | Restart from the last completed step. |
| `--output` | `output` | Output directory. |

## Outputs

Written into `--output`:

- `transition.xtc` — the computed conformational transition trajectory (converted from GOdMD's native `.mdcrd`).
- `transition.pdb` — reference topology/structure matching `transition.xtc` (GOdMD's coarse-grain atom set, not the full-atom input).
- `godmd.log`, `godmd_energy.out` — GOdMD run log and energy trace.
- `manifest.yaml`, `config.yml`, `log.out`, and per-step directories for inspection.

## Limitations

- **Coarse grain approximation.** The simulation uses discrete Molecular Dynamics with GO-like potentials representing the end states of the transition. What happens in between is not guaranteed to be physically meaningful, and the simulation is not a substitute for a full atomistic MD simulation. It only provides a rough guess of how the transition may occur, and is useful for generating starting points for further atomistic simulations.
- **Requires a residue-level mapping between origin and target.** The pairwise alignment step (EMBOSS `water`) needs origin and target sequences to be alignable; large insertions/deletions or very low sequence identity between the two states can produce a poor or unusable mapping.
- **Waters and ligands are always stripped.** GOdMD only operates on the protein; cofactors present in one state and absent in the other are not modeled.

## References

**Exploration of conformational transition pathways from coarse-grained simulations.**<br>
*Sfriso P, Hospital A, Emperador A, Orozco M.*<br>
*Bioinformatics, 129(16):1980-6.*<br>
*Available at: https://doi.org/10.1093/bioinformatics/btt324*
