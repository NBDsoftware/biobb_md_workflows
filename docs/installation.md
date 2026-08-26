# Installation

Requirements: `git`, `conda`. Creates the `biobb_md` conda environment and install the package:

```bash
git clone https://github.com/NBDsoftware/biobb_md_workflows.git
cd biobb_md_workflows
export KEY_MODELLER="HERE YOUR MODELLER KEY"   # only for academic use
conda env create -f environment.yml
conda activate biobb_md
```

This exposes the four workflow commands (`protein_preparation`, `ligand_parameterization`,
`md_gromacs`, `traj_postprocessing`).

## Releasing (maintainers)

Every version string in the repo (`pyproject.toml`, the three conda env files, `CITATION.cff`,
`docs/conf.py`) is generated — do not edit them by hand. Instead:

1. Describe the changes under `## [Unreleased]` in `CHANGELOG.md` as you work.
2. **Actions ▸ Release (1/2) prepare ▸ Run workflow**, branch `master`, version e.g. `0.2.0`
   (no leading `v`). This rewrites the version everywhere, renames `[Unreleased]` to the new
   version, and pushes a `release/0.2.0` branch. The job summary links you to the PR.
3. Open that PR (keep the title `release: 0.2.0`) and merge it.
4. **Release (2/2) publish** then fires automatically: it tags the release commit and publishes
   the GitHub Release, using that CHANGELOG section as the notes.

It is split in two because `master` requires pull requests, so the version change cannot be pushed
to it directly.
