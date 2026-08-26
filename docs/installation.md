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
2. Go to **Actions ▸ Release ▸ Run workflow**, pick the branch (normally `master`), and enter the
   new version (e.g. `0.2.0`, no leading `v`).

The workflow rewrites the version everywhere, renames `[Unreleased]` to the new version, commits,
tags, and publishes the GitHub Release using that CHANGELOG section as the release notes. Preview
what it would do locally with `python scripts/bump_version.py 0.2.0 --dry-run`.
