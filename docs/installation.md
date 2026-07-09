# Installation

Requirements: `git`, `conda`.

```bash
git clone https://github.com/NBDsoftware/biobb_md_workflows.git
cd biobb_md_workflows
export KEY_MODELLER="HERE YOUR MODELLER KEY"   # only for academic use
conda env create -f environment.yml
conda activate biobb_md
pip install .                                  # Use pip install -e . to develop
```

This creates the `biobb_md` conda environment and installs the package, exposing
the four workflow commands (`protein_preparation`, `ligand_parameterization`,
`md_gromacs`, `traj_postprocessing`).
