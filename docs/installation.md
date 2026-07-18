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
