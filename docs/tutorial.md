# Tutorial

A hands-on notebook that runs a complete **protein–ligand molecular dynamics** simulation with
the workflows in this repository — in just **three commands**. The system is the T4 lysozyme
L99A/M102Q mutant (PDB [`3HTB`](https://www.rcsb.org/structure/3HTB)) in complex with
2-propylphenol (ligand [`JZ4`](https://www.rcsb.org/ligand/JZ4)), the classic
[GROMACS protein–ligand example](http://www.mdtutorials.com/gmx/complex/index.html).

| Step | Command | What it does |
|------|---------|--------------|
| 1 | `protein_preparation` | Clean, fix and protonate the protein |
| 2 | `ligand_parameterization` | Build a GROMACS topology + coordinates for JZ4 |
| 3 | `md_gromacs` | Setup → energy minimization → NVT/NPT equilibration → production MD → analysis |

The notebook then visualizes the trajectory and plots the analyses (RMSD, radius of gyration,
RMSF, ligand RMSD) produced by the MD workflow.

The tutorial can be run **either on Google Colab (no local install) or on a local Jupyter**.

## Run on Google Colab

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/NBDsoftware/biobb_md_workflows/blob/master/notebooks/notebook_tutorial.ipynb)

Click the badge above and run the first setup cell. It builds the `biobb_md_tutorial` conda
environment from `notebooks/colab_environment.yml` (GROMACS, AmberTools, acpype, the BioBB
packages, this repo's workflows, and the visualization libraries) — nothing to install on your
machine. The conda solve is the slow part; run it once and wait.

## Run locally

Create the environment from `notebooks/local_environment.yml`, then launch Jupyter from it:

```bash
git clone https://github.com/NBDsoftware/biobb_md_workflows.git
cd biobb_md_workflows
conda env create -f notebooks/local_environment.yml
conda activate biobb_md_tutorial
jupyter notebook notebooks/notebook_tutorial.ipynb
```

The notebook lives at
[`notebooks/notebook_tutorial.ipynb`](https://github.com/NBDsoftware/biobb_md_workflows/blob/master/notebooks/notebook_tutorial.ipynb)
in the repository.
