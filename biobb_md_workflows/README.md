# BioBB MD Workflows

BioBB workflows are ready-to-use pipelines built using BioExcel Building Blocks (BioBB) to perform common tasks in biomolecular simulations and modeling. Each workflow is a CLI entrypoint of the `biobb_md_workflows` package.

Different workflows can be combined to create more complex pipelines. This repo covers the GROMACS/AmberTools-based workflows; virtual-screening workflows (`cavity_analysis`, `vs_autodock`) live in [biobb_vs_workflows](https://github.com/NBDsoftware/biobb_vs_workflows).

## Single protein pipeline

![alt text](../img/Single_protein_scheme.png?raw=true)

## Workflows

**protein_preparation**: prepares a protein structure from a PDB file for further simulations or modeling. Includes adding missing residues/atoms and selecting the protonation state of residues at a given pH.

**ligand_parameterization**: parameterizes a small molecule using Antechamber and the GAFF force field. Outputs GROMACS or AMBER topology files.

**md_gromacs**: prepares and launches a GROMACS MD simulation starting from a prepared PDB structure file.

**traj_postprocessing**: post-processes a GROMACS MD trajectory: strips solvent, centers, images and fits.

## Installation

Requirements: `git`, `conda`

```bash
git clone https://github.com/NBDsoftware/biobb_md_workflows.git
cd biobb_md_workflows
```

> To use MODELLER in `protein_preparation`, export your key before creating the environment:
> ```bash
> export KEY_MODELLER="YOUR_MODELLER_KEY"
> ```

```bash
conda env create -f environment.yml
conda activate biobb_md
pip install -e .
```

## Usage

Once installed, each workflow is available as a CLI command:

| Command | Workflow |
|---------|----------|
| `md_gromacs` | md_gromacs |
| `ligand_parameterization` | ligand_parameterization |
| `protein_preparation` | protein_preparation |
| `traj_postprocessing` | traj_postprocessing |

```bash
md_gromacs --help
protein_preparation --input_pdb protein.pdb --output output/
```
