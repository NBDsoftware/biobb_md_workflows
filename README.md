# BioBB workflows for MD simulations with GROMACS

These workflows use BioExcel Building Blocks (BioBB) to integrate different popular biomolecular simulation tools. This software has been developed for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU Horizon Europe [101093290](https://cordis.europa.eu/project/id/101093290)).

This repo covers the protein preparation, ligand parameterization, MD with GROMACS and trajectory post processing workflows.  Virtual-screening workflows (cavity_analysis, vs_autodock) live in [biobb_vs_workflows](https://github.com/NBDsoftware/biobb_vs_workflows).

## Workflows

**protein_preparation**: prepares a protein structure for molecular dynamics simulations with GROMACS. It includes the following steps: structure cleaning, protonation, missing atom reconstruction, and topology generation.

**ligand_parameterization**: prepares a ligand structure for molecular dynamics simulations with GROMACS. It reads in custom ligand parameters (e.g. those from the [Manchester University database](http://amber.manchester.ac.uk/)) and generates the topology and coordinate files for GROMACS. If no custom parameters are provided, it will generate them using the antechamber and GAFF.

**md_gromacs**: runs a molecular dynamics simulation with GROMACS. It includes the following steps: energy minimization, equilibration, production run and trajectory post-processing. The workflow can be run with either a single or multiple replicas.

**traj_postprocessing**: post-processes the raw trajectory of a molecular dynamics simulation with the recomended centering and imaging steps of the GROMACS manual. Note that results should be checked for each trajectory, as the centering and imaging steps may not be appropriate for all systems.

## Installation

Requirements: `git`, `conda`

```bash
git clone https://github.com/NBDsoftware/biobb_md_workflows.git
cd biobb_md_workflows
conda env create -f environment.yml
conda activate biobb_md
```

## Usage

Once installed, each workflow is available as a CLI command:

| Command | Workflow |
|---------|----------|
| `protein_preparation` | [Protein preparation](biobb_md_workflows/protein_preparation/README.md) |
| `ligand_parameterization` | [Ligand parameterization](biobb_md_workflows/ligand_parameterization/README.md) |
| `md_gromacs` | [Molecular dynamics with GROMACS](biobb_md_workflows/md_gromacs/README.md) |
| `traj_postprocessing` | [Trajectory post-processing](biobb_md_workflows/traj_postprocessing/README.md) |

```bash
protein_preparation --help
```

## Licensing

This project is offered under a dual-license model, intended to make the software freely available for academic and non-commercial use while preventing its use for profit.

### 1. Academic and Non-Commercial Use

For academic, research, and other non-commercial purposes, this software is licensed under the **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**.

Under this license, you are free to:
*   **Share** — copy and redistribute the material in any medium or format.
*   **Adapt** — remix, transform, and build upon the material.

As long as you follow the license terms:
*   **Attribution** — You must give appropriate credit.
*   **Non-Commercial** — You may not use the material for commercial purposes.
*   **Share-Alike** — If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

A full copy of the license is available in the [LICENSE](LICENSE) file in this repository.

### 2. Commercial Use

**Use of this software for commercial purposes is not permitted under the CC BY-NC-SA 4.0 license.**

If you wish to use this software in a commercial product, for-profit service, or any other commercial context, you must obtain a separate commercial license.

Please contact **it@nostrumbiodiscovery.com** to inquire about purchasing a commercial license.

![](https://bioexcel.eu/wp-content/uploads/2019/04/Bioexcell_logo_1080px_transp.png "Bioexcel")