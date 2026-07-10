# BioBB workflows for MD simulations with GROMACS

These workflows use BioExcel Building Blocks (BioBB) to integrate different popular biomolecular simulation tools. This software has been developed for the [European BioExcel](http://bioexcel.eu/), funded by the European Commission (EU Horizon Europe [101093290](https://cordis.europa.eu/project/id/101093290)).

This repo covers the protein preparation, ligand parameterization, MD with GROMACS and trajectory post processing workflows.  Virtual-screening workflows (cavity_analysis, vs_autodock) live in [biobb_vs_workflows](https://github.com/NBDsoftware/biobb_vs_workflows).

## Tutorial

New here? The [Colab tutorial](notebooks/biobb_md_workflows_colab_tutorial.ipynb) runs a full protein–ligand MD (PDB `3HTB` + ligand `JZ4`) end to end using the workflow commands — no local install required:

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/NBDsoftware/biobb_md_workflows/blob/master/notebooks/biobb_md_workflows_colab_tutorial.ipynb)

## Documentation

You can find workflow descriptions, installation, usage, and known limitations in the [Documentation](https://nbdsoftware.github.io/biobb_md_workflows/).

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