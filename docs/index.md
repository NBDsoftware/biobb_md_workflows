# BioBB MD Workflows

Ready-to-use command-line pipelines for molecular dynamics (MD) simulations with
GROMACS, built on top of [BioExcel Building Blocks (BioBB)](https://mmb.irbbarcelona.org/biobb/).
Developed for the [European BioExcel](http://bioexcel.eu/) project, funded by the
European Commission (EU Horizon Europe
[101093290](https://cordis.europa.eu/project/id/101093290)).

This repository covers protein preparation, ligand parameterization, MD with
GROMACS, and trajectory post-processing. Related pipelines live in separate
repositories: virtual screening in
[biobb_vs_workflows](https://github.com/NBDsoftware/biobb_vs_workflows).

## Workflows

| Command | Purpose |
|---------|---------|
| [`protein_preparation`](workflows/protein_preparation.md) | Clean, protonate, and fix a raw PDB structure for MD. |
| [`ligand_parameterization`](workflows/ligand_parameterization.md) | Generate GROMACS/AMBER topology + coordinates for ligands and cofactors. |
| [`md_gromacs`](workflows/md_gromacs.md) | Full GROMACS MD: setup → equilibration → production → post-processing. |
| [`traj_postprocessing`](workflows/traj_postprocessing.md) | Strip solvent, center, image, and fit an existing trajectory. |

```{toctree}
:maxdepth: 2
:hidden:

installation
workflows/protein_preparation
workflows/ligand_parameterization
workflows/md_gromacs
workflows/traj_postprocessing
```

## Licensing

Offered under a dual-license model: free for academic and non-commercial use under
**CC BY-NC-SA 4.0**; a separate commercial license is required for for-profit use
(contact `it@nostrumbiodiscovery.com`). See the `LICENSE` file in the repository.