# Changelog

All notable changes to this project are documented here.

## [0.1.0] - 2026-07-23

### ⚠️ Breaking changes
- **ligand_parameterization**: renamed `--custom_parameters` to `--ligand_parameters`
  (update any scripts/configs that used the old flag).
- **Python 3.12 required** (was 3.11), to match Google Colab.
- **md_gromacs**: `.trr` trajectories are no longer written (minimization, NVT, NPT,
  production) — only `.xtc` is produced.

### Added
- **md_gromacs**: per-ligand RMSD (pose stability) — heavy-atom ligand RMSD on the fitted
  trajectory (`nofit`), one `md_rmsd_ligand_<LIG>.xmgr` per ligand.
- **All workflows**: machine-readable `manifest.yaml` output (stable schema for external
  consumers), listing the produced structures/topologies/trajectories/analysis files and
  the stage reached.
- **protein_preparation, ligand_parameterization**: `--debug` flag to keep temporary files.
- **ligand_parameterization**: early warning when a ligand has an odd electron count and no
  charge is set (prevents an acpype/SQM failure due to asking an odd electron count calculation 
  with spin multiplicity = 1).
- **All workflows**: log the package version at startup.
- **Tutorial**: end-to-end protein–ligand MD notebook (PDB `3HTB` + ligand `JZ4`),
  runnable on Google Colab or a local Jupyter, plus a documentation page.

### Fixed
- **md_gromacs**: RMSF/RMSD/Rgyr are now computed on the cleaned (dry, imaged, centered,
  fitted) trajectory, removing PBC-induced artifacts such as spurious RMSF spikes on
  surface residues.
- **md_gromacs**: correct RMSD reference structures for every input mode
  (`input_pdb`, `input_gro_top`, `restart_simulation`).
- **md_gromacs**: restart no longer fails looking for `.trr` files that were never written.
- Default installation is non-editable and installs from a pinned git tag.

### Dependencies
- Pinned the `biobb_gromacs` and `biobb_analysis` NBD forks to the `nbd-5.3.0` and `nbd-5.2.1` tags for
  reproducible installs.
- Pinned `mdanalysis==2.10.0` and `biopython==1.87`; added `pyyaml`.

### Documentation
- New tutorial page and README section (with Colab badge), project logos, a FAIR4RS
  assessment (`docs/fair4rs_assessment.md`), a `CITATION.cff`, and Sphinx warning fixes.

[0.1.0]: https://github.com/NBDsoftware/biobb_md_workflows/releases/tag/0.1.0
