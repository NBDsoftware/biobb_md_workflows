# FAIR4RS & Robustness Assessment

A review of how well `biobb_md_workflows` meets **FAIR4RS** and
data-pipeline-engineering principles, plus a phased roadmap to close the gaps.
Written for a small team with limited engineering bandwidth: recommendations are
ordered so the cheapest, highest-trust changes come first.

## Where the project stands today

Legend: 🟢 real · 🟡 partial/aspirational · 🔴 missing or undermined.

### FAIR4RS

| Principle | Reality in the repo | Status |
|---|---|---|
| **Findable** | Public GitHub repo, git tags (`0.0.1`–`0.0.5`), GitHub Pages docs, Colab tutorial. No DOI/Zenodo, no `CITATION.cff`, no `__version__` in code. Still needs to link to iternal docs properly. | 🟡 |
| **Accessible** | Repo + tags + `LICENSE` present. Bus-factor risk: two org forks pulled from moving branches. License is CC BY-NC-SA — a *content* license, awkward for software (NC also restricts reuse). | 🟡 |
| **Interoperable** | Genuinely the strongest axis: standard formats in/out (PDB/GRO/TOP/XTC), YAML config, thin argparse CLIs, importable workflow functions. | 🟢 |
| **Reusable** | `--help` + docs exist, but reusability is undermined by the reproducibility gaps below and the NC license. Forks float on default branches. | 🟡 |

### Engineering / pipeline robustness

| Principle | Reality in the repo | Status |
|---|---|---|
| **Reproducibility** | Not guaranteed even within one version: `biobb_gromacs` & `biobb_analysis` are `git+https` with no commit/tag (default-branch HEAD, a moving target); `mdanalysis`, `biopython`, `gromacs`, `ambertools`, `acpype` unpinned; no lockfile. | 🔴 |
| **Provenance** | `config.yml` + absolute input path are written per run; `log.out` and per-step logs exist. But software version, git commit, full CLI command, and input-file checksums are recorded nowhere. | 🟡 |
| **Environment portability** | Conda env files exist, but with unpinned + moving git deps; no container; Python is `==3.12` in the package yet `3.11` in the docs env. | 🟡 |
| **Modularity** | Genuinely modular at the BioBB-step level, but the payoff (isolated testing) is unrealized — no pytest, no CI tests. `create_config_file`/config scaffolding is copy-pasted across all four workflows (only `common.to_yaml` is shared). | 🟡 |
| **Validation** | Inconsistent: only `md_gromacs` has a real `check_inputs()`; no argparse `choices`, no `--ph` bounds, no check that a step's output exists before the next consumes it. Several bad-arg paths log an error then silently `return`. | 🔴 |

### Principles worth adding

- **Honest failure reporting / observability (🔴 — highest-priority defect).** On a
  caught mid-pipeline failure, `md_gromacs` and `traj_postprocessing` still log
  `"Execution successful"` and exit 0; `traj_postprocessing` then renames outputs
  that may not exist. Core setup/equilibration/production and *all* of ligand/protein
  prep run with no error handling at all. Silent success-on-failure is the worst mode
  to debug — it defeats validation no matter how many input checks are added.
- **Bus factor / dependency risk** — the concrete form of the "Reusable" goal; name it explicitly.
- **Testing & CI** — the only CI today deploys the docs. No lint, build check, or smoke test.

## Roadmap (cheapest → most involved)

Effort tags are rough: **S** = hours · **M** = a day or two · **L** = a week+.

### Phase 1 — Stop reporting false success *(S, do first)*
- Gate the final `"Execution successful"` message on real completion; on a caught
  failure, log the error and `sys.exit(1)`.
- Before each step consumes a prior step's output, assert the expected file exists and fail loudly with the step name. Check if this is already taken into account by the biobbs. If it is, make sure the message is understandable. 

### Phase 2 — Reproducibility floor *(M)*
- Pin the two forks to a tag or commit SHA everywhere (`pyproject.toml`, all
  `environment.yml`, `notebooks/*_environment.yml`). Tag the forks so the pins are readable.
- Pin `mdanalysis`, `biopython`, `gromacs`, `ambertools`, `acpype`, `openbabel`.
- Add a lockfile (`conda-lock`) generated per release; document `conda-lock install`.
- Reconcile Python 3.11 (docs) vs 3.12 (package).

### Phase 3 — Provenance manifest *(M, high value)*
- One shared helper writing `output/run_manifest.json` from every workflow: package
  version, git commit (if resolvable), full `sys.argv`, SHA-256 of each input file, UTC timestamp, resolved `config.yml` path, and a `conda env export` / `pip freeze` snapshot.
- Add `__version__` (via `importlib.metadata`); log version + commit in the `log.out` header.

### Phase 4 — Validation hardening *(M)*
- argparse `choices=` for the `format`/`output_format`/protonation-tool enums; `--ph` bounds `[0, 14]`; convert "log-error-then-`return`" silent exits into `raise`.
- Generalize the `md_gromacs.check_inputs()` pattern to the other three workflows.

### Phase 5 — Testing & automation *(L)*
- Minimal GitHub Actions CI: `ruff` lint + `pip install .` import-smoke + package build.
- `pre-commit` with `ruff` for local guardrails.
- One tiny end-to-end smoke test on a small system (reuse `tests/*/input.yml`);
  parameterize the hard-coded env path in `run.sl`.
- Dependency automation (Dependabot/Renovate) watching the two forks.
- De-duplicate the copy-pasted `create_config_file`/config scaffolding into `common/`.

### Phase 6 — FAIR4RS metadata & portability *(L, stretch)*
- `CITATION.cff` + `CHANGELOG.md`; mint a Zenodo DOI on each release (a persistent
  identifier is what makes the software truly Findable).
- License review *(business decision)*: CC BY-NC-SA is a content license — decide
  between a properly-worded source-available/commercial dual model or an OSI license. Flag to management; do not change unilaterally.
- Build an Apptainer/Singularity image per release for one-command portable execution.
