# FAIR4RS & Robustness Assessment

A review of how well `biobb_md_workflows` meets **FAIR4RS** and
data-pipeline-engineering principles, plus a phased roadmap to close the gaps.

## Where the project stands today

Legend: 🟢 real / settled · 🟡 partial/aspirational · 🔴 missing or undermined.

### FAIR4RS

| Principle | Reality in the repo | Status |
|---|---|---|
| **Findable** | Public GitHub repo, git tags (`0.0.1`–`0.0.5`), GitHub Pages docs, Colab tutorial. `__version__` now exposed (via `importlib.metadata`) and logged at the start of every run. `CITATION.cff` added. Zenodo DOI planned. | 🟡 |
| **Accessible** | Repo + tags + `LICENSE` present.  License is **CC BY-NC-SA** — academic-open / commercially-restricted choice. | 🟢 |
| **Interoperable** | Standard formats in/out (PDB/GRO/TOP/XTC), YAML config, thin argparse CLIs, importable workflow functions. | 🟢 |
| **Reusable** | `--help` + docs exist; reproducibility in place (pinned deps + tagged forks). | 🟢 |

### Engineering / pipeline robustness

| Principle | Reality in the repo | Status |
|---|---|---|
| **Reproducibility** | Floor now in place: all forks are pinned (`biobb_gromacs`, `biobb_analysis`) to fork-namespaced tags (`@nbd-5.2.1`); `mdanalysis`, `biopython`, `ambertools`, `acpype`, `openbabel`, `propka` are version-pinned. | 🟡 |
| **Provenance** | `config.yml` + absolute input path are written per run; `log.out` and per-step logs exist. Package version is now logged in the `log.out` header. Full CLI command, and input-file checksums are still recorded nowhere (full `run_manifest.json` remains future work). | 🟡 |
| **Environment portability** | Conda env files exist; deps + forks now pinned. No container yet. | 🟡 |
| **Modularity** | Genuinely modular at the BioBB-step level. `create_config_file`/config scaffolding is copy-pasted across all four workflows (only `common.to_yaml` is shared). **De-duplication is deferred — delicate** (each workflow's config differs subtly). | 🟡 |
| **Validation** | Step outputs are validated **internally by the biobbs** — each building block checks its own outputs, so the workflow layer does not need to re-assert file existence between steps. Remaining gaps are at the CLI boundary: no argparse `choices` for enums, no `--ph` bounds, and several bad-arg paths log an error then silently `return`. Only `md_gromacs` has a real `check_inputs()`. | 🟡 |


## Roadmap (cheapest → most involved)

Effort tags are rough: **S** = hours · **M** = a day or two · **L** = a week+.

### Phase 1 — Stop reporting false success *(S, do first)*
- Gate the final `"Execution successful"` message on real completion; on a caught
  failure, log the error and `sys.exit(1)`.
- (Step-output existence between steps is already handled internally by the biobbs;
  where their failure message is unclear, surface a clearer workflow-level message.)

### Phase 2 — Provenance *(M, high value)*

- Future: one shared helper writing `output/run_manifest.json` — git commit (if
  resolvable), full `sys.argv`, SHA-256 of each input file, UTC timestamp, resolved
  `config.yml` path, and a `conda env export` / `pip freeze` snapshot.

### Phase 4 — Validation hardening *(M)*
- argparse `choices=` for the `format`/`output_format`/protonation-tool enums; `--ph`
  bounds `[0, 14]`; convert "log-error-then-`return`" silent exits into `raise`.
- Generalize the `md_gromacs.check_inputs()` pattern to the other three workflows.

### Phase 5 — Testing & automation *(L)*
- Minimal GitHub Actions CI: `ruff` lint + `pip install .` import-smoke + package build.
- `pre-commit` with `ruff` for local guardrails.
- One tiny end-to-end smoke test on a small system (reuse `tests/*/input.yml`);
  parameterize the hard-coded env path in `run.sl`.
- Dependency automation (Dependabot/Renovate) watching the two forks.
- De-duplicate the copy-pasted `create_config_file`/config scaffolding into `common/`
  (deferred — delicate).

### Phase 6 — FAIR4RS metadata & portability *(L, stretch)*

- Mint a **Zenodo DOI** for easy academic citation: (1) an NBDsoftware org admin
  authorizes Zenodo for the org, (2) enable the repo in Zenodo, (3) cut a GitHub Release
  → DOI auto-minted, (4) paste the concept-DOI into `CITATION.cff` and a README badge.
- `CHANGELOG.md`.
- License: **settled** — Explore other options more appropiate for code.
- Build an Apptainer/Singularity image per release for one-command portable execution.
