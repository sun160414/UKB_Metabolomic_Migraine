# Verification of study safeguards and Figure 1D completeness

Verified on 2026-10-29 (UTC). All inputs used here were artificial. No individual-level UK Biobank data or original study results were used or uploaded.

## Tested source and scope

The baseline is repository commit [`b83921d2e54e65fba8fbfe309ae5fd4ac5b4cc0c`](https://github.com/sun160414/UKB_Metabolomic_Migraine/commit/b83921d2e54e65fba8fbfe309ae5fd4ac5b4cc0c). The revised files are identified by [source_fingerprints.csv](source_fingerprints.csv), including MD5, SHA-256 and Git blob hashes. This avoids a self-referential commit identifier inside the record. [source_md5.csv](source_md5.csv) is the runner-generated fingerprint table from the final run; [baseline_source_md5.csv](baseline_source_md5.csv) identifies the baseline demo scripts.

The root study scripts were checked with artificial regression fixtures for explicit sex coding, preservation of numeric factor labels, exclusion of prevalent cases, strict incident-status validation and supported LOESS predictions. The full original UK Biobank analysis was **not** rerun. Its participant data and locally prepared inputs were unavailable for this verification.

The standalone demonstration was run end-to-end for the baseline and twice for the final revised source. Each run used 5,000 artificial participants, 20 metabolites, seed `20260801` and the unchanged primary configuration. All runs used the same R/package environment.

## Results

| Check | Result |
|---|---|
| Study-script safeguards | 26/26 passed |
| Figure 1D completeness regression checks | 9/9 passed |
| Existing balance regression checks | 53/53 passed |
| Pipeline execution checks | 8/8 passed in each final revised run |
| Baseline/repeat comparisons, independent completeness calculation and PNG checks | 268/268 passed |
| Figure decoding | All 10 baseline PNGs and all 11 PNGs in each final revised run fully decoded |

The revised default run contained 5,000 eligible participants, 322 incident cases, 322 matched cases and 3,220 matched controls. Twenty metabolites were tested; 16 were selected and clustered, with no selected curves excluded from clustering. Figure 1D reports **20/20 metabolites complete in every retained matched set (100%)**.

The original numerical outputs were identical after parsing the CSV/RDS data. This covers generated inputs, Cox results, matching allocations, residuals, matched-set Z-scores, LOESS predictions, clustering, sensitivity analyses and run summaries. The revised runs were also identical to each other. Figure 1D now uses the manuscript's aggregate definition; the per-metabolite completeness view is retained as a separate diagnostic.

The maximum post-match absolute SMD remained **0.111126042215701** for age. One of 13 covariates exceeded the 0.10 diagnostic threshold; no covariate/indicator diagnostic was unevaluable. The original warning is retained in [warnings.log](warnings.log). Completion is not a claim of perfect balance or biological validity.

## Comparison methods and exclusions

[compare_runs.R](compare_runs.R) compares every original CSV and RDS data output with `identical()` after loading it. File inventories permit only the added Figure 1D aggregate source table in CSV/RDS form. The baseline/revised comparison excludes `analysis/source_md5.csv`, which necessarily changes with the edited source; the repeat comparison includes it.

The serialized MatchIt object contains formula environments that differ across fresh R sessions. Its six substantive components are compared separately: `match.matrix`, `subclass`, `weights`, `discarded`, `treat` and `X`. Other RDS objects are compared in full. There are no tolerances that could conceal numerical differences.

Figure 1D is independently recalculated from the finite matched-set Z-scores, without calling the new completeness helper. Plot inventories are checked, and each PNG must be nonempty with the correct signature. A separate Pillow verification fully decoded every final image; [figure_file_checks.csv](figure_file_checks.csv) records dimensions, sizes and SHA-256 hashes. Exported panels were also visually inspected, including the aggregate Figure 1D, its detailed diagnostic, balance diagnostics and sensitivity panels.

During a preliminary parallel run before the export guard was added, one `figure2_cluster_curves.png` file was empty despite a completed status. Its numerical outputs passed comparison, and a repeat produced the image normally; the cause was not established. The final source now stops on missing, empty or invalid-signature PNG exports. The two final revised runs recorded here both passed the export guard and full decoding. The preliminary run is not counted as a successful final validation.

## Environment

R 4.3.3 on Ubuntu 24.04.3 LTS, x86_64 Linux:

| Package | Version |
|---|---|
| MatchIt | 4.5.5 |
| survival | 3.5-8 |
| ggplot2 | 3.4.4 |
| mclust | 6.0.1 |
| cluster | 2.1.6 |
| dplyr | 1.1.4 |
| tidyr | 1.3.1 |
| tibble | 3.2.1 |
| purrr | 1.0.2 |

See [sessionInfo.txt](sessionInfo.txt) for the final pipeline session. R and dependencies were extracted from signed Ubuntu snapshot packages into an isolated local runtime; no system package installation was required. This is a new Linux verification, separate from the earlier macOS/R 4.4.2 balance-synchronization record. The two environments are not treated as an identical software lockfile.

## Reproduce

From the revised repository root, after installing the documented dependencies:

```bash
Rscript --vanilla tests/test_study_safeguards.R
Rscript --vanilla Protocol_demo/tests/test_completeness.R
Rscript --vanilla Protocol_demo/tests/test_balance.R
Rscript --vanilla Protocol_demo/run_synthetic_example.R demo_run_01 5000 20 20260801
Rscript --vanilla Protocol_demo/run_synthetic_example.R demo_run_02 5000 20 20260801
```

In a separate checkout of the baseline commit, run its unchanged driver in the same R environment, writing to a new `baseline_run` directory. Then compare all three output directories:

```bash
Rscript --vanilla Protocol_demo/validation/2026-10-29_repository_fixes/compare_runs.R \
  /path/to/baseline_run demo_run_01 demo_run_02 verification_checks.csv
```

The baseline path above is an explicit placeholder. Use new or empty output directories for every run. The checked-in logs replace the temporary absolute output directory with `<baseline>`, `<run_1>` or `<run_2>`; other output is retained. Historical records under `2026-10-29_balance_sync/` are unchanged.

## Retained records

- `test_summary.csv` and the three `*_tests.log` files: 88 focused checks.
- `baseline_execution.log`, `run_1_execution.log` and `run_2_execution.log`: full execution messages.
- `execution_status.txt`, `execution_checks.csv`, `run_summary.csv`, `warnings.log` and `sessionInfo.txt`: final revised run status and environment.
- `comparison_checks.csv` and `comparison.log`: all 268 comparison/export checks and their summary.
- `figure1D_z_completeness_source.csv` and `covariate_balance_plot.csv`: aggregate diagnostic results.
- `figure_file_checks.csv`: fully decoded image inventory and hashes.
- `source_fingerprints.csv`, `source_md5.csv` and `baseline_source_md5.csv`: exact tested source contents.

These compact records omit participant-level synthetic intermediates. The scripts regenerate the full local outputs when run with the recorded configuration.
