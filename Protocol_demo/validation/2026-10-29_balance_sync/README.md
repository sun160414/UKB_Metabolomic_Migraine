# Balance-diagnostic synchronization: execution record

Verified on **2026-10-29** using fully synthetic data only. These records apply to
the tested runner/generator contents, not automatically to future changes.

- Tested source revision: [`2a40c0e87cfdca0b402a2742aad6755bf3474137`](https://github.com/sun160414/UKB_Metabolomic_Migraine/commit/2a40c0e87cfdca0b402a2742aad6755bf3474137).
- Baseline revision: [`9f7c8125e35257a3b8743fef3dae209a4f17e3af`](https://github.com/sun160414/UKB_Metabolomic_Migraine/commit/9f7c8125e35257a3b8743fef3dae209a4f17e3af).
- Configuration: 5,000 synthetic participants, 20 metabolites, seed `20260801`;
  all other settings are `protocol_config()` defaults.
- Environment: R 4.4.2, MatchIt 4.7.2, survival 3.7-0, ggplot2 3.5.2,
  mclust 6.1.2 and cluster 2.1.6. Full initial/final sessions are retained.

## Results

Two fresh end-to-end runs completed and each passed all eight execution checks.
They retained 322 cases and 3,220 controls; 16 of 20 metabolites were selected and
clustered. The 53 balance regression checks also passed.

All **227** recorded comparison checks passed using strict R `identical()` on
parsed data or explicitly identified components:

- Baseline versus revised run: 47 unchanged CSV tables and 51 unchanged RDS
  objects; matched participants, residuals, Cox results, Z-scores, LOESS,
  clustering and sensitivity outputs were unchanged.
- Revised run versus repeat: 52 CSV tables and 56 RDS objects were identical.
- MatchIt objects contain R environments, so their matching matrix, subclasses,
  weights, discarded flags, treatment vector and covariate matrix were compared
  separately. All were identical in both comparisons.
- File inventories, original nonbalance summary metrics, matched-set integrity,
  preservation of the legacy MatchIt diagnostics, and independently recomputed
  SMDs and per-covariate maxima also passed.

The new balance tables contain **24** numeric/indicator-level results and **13**
covariate summaries. Synthetic education has three factor levels. Every SMD was
evaluable. Age had absolute SMD **0.111126**, above the diagnostic reference 0.10;
the warning was retained, with no automatic rematching. Computational completion
does not imply satisfactory balance for every covariate or biological validation.

All ten exported figures were visually inspected for titles, axes, legends and
synthetic-data labels. Figure 1A has one point per covariate. Ward.D2 appears only
in the cluster-agreement panel, not the curve-correlation panel. The regression
fixture with an undefined SMD was also inspected: its covariate remains visible
with a `Not evaluable` label and no false zero point.

R emitted the startup warning `Setting LC_CTYPE failed, using "C"`; it did not
prevent execution. This startup message occurs before the pipeline's warning
recorder. The pipeline's actual diagnostic warning is preserved in `warnings.log`.

## Retained evidence

`sessionInfo_initial.txt`, `sessionInfo.txt`, `source_md5.csv`,
`execution_status.txt`, `warnings.log`, `execution_checks.csv` and `run_summary.csv`
come from the first completed revised run. `test_balance_output.txt` records a
fresh successful execution of the 53 regression checks. `comparison_checks.csv`
records each comparison and its method; `compare_runs.R` regenerates those checks.
The three `covariate_balance*.csv` files are synthetic aggregate diagnostics.
No real participant data or personal filesystem paths are included.

The SMD definition and plot summary intentionally changed. Accordingly, baseline
comparisons exclude the old balance table, whole QC container and balance part
of the summary, then check preserved legacy values and unchanged QC components
separately. Serialized bytes, R environment identity, session text and PNG bytes
are not used as tests of numerical reproducibility.

## Reproduce

At the tested source revision, use new output directories:

```bash
Rscript --vanilla Protocol_demo/tests/test_balance.R
Rscript --vanilla Protocol_demo/run_synthetic_example.R run_1 5000 20 20260801
Rscript --vanilla Protocol_demo/run_synthetic_example.R run_2 5000 20 20260801
```

In a separate checkout of the baseline revision, run its unmodified demo with
the same arguments into another new directory, `baseline_run`. From this
validation directory, supply the three run-root paths (each containing `input/`
and `analysis/`) and a new output CSV path:

```bash
Rscript --vanilla compare_runs.R /path/to/run_1 /path/to/run_2 /path/to/baseline_run new_comparison_checks.csv
```

The comparison script is specific to this before/after change and default
configuration. The same environment was used for both runs and the baseline;
these results do not establish equality across other R/package versions.
