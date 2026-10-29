# Standalone synthetic demonstration of matched time-to-diagnosis gradients

[Back to repository overview](../README.md)

A standalone companion for **Protocol for matched time-to-diagnosis gradient analysis of baseline plasma metabolomics in prospective cohorts**.

**All default inputs and outputs are fully synthetic. No UK Biobank participant records, fitted UKB parameters, private intermediate files, or original study result tables are used.** The example is a software demonstration, not a new biological analysis, independent validation, or a reproduction of the published migraine results.

## Execution status

See the [versioned validation records](validation/2026-10-29_balance_sync/) for the tested code revision, configuration, execution checks, warnings and comparison results. A successful run applies to the recorded revision and configuration, not automatically to later code changes. Successful execution confirms computational functionality but does not constitute independent biological validation.

This companion is independent of the original study scripts in the repository root. It generates its own inputs and recomputes Cox results locally; it does not modify or source the study-specific analysis files.

## 1. Files

Keep the two R scripts together in the repository's **`Protocol_demo/`** folder. The directory name starts with an uppercase **P**; the R filenames start with lowercase letters.

```text
Protocol_demo/
├── generate_synthetic_data.R
├── run_synthetic_example.R
├── README.md
├── tests/test_balance.R
└── validation/2026-10-29_balance_sync/
```

- `generate_synthetic_data.R`: generates a participant-level cohort, a raw metabolite matrix, a preprocessed matrix, a metabolite dictionary, and preprocessing audit files. It requires base R only. It does **not** fit Cox models or produce final analysis figures.
- `run_synthetic_example.R`: calls the generator, fits Cox models from the newly generated data, performs matching and the downstream protocol, and writes analysis tables, RDS checkpoints, PNG figures, configuration, warnings and execution checks. It can also be sourced to expose functions without executing the default example.
- `README.md`: installation, run instructions, inputs, outputs, parameter changes and limitations.
- `tests/test_balance.R`: focused regression checks for SMD calculation, categorical aggregation and unevaluable diagnostics.
- `validation/`: compact, version-specific verification records; not participant-level study data.

No other files from the study repository need to be sourced. In particular, this example does not depend on `NMR.csv`, `Map.txt`, `dat`, `Data_Met_ascii.xlsx`, or an existing R workspace.

## 2. Dependencies

Use R 4.1 or later. Install the required packages once in R:

```r
install.packages(c("MatchIt", "survival", "ggplot2", "mclust", "cluster"))
```

The plotting code uses `linewidth`, so use **ggplot2 3.4.0 or later**. Matching uses the documented `matchit()` and `match.data()` APIs. The verification environment is R 4.4.2 with MatchIt 4.7.2, survival 3.7-0, ggplot2 3.5.2, mclust 6.1.2 and cluster 2.1.6; see the versioned records above for the associated source revision and execution status. These versions describe the verification environment, not an enforced lockfile or a guarantee for other versions. The script reports missing packages rather than installing software automatically.

The standalone demo uses `survival` to fit Cox models, `mclust` for adjusted Rand indices, and `cluster` for silhouette diagnostics. The manuscript code boxes also directly use `data.table` 1.17.0, `readxl` 1.4.5, `patchwork` 1.3.0 and `scales` 1.3.0, as recorded in the original analysis sessions. Those four are not additional direct dependencies of this demo: it uses base R for data handling and exports separate figure panels. Packages such as `scales` may still be loaded indirectly by ggplot2.

Each run saves `analysis/sessionInfo_initial.txt` at entry and `analysis/sessionInfo.txt` at exit, including after an analysis error. Inspect `execution_status.txt` alongside either file; a session record alone does not prove completion. `analysis/source_md5.csv` identifies the runner and generator contents used by the run.

## 3. Run the entire example

### Command line

From the repository root:

```bash
Rscript --vanilla Protocol_demo/run_synthetic_example.R
```

This creates `Protocol_demo/synthetic_demo_output/`. The runner locates the generator relative to its own file, not through a private absolute path. To choose the output directory, sample size, metabolite count and seed:

```bash
Rscript --vanilla Protocol_demo/run_synthetic_example.R demo_run_01 5000 20 20260801
```

### R or RStudio

From the repository root in a fresh session:

```r
source("Protocol_demo/run_synthetic_example.R")
result <- run_synthetic_example(out_dir = "demo_run_01")
result$summary
result$checks
```

Sourcing the runner defines its functions without automatically performing the analysis. Use one of the two methods above, not both with the same output path.

**Use an empty/new output directory for every run.** The analysis driver refuses to write into a nonempty directory, preventing accidental mixtures of old and new results. It does not delete old runs. The standalone generator itself writes its named outputs, so likewise give it a new folder.

## 4. What the default example does

The default settings are 5,000 artificial participants, 20 generic metabolite measurements and seed `20260801`. Covariate distributions, correlated latent metabolite blocks, event-time coefficients, censoring, missingness and nonpositive test values are arbitrary simulation choices, not estimates from UKB. The endpoint is generic, not a calibrated migraine phenotype. The realized number of cases and selected metabolites must be read from the actual run; no particular clustering or sensitivity result is guaranteed.

The sequence is:

1. Generate baseline covariates, dates, a generic future event, correlated positive metabolite values, missing cells and a small number of nonpositive test values.
2. Impute each metabolite using its observed median in the preprocessing dataset, replace nonpositive values with half its smallest positive value, natural-log transform, then standardize using the sample SD. All-missing, no-positive-value and zero-variance columns cause explicit errors rather than becoming all-zero columns.
3. Construct the eligible prospective cohort from dates; require complete model/design covariates and valid preprocessed metabolites. Invalid observation intervals and prevalent diagnoses are excluded and recorded. Missing covariates are not silently imputed by the analysis driver.
4. Fit one covariate-adjusted Cox model per metabolite on the full eligible synthetic cohort, using time from baseline to event/censoring and Efron handling of ties. Calculate Bonferroni and BH-adjusted P values from these newly fitted models. No old results table is reused.
5. Match incident cases to participants with no recorded event in their own available follow-up, using 1:10 nearest-neighbour Mahalanobis matching, no replacement and exact matching on sex. Data order is made explicit. Verify one case, the specified control count and exact-match agreement per set; save overall balance diagnostics.
6. Fit covariate-adjusted linear models in the matched data. For each metabolite and matched set, compute `(case residual - control residual mean) / control residual sample SD`, retaining diagnostic numerator/denominator values. Invalid Z-scores are left missing.
7. Fit descriptive LOESS curves for all metabolites using the fixed primary settings, then select metabolites with same-cohort Cox Bonferroni-adjusted P < 0.05 for visualization and clustering. The prediction grid is rounded **inward**, never extended beyond observed time support.
8. Cluster complete selected prediction vectors using Euclidean distance, `ward.D` linkage and three clusters. Raw cluster labels are preserved; display labels are deterministically ordered by mean profile. Labels are neutral, not biological subtype names.
9. Save diagnostic plots, a heatmap, cluster curves and their source tables. Run alternative spans (0.50 and 1.00), exclusion of diagnoses within one year, 1st/99th-percentile per-metabolite Z-score winsorization, and `ward.D2` clustering. Ward.D2 changes clustering only and therefore is **not** shown in the curve-correlation panel.
10. Optionally report mean silhouette widths and neighboring cluster-count memberships/merge–split tables, without selecting a new primary cluster count.

### Primary settings and deliberate safeguards

- **Control selection remains the original observed-noncase design**, not case-date-specific risk-set sampling. Future cases do not serve as controls in this implementation. No check that every control remains under follow-up until its paired case's date is imposed.
- Covariates are explicitly supplied separately for matching, exact matching, residualization and Cox regression. The default covariate names follow the protocol example; synthetic education uses a different encoding, described below. Medication use, fasting, diet and physical activity are not added silently.
- A post-match absolute SMD threshold of 0.10 produces a diagnostic warning, not an automatic declaration that the analysis is valid or invalid. It does not trigger automatic rematching or exclusion. Exact-matching variables are included in the balance diagnostics.
- Minimum two valid control residuals and control SD greater than `1e-12` are required for a Z-score. This is a numerical safeguard, **not** proof of denominator stability. Additional control-SD quantiles/plots are diagnostic only and do not alter primary Z-scores.
- Primary LOESS uses span 0.75, degree 2, Gaussian family, `surface = "interpolate"`, grid interval 0.20 years, at least 50 valid sets and 10 distinct time values per metabolite. Unsupported predictions remain missing.
- Curves with missing grid predictions are excluded from clustering and recorded. A run with too few significant or complete curves stops; the script does not loosen the P-value threshold, choose a different seed, or force a top-N list to obtain a figure.
- Winsorization changes sensitivity-analysis inputs only; primary values remain intact.
- Alternative scenarios are compared over common supported time points. In the early-diagnosis exclusion, cluster agreement also reflects the restricted time window. RMSE is exported in addition to the protocol's curve correlations; it is an additional descriptive diagnostic, not an inferential test.
- The extra silhouette and cluster-count comparison outputs implement a reviewer-suggested diagnostic. Set `run_cluster_count_diagnostics = FALSE` to omit them. They do not automatically optimize or replace the prespecified primary `k`.

### Figure 1A balance definition and source tables

To match the manuscript figure, post-match SMDs use all matched cases and controls **without weighting**, under the implemented fixed-ratio, no-replacement design:

```text
SMD = (mean_case - mean_control) / sqrt((var_case + var_control) / 2)
```

Both variances are sample variances calculated after matching. Numeric covariates contribute one SMD; categorical covariates contribute one SMD for each observed 0/1 indicator level. The plot shows one point per evaluable covariate, using the maximum absolute level-specific SMD for a categorical covariate. This is a display summary, not a replacement for the full level-specific results.

- `covariate_balance.csv` and its RDS copy retain the complete covariate/indicator-level results, group statistics and diagnostic status.
- `covariate_balance_plot.csv` and its RDS copy contain one summary row per covariate and are the Figure 1A plotting source.
- `covariate_balance_matchit.csv` and its RDS copy retain the separate MatchIt summary, including `abs_smd_matchit` and `abs_smd_before_matchit`. These use MatchIt's own standardization and are not the values plotted in Figure 1A; do not combine their before-match values with the new matched-group SMDs as if the denominators were identical.

Missing/nonfinite values, insufficient group sizes, or unequal constant groups make an SMD unevaluable. Identical constant groups are recorded as zero. If any indicator level is unevaluable, the covariate's plotting summary is `NA`, not the maximum of only its remaining levels. The plot retains its row with **Not evaluable** instead of placing a point at zero or dropping it silently. Unevaluable diagnostics are counted separately; `max_postmatch_absolute_smd` is `NA` if any covariate is unevaluable. A low maximum among a subset must not be interpreted as complete balance.

### Education encoding

The synthetic generator deliberately creates `education` as the arbitrary categories `Low`, `Medium` and `High`; the demo therefore lists it in `cfg$factor_covariates`. In the original worked example, `education` corresponds to the derived **numeric** `Qualification` covariate. These encodings are not equivalent. For a dataset with the original numeric education variable, retain its numeric values and remove `education` from `cfg$factor_covariates`. Do not simply remove that setting while using the synthetic text categories, or recode `Low/Medium/High` as `1/2/3` and claim equivalence to the original variable.

## 5. Expected output structure

The following are **expected files on successful execution**. The compact versioned verification records are separate from these full, locally generated outputs:

```text
synthetic_demo_output/
├── SYNTHETIC_ONLY.txt
├── input/
│   ├── cohort_raw.csv
│   ├── cohort.csv
│   ├── cohort_raw.rds
│   ├── cohort.rds
│   ├── metabolite_dictionary.csv
│   ├── metabolite_preprocessing_qc.csv
│   ├── validation_checks_R.csv
│   └── simulation_parameters.rds
└── analysis/
    ├── analysis_configuration.R
    ├── sessionInfo_initial.txt
    ├── sessionInfo.txt
    ├── source_md5.csv
    ├── warnings.log
    ├── execution_status.txt
    ├── data_processed/   # RDS checkpoints and table copies
    ├── tables/           # CSV results and figure-source tables
    └── figures/          # PNG demonstration figures
```

Important tables include `cohort_flow.csv`, `cohort_exclusions.csv`, `synthetic_cox_results.csv`, `matched_set_integrity.csv`, the three covariate-balance tables described above, `case_metabolite_z_scores_long.csv`, `z_score_completeness.csv`, `control_sd_diagnostics.csv`, `loess_fit_audit.csv`, `selected_metabolites.csv`, `cluster_membership.csv`, `cluster_summary.csv`, `clustering_excluded_metabolites.csv`, `sensitivity_curve_similarity.csv`, `sensitivity_cluster_agreement.csv`, `execution_checks.csv` and `run_summary.csv`. Tables written by the analysis driver have matching RDS copies in `data_processed/`.

Demonstration figure panels are exported separately as PNG files rather than replacing the published composite figures. The main families are `figure1A_*` to `figure1D_*`, `figure2_*`, `figure3A_*` and `figure3B_*`. A denominator plot and, by default, a cluster-count diagnostic are also generated. Figures are prominently labeled **SYNTHETIC DEMONSTRATION**.

A `COMPLETED` status means that the code reached its end and passed its programmed execution checks. It does **not** mean all covariates are balanced, there is no selection bias, or the method has been biologically validated. Inspect `warnings.log` even when execution completes. After an error, use a new output folder for the next attempt.

## 6. Change parameters without editing every model

```r
source("Protocol_demo/run_synthetic_example.R")
cfg <- protocol_config()
cfg$n <- 6000L
cfg$n_metabolites <- 24L
cfg$seed <- 20260802L
cfg$balance_threshold <- 0.10
cfg$match_ratio <- 10L
cfg$loess_span <- 0.75
cfg$n_clusters <- 3L
cfg$make_plots <- TRUE
cfg$run_cluster_count_diagnostics <- TRUE
run_synthetic_example(out_dir = "demo_run_02", config = cfg)
```

Other parameters include the separate covariate vectors, minimum valid controls/SD, time-grid interval, minimum LOESS sample support, selection threshold, sensitivity spans, exclusion interval, winsorization percentiles, linkage method and diagnostic reference values. Prespecify changes and report them; do not optimize settings merely to obtain a desired synthetic result.

Changed simulation settings can produce insufficient cases, inadequate exact-stratum control counts or too few selected metabolites. The generator/driver should then fail explicitly. That is an informative failure, not a reason to silently relax the primary method. Cox models are recomputed automatically after every new simulation.

## 7. Use the generator only

```r
source("Protocol_demo/generate_synthetic_data.R")
x <- generate_synthetic_data("new_synthetic_input", n = 5000,
                             n_metabolites = 20, seed = 20260801)
```

Or:

```bash
Rscript --vanilla Protocol_demo/generate_synthetic_data.R new_synthetic_input 5000 20 20260801
```

`cohort_raw.csv` contains the deliberately introduced missing/nonpositive values. `cohort.csv` is already imputed, log-transformed and standardized. **Do not log-transform the processed file again**; its negative values are valid standardized measurements. The two files represent the same artificial participants at different preprocessing stages, not independent cohorts.

## 8. Adapt the analysis interface to another approved dataset

The computational functions take ordinary data frames and explicit arguments. They do not fetch UKB or depend on UKB field numbers. This software interface is not evidence of cross-cohort or cross-platform validation.

Your preprocessed table should have one row per participant, valid dates, the prespecified covariates and an explicit set of numeric metabolite columns. The dictionary must contain unique `Metabolite_raw` identifiers and unique `Metabolite_name` labels. Rename nonstandard model-variable names once at import.

```r
source("Protocol_demo/run_synthetic_example.R")
cfg <- protocol_config()
cfg$columns <- c(id = "id", baseline = "sample_date",
                 diagnosis = "event_date", followup_end = "observation_end")
# Example only: choose covariates based on your actual study, not this short list.
cfg$matching_covariates <- c("age", "bmi")
cfg$exact_covariates <- "sex"
cfg$adjustment_covariates <- c("age", "sex", "bmi")
cfg$cox_covariates <- c("age", "sex", "bmi")
cfg$factor_covariates <- "sex"

# my_preprocessed_cohort and my_dictionary must be supplied by you.
# Do not place restricted participant-level outputs in a public repository.
# run_gradient_pipeline(my_preprocessed_cohort,
#   metabolite_cols = c("marker_01", "marker_02", "marker_03", "marker_04"),
#   dictionary = my_dictionary, out_dir = "local_approved_analysis",
#   config = cfg, synthetic = FALSE)
```

`followup_end` means the end of available outcome observation, not the event date substituted for every case. The function derives case status and time-to-event/censoring from dates. It does not trust an existing disease flag or infer event timing from row order. Missing/unparsable diagnosis dates must be resolved before use, so they are not accidentally reclassified as controls. Platform-specific measurement quality control remains a prerequisite; no purported UKB QC flag is introduced in the generic example.

## 9. Scientific interpretation and local verification

One baseline measurement per participant supports **between-person descriptive gradients**, not within-person longitudinal change. Same-cohort Cox selection and curve characterization are not independent; no sample split or cross-validation of that selection is implemented. Neither a successful run nor the synthetic generation model demonstrates that real-study gradients are unbiased. The observed-noncase control strategy, possible medication/fasting confounding, small-set standardization and strong metabolite correlations retain their limitations.

Before uploading a tested release, locally:

1. Run the regression checks, then the command-line example with `--vanilla` after installing dependencies:

   ```bash
   Rscript --vanilla Protocol_demo/tests/test_balance.R
   Rscript --vanilla Protocol_demo/run_synthetic_example.R demo_run_01 5000 20 20260801
   ```

2. Check the execution status, all execution-check rows, unmatched cases, balance diagnostics, excluded curves and warnings.
3. Inspect every exported figure and its source table. In particular, verify one summary row per matching/exact covariate in Figure 1A, preservation of the full indicator-level results, and absence of Ward.D2 from the curve-correlation panel.
4. Rerun the same configuration into a new folder and compare numerical tables (rather than RDS binary hashes or timestamp/session text).
5. Retain the tested R/package versions, code commit and source fingerprints. Do not claim an end-to-end test until these steps have actually succeeded. New code requires new verification records; earlier successful logs do not validate a later revision. Remove personal paths/usernames from public logs and never upload restricted participant-level inputs or intermediates.

## 10. Implementation references

The functions implement the supplied protocol's primary settings. The substantive control-selection and same-cohort selection decisions have not been changed. The new driver supplies the missing data connection, Cox calculation, parameter interfaces, plotting/export steps, checks and explicit failure handling.

- Associated study repository: https://github.com/sun160414/UKB_Metabolomic_Migraine
- MatchIt nearest-neighbour matching documentation: https://kosukeimai.github.io/MatchIt/reference/method_nearest.html
- MatchIt balance-summary documentation: https://kosukeimai.github.io/MatchIt/reference/summary.matchit.html
- R LOESS prediction documentation: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/predict.loess.html
- R Cox model documentation: https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html

Installation and run commands assume the `Protocol_demo/` layout shown above. The original study scripts remain in the repository root; do not source them to run this demonstration.
