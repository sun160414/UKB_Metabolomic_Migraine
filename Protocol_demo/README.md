# Standalone synthetic demonstration of matched time-to-diagnosis gradients

[Back to repository overview](../README.md)

A standalone companion for **Protocol for matched time-to-diagnosis gradient analysis of baseline plasma metabolomics in prospective cohorts**.

**All default inputs and outputs are fully synthetic. No UK Biobank participant records, fitted UKB parameters, private intermediate files, or original study result tables are used.** The example is a software demonstration, not a new biological analysis, independent validation, or a reproduction of the published migraine results.

## Execution status

The scripts have received static inspection, but **no completed R execution is documented in this update**. Static inspection is not an R parser check or a successful end-to-end test. Run the default example in a clean R session, inspect all warnings, tables and figures, and retain the generated session information before reporting successful execution.

This companion is independent of the original study scripts in the repository root. It generates its own inputs and recomputes Cox results locally; it does not modify or source the study-specific analysis files.

## 1. Files

Keep the three files together in the repository's **`Protocol_demo/`** folder. The directory name starts with an uppercase **P**; the R filenames start with lowercase letters.

```text
Protocol_demo/
├── generate_synthetic_data.R
├── run_synthetic_example.R
└── README.md
```

- `generate_synthetic_data.R`: generates a participant-level cohort, a raw metabolite matrix, a preprocessed matrix, a metabolite dictionary, and preprocessing audit files. It requires base R only. It does **not** fit Cox models or produce final analysis figures.
- `run_synthetic_example.R`: calls the generator, fits Cox models from the newly generated data, performs matching and the downstream protocol, and writes analysis tables, RDS checkpoints, PNG figures, configuration, warnings and execution checks. It can also be sourced to expose functions without executing the default example.
- `README.md`: installation, run instructions, inputs, outputs, parameter changes and limitations.

No other files from the study repository need to be sourced. In particular, this example does not depend on `NMR.csv`, `Map.txt`, `dat`, `Data_Met_ascii.xlsx`, or an existing R workspace.

## 2. Dependencies

Use R 4.1 or later. Install the required packages once in R:

```r
install.packages(c("MatchIt", "survival", "ggplot2", "mclust", "cluster"))
```

The plotting code uses `linewidth`, so use **ggplot2 3.4.0 or later**. Matching uses the documented `matchit()` and `match.data()` APIs. Exact package versions in a successful local run are recorded in `analysis/sessionInfo.txt`; this bundle is not accompanied by a tested lockfile. The script reports missing packages rather than installing software automatically.

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
- Covariates are explicitly supplied separately for matching, exact matching, residualization and Cox regression. The default list follows the protocol example. Medication use, fasting, diet and physical activity are not added silently.
- A post-match absolute SMD threshold of 0.10 produces a diagnostic warning, not an automatic declaration that the analysis is valid or invalid. Exact-matching factors and factor levels are included in MatchIt balance summaries.
- Minimum two valid control residuals and control SD greater than `1e-12` are required for a Z-score. This is a numerical safeguard, **not** proof of denominator stability. Additional control-SD quantiles/plots are diagnostic only and do not alter primary Z-scores.
- Primary LOESS uses span 0.75, degree 2, Gaussian family, `surface = "interpolate"`, grid interval 0.20 years, at least 50 valid sets and 10 distinct time values per metabolite. Unsupported predictions remain missing.
- Curves with missing grid predictions are excluded from clustering and recorded. A run with too few significant or complete curves stops; the script does not loosen the P-value threshold, choose a different seed, or force a top-N list to obtain a figure.
- Winsorization changes sensitivity-analysis inputs only; primary values remain intact.
- Alternative scenarios are compared over common supported time points. In the early-diagnosis exclusion, cluster agreement also reflects the restricted time window. RMSE is exported in addition to the protocol's curve correlations; it is an additional descriptive diagnostic, not an inferential test.
- The extra silhouette and cluster-count comparison outputs implement a reviewer-suggested diagnostic. Set `run_cluster_count_diagnostics = FALSE` to omit them. They do not automatically optimize or replace the prespecified primary `k`.

## 5. Expected output structure

The following are **expected files on successful execution**, not precomputed or verified R outputs supplied with this bundle:

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
    ├── sessionInfo.txt
    ├── warnings.log
    ├── execution_status.txt
    ├── data_processed/   # RDS checkpoints and table copies
    ├── tables/           # CSV results and figure-source tables
    └── figures/          # PNG demonstration figures
```

Important tables include `cohort_flow.csv`, `cohort_exclusions.csv`, `synthetic_cox_results.csv`, `matched_set_integrity.csv`, `covariate_balance.csv`, `case_metabolite_z_scores_long.csv`, `z_score_completeness.csv`, `control_sd_diagnostics.csv`, `loess_fit_audit.csv`, `selected_metabolites.csv`, `cluster_membership.csv`, `cluster_summary.csv`, `clustering_excluded_metabolites.csv`, `sensitivity_curve_similarity.csv`, `sensitivity_cluster_agreement.csv`, `execution_checks.csv` and `run_summary.csv`.

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

1. Run the command-line example with `--vanilla` after installing dependencies.
2. Check the execution status, all execution-check rows, unmatched cases, balance diagnostics, excluded curves and warnings.
3. Inspect every exported figure and its source table. In particular, check that Ward.D2 is absent from the curve-correlation panel.
4. Rerun the same configuration into a new folder and compare numerical tables (rather than RDS binary hashes or timestamp/session text).
5. Retain the tested R/package versions and code commit. Do not claim an end-to-end test until these steps have actually succeeded.

## 10. Implementation references

The functions implement the supplied protocol's primary settings. The substantive control-selection and same-cohort selection decisions have not been changed. The new driver supplies the missing data connection, Cox calculation, parameter interfaces, plotting/export steps, checks and explicit failure handling.

- Associated study repository: https://github.com/sun160414/UKB_Metabolomic_Migraine
- MatchIt nearest-neighbour matching documentation: https://kosukeimai.github.io/MatchIt/reference/method_nearest.html
- MatchIt balance-summary documentation: https://kosukeimai.github.io/MatchIt/reference/summary.matchit.html
- R LOESS prediction documentation: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/predict.loess.html
- R Cox model documentation: https://stat.ethz.ch/R-manual/R-devel/library/survival/html/coxph.html

Installation and run commands assume the `Protocol_demo/` layout shown above. The original study scripts remain in the repository root; do not source them to run this demonstration.
