# Plasma Metabolomic Signatures of Migraine

<img width="6625" height="5852" alt="Fig1" src="https://github.com/user-attachments/assets/88834f66-8275-436b-90f7-3c935c6d4113" />

This repository contains analysis code for:

**Plasma metabolomic signatures of migraine in 479,760 adults**

The project investigates associations between baseline plasma nuclear magnetic resonance (NMR) metabolomic measures and hospital-diagnosed migraine in UK Biobank, with downstream analyses integrating matched time-to-diagnosis gradients, feature prioritization, brain MRI phenotypes, polygenic risk scores, sleep/affective traits, and exploratory structural equation modeling.

The repository also includes a separate synthetic-data companion for **Protocol for matched time-to-diagnosis gradient analysis of baseline plasma metabolomics in prospective cohorts**.

> **Choose the appropriate entry point.** To try the protocol without UK Biobank data, start with the [synthetic demonstration](Protocol_demo/README.md). The study-specific scripts in the repository root require locally prepared, authorized research data and are not the synthetic-demo entry point.

> **Data note.** No individual-level UK Biobank data are distributed here. The demonstration generates artificial records locally and does not use real participants, fitted UKB parameters, or the original study's results tables. Reanalysis of the original study requires authorized data access and local input preparation.

## Quick start: STAR Protocols synthetic demonstration

### 1. Install the demonstration dependencies

Use R 4.1 or later and ggplot2 3.4.0 or later. In R, install:

```r
install.packages(c("MatchIt", "survival", "ggplot2", "mclust", "cluster"))
```

These are the dependencies of the **synthetic demonstration**. The original study scripts have additional dependencies listed below.

### 2. Run from the repository root

After downloading or cloning the repository, run in a terminal:

```bash
Rscript --vanilla Protocol_demo/run_synthetic_example.R
```

The runner locates its companion generator in the same directory. The default output folder is `Protocol_demo/synthetic_demo_output/`; it must be new or empty. The directory name starts with an uppercase **P**; use the exact capitalization shown here.

Alternatively, from R or RStudio with the repository root as the working directory:

```r
source("Protocol_demo/run_synthetic_example.R")
result <- run_synthetic_example(out_dir = "demo_run_01")
result$summary
result$checks
```

Sourcing the runner only defines functions; the explicit `run_synthetic_example()` call starts the analysis. Use a new output directory for each run.

### 3. What the demonstration is designed to produce

The default configuration generates **5,000 artificial participants and 20 generic metabolite measurements**, using seed `20260801`. The script connects data generation and preprocessing to cohort construction, Cox association testing, matching, balance assessment, residualization, matched-set Z-scores, LOESS, descriptive clustering, sensitivity analyses, figures, and output tables.

It calculates Cox results from the generated cohort and does not require `NMR.csv`, `Map.txt`, `Data_Met_ascii.xlsx`, a preloaded `dat` object, or any of the original study scripts.

After a successful run, inspect:

```text
Protocol_demo/synthetic_demo_output/
├── SYNTHETIC_ONLY.txt
├── input/                       # Synthetic data, dictionary, preprocessing checks
└── analysis/
    ├── analysis_configuration.R
    ├── sessionInfo.txt
    ├── warnings.log
    ├── execution_status.txt
    ├── data_processed/          # RDS checkpoints
    ├── tables/                  # CSV results and figure-source tables
    └── figures/                 # Labeled synthetic demonstration PNGs
```

**Execution status:** the supplied scripts have received static inspection, but no completed R execution is documented in this update. Static inspection is not an R parser check or an end-to-end test. Run the example in a clean R session, inspect warnings and outputs, and retain the session information before reporting successful execution. A completed software run does not establish biological validity or eliminate confounding or selection bias.

See [the demonstration README](Protocol_demo/README.md) for parameter settings, output details, input-field mapping, scientific limitations, and verification steps.

## Repository structure

```text
.
├── README.md
├── 1.Cox+linear.R
├── 2.Trajectories+clusters.R
├── 3.PRS
├── 4.SEM.R
└── Protocol_demo/
    ├── generate_synthetic_data.R
    ├── run_synthetic_example.R
    └── README.md
```

The original study files remain in the **repository root**, not in a `scripts/` subfolder. The synthetic demonstration is separate and does not modify or source those files.

## Original study: analysis overview

### 1. Metabolite association analyses

`1.Cox+linear.R` contains the study's association-analysis workflow:

- Cox regression for incident hospital-diagnosed migraine, including age- and sex-stratified analyses.
- Cross-sectional linear regression for prevalent migraine.
- Bonferroni and Benjamini–Hochberg multiple-testing correction and volcano plots.

The original workflow documentation identifies the following main outputs:

```text
results/cox_overall.csv
results/cox_age_lt55.csv
results/cox_age_ge55.csv
results/cox_male.csv
results/cox_female.csv
results/linear_prevalent.csv
figures/Fig1A_cox_overall_volcano.pdf
figures/Fig1B_cox_age_lt55_volcano.pdf
figures/Fig1C_cox_age_ge55_volcano.pdf
figures/Fig1D_cox_male_volcano.pdf
figures/Fig1E_cox_female_volcano.pdf
figures/Fig1F_linear_prevalent_volcano.pdf
```

These are study-specific outputs, not files distributed as results of the synthetic demonstration.

### 2. Matched time-to-diagnosis gradient analysis

`2.Trajectories+clusters.R` contains the matched case-control gradient workflow. Incident cases are matched to controls using nearest-neighbor Mahalanobis matching with exact matching on sex. Metabolite values are residualized for covariates, and case-control differences are expressed as matched-set standardized Z-scores. LOESS and hierarchical clustering summarize descriptive baseline differences aligned to future diagnosis time.

The curves describe **between-person baseline gradients**, not within-person longitudinal trajectories.

The original workflow documentation identifies the following main outputs:

```text
results/match_incident_migraine_metabolome.csv
results/resid_ukb_metabolome.csv
results/matched_case_control_z_scores.csv
results/loess_trajectories_significant_long.csv
results/loess_trajectories_significant_wide.csv
results/cluster_metabolite_list_long.csv
results/cluster_metabolite_collapsed.csv
results/cluster_metabolite_list_wide.csv
figures/Fig2A_heatmap_metabolome.pdf
figures/Fig2B_trajectories_metabolome.pdf
```

### 3. PRS analysis

`3.PRS` is the original study's PRS file. Review its input requirements and commands before running it. PRS analysis is outside the scope of the standalone matched-gradient demonstration.

### 4. Exploratory SEM

`4.SEM.R` contains exploratory structural equation models relating sleep/affective traits, a latent factor based on 12 LASSO-prioritized metabolites, and incident hospital-diagnosed migraine. The documented exposure traits include depression, anxiety, insomnia, and sleep duration. These models describe statistical interrelationships and are not evidence of formal causal mediation.

The original workflow documentation identifies:

```text
results/sem_latent_metabolite_factor_loadings.csv
results/sem_model_based_association_parameters.csv
results/sem_fit_indices.csv
results/sem_latent_factor_results.rds
```

## Original study: expected inputs

This section applies to the **study-specific root scripts**, not to the canonical field names used by `Protocol_demo/`.

The participant-level study table has one row per participant and includes:

| Variable | Description |
|---|---|
| `eid` | Participant identifier |
| `prevalent_migraine` | Diagnosis before or at baseline |
| `incident_migraine` | Incident diagnosis after baseline |
| `followup_years` | Time from baseline to event or censoring |
| `migraine_years` | Diagnosis or censoring interval used by the study gradient script |
| `age`, `sex`, `ethn` | Baseline demographic variables |
| `Qualification` | Educational attainment |
| `bmi`, `Socioeconomic` | BMI and deprivation indicator |
| `Smoking_status`, `Alcohol_consumption` | Baseline smoking and alcohol categories |
| `screen time (TV)`, `screen time (computer)` | Daily screen-time variables |
| `sleep duration` | Average sleep duration |
| `diabetes_status`, `CVD_status` | Baseline comorbidities |
| `Meta_*` | Preprocessed metabolite variables, such as `Meta_1` |

The metabolite map `data/metabolite_name_map.csv` contains:

| Column | Description |
|---|---|
| `Meta` | Internal metabolite identifier, such as `Meta_1` |
| `Original_Metabolite` | Human-readable metabolite name |

Consult the configuration and formulas in each study script for its exact input requirements. Restricted source data and private intermediate objects must be prepared locally and must not be uploaded to this repository.

## Original study: dependencies and paths

The original README lists these R dependencies for the study scripts:

```r
install.packages(c(
  "data.table", "dplyr", "tidyr", "tibble", "purrr", "stringr",
  "survival", "ggplot2", "ggrepel", "MatchIt", "patchwork",
  "circlize", "lavaan", "fastDummies"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")
```

After preparing authorized inputs and checking each script's paths, the root-relative R entry points are:

```r
# Study-specific commands: not a standalone synthetic example.
# Prepare the required dat object and file inputs before sourcing these scripts.
source("1.Cox+linear.R")

# Prepare data/mydata_base.csv, the dictionary, and other required inputs first.
source("2.Trajectories+clusters.R")

# Prepare the required SEM inputs before this step.
dat_base <- dat
source("4.SEM.R")
```

Do not prefix these paths with `scripts/` unless the files are deliberately moved to that directory. For the fully synthetic example, use only the quick-start commands at the top of this README.

## Interpretation and scope

- The original study concerns hospital-diagnosed migraine and should not be interpreted as representing all community-managed migraine.
- Many NMR lipid and lipoprotein measures are correlated; individual signals and descriptive clusters are not independent biological discoveries.
- The matched-gradient workflow uses one baseline measurement per participant and cannot estimate within-person change.
- The demonstration retains observed-noncase controls rather than case-date-specific risk-set sampling, and same-cohort Cox selection rather than sample-split validation. These choices are explicit in its documentation.
- The study's internally prioritized metabolites are not an externally validated prediction panel, and the SEM analyses do not establish causal mediation.
- Synthetic execution is a software check only. It does not validate the original results or demonstrate generalizability across cohorts, omics platforms, or outcomes.

## Citation

For the associated study:

> Hong Y, Chen F, Wang Y, and Huang X. F. (2026). Plasma metabolomic signatures of migraine in 479,760 adults. *iScience*, 29(8), 117031. https://doi.org/10.1016/j.isci.2026.117031

The companion protocol is titled **Protocol for matched time-to-diagnosis gradient analysis of baseline plasma metabolomics in prospective cohorts**. Cite the final protocol publication when its bibliographic details are available; the synthetic example does not constitute a separate biological study.

## Contact

For questions about the original study or protocol implementation, contact the corresponding author listed in the associated manuscript.
