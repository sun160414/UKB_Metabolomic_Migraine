# Plasma Metabolomic Signatures of Migraine

<img width="7300" height="7049" alt="Fig1" src="https://github.com/user-attachments/assets/8d3f0194-46bb-406a-ae68-b8ba40cba8ae" />

This repository contains analysis code for the project:

**Plasma metabolomic signatures of migraine in 479,760 adults**

The project investigates associations between baseline plasma nuclear magnetic resonance (NMR) metabolomic measures and hospital-diagnosed migraine in UK Biobank, with downstream analyses integrating matched time-to-diagnosis gradients, feature prioritization, brain MRI phenotypes, polygenic risk scores, sleep/affective traits, and exploratory structural equation modeling.

> **Data note**
> UK Biobank individual-level data are not included in this repository. Users must obtain access to UK Biobank data through the official application process and prepare the required input files locally.

---

## Repository structure

```text
.
├── README.md
├── scripts/
│   ├── 1.Cox+linear.R
│   ├── 2.Trajectories+clusters.R
│   └── 3.PRS
│   └── 4.SEM.R
├── data/
│   ├── analysis_dataset.rds              # not provided
│   ├── mydata_base.csv                   # not provided
│   └── metabolite_name_map.csv           # not provided
├── results/
└── figures/
```

---

## Analysis overview

### 1. Metabolite association analyses

`scripts/01_cox_linear_github.R` performs:

* Cox proportional hazards regression for incident hospital-diagnosed migraine.
* Age-stratified Cox analyses: `<55 years` and `≥55 years`.
* Sex-stratified Cox analyses: male and female participants.
* Cross-sectional linear regression for prevalent migraine.
* Multiple testing correction using Bonferroni and Benjamini-Hochberg FDR.
* Volcano plots for the main and stratified analyses.

Main outputs:

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

### 2. Matched time-to-diagnosis gradient analysis

`scripts/02_trajectories_clusters_github.R` performs a matched case-control time-to-diagnosis gradient analysis. Incident migraine cases are matched to controls using nearest-neighbor Mahalanobis matching with exact matching on sex. Metabolite values are residualized for covariates, and case-control differences are expressed as standardized Z-scores within each matched subclass. LOESS smoothing and hierarchical clustering are then used to summarize descriptive baseline metabolite differences aligned to future diagnosis time.

This analysis should be interpreted as a **between-person descriptive gradient analysis**, not as a within-person longitudinal trajectory analysis, because metabolites were measured once at baseline.

Main outputs:

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

### 3. Exploratory SEM analysis

`scripts/03_sem_latent_factor_github.R` fits exploratory structural equation models linking selected sleep/affective traits, a latent metabolite factor based on 12 LASSO-prioritized metabolites, and incident hospital-diagnosed migraine.

The default exposures are:

```text
Depression
Anxiety
Insomnia
Sleep duration
```

The SEM results are intended to summarize statistical interrelationships and should not be interpreted as formal causal mediation evidence.

Main outputs:

```text
results/sem_latent_metabolite_factor_loadings.csv
results/sem_model_based_association_parameters.csv
results/sem_fit_indices.csv
results/sem_latent_factor_results.rds
```

---

## Expected input data

### Main analysis dataset

The main analysis dataset should contain one row per participant and include at least the following variables:

| Variable                 | Description                                                                                             |
| ------------------------ | ------------------------------------------------------------------------------------------------------- |
| `eid`                    | Participant identifier                                                                                  |
| `prevalent_migraine`     | Migraine diagnosis before or at baseline                                                                |
| `incident_migraine`      | Incident migraine diagnosis after baseline                                                              |
| `followup_years`         | Follow-up time from baseline to event or censoring                                                      |
| `migraine_years`         | Time from baseline to incident migraine diagnosis or censoring; required by the matched gradient script |
| `age`                    | Age at baseline                                                                                         |
| `sex`                    | Sex                                                                                                     |
| `ethn`                   | Ethnicity category                                                                                      |
| `Qualification`          | Educational attainment variable                                                                         |
| `bmi`                    | Body mass index                                                                                         |
| `Socioeconomic`          | Townsend deprivation index or equivalent socioeconomic indicator                                        |
| `Smoking_status`         | Smoking status                                                                                          |
| `Alcohol_consumption`    | Alcohol consumption status                                                                              |
| `screen time (TV)`       | Daily TV screen time                                                                                    |
| `screen time (computer)` | Daily computer screen time                                                                              |
| `sleep duration`         | Average sleep duration                                                                                  |
| `diabetes_status`        | Baseline diabetes status                                                                                |
| `CVD_status`             | Baseline cardiovascular disease status                                                                  |
| `Meta_*`                 | Standardized metabolite variables, such as `Meta_1`, `Meta_2`, ...                                      |

### Metabolite name map

`data/metabolite_name_map.csv` should include:

| Column                | Description                                  |
| --------------------- | -------------------------------------------- |
| `Meta`                | Internal metabolite ID, for example `Meta_1` |
| `Original_Metabolite` | Human-readable metabolite name               |

---

## Installation

The scripts were developed in R and require the following packages:

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

---

## Example workflow

### Step 1: Prepare data

```r
library(data.table)
library(dplyr)

# Individual-level UK Biobank data are not provided.
```

### Step 2: Run association analyses

```r
source("scripts/1.Cox+linear.R")
```

### Step 3: Prepare input for matched gradient analysis

```r
mydata_base <- dat %>%
  filter(prevalent_migraine == 0) %>%
  mutate(
    migraine_status = as.integer(incident_migraine == 1),
    migraine_years = as.numeric(followup_years)
  )

fwrite(mydata_base, "data/mydata_base.csv")
```

### Step 4: Run matched time-to-diagnosis gradient analysis

```r
source("scripts/2.Trajectories+clusters.R")
```

### Step 5: Run PRS

```
source("3.PRS")
```

---

### Step 6: Run exploratory SEM analysis

```r
dat_base <- dat
source("scripts/03_sem_latent_factor_github.R")
```

---

## Important interpretation notes

1. Migraine was defined using linked hospital inpatient records and ICD-10 code `G43`; therefore, the outcome should be interpreted as **hospital-diagnosed migraine**, not all migraine cases in the community.
2. Many NMR-derived lipoprotein metabolites are highly correlated. Results should be interpreted at the lipid/lipoprotein module level rather than as fully independent single-metabolite effects.
3. The matched time-to-diagnosis analysis uses baseline metabolite measurements only. It describes between-person gradients aligned to future diagnosis time and does not estimate within-person longitudinal metabolic change.
4. The LASSO-prioritized metabolites are internally selected features and should not be treated as an externally validated clinical prediction panel without independent validation.
5. The SEM analysis is exploratory and should not be interpreted as formal causal mediation.

---

## Citation

If you use this code, please cite the associated manuscript:

> Hong Y, Chen F, Wang Y, Huang X-F. *Plasma metabolomic signatures of migraine in 479,760 adults*.

---

## Contact

For questions about the analysis code, please contact the corresponding author listed in the associated manuscript.
