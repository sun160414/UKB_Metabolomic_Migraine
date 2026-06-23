#!/usr/bin/env Rscript

# ============================================================
# SEM analysis: latent metabolite factor and incident migraine
# ============================================================
# Purpose:
#   Fit exploratory structural equation models linking selected
#   sleep/affective traits, a latent metabolite factor, and
#   incident migraine.
#
# Required input object:
#   dat_base: individual-level analysis dataset containing:
#     - incident_migraine
#     - prevalent_depression, prevalent_anxiety, insomnia, sleep duration
#     - 12 LASSO-prioritized metabolite variables
#     - covariates listed below
#
# Optional input object:
#   meta_reverse_map: named vector mapping Met IDs to metabolite names
#
# Note:
#   These SEMs are exploratory statistical interrelationship models.
#   They should not be interpreted as formal causal mediation analyses.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(lavaan)
  library(fastDummies)
  library(stringr)
  library(purrr)
})

# -----------------------------
# User parameters
# -----------------------------
result_dir <- "results"
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

lasso_meta_cols <- c(
  "Met38", "Met57", "Met60", "Met82",
  "Met88", "Met111", "Met120", "Met125",
  "Met146", "Met201", "Met245", "Met246"
)

exposures <- c(
  dep_status   = "Depression",
  anx_status   = "Anxiety",
  insom_status = "Insomnia",
  sleep_hours  = "Sleep duration"
)

base_covars <- c(
  "age", "sex", "ethn", "Qualification", "bmi", "Socioeconomic",
  "Smoking_status", "Alcohol_consumption"
)

# -----------------------------
# Helper functions
# -----------------------------
check_required_vars <- function(dat, vars) {
  missing_vars <- setdiff(vars, names(dat))
  if (length(missing_vars) > 0) {
    stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
  }
}

make_metabolite_map <- function(meta_cols) {
  if (exists("meta_reverse_map")) {
    data.frame(
      Metabolite_ID = names(meta_reverse_map),
      Metabolite_name = unname(meta_reverse_map),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      Metabolite_ID = meta_cols,
      Metabolite_name = meta_cols,
      stringsAsFactors = FALSE
    )
  }
}

prepare_sem_data <- function(dat, meta_cols) {
  dat <- dat %>%
    mutate(
      migraine_incident = as.integer(incident_migraine == 1),
      dep_status   = as.integer(prevalent_depression == 1),
      anx_status   = as.integer(prevalent_anxiety == 1),
      insom_status = as.numeric(insomnia),
      sleep_hours  = as.numeric(`sleep duration`),
      age = as.numeric(age),
      bmi = as.numeric(bmi),
      Socioeconomic = as.numeric(Socioeconomic),
      Qualification = as.numeric(Qualification),
      sex = factor(sex),
      ethn = factor(ethn),
      Smoking_status = factor(Smoking_status),
      Alcohol_consumption = factor(Alcohol_consumption)
    ) %>%
    as.data.frame()

  dat[, meta_cols] <- lapply(dat[, meta_cols, drop = FALSE], function(x) {
    as.numeric(scale(as.numeric(x)))
  })

  dat <- fastDummies::dummy_cols(
    dat,
    select_columns = c("sex", "ethn", "Smoking_status", "Alcohol_consumption"),
    remove_first_dummy = TRUE,
    remove_selected_columns = TRUE,
    ignore_na = TRUE
  )

  dat
}

get_sem_covars <- function(dat) {
  c(
    "age", "Qualification", "bmi", "Socioeconomic",
    grep("^(sex_|ethn_|Smoking_status_|Alcohol_consumption_)",
         names(dat), value = TRUE)
  ) %>%
    unique() %>%
    intersect(names(dat))
}

remove_invalid_covars <- function(df, covars) {
  covars[vapply(covars, function(v) {
    x <- df[[v]]
    length(unique(x[!is.na(x)])) > 1
  }, logical(1))]
}

build_sem_model <- function(exposure, meta_cols, outcome, covars) {
  measurement_part <- paste0(
    "MetFactor =~ ",
    paste(meta_cols, collapse = " + ")
  )

  covar_rhs <- if (length(covars) > 0) {
    paste0(" + ", paste(covars, collapse = " + "))
  } else {
    ""
  }

  structural_part <- paste0(
    "\n",
    "MetFactor ~ a*", exposure, covar_rhs, "\n",
    outcome, " ~ b*MetFactor + cprime*", exposure, covar_rhs, "\n\n",
    "IE := a*b\n",
    "DE := cprime\n",
    "TE := cprime + (a*b)\n"
  )

  paste(measurement_part, structural_part, sep = "\n")
}

run_sem_latent <- function(dat, exposure, meta_cols,
                           outcome = "migraine_incident",
                           covars,
                           estimator = "WLSMV") {
  use_vars <- unique(c(exposure, meta_cols, outcome, covars))

  df <- dat %>%
    select(all_of(use_vars)) %>%
    filter(complete.cases(.)) %>%
    as.data.frame()

  if (nrow(df) < 100 || length(unique(df[[outcome]])) < 2) {
    warning("Skipping ", exposure, ": insufficient complete cases or outcome variation.")
    return(NULL)
  }

  df[[exposure]] <- as.numeric(scale(as.numeric(df[[exposure]])))
  valid_covars <- remove_invalid_covars(df, covars)

  model <- build_sem_model(
    exposure = exposure,
    meta_cols = meta_cols,
    outcome = outcome,
    covars = valid_covars
  )

  fit <- lavaan::sem(
    model,
    data = df,
    ordered = outcome,
    estimator = estimator,
    std.lv = TRUE
  )

  pe <- lavaan::parameterEstimates(fit, standardized = TRUE) %>%
    as.data.frame()

  loadings <- pe %>%
    filter(op == "=~", lhs == "MetFactor") %>%
    transmute(
      Exposure = exposure,
      Metabolite_ID = rhs,
      Loading = est,
      SE = se,
      Z = z,
      P_value = pvalue,
      Std_loading = std.all,
      N = nrow(df)
    )

  effects <- pe %>%
    filter(
      (op == "~" & lhs %in% c("MetFactor", outcome)) |
        (op == ":=" & lhs %in% c("IE", "DE", "TE"))
    ) %>%
    mutate(
      Exposure = exposure,
      N = nrow(df)
    )

  fit_indices <- data.frame(
    Exposure = exposure,
    N = nrow(df),
    Chisq = lavaan::fitMeasures(fit, "chisq"),
    DF = lavaan::fitMeasures(fit, "df"),
    P_value = lavaan::fitMeasures(fit, "pvalue"),
    CFI = lavaan::fitMeasures(fit, "cfi"),
    TLI = lavaan::fitMeasures(fit, "tli"),
    RMSEA = lavaan::fitMeasures(fit, "rmsea"),
    SRMR = lavaan::fitMeasures(fit, "srmr")
  )

  list(
    exposure = exposure,
    model = model,
    fit = fit,
    loadings = loadings,
    effects = effects,
    fit_indices = fit_indices
  )
}

format_sem_tables <- function(sem_results, exposure_labels, meta_name_df) {
  valid_results <- compact(sem_results)

  loading_table <- bind_rows(map(valid_results, "loadings")) %>%
    left_join(meta_name_df, by = "Metabolite_ID") %>%
    mutate(
      Exposure = recode(Exposure, !!!as.list(exposure_labels)),
      Loading_95CI = sprintf(
        "%.4f (%.4f, %.4f)",
        Loading,
        Loading - 1.96 * SE,
        Loading + 1.96 * SE
      ),
      Direction = case_when(
        Std_loading > 0 ~ "Positive",
        Std_loading < 0 ~ "Negative",
        TRUE ~ "Neutral"
      )
    ) %>%
    select(
      Exposure, Metabolite_ID, Metabolite_name,
      Loading, SE, Z, P_value, Std_loading,
      Loading_95CI, Direction, N
    )

  effect_table <- bind_rows(map(valid_results, "effects")) %>%
    mutate(Exposure = recode(Exposure, !!!as.list(exposure_labels)))

  fit_table <- bind_rows(map(valid_results, "fit_indices")) %>%
    mutate(Exposure = recode(Exposure, !!!as.list(exposure_labels)))

  list(
    loadings = loading_table,
    effects = effect_table,
    fit_indices = fit_table
  )
}

# -----------------------------
# Main analysis
# -----------------------------
required_vars <- c(
  "incident_migraine", "prevalent_depression", "prevalent_anxiety",
  "insomnia", "sleep duration", lasso_meta_cols, base_covars
)
check_required_vars(dat_base, required_vars)

meta_name_df <- make_metabolite_map(lasso_meta_cols)
dat_sem <- prepare_sem_data(dat_base, lasso_meta_cols)
covars_sem <- get_sem_covars(dat_sem)

sem_results <- lapply(names(exposures), function(expo) {
  tryCatch(
    run_sem_latent(
      dat = dat_sem,
      exposure = expo,
      meta_cols = lasso_meta_cols,
      outcome = "migraine_incident",
      covars = covars_sem,
      estimator = "WLSMV"
    ),
    error = function(e) {
      warning("SEM failed for ", expo, ": ", conditionMessage(e))
      NULL
    }
  )
})
names(sem_results) <- names(exposures)

sem_tables <- format_sem_tables(sem_results, exposures, meta_name_df)

# -----------------------------
# Export results
# -----------------------------
fwrite(
  sem_tables$loadings,
  file.path(result_dir, "sem_latent_metabolite_factor_loadings.csv")
)

fwrite(
  sem_tables$effects,
  file.path(result_dir, "sem_model_based_association_parameters.csv")
)

fwrite(
  sem_tables$fit_indices,
  file.path(result_dir, "sem_fit_indices.csv")
)

saveRDS(
  sem_results,
  file.path(result_dir, "sem_latent_factor_results.rds")
)

message("SEM analysis completed. Results saved to: ", normalizePath(result_dir))
