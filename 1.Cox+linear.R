#!/usr/bin/env Rscript

# 01_cox_linear.R
# Association analyses between plasma metabolites and migraine
#
# Required objects:
#   dat: data frame containing phenotype, covariates, and metabolite columns named Meta_*
#   metabolome_name_map: data frame with columns Meta and Original_Metabolite
#
# Example:
#   dat <- readRDS("data/analysis_dataset.rds")
#   metabolome_name_map <- data.table::fread("data/metabolite_name_map.csv")
#   source("scripts/01_cox_linear.R")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(purrr)
  library(survival)
  library(ggplot2)
  library(ggrepel)
})

# ---------------------------
# Configuration
# ---------------------------

output_dir <- "results"
figure_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

covariates <- c(
  "age", "sex", "ethn", "Qualification", "bmi", "Socioeconomic",
  "Smoking_status", "Alcohol_consumption",
  "screen time (TV)", "screen time (computer)", "sleep duration",
  "diabetes_status", "CVD_status"
)

# ---------------------------
# Helper functions
# ---------------------------

check_required_columns <- function(data, cols, object_name = "dat") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}

quote_if_needed <- function(x) {
  ifelse(make.names(x) != x, paste0("`", x, "`"), x)
}

prepare_covariates <- function(data) {
  data %>%
    mutate(
      age = as.numeric(age),
      sex = factor(sex),
      ethn = factor(ethn),
      Qualification = as.numeric(Qualification),
      bmi = as.numeric(bmi),
      Socioeconomic = as.numeric(Socioeconomic),
      Smoking_status = factor(Smoking_status),
      Alcohol_consumption = factor(Alcohol_consumption),
      `screen time (TV)` = as.numeric(`screen time (TV)`),
      `screen time (computer)` = as.numeric(`screen time (computer)`),
      `sleep duration` = as.numeric(`sleep duration`),
      diabetes_status = factor(diabetes_status),
      CVD_status = factor(CVD_status)
    )
}

valid_covariates <- function(data, covars) {
  keep <- vapply(covars, function(v) {
    x <- data[[v]]
    if (is.factor(x) || is.character(x)) {
      nlevels(droplevels(factor(x))) >= 2
    } else {
      length(unique(x[!is.na(x)])) >= 2
    }
  }, logical(1))
  covars[keep]
}

safe_scale <- function(x) {
  x <- as.numeric(x)
  sx <- stats::sd(x, na.rm = TRUE)
  if (is.na(sx) || sx == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / sx
}

format_results <- function(res, p_col, metabolite_map = NULL) {
  out <- res %>%
    mutate(
      p_bonf = p.adjust(.data[[p_col]], method = "bonferroni"),
      p_fdr = p.adjust(.data[[p_col]], method = "BH"),
      Metabolite_raw = trimws(Metabolite)
    )

  if (!is.null(metabolite_map)) {
    reverse_map <- setNames(
      metabolite_map$Original_Metabolite,
      metabolite_map$Meta
    )
    out <- out %>%
      mutate(
        Metabolite_name = unname(reverse_map[Metabolite_raw]),
        Metabolite_label = ifelse(
          is.na(Metabolite_name),
          Metabolite_raw,
          Metabolite_name
        )
      )
  } else {
    out <- out %>%
      mutate(
        Metabolite_name = NA_character_,
        Metabolite_label = Metabolite_raw
      )
  }

  out
}

write_result <- function(x, file_name) {
  fwrite(x, file.path(output_dir, file_name))
}

# ---------------------------
# Cox regression
# ---------------------------

run_cox_scan <- function(data, metabolites, covars, min_n = 50, min_events = 5) {
  map_dfr(metabolites, function(met) {
    model_df <- data %>%
      transmute(
        followup_years = as.numeric(followup_years),
        event = as.integer(incident_migraine == 1),
        Metabolite_value = as.numeric(.data[[met]]),
        across(all_of(covars))
      ) %>%
      mutate(Metabolite_z = safe_scale(Metabolite_value)) %>%
      filter(
        !is.na(followup_years),
        followup_years > 0,
        !is.na(event),
        complete.cases(.)
      )

    n_events <- sum(model_df$event == 1, na.rm = TRUE)

    if (nrow(model_df) < min_n || n_events < min_events) {
      return(tibble(
        Metabolite = met,
        n = nrow(model_df),
        nevent = n_events,
        beta = NA_real_,
        se = NA_real_,
        HR = NA_real_,
        HR_lower = NA_real_,
        HR_upper = NA_real_,
        p = NA_real_,
        error = NA_character_
      ))
    }

    covars_use <- valid_covariates(model_df, covars)
    rhs <- paste(c("Metabolite_z", quote_if_needed(covars_use)), collapse = " + ")
    form <- as.formula(paste0("Surv(followup_years, event) ~ ", rhs))

    fit <- tryCatch(
      coxph(form, data = model_df),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(tibble(
        Metabolite = met,
        n = nrow(model_df),
        nevent = n_events,
        beta = NA_real_,
        se = NA_real_,
        HR = NA_real_,
        HR_lower = NA_real_,
        HR_upper = NA_real_,
        p = NA_real_,
        error = conditionMessage(fit)
      ))
    }

    fit_sum <- summary(fit)
    coef_row <- fit_sum$coefficients["Metabolite_z", ]
    ci_row <- fit_sum$conf.int["Metabolite_z", ]

    tibble(
      Metabolite = met,
      n = fit_sum$n,
      nevent = fit_sum$nevent,
      beta = unname(coef_row["coef"]),
      se = unname(coef_row["se(coef)"]),
      HR = unname(ci_row["exp(coef)"]),
      HR_lower = unname(ci_row["lower .95"]),
      HR_upper = unname(ci_row["upper .95"]),
      p = unname(coef_row["Pr(>|z|)"]),
      error = NA_character_
    )
  })
}

# ---------------------------
# Linear regression
# ---------------------------

run_lm_scan <- function(data, metabolites, covars, min_n = 50) {
  map_dfr(metabolites, function(met) {
    model_df <- data %>%
      transmute(
        prevalent_migraine = as.integer(prevalent_migraine == 1),
        Metabolite_value = as.numeric(.data[[met]]),
        across(all_of(covars))
      ) %>%
      mutate(Metabolite_z = safe_scale(Metabolite_value)) %>%
      filter(!is.na(prevalent_migraine), complete.cases(.))

    if (nrow(model_df) < min_n) {
      return(tibble(
        Metabolite = met,
        n = nrow(model_df),
        beta = NA_real_,
        se = NA_real_,
        beta_lower = NA_real_,
        beta_upper = NA_real_,
        p = NA_real_,
        error = NA_character_
      ))
    }

    covars_use <- valid_covariates(model_df, covars)
    rhs <- paste(c("prevalent_migraine", quote_if_needed(covars_use)), collapse = " + ")
    form <- as.formula(paste0("Metabolite_z ~ ", rhs))

    fit <- tryCatch(
      lm(form, data = model_df),
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(tibble(
        Metabolite = met,
        n = nrow(model_df),
        beta = NA_real_,
        se = NA_real_,
        beta_lower = NA_real_,
        beta_upper = NA_real_,
        p = NA_real_,
        error = conditionMessage(fit)
      ))
    }

    coef_tab <- summary(fit)$coefficients
    beta <- unname(coef_tab["prevalent_migraine", "Estimate"])
    se <- unname(coef_tab["prevalent_migraine", "Std. Error"])

    tibble(
      Metabolite = met,
      n = nrow(model_df),
      beta = beta,
      se = se,
      beta_lower = beta - 1.96 * se,
      beta_upper = beta + 1.96 * se,
      p = unname(coef_tab["prevalent_migraine", "Pr(>|t|)"]),
      error = NA_character_
    )
  })
}

# ---------------------------
# Volcano plots
# ---------------------------

make_volcano <- function(
  data,
  effect_col,
  p_col = "p",
  effect_ref,
  x_lab,
  file_name,
  label_n = 10
) {
  plot_df <- data %>%
    mutate(
      neglog10p = -log10(pmax(.data[[p_col]], .Machine$double.xmin)),
      direction = case_when(
        p_bonf < 0.05 & .data[[effect_col]] > effect_ref ~ "Positive",
        p_bonf < 0.05 & .data[[effect_col]] < effect_ref ~ "Negative",
        TRUE ~ "NS"
      )
    )

  label_df <- bind_rows(
    plot_df %>% filter(direction == "Positive") %>% arrange(p_bonf) %>% slice_head(n = label_n),
    plot_df %>% filter(direction == "Negative") %>% arrange(p_bonf) %>% slice_head(n = label_n)
  )

  p <- ggplot(plot_df, aes(x = .data[[effect_col]], y = neglog10p)) +
    geom_hline(
      yintercept = -log10(0.05 / sum(!is.na(plot_df[[p_col]]))),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_vline(xintercept = effect_ref, linetype = "dashed", linewidth = 0.4) +
    geom_point(aes(color = direction), size = 2.2, alpha = 0.75) +
    scale_color_manual(
      values = c(Positive = "#a73336", Negative = "#333aab", NS = "#6b6b6b")
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = Metabolite_label),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank()
    ) +
    labs(x = x_lab, y = "-log10(P value)")

  ggsave(
    filename = file.path(figure_dir, file_name),
    plot = p,
    width = 5,
    height = 5
  )

  p
}

# ---------------------------
# Main analysis
# ---------------------------

check_required_columns(
  dat,
  c(
    "prevalent_migraine", "incident_migraine", "followup_years",
    covariates
  )
)

if (!exists("metabolome_name_map")) {
  metabolome_name_map <- NULL
} else {
  check_required_columns(
    metabolome_name_map,
    c("Meta", "Original_Metabolite"),
    object_name = "metabolome_name_map"
  )
}

dat <- prepare_covariates(dat)
metabolites <- grep("^Meta_", names(dat), value = TRUE)

if (length(metabolites) == 0) {
  stop("No metabolite columns found. Expected columns starting with 'Meta_'.", call. = FALSE)
}

message("Detected ", length(metabolites), " metabolite columns.")

# Prospective incident migraine dataset
dat_incident <- dat %>%
  filter(prevalent_migraine == 0)

# Overall Cox analysis
cox_overall <- run_cox_scan(dat_incident, metabolites, covariates) %>%
  format_results(p_col = "p", metabolite_map = metabolome_name_map)
write_result(cox_overall, "cox_overall.csv")

# Stratified Cox analyses
cox_strata <- list(
  age_lt55 = list(data = dat_incident %>% filter(age < 55), covars = setdiff(covariates, "age")),
  age_ge55 = list(data = dat_incident %>% filter(age >= 55), covars = setdiff(covariates, "age")),
  male = list(data = dat_incident %>% filter(sex == levels(droplevels(dat_incident$sex))[1]), covars = setdiff(covariates, "sex")),
  female = list(data = dat_incident %>% filter(sex == levels(droplevels(dat_incident$sex))[2]), covars = setdiff(covariates, "sex"))
)

cox_strata_results <- imap(cox_strata, function(x, nm) {
  res <- run_cox_scan(x$data, metabolites, x$covars) %>%
    format_results(p_col = "p", metabolite_map = metabolome_name_map)
  write_result(res, paste0("cox_", nm, ".csv"))
  res
})

# Cross-sectional prevalent migraine analysis
lm_prevalent <- run_lm_scan(dat, metabolites, covariates) %>%
  format_results(p_col = "p", metabolite_map = metabolome_name_map)
write_result(lm_prevalent, "linear_prevalent.csv")

# Volcano plots
make_volcano(
  cox_overall,
  effect_col = "HR",
  effect_ref = 1,
  x_lab = "Hazard ratio",
  file_name = "Fig1A_cox_overall_volcano.pdf"
)

make_volcano(
  cox_strata_results$age_lt55,
  effect_col = "HR",
  effect_ref = 1,
  x_lab = "Hazard ratio",
  file_name = "Fig1B_cox_age_lt55_volcano.pdf"
)

make_volcano(
  cox_strata_results$age_ge55,
  effect_col = "HR",
  effect_ref = 1,
  x_lab = "Hazard ratio",
  file_name = "Fig1C_cox_age_ge55_volcano.pdf"
)

make_volcano(
  cox_strata_results$male,
  effect_col = "HR",
  effect_ref = 1,
  x_lab = "Hazard ratio",
  file_name = "Fig1D_cox_male_volcano.pdf"
)

make_volcano(
  cox_strata_results$female,
  effect_col = "HR",
  effect_ref = 1,
  x_lab = "Hazard ratio",
  file_name = "Fig1E_cox_female_volcano.pdf"
)

make_volcano(
  lm_prevalent,
  effect_col = "beta",
  effect_ref = 0,
  x_lab = "Beta",
  file_name = "Fig1F_linear_prevalent_volcano.pdf"
)

message("Analysis finished.")
message("Results: ", normalizePath(output_dir))
message("Figures: ", normalizePath(figure_dir))
