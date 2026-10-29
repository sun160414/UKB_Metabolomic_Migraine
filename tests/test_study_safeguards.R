#!/usr/bin/env Rscript
# Focused regression checks for the original study scripts.
# Run from the repository root: Rscript --vanilla tests/test_study_safeguards.R
# These checks use artificial data and load function definitions only. They do
# not execute the study scripts' main analyses or require UK Biobank inputs.

required <- c("dplyr", "tidyr", "tibble", "purrr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install test dependencies: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
})

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this test with Rscript --vanilla.")
test_file <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
repo_root <- dirname(dirname(test_file))

load_functions <- function(path, selected) {
  env <- new.env(parent = globalenv())
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% selected &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
      eval(expr, envir = env)
    }
  }
  if (!all(vapply(selected, exists, logical(1), envir = env, inherits = FALSE))) {
    stop("Could not find all requested function definitions in ", path)
  }
  env
}

cox <- load_functions(file.path(repo_root, "1.Cox+linear.R"),
  c("split_sex_strata", "prepare_covariates"))
gradient <- load_functions(file.path(repo_root, "2.Trajectories+clusters.R"),
  c("standardize_covariates", "parse_binary_status", "prepare_matching_data",
    "fit_loess_grid", "cluster_loess_matrix"))

n_checks <- 0L
check <- function(label, condition) {
  if (!isTRUE(condition)) stop("FAIL: ", label, call. = FALSE)
  n_checks <<- n_checks + 1L
  cat("PASS: ", label, "\n", sep = "")
}
expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  inherits(error, "error") && grepl(pattern, conditionMessage(error))
}

# Numeric factor labels must survive covariate preparation. Already-numeric
# values must not round-trip through character strings and lose precision.
precision_values <- c(pi, 1 + .Machine$double.eps)
covariate_data <- data.frame(
  age = factor(c("55", "68"), levels = c("68", "55")),
  sex = c(0L, 1L), ethn = c(1L, 2L),
  Qualification = factor(c("10", "20"), levels = c("20", "10")),
  bmi = precision_values, Socioeconomic = c(-1, 2),
  Smoking_status = c(0L, 1L), Alcohol_consumption = c(0L, 1L),
  `screen time (TV)` = c(2, 3), `screen time (computer)` = c(1, 2),
  `sleep duration` = c(7, 8), diabetes_status = c(0L, 1L),
  CVD_status = c(0L, 1L), check.names = FALSE
)
prepared_covariates <- cox$prepare_covariates(covariate_data)
check("Cox covariate preparation preserves numeric factor labels",
      identical(prepared_covariates$age, c(55, 68)) &&
      identical(prepared_covariates$Qualification, c(10, 20)))
gradient_covariates <- gradient$standardize_covariates(data.frame(
  age = covariate_data$age,
  migraine_years = factor(c("2.25", "10.5"), levels = c("10.5", "2.25")),
  bmi = precision_values
))
check("Gradient covariate preparation preserves age and interval factor labels",
      identical(gradient_covariates$age, c(55, 68)) &&
      identical(gradient_covariates$migraine_years, c(2.25, 10.5)))
check("Both covariate preparers preserve already-numeric precision exactly",
      identical(prepared_covariates$bmi, precision_values) &&
      identical(gradient_covariates$bmi, precision_values))

# Sex-stratum membership must depend on the declared input values, not factor
# ordering, alphabetical labels or the order in which participants appear.
sex_data <- data.frame(
  id = c("f1", "m1", "f2", "m2"),
  sex = factor(c("0", "1", "0", "1"), levels = c("1", "0"))
)
sex_groups <- cox$split_sex_strata(sex_data, c(female = "0", male = "1"))
check("UKB sex codes remain correct with reversed factor levels",
      identical(sex_groups$female$id, c("f1", "f2")) &&
      identical(sex_groups$male$id, c("m1", "m2")))
text_sex <- data.frame(id = c("m", "f"), sex = factor(c("Man", "Woman")))
text_groups <- cox$split_sex_strata(text_sex, c(female = "Woman", male = "Man"))
check("Explicit text-label mapping is honored",
      identical(text_groups$male$id, "m") && identical(text_groups$female$id, "f"))
check("Unmapped sex labels are rejected",
      expect_error(cox$split_sex_strata(text_sex, c(female = "0", male = "1")),
                   "Unmapped sex"))
check("A missing mapped sex group is rejected",
      expect_error(cox$split_sex_strata(sex_data[c(1, 3), ], c(female = "0", male = "1")),
                   "No participants in the mapped sex stratum"))
check("Ambiguous sex-code mappings are rejected",
      expect_error(cox$split_sex_strata(sex_data, c(female = "0", male = "0")),
                   "must be nonmissing, nonempty and distinct"))

# A prevalent participant with a positive observation interval was eligible as a
# control in the previous matching entry point. It must now be removed.
cohort <- data.frame(
  eid = c("case", "prevalent", "control1", "control2"),
  prevalent_migraine = factor(c("0", "1", "0", "0"), levels = c("1", "0")),
  incident_migraine = factor(c("1", "0", "0", "0"), levels = c("1", "0")),
  migraine_years = c(2, 3, 4, 5)
)
prepared <- gradient$prepare_matching_data(cohort, character())
check("Prevalent cases are excluded and binary factor labels are preserved",
      identical(prepared$eid, c("case", "control1", "control2")) &&
      identical(prepared$migraine_status, c(1L, 0L, 0L)))
legacy <- cohort
legacy$migraine_status <- legacy$incident_migraine
legacy$incident_migraine <- NULL
check("The incident-only legacy flag follows the same prevalence exclusion",
      identical(gradient$prepare_matching_data(legacy, character())$eid, prepared$eid))
without_prevalence <- cohort
without_prevalence$prevalent_migraine <- NULL
check("The baseline prevalence field is required",
      expect_error(gradient$prepare_matching_data(without_prevalence, character()),
                   "must contain prevalent_migraine"))
invalid_incident <- cohort
invalid_incident$incident_migraine <- c(1L, 0L, 2L, 0L)
check("Invalid incident codes cannot silently become controls",
      expect_error(gradient$prepare_matching_data(invalid_incident, character()),
                   "incident_migraine must contain only"))
invalid_prevalent <- cohort
invalid_prevalent$prevalent_migraine <- c(0L, 1L, 2L, 0L)
check("Invalid prevalent codes are rejected",
      expect_error(gradient$prepare_matching_data(invalid_prevalent, character()),
                   "prevalent_migraine must contain only"))
contradictory <- cohort
contradictory$incident_migraine <- c(1L, 1L, 0L, 0L)
check("A case cannot be simultaneously prevalent and incident",
      expect_error(gradient$prepare_matching_data(contradictory, character()),
                   "both prevalent and incident"))
disagreeing <- cohort
disagreeing$migraine_status <- c(1L, 0L, 1L, 0L)
check("Two supplied incident flags must agree",
      expect_error(gradient$prepare_matching_data(disagreeing, character()),
                   "incident_migraine and migraine_status disagree"))
missing_flags <- cohort
missing_flags$prevalent_migraine <- c(0L, 1L, NA_integer_, 0L)
missing_flags$incident_migraine <- c(1L, 0L, 0L, NA_integer_)
status_warnings <- character()
missing_result <- withCallingHandlers(
  gradient$prepare_matching_data(missing_flags, character()),
  warning = function(w) {
    status_warnings <<- c(status_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
check("Missing outcome flags are explicitly reported and excluded",
      identical(missing_result$eid, "case") &&
      any(grepl("Excluded 2 participants with missing", status_warnings)))

# Noninteger observed endpoints reproduce the former outward-grid failure.
year <- seq(-14.83, -0.11, length.out = 120L)
z <- data.frame(Year = year,
                Meta_1 = sin(year / 3) + 0.02 * cos(year),
                Meta_2 = cos(year / 4), Meta_3 = 0.03 * year - 0.2)
metabolites <- c("Meta_1", "Meta_2", "Meta_3")
outward <- seq(floor(min(year)), ceiling(max(year)), by = 0.2)
old_prediction <- predict(loess(Meta_1 ~ Year, data = z, span = 0.75),
                          newdata = data.frame(Year = outward))
check("The regression fixture reproduces unsupported outward-grid predictions",
      anyNA(old_prediction))
pred <- gradient$fit_loess_grid(z, metabolites)
check("The revised common grid stays within observed support",
      min(pred$Year) >= min(year) && max(pred$Year) <= max(year) &&
      length(unique(pred$Year)) > 2L)
check("All fully supported artificial curves remain finite",
      all(is.finite(pred$Estimate_loess)))
clusters <- gradient$cluster_loess_matrix(pred, k = 3L)
check("Supported curves survive to the clustering step",
      nrow(clusters$matrix) == 3L && all(is.finite(clusters$matrix)))
nonfinite <- rbind(z, data.frame(Year = c(Inf, -Inf, NA_real_, -99),
                               Meta_1 = c(1, 1, 1, Inf),
                               Meta_2 = c(1, 1, 1, NaN),
                               Meta_3 = c(1, 1, 1, NA_real_)))
finite_pred <- gradient$fit_loess_grid(nonfinite, metabolites)
check("Nonfinite observations cannot extend the prediction grid",
      isTRUE(all.equal(pred, finite_pred, check.attributes = FALSE)))
check("An all-invalid LOESS input fails clearly",
      expect_error(gradient$fit_loess_grid(
        data.frame(Year = c(NA_real_, Inf), Meta_1 = c(1, 2)), "Meta_1"),
        "No finite matched-set Z-scores"))
check("A nonpositive grid step is rejected",
      expect_error(gradient$fit_loess_grid(z, metabolites, step = 0),
                   "finite positive number"))
narrow <- z
narrow$Year <- seq(-0.11, -0.01, length.out = nrow(narrow))
check("Insufficient support for the requested grid fails clearly",
      expect_error(gradient$fit_loess_grid(narrow, metabolites, step = 0.2),
                   "Insufficient observed time range"))
partial <- z
partial$Meta_3[partial$Year > -1] <- NA_real_
partial_pred <- gradient$fit_loess_grid(partial, metabolites)
unsupported <- partial_pred$Metabolite_raw == "Meta_3" &
  partial_pred$Year > max(partial$Year[is.finite(partial$Meta_3)])
check("A metabolite with shorter support is not extrapolated",
      any(unsupported) && all(is.na(partial_pred$Estimate_loess[unsupported])))
nonfinite_curve <- pred
nonfinite_curve$Estimate_loess[which(nonfinite_curve$Metabolite_raw == "Meta_1")[1]] <- Inf
finite_clusters <- gradient$cluster_loess_matrix(nonfinite_curve, k = 2L)
check("Nonfinite prediction vectors cannot enter clustering",
      nrow(finite_clusters$matrix) == 2L &&
      !"Meta_1" %in% rownames(finite_clusters$matrix))

cat("All ", n_checks, " study-safeguard checks passed.\n", sep = "")
