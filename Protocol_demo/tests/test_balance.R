#!/usr/bin/env Rscript
# Run from any working directory:
# Rscript --vanilla Protocol_demo/tests/test_balance.R
# Uses base R assertions and the demo's existing ggplot2 dependency only.

test_file <- local({
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  if (length(paths)) return(normalizePath(tail(paths, 1L), mustWork = TRUE))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(arg)) stop("Run this file with Rscript or source().")
  normalizePath(sub("^--file=", "", arg[1L]), mustWork = TRUE)
})
demo_dir <- dirname(dirname(test_file))
runner <- file.path(demo_dir, "run_synthetic_example.R")
default_output <- file.path(demo_dir, "synthetic_demo_output")
output_existed <- dir.exists(default_output)
env <- new.env(parent = globalenv())
source(runner, local = env)

n_checks <- 0L
expect_true <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAILED: ", label, call. = FALSE)
  n_checks <<- n_checks + 1L
  cat(sprintf("ok %02d - %s\n", n_checks, label))
}
expect_equal <- function(actual, expected, label, tolerance = 1e-12) {
  expect_true(isTRUE(all.equal(actual, expected, tolerance = tolerance,
                               check.attributes = FALSE)), label)
}
expect_error <- function(expr, label) {
  result <- tryCatch({ force(expr); NULL }, error = identity)
  expect_true(inherits(result, "error"), label)
}
balance <- env$matched_covariate_balance
strict_max <- env$strict_balance_max
plot_balance <- env$plot_covariate_balance
expect_true(is.function(balance) && is.function(strict_max) &&
              is.function(plot_balance), "balance helpers are available after sourcing")
expect_true(identical(dir.exists(default_output), output_existed),
            "sourcing the runner does not execute the default pipeline")
expect_true("education" %in% env$protocol_config()$factor_covariates,
            "default synthetic education remains explicitly categorical")

# Independent calculation: the SD denominator is the unweighted mean of
# the two matched groups' sample variances, not MatchIt's default denominator.
manual_smd <- function(x, status) {
  a <- x[status == 1L]
  b <- x[status == 0L]
  (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
}
d <- data.frame(
  status = c(rep(1L, 4L), rep(0L, 8L)),
  age = c(51, 55, 58, 62, 45, 47, 50, 54, 56, 59, 61, 64),
  education = c(11, 13, 14, 17, 8, 10, 11, 12, 13, 14, 15, 16),
  group = factor(c("A", "A", "A", "B", "A", "A", "B", "B",
                   "B", "C", "C", "C"), levels = c("A", "B", "C", "unused")),
  sex = factor(c("F", "M", "F", "M", "F", "M", "F", "M",
                 "F", "M", "F", "M"))
)
variables <- c("age", "education", "group", "sex")
b <- balance(d, variables)
required_levels <- c("variable", "variable_type", "level", "n_case", "n_control",
                     "mean_case", "mean_control", "var_case", "var_control",
                     "denominator", "smd", "abs_smd", "evaluable",
                     "diagnostic_status", "above_threshold")
required_summary <- c("variable", "variable_type", "n_levels", "n_unevaluable",
                      "abs_smd", "evaluable", "above_threshold", "diagnostic_status")
expect_true(all(required_levels %in% names(b$levels)) &&
              all(required_summary %in% names(b$summary)),
            "complete level and plotting-summary schemas are exported")
expect_true(nrow(b$summary) == length(variables) &&
              setequal(b$summary$variable, variables) && !anyDuplicated(b$summary$variable),
            "the plotting summary contains exactly one row per covariate")
expect_true(nrow(b$levels) == 7L && !"unused" %in% b$levels$level,
            "only observed categorical indicator levels are included")
expect_true(all(b$levels$n_case == 4L) && all(b$levels$n_control == 8L),
            "group sizes are recorded without treating the unequal ratio as weights")

for (v in c("age", "education")) {
  row <- b$levels[b$levels$variable == v, , drop = FALSE]
  expect_true(nrow(row) == 1L && row$variable_type == "numeric" && is.na(row$level),
              paste(v, "is retained as one numeric covariate"))
  expect_equal(row$smd, manual_smd(d[[v]], d$status), paste(v, "SMD matches the explicit formula"))
  expect_equal(row$denominator,
               sqrt((var(d[[v]][d$status == 1L]) + var(d[[v]][d$status == 0L])) / 2),
               paste(v, "records the matched-group SD denominator"))
}
for (level in c("A", "B", "C")) {
  row <- b$levels[b$levels$variable == "group" & b$levels$level == level &
                    !is.na(b$levels$level), , drop = FALSE]
  expect_equal(row$smd, manual_smd(as.numeric(d$group == level), d$status),
               paste("factor indicator", level, "uses the same sample-variance formula"))
}
expected_group <- max(abs(vapply(c("A", "B", "C"), function(level)
  manual_smd(as.numeric(d$group == level), d$status), numeric(1))))
expect_equal(b$summary$abs_smd[b$summary$variable == "group"], expected_group,
             "categorical plotting value is the maximum absolute indicator SMD")
expect_true(all(b$summary$evaluable) && all(b$levels$evaluable),
            "all valid fixture diagnostics are evaluable")
expect_equal(b$summary$abs_smd[b$summary$variable == "sex"], 0,
             "an exactly balanced factor has a zero summary SMD")
boundary <- balance(d, "age", threshold = abs(manual_smd(d$age, d$status)))
expect_true(boundary$levels$above_threshold && boundary$summary$above_threshold,
            "a value equal to the diagnostic threshold is flagged")

# Relabeling a factor must not change the numeric results or use integer codes.
renamed <- d
renamed$group <- factor(c("third", "third", "third", "first", "third", "third",
                          "first", "first", "first", "second", "second", "second"),
                        levels = c("first", "second", "third"))
renamed_balance <- balance(renamed, "group")
expect_equal(renamed_balance$summary$abs_smd, expected_group,
             "factor relabeling leaves the categorical maximum unchanged")

# Undefined diagnostics must remain NA, with a reason and no implied pass.
constant <- data.frame(status = c(1L, 1L, 0L, 0L), same = rep(5, 4),
                       different = c(1, 1, 0, 0))
cb <- balance(constant, c("same", "different"))
expect_equal(cb$levels$smd[cb$levels$variable == "same"], 0,
             "identical constants in both groups yield zero")
bad_row <- cb$levels[cb$levels$variable == "different", , drop = FALSE]
expect_true(is.na(bad_row$smd) && !bad_row$evaluable && is.na(bad_row$above_threshold) &&
              !is.na(bad_row$diagnostic_status) && nzchar(bad_row$diagnostic_status),
            "different constants are undefined rather than zero or a balance pass")
for (invalid in c(NA_real_, Inf, -Inf, NaN)) {
  bad <- d
  bad$age[1L] <- invalid
  rb <- balance(bad, "age")
  expect_true(is.na(rb$levels$smd) && !rb$levels$evaluable &&
                is.na(rb$summary$abs_smd) && is.na(rb$summary$above_threshold),
              paste("nonfinite input", format(invalid), "is not silently omitted"))
}
few <- balance(d[c(1L, 5L, 6L), ], "age")
expect_true(!few$levels$evaluable && is.na(few$summary$abs_smd),
            "fewer than two observations in one group is unevaluable")
missing_factor <- data.frame(status = c(1L, 1L, 0L, 0L),
                             group = factor(rep(NA_character_, 4L)))
mb <- balance(missing_factor, "group")
expect_true(nrow(mb$levels) >= 1L && nrow(mb$summary) == 1L &&
              !mb$summary$evaluable && is.na(mb$summary$abs_smd),
            "an all-missing factor is retained as an unevaluable diagnostic")
partial <- data.frame(status = c(1L, 1L, 0L, 0L, 0L, 0L),
                      group = factor(c("A", "A", "B", "C", "B", "C")))
pb <- balance(partial, "group")
expect_true(any(pb$levels$evaluable) && any(!pb$levels$evaluable),
            "the partial-factor fixture includes both finite and undefined levels")
expect_true(is.na(pb$summary$abs_smd) && !pb$summary$evaluable &&
              pb$summary$n_unevaluable == sum(!pb$levels$evaluable) &&
              is.na(pb$summary$above_threshold),
            "an undefined factor level prevents a misleading partial maximum")

expect_equal(strict_max(c(0, .1, .4)), .4, "the complete finite maximum is computed")
for (x in list(numeric(), NA_real_, c(.1, NA_real_), c(.1, Inf), c(.1, NaN))) {
  expect_true(is.na(strict_max(x)), "incomplete or empty overall maxima remain NA")
}
for (threshold in list(0, -1, NA_real_, Inf, numeric(), c(.1, .2), "invalid")) {
  expect_error(balance(d, "age", threshold), "invalid thresholds are rejected")
}
expect_error(balance(d, "absent"), "unknown covariates are rejected")
bad <- d
bad$status[1L] <- 2L
expect_error(balance(bad, "age"), "nonbinary status is rejected")
bad$status[1L] <- NA_integer_
expect_error(balance(bad, "age"), "missing status is rejected")
bad <- d
bad$group <- as.character(bad$group)
expect_error(balance(bad, "group"), "unprepared character covariates are rejected")

# Plot checks inspect ggplot objects, so the tests create no figure files.
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install the demo dependency ggplot2 first.")
plot_data <- rbind(b$summary, cb$summary[cb$summary$variable == "different", , drop = FALSE])
cfg <- env$protocol_config()
cfg$plot_title_prefix <- "SYNTHETIC TEST:"
p <- plot_balance(plot_data, cfg)
built <- ggplot2::ggplot_build(p)
point_layers <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))
text_layers <- which(vapply(p$layers, function(layer) inherits(layer$geom, "GeomText"), logical(1)))
expect_true(length(point_layers) == 1L &&
              nrow(built$data[[point_layers]]) == sum(plot_data$evaluable),
            "only evaluable covariates receive points")
expect_true(setequal(as.character(built$layout$panel_scales_y[[1L]]$get_limits()),
                     as.character(plot_data$variable)),
            "unevaluable covariates remain on the plot axis")
labels <- unlist(lapply(text_layers, function(i) built$data[[i]]$label), use.names = FALSE)
expect_true(any(labels == "Not evaluable"),
            "undefined SMDs receive an explicit plot annotation")
all_missing <- cb$summary[cb$summary$variable == "different", , drop = FALSE]
p_na <- plot_balance(all_missing, cfg)
build_na <- ggplot2::ggplot_build(p_na)
points_na <- which(vapply(p_na$layers, function(layer) inherits(layer$geom, "GeomPoint"), logical(1)))
expect_true(length(points_na) == 1L && nrow(build_na$data[[points_na]]) == 0L,
            "an entirely unevaluable plot contains no false zero point")
expect_true(setequal(as.character(build_na$layout$panel_scales_y[[1L]]$get_limits()),
                     as.character(all_missing$variable)),
            "an entirely unevaluable plot still retains its covariate row")

cat(sprintf("PASS: %d balance regression checks.\n", n_checks))
