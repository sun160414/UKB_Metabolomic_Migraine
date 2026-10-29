#!/usr/bin/env Rscript
# Rscript --vanilla Protocol_demo/tests/test_completeness.R
# Base R only. Sourcing the runner defines functions without running the demo.

test_file <- local({
  paths <- unlist(lapply(sys.frames(), function(x) x$ofile), use.names = FALSE)
  if (length(paths)) return(normalizePath(tail(paths, 1L), mustWork = TRUE))
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (!length(arg)) stop("Run this file with Rscript or source().")
  normalizePath(sub("^--file=", "", arg[1L]), mustWork = TRUE)
})
demo_dir <- dirname(dirname(test_file))
env <- new.env(parent = globalenv())
source(file.path(demo_dir, "run_synthetic_example.R"), local = env)
summarize <- env$summarize_z_completeness

n_checks <- 0L
expect_true <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAILED: ", label, call. = FALSE)
  n_checks <<- n_checks + 1L
  cat(sprintf("ok %02d - %s\n", n_checks, label))
}
expect_error <- function(expr, label) {
  result <- tryCatch({ force(expr); NULL }, error = identity)
  expect_true(inherits(result, "error"), label)
}
fixture <- function(fraction) data.frame(
  Metabolite_raw = paste0("Met", seq_along(fraction)), valid_fraction = fraction)

all_complete <- summarize(fixture(c(1, 1, 1)))
expect_true(all_complete$n_metabolites == 3L &&
  all_complete$n_complete_metabolites == 3L && all_complete$percent_complete == 100,
  "all-complete metabolites produce 3/3 and 100 percent")
expect_true(all_complete$n_unevaluable_metabolites == 0L &&
  all_complete$diagnostic_status == "computed",
  "a fully known summary is evaluable")

partial <- summarize(fixture(c(1, 0.5, 0, 1)))
expect_true(partial$n_metabolites == 4L && partial$n_complete_metabolites == 2L &&
  partial$percent_complete == 50,
  "only metabolites valid in every matched set enter the numerator")
expect_true(partial$percent_complete != 62.5,
  "Figure 1D is not the average per-metabolite valid fraction")
none_complete <- summarize(fixture(c(0.5, 0.9)))
expect_true(none_complete$n_complete_metabolites == 0L &&
  none_complete$percent_complete == 0,
  "partially observed metabolites do not count as complete")

unknown <- summarize(fixture(c(1, NA_real_)))
expect_true(unknown$n_metabolites == 2L && unknown$n_complete_metabolites == 1L &&
  unknown$n_unevaluable_metabolites == 1L && is.na(unknown$percent_complete),
  "a missing fraction is retained and cannot turn one known complete metabolite into 100 percent")
nonfinite <- summarize(fixture(c(1, Inf)))
expect_true(is.na(nonfinite$percent_complete) &&
  nonfinite$diagnostic_status == "unevaluable_completeness",
  "a nonfinite fraction makes the aggregate explicitly unevaluable")
expect_error(summarize(data.frame(Metabolite_raw = character(), valid_fraction = numeric())),
  "an empty metabolite collection fails explicitly")
expect_error(summarize(fixture(c(1, 1.1))),
  "out-of-range completeness fractions are rejected")

cat(sprintf("PASS: %d completeness checks.\n", n_checks))
