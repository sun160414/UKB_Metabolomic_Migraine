#!/usr/bin/env Rscript
# Read-only numerical regression/reproducibility checks, apart from output_csv.
# Rscript --vanilla compare_runs.R run_1 run_2 baseline output_csv
# All paths in the resulting CSV are relative; no participant data are exported.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Expected: run_1 run_2 baseline output_csv")
run_1 <- args[1]; run_2 <- args[2]; baseline <- args[3]; output_csv <- args[4]
if (!all(dir.exists(c(run_1, run_2, baseline)))) stop("A run directory is missing")
records <- list()
record <- function(comparison, relative_file, method, passed) {
  records[[length(records) + 1L]] <<- data.frame(
    comparison = comparison, relative_file = relative_file,
    comparison_method = method, passed = isTRUE(passed), stringsAsFactors = FALSE)
}
read_data <- function(path) {
  if (grepl("[.]rds$", path)) readRDS(path) else
    read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}
compare_files <- function(left, right, comparison, pattern, exclude = character(),
                          expected_new = character()) {
  a <- list.files(left, pattern = pattern, recursive = TRUE)
  b <- list.files(right, pattern = pattern, recursive = TRUE)
  extension <- if (pattern == "[.]csv$") "csv" else "rds"
  record(comparison, paste0("*.", extension), "expected_file_inventory",
         setequal(b, union(a, expected_new)))
  for (f in setdiff(a, exclude)) {
    ok <- f %in% b && identical(read_data(file.path(left, f)),
                               read_data(file.path(right, f)))
    record(comparison, f, paste0("identical(parsed_", extension, ")"), ok)
  }
}
compare_files(baseline, run_1, "baseline_vs_run_1", "[.]csv$",
  exclude = c("analysis/tables/covariate_balance.csv", "analysis/tables/run_summary.csv"),
  expected_new = c("analysis/source_md5.csv", "analysis/tables/covariate_balance_matchit.csv",
                   "analysis/tables/covariate_balance_plot.csv"))
compare_files(baseline, run_1, "baseline_vs_run_1", "[.]rds$",
  exclude = c("analysis/data_processed/covariate_balance.rds",
              "analysis/data_processed/matched_quality_control.rds",
              "analysis/data_processed/run_summary.rds",
              "analysis/data_processed/matchit_object.rds"),
  expected_new = c("analysis/data_processed/covariate_balance_matchit.rds",
                   "analysis/data_processed/covariate_balance_plot.rds"))
compare_files(run_1, run_2, "run_1_vs_run_2", "[.]csv$")
compare_files(run_1, run_2, "run_1_vs_run_2", "[.]rds$",
  exclude = "analysis/data_processed/matchit_object.rds")

# matchit objects contain formula environments; compare their substantive data
# components explicitly instead of serialized bytes or R environment identity.
matchit_file <- "analysis/data_processed/matchit_object.rds"
components <- c("match.matrix", "subclass", "weights", "discarded", "treat", "X")
for (comparison in c("baseline_vs_run_1", "run_1_vs_run_2")) {
  left <- if (comparison == "baseline_vs_run_1") baseline else run_1
  right <- if (comparison == "baseline_vs_run_1") run_1 else run_2
  a <- readRDS(file.path(left, matchit_file)); b <- readRDS(file.path(right, matchit_file))
  for (component in components) {
    record(comparison, paste0(matchit_file, "#", component),
           "identical(MatchIt_data_component)",
           component %in% names(a) && component %in% names(b) &&
             identical(a[[component]], b[[component]]))
  }
}

# Verify that the old MatchIt values are preserved despite the column renaming.
a <- readRDS(file.path(baseline, "analysis/data_processed/covariate_balance.rds"))
b <- readRDS(file.path(run_1, "analysis/data_processed/covariate_balance_matchit.rds"))
names(b)[names(b) == "matchit_variable"] <- "variable"
names(b)[names(b) == "abs_smd_matchit"] <- "abs_smd"
names(b)[names(b) == "abs_smd_before_matchit"] <- "abs_smd_before"
record("baseline_vs_run_1", "analysis/data_processed/covariate_balance_matchit.rds",
       "identical(legacy_MatchIt_columns_after_renaming)",
       all(names(b) %in% names(a)) && identical(a[, names(b), drop = FALSE], b))

# The former seven science/count metrics are unchanged; the maximum SMD and
# additional balance diagnostics are expected to change under the new formula.
summary_file <- "analysis/data_processed/run_summary.rds"
a <- readRDS(file.path(baseline, summary_file)); b <- readRDS(file.path(run_1, summary_file))
a <- a[a$metric != "max_postmatch_absolute_smd", , drop = FALSE]
b <- b[match(a$metric, b$metric), , drop = FALSE]
rownames(a) <- rownames(b) <- NULL
record("baseline_vs_run_1", paste0(summary_file, "#original_nonbalance_metrics"),
       "identical(original_nonbalance_metrics)", identical(a, b))

# Matched-set integrity is unchanged although its containing QC object now also
# includes the new balance tables.
qc_file <- "analysis/data_processed/matched_quality_control.rds"
a <- readRDS(file.path(baseline, qc_file)); b <- readRDS(file.path(run_1, qc_file))
record("baseline_vs_run_1", paste0(qc_file, "#set_qc"),
       "identical(matched_set_integrity)", identical(a$set_qc, b$set_qc))

# Recompute the new SMD from the matched participant data without calling the
# runner's balance helper; the two independent formulas should be identical.
md <- readRDS(file.path(run_1, "analysis/data_processed/matched_dataset.rds"))
levels <- readRDS(file.path(run_1, "analysis/data_processed/covariate_balance.rds"))
summary <- readRDS(file.path(run_1, "analysis/data_processed/covariate_balance_plot.rds"))
manual <- vapply(seq_len(nrow(levels)), function(i) {
  x <- md[[levels$variable[i]]]
  if (levels$variable_type[i] == "categorical") x <- as.numeric(x == levels$level[i])
  a <- x[md$status == 1L]; b <- x[md$status == 0L]
  (mean(a) - mean(b)) / sqrt((var(a) + var(b)) / 2)
}, numeric(1))
record("run_1_independent_formula", "analysis/data_processed/covariate_balance.rds#smd",
       "identical(independent_sample_variance_SMD)", identical(manual, levels$smd))
manual_max <- vapply(summary$variable, function(v)
  max(abs(manual[levels$variable == v])), numeric(1))
record("run_1_independent_formula", "analysis/data_processed/covariate_balance_plot.rds#abs_smd",
       "identical(independent_per_covariate_maximum)",
       identical(unname(manual_max), summary$abs_smd))

checks <- do.call(rbind, records)
write.csv(checks, output_csv, row.names = FALSE)
print(with(checks, table(comparison, comparison_method, passed)))
cat("TOTAL:", nrow(checks), "PASS:", sum(checks$passed), "FAIL:", sum(!checks$passed), "\n")
if (!all(checks$passed)) {
  print(checks[!checks$passed, ])
  quit(save = "no", status = 1L)
}
