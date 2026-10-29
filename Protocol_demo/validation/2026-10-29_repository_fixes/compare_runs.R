#!/usr/bin/env Rscript
# Read-only data comparisons, apart from output_csv. No participant records are
# included in the resulting report. Run all three analyses in the same R setup.
# Rscript --vanilla compare_runs.R baseline revised repeated output_csv
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("Expected: baseline revised repeated output_csv")
baseline <- args[1L]; revised <- args[2L]; repeated <- args[3L]
if (!all(dir.exists(c(baseline, revised, repeated)))) stop("A run directory is missing")
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
compare_files <- function(left, right, comparison, extension, exclude = character(),
                          expected_new = character()) {
  pattern <- paste0("[.]", extension, "$")
  a <- list.files(left, pattern = pattern, recursive = TRUE)
  b <- list.files(right, pattern = pattern, recursive = TRUE)
  record(comparison, paste0("*.", extension), "expected_file_inventory",
         setequal(b, union(a, expected_new)))
  for (f in setdiff(a, exclude)) {
    record(comparison, f, paste0("identical(parsed_", extension, ")"),
           f %in% b && identical(read_data(file.path(left, f)),
                                 read_data(file.path(right, f))))
  }
}
matchit_file <- "analysis/data_processed/matchit_object.rds"
compare_files(baseline, revised, "baseline_vs_revised", "csv",
  exclude = "analysis/source_md5.csv",
  expected_new = "analysis/tables/figure1D_z_completeness_source.csv")
compare_files(baseline, revised, "baseline_vs_revised", "rds",
  exclude = matchit_file,
  expected_new = "analysis/data_processed/figure1D_z_completeness_source.rds")
compare_files(revised, repeated, "revised_vs_repeated", "csv")
compare_files(revised, repeated, "revised_vs_repeated", "rds", exclude = matchit_file)

# Formula environments differ between fresh R sessions. Check the complete
# matching allocation and substantive MatchIt input/output data explicitly.
components <- c("match.matrix", "subclass", "weights", "discarded", "treat", "X")
for (comparison in c("baseline_vs_revised", "revised_vs_repeated")) {
  left <- if (comparison == "baseline_vs_revised") baseline else revised
  right <- if (comparison == "baseline_vs_revised") revised else repeated
  a <- readRDS(file.path(left, matchit_file))
  b <- readRDS(file.path(right, matchit_file))
  for (component in components) {
    record(comparison, paste0(matchit_file, "#", component),
           "identical(MatchIt_data_component)",
           component %in% names(a) && component %in% names(b) &&
             identical(a[[component]], b[[component]]))
  }
}

# Calculate the aggregate independently from matched-set Z-scores. This does
# not call the revised completeness helper or trust the exported fractions.
z <- readRDS(file.path(revised, "analysis/data_processed/case_metabolite_z_scores_long.rds"))
by_metabolite <- split(z$Estimate, z$Metabolite_raw)
complete <- vapply(by_metabolite, function(x) all(is.finite(x)), logical(1))
overall <- readRDS(file.path(revised, "analysis/data_processed/figure1D_z_completeness_source.rds"))
record("revised_independent_formula", "analysis/data_processed/figure1D_z_completeness_source.rds",
       "aggregate_recomputed_from_matched_set_Z_scores",
       nrow(overall) == 1L && overall$n_metabolites == length(complete) &&
         overall$n_complete_metabolites == sum(complete) &&
         overall$n_unevaluable_metabolites == 0L &&
         overall$percent_complete == 100 * sum(complete) / length(complete) &&
         overall$diagnostic_status == "computed")

# Existing figure outputs remain; only the aggregate Figure 1D is replaced,
# while the old per-metabolite view gains its own diagnostic filename.
figure_files <- function(path) list.files(file.path(path, "analysis", "figures"), pattern = "[.]png$")
record("baseline_vs_revised", "analysis/figures/*.png", "expected_figure_inventory",
       setequal(figure_files(revised), union(figure_files(baseline), "diagnostic_z_completeness.png")))
record("revised_vs_repeated", "analysis/figures/*.png", "expected_figure_inventory",
       setequal(figure_files(revised), figure_files(repeated)))
png_signature <- as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
for (run_name in c("baseline", "revised", "repeated")) {
  run_path <- switch(run_name, baseline = baseline, revised = revised, repeated = repeated)
  for (f in figure_files(run_path)) {
    path <- file.path(run_path, "analysis", "figures", f)
    record(paste0(run_name, "_figure_export"), paste0("analysis/figures/", f),
           "nonempty_PNG_with_signature",
           file.info(path)$size > 8L &&
             identical(readBin(path, what = "raw", n = 8L), png_signature))
  }
}
checks <- do.call(rbind, records)
write.csv(checks, args[4L], row.names = FALSE)
print(with(checks, table(comparison, passed)))
cat("TOTAL:", nrow(checks), "PASS:", sum(checks$passed), "FAIL:", sum(!checks$passed), "\n")
if (!all(checks$passed)) {
  print(checks[!checks$passed, ])
  quit(save = "no", status = 1L)
}
