#!/usr/bin/env Rscript
# Standalone synthetic demonstration of the matched time-to-diagnosis protocol.
# Keep this file beside generate_synthetic_data.R. No private inputs are required.
# CLI: Rscript --vanilla run_synthetic_example.R [output_dir] [n] [p] [seed]
# Source without running: source("run_synthetic_example.R")
# Then run: run_synthetic_example(out_dir = "demo_output")
# Dependencies: MatchIt, survival, ggplot2, mclust, cluster.
# IMPORTANT: read the testing status and scientific limitations in README.md.

# Remember this file's location both under Rscript and under source().
.protocol_demo_dir <- local({
  frames <- sys.frames()
  paths <- lapply(frames, function(x) x$ofile)
  paths <- unlist(paths[!vapply(paths, is.null, logical(1))], use.names = FALSE)
  if (length(paths)) {
    dirname(normalizePath(tail(paths, 1L), mustWork = TRUE))
  } else {
    arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
    if (length(arg)) dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)) else getwd()
  }
})

# ------------------------- Configuration -------------------------------------
# Parameters are centralized here; the analysis functions accept arguments.
# The synthetic endpoint is generic and is NOT a calibrated migraine phenotype.
protocol_config <- function() {
  baseline_covariates <- c(
    "age", "sex", "ethnicity", "education", "deprivation", "bmi",
    "smoking", "alcohol", "tv_time", "computer_time", "sleep", "diabetes", "cvd"
  )
  list(
    n = 5000L, n_metabolites = 20L, seed = 20260801L,
    plot_title_prefix = "SYNTHETIC DEMONSTRATION:",
    columns = c(id = "participant_id", baseline = "baseline_date",
                diagnosis = "diagnosis_date", followup_end = "followup_end"),
    matching_covariates = setdiff(baseline_covariates, "sex"),
    exact_covariates = "sex",
    adjustment_covariates = baseline_covariates,
    cox_covariates = baseline_covariates,
    factor_covariates = c("sex", "ethnicity", "education", "smoking",
                          "alcohol", "diabetes", "cvd"),
    match_ratio = 10L, matching_order = "data", balance_threshold = 0.10,
    min_valid_controls = 2L, min_control_sd = 1e-12,
    denominator_diagnostic_quantile = 0.01,
    loess_span = 0.75, loess_degree = 2L, grid_interval = 0.20,
    min_loess_sets = 50L, min_unique_times = 10L,
    selection_alpha = 0.05, linkage = "ward.D", n_clusters = 3L,
    sensitivity_spans = c(0.50, 1.00), exclude_early_years = 1,
    winsor_probs = c(0.01, 0.99), alternative_linkage = "ward.D2",
    curve_reference = 0.90, ari_reference = 0.80,
    heatmap_display_limit = 1, max_z_plot = 100000L,
    # Additional descriptive diagnostics requested by Reviewer 2; no model refit.
    run_cluster_count_diagnostics = TRUE, silhouette_k = 2:6,
    make_plots = TRUE, plot_dpi = 150L
  )
}

check_packages <- function() {
  required <- c("MatchIt", "survival", "ggplot2", "mclust", "cluster")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop(
    "Install missing packages first:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing), collapse = ", "), "))", call. = FALSE
  )
  if (utils::packageVersion("ggplot2") < "3.4.0")
    stop("Update ggplot2 to version 3.4.0 or later (linewidth is used in plots)")
  invisible(required)
}

check_config <- function(cfg) {
  ints <- c("n", "n_metabolites", "seed", "match_ratio", "min_valid_controls",
            "min_loess_sets", "min_unique_times", "n_clusters", "max_z_plot", "plot_dpi")
  for (nm in ints) {
    x <- cfg[[nm]]
    if (length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x))
      stop(nm, " must be one positive integer")
  }
  if (cfg$min_valid_controls < 2L || cfg$match_ratio < cfg$min_valid_controls)
    stop("Need match_ratio >= min_valid_controls >= 2")
  if (cfg$n_clusters < 2L) stop("n_clusters must be at least 2")
  for (nm in c("loess_span", "grid_interval", "min_control_sd", "heatmap_display_limit")) {
    x <- cfg[[nm]]
    if (length(x) != 1L || !is.finite(x) || x <= 0) stop(nm, " must be positive")
  }
  if (!cfg$loess_degree %in% 1:2) stop("loess_degree must be 1 or 2")
  if (cfg$selection_alpha <= 0 || cfg$selection_alpha >= 1) stop("Invalid selection_alpha")
  if (cfg$balance_threshold <= 0) stop("balance_threshold must be positive")
  if (length(cfg$winsor_probs) != 2L || any(!is.finite(cfg$winsor_probs)) ||
      cfg$winsor_probs[1] < 0 || cfg$winsor_probs[2] > 1 ||
      cfg$winsor_probs[1] >= cfg$winsor_probs[2]) stop("Invalid winsor_probs")
  if (cfg$denominator_diagnostic_quantile <= 0 || cfg$denominator_diagnostic_quantile >= 1)
    stop("Invalid denominator diagnostic quantile")
  if (cfg$exclude_early_years < 0 || any(!is.finite(cfg$sensitivity_spans)) ||
      any(cfg$sensitivity_spans <= 0)) stop("Invalid sensitivity settings")
  invisible(cfg)
}

write_table <- function(x, name, dirs) {
  write.csv(x, file.path(dirs$tables, paste0(name, ".csv")), row.names = FALSE, na = "")
  saveRDS(x, file.path(dirs$data, paste0(name, ".rds")))
  invisible(x)
}

make_directories <- function(out_dir) {
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    stop("Output directory is not empty. Choose a new directory: ", out_dir)
  dirs <- list(root = out_dir, data = file.path(out_dir, "data_processed"),
               tables = file.path(out_dir, "tables"), figures = file.path(out_dir, "figures"))
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

safe_mean <- function(x) if (any(is.finite(x))) mean(x[is.finite(x)]) else NA_real_
safe_cor <- function(x, y) {
  keep <- is.finite(x) & is.finite(y)
  if (sum(keep) < 3L || sd(x[keep]) == 0 || sd(y[keep]) == 0) return(NA_real_)
  cor(x[keep], y[keep])
}

# ------------------------- Cohort construction ------------------------------
# The control definition is deliberately the original noncase-through-follow-up
# design. This function DOES NOT construct case-date-specific risk sets.
prepare_cohort <- function(data, metabolite_cols, cfg) {
  d <- as.data.frame(data, stringsAsFactors = FALSE)
  cols <- cfg$columns
  if (!identical(sort(names(cols)), sort(c("id", "baseline", "diagnosis", "followup_end"))))
    stop("columns must map id, baseline, diagnosis, and followup_end")
  if (anyDuplicated(unname(cols))) stop("Each date/ID mapping must name a different input column")
  covars <- unique(c(cfg$matching_covariates, cfg$exact_covariates,
                    cfg$adjustment_covariates, cfg$cox_covariates))
  required <- unique(c(unname(cols), covars, metabolite_cols))
  absent <- setdiff(required, names(d))
  if (length(absent)) stop("Missing input columns: ", paste(absent, collapse = ", "))
  if (!length(metabolite_cols) || anyDuplicated(metabolite_cols)) stop("Supply unique metabolite_cols")
  if (length(intersect(metabolite_cols, c(covars, unname(cols)))))
    stop("Metabolite columns must not overlap covariates or date/ID fields")
  # Formula variables should have standard R names. Rename nonstandard fields once
  # at import; do not edit each downstream model separately.
  if (any(make.names(c(covars, metabolite_cols)) != c(covars, metabolite_cols)))
    stop("Map model variables to syntactically valid, unique R column names")
  if (anyDuplicated(names(d))) stop("Duplicate column names")
  d$participant_id <- as.character(d[[cols[["id"]]]])
  if (anyNA(d$participant_id) || any(!nzchar(d$participant_id)) || anyDuplicated(d$participant_id))
    stop("Participant identifiers must be nonmissing and unique")
  parse_date <- function(x, label) {
    if (inherits(x, "Date")) return(x)
    s <- as.character(x)
    s[trimws(s) == ""] <- NA_character_
    v <- as.Date(s, format = "%Y-%m-%d")
    if (any(!is.na(s) & is.na(v))) stop(label, ": use valid ISO YYYY-MM-DD dates")
    v
  }
  d$baseline_date <- parse_date(d[[cols[["baseline"]]]], "baseline")
  d$diagnosis_date <- parse_date(d[[cols[["diagnosis"]]]], "diagnosis")
  d$followup_end <- parse_date(d[[cols[["followup_end"]]]], "followup_end")
  reason <- rep(NA_character_, nrow(d))
  bad_dates <- is.na(d$baseline_date) | is.na(d$followup_end) |
    (!is.na(d$baseline_date) & !is.na(d$followup_end) & d$baseline_date >= d$followup_end)
  reason[bad_dates] <- "missing_or_nonpositive_observation_interval"
  prevalent <- !is.na(d$diagnosis_date) & !is.na(d$baseline_date) &
    d$diagnosis_date <= d$baseline_date
  reason[is.na(reason) & prevalent] <- "diagnosis_on_or_before_baseline"
  for (v in covars) {
    if (v %in% cfg$factor_covariates) {
      x <- as.character(d[[v]])
      x[!is.na(x) & trimws(x) == ""] <- NA_character_
      d[[v]] <- factor(x)
    } else {
      x <- d[[v]]
      if (!is.numeric(x)) {
        original <- as.character(x)
        x <- suppressWarnings(as.numeric(original))
        if (any(!is.na(original) & nzchar(trimws(original)) & is.na(x)))
          stop(v, ": invalid numeric covariate")
      }
      x[!is.finite(x)] <- NA_real_
      d[[v]] <- x
    }
  }
  complete <- if (length(covars)) complete.cases(d[, covars, drop = FALSE]) else rep(TRUE, nrow(d))
  reason[is.na(reason) & !complete] <- "incomplete_design_or_adjustment_covariates"
  excluded <- data.frame(participant_id = d$participant_id[!is.na(reason)],
                         reason = reason[!is.na(reason)])
  n_input <- nrow(d)
  d <- droplevels(d[is.na(reason), , drop = FALSE])
  if (!nrow(d)) stop("No eligible participants")
  d$status <- as.integer(!is.na(d$diagnosis_date) & d$diagnosis_date <= d$followup_end)
  d$years_to_diagnosis <- NA_real_
  is_case <- d$status == 1L
  d$years_to_diagnosis[is_case] <- as.numeric(d$diagnosis_date[is_case] - d$baseline_date[is_case]) / 365.25
  d$time_at_risk_years <- as.numeric(d$followup_end - d$baseline_date) / 365.25
  d$time_at_risk_years[is_case] <- d$years_to_diagnosis[is_case]
  if (!all(c(0L, 1L) %in% d$status)) stop("Both cases and controls are required")
  for (v in covars) if (length(unique(d[[v]])) < 2L)
    stop(v, ": only one observed value; revise the prespecified covariate set")
  for (m in metabolite_cols) {
    if (!is.numeric(d[[m]]) || any(!is.finite(d[[m]])) || sd(d[[m]]) == 0)
      stop(m, ": expected finite, nonconstant preprocessed measurements")
  }
  # Stable matching order and factor reference levels are retained in the RDS.
  d <- d[order(d$participant_id), , drop = FALSE]
  rownames(d) <- NULL
  flow <- data.frame(stage = c("input", "excluded", "eligible", "cases", "noncase_controls"),
                     n = c(n_input, nrow(excluded), nrow(d), sum(is_case), sum(!is_case)))
  list(data = d, excluded = excluded, flow = flow)
}

# ------------------------- Cox selection ------------------------------------
# Fit on the whole eligible synthetic cohort, not on the matched subsample.
# This intentionally retains the same-cohort post-selection design of the example.
fit_cox_results <- function(data, metabolite_cols, dictionary, covariates) {
  rows <- lapply(metabolite_cols, function(m) {
    f <- as.formula(paste("survival::Surv(time_at_risk_years, status) ~",
                         paste(c(m, covariates), collapse = " + ")))
    fit <- tryCatch(survival::coxph(f, data = data, ties = "efron", na.action = na.fail),
                    error = function(e) e)
    if (inherits(fit, "error")) return(data.frame(
      Metabolite_raw = m, HR = NA_real_, lower_95 = NA_real_, upper_95 = NA_real_,
      p_value = NA_real_, n = nrow(data), events = sum(data$status),
      fit_status = conditionMessage(fit)))
    sm <- summary(fit)
    b <- unname(coef(fit)[m]); se <- sqrt(unname(vcov(fit)[m, m]))
    data.frame(Metabolite_raw = m, HR = exp(b), lower_95 = exp(b - qnorm(.975) * se),
               upper_95 = exp(b + qnorm(.975) * se),
               p_value = unname(sm$coefficients[m, "Pr(>|z|)"]),
               n = fit$n, events = fit$nevent,
               fit_status = if (is.finite(b) && is.finite(se)) "fitted" else "nonfinite_coefficient")
  })
  ans <- do.call(rbind, rows)
  ans$p_bonf <- p.adjust(ans$p_value, method = "bonferroni", n = length(metabolite_cols))
  ans$p_fdr <- p.adjust(ans$p_value, method = "BH", n = length(metabolite_cols))
  ix <- match(ans$Metabolite_raw, dictionary$Metabolite_raw)
  ans$Metabolite_name <- dictionary$Metabolite_name[ix]
  ans$Metabolite <- if ("Metabolite" %in% names(dictionary)) dictionary$Metabolite[ix] else ans$Metabolite_raw
  ans
}

# ------------------------- Matching and balance -----------------------------
match_and_assess <- function(data, cfg) {
  exact <- if (length(cfg$exact_covariates)) reformulate(cfg$exact_covariates) else NULL
  m <- MatchIt::matchit(
    reformulate(cfg$matching_covariates, response = "status"), data = data,
    method = "nearest", distance = "mahalanobis", ratio = cfg$match_ratio,
    replace = FALSE, exact = exact, m.order = cfg$matching_order
  )
  # Provide data explicitly, so no caller-local object is required at extraction.
  md <- as.data.frame(MatchIt::match.data(m, data = data))
  if (!nrow(md) || !"subclass" %in% names(md)) stop("No matched subclasses returned")
  md$subclass <- as.character(md$subclass)
  md <- md[order(md$subclass, -md$status, md$participant_id), , drop = FALSE]
  groups <- split(seq_len(nrow(md)), md$subclass)
  qc <- do.call(rbind, lapply(names(groups), function(s) {
    x <- md[groups[[s]], , drop = FALSE]
    exact_ok <- all(vapply(cfg$exact_covariates, function(v) length(unique(x[[v]])) == 1L, logical(1)))
    data.frame(subclass = s, n_total = nrow(x), n_case = sum(x$status == 1L),
               n_control = sum(x$status == 0L), exact_match_ok = exact_ok)
  }))
  bad <- qc$n_case != 1L | qc$n_control != cfg$match_ratio | !qc$exact_match_ok
  # Do not silently drop partial sets or pad controls to obtain the requested ratio.
  if (any(bad)) stop("Incomplete/invalid matched sets: ", paste(head(qc$subclass[bad], 10), collapse = ", "),
                     ". No automatic change to ratio or eligibility is made.")
  if (anyDuplicated(md$participant_id)) stop("Unexpected duplicate participants in no-replacement matching")
  sm <- summary(m, standardize = TRUE, pair.dist = FALSE, un = TRUE)
  matched <- sm$sum.matched
  smd_col <- grep("Std.*Mean.*Diff", colnames(matched), value = TRUE)
  if (length(smd_col) != 1L) stop("Cannot identify MatchIt standardized mean difference column")
  bal <- data.frame(variable = rownames(matched), matched, check.names = FALSE, row.names = NULL)
  bal$abs_smd <- abs(matched[, smd_col])
  bal$abs_smd_before <- abs(sm$sum.all[rownames(matched), smd_col])
  bal$above_threshold <- is.finite(bal$abs_smd) & bal$abs_smd >= cfg$balance_threshold
  if (any(bal$above_threshold)) warning("Some post-match SMDs exceed the diagnostic threshold; inspect covariate_balance.csv.")
  if (any(!is.finite(bal$abs_smd))) warning("Some balance diagnostics are nonfinite; inspect covariate_balance.csv.")
  unused_cases <- setdiff(data$participant_id[data$status == 1L], md$participant_id[md$status == 1L])
  list(matched = md, matchit = m, integrity = qc, balance = bal,
       unmatched_cases = data.frame(participant_id = unused_cases))
}

# ------------------------- Residuals and matched Z --------------------------
residualize_metabolites <- function(matched, metabolite_cols, covariates) {
  out <- matched[, c("participant_id", "subclass", "status", "years_to_diagnosis"), drop = FALSE]
  for (m in metabolite_cols) {
    fit <- lm(reformulate(covariates, response = m), data = matched, na.action = na.exclude)
    if (fit$rank < ncol(model.matrix(fit))) warning("Rank-deficient residualization model for ", m)
    out[[paste0(m, "_res")]] <- as.numeric(residuals(fit))
  }
  out
}

calculate_z_scores <- function(residuals, metabolite_cols, min_controls = 2L, min_sd = 1e-12) {
  groups <- split(seq_len(nrow(residuals)), residuals$subclass)
  ans <- lapply(metabolite_cols, function(m) {
    do.call(rbind, lapply(names(groups), function(s) {
      d <- residuals[groups[[s]], , drop = FALSE]
      x <- d[[paste0(m, "_res")]]
      ci <- which(d$status == 1L)
      controls <- x[d$status == 0L & is.finite(x)]
      csd <- sd(controls)
      case <- if (length(ci) == 1L) x[ci] else NA_real_
      reason <- if (length(ci) != 1L) "not_one_case" else if (!is.finite(case)) "missing_case" else
        if (length(controls) < min_controls) "too_few_controls" else
          if (!is.finite(csd) || csd <= min_sd) "zero_or_near_zero_sd" else "valid"
      data.frame(subclass = s, Metabolite_raw = m,
                 Year = if (length(ci) == 1L) -d$years_to_diagnosis[ci] else NA_real_,
                 case_residual = case, control_mean = if (length(controls)) mean(controls) else NA_real_,
                 control_sd = csd, n_valid_controls = length(controls),
                 Estimate = if (reason == "valid") (case - mean(controls)) / csd else NA_real_,
                 z_status = reason)
    }))
  })
  do.call(rbind, ans)
}

# ------------------------- LOESS --------------------------------------------
make_prediction_grid <- function(z, interval) {
  y <- z$Year[is.finite(z$Year) & is.finite(z$Estimate)]
  if (!length(y)) stop("No finite matched-set Z-scores")
  lower <- ceiling(min(y) / interval) * interval
  upper <- floor(max(y) / interval) * interval
  if (lower >= upper) stop("Insufficient time support for the requested grid")
  seq(lower, upper, by = interval)
}

fit_gradients <- function(z, metabolites, grid, span, degree = 2L,
                          min_sets = 50L, min_times = 10L, estimate_col = "Estimate") {
  audit <- vector("list", length(metabolites))
  predictions <- vector("list", length(metabolites))
  for (j in seq_along(metabolites)) {
    m <- metabolites[j]
    keep <- z$Metabolite_raw == m & is.finite(z$Year) & is.finite(z[[estimate_col]])
    d <- data.frame(Year = z$Year[keep], Estimate = z[[estimate_col]][keep])
    p <- rep(NA_real_, length(grid)); fit_status <- "insufficient_data"
    if (nrow(d) >= min_sets && length(unique(d$Year)) >= min_times) {
      result <- tryCatch({
        fit <- loess(Estimate ~ Year, data = d, span = span, degree = degree,
                     family = "gaussian", control = loess.control(surface = "interpolate"))
        as.numeric(predict(fit, newdata = data.frame(Year = grid)))
      }, error = function(e) e)
      if (inherits(result, "error")) {
        fit_status <- paste0("error: ", conditionMessage(result))
      } else {
        p <- result
        # Do not extrapolate at an individual metabolite's unsupported edges.
        p[grid < min(d$Year) | grid > max(d$Year)] <- NA_real_
        fit_status <- if (all(is.finite(p))) "fitted" else "fitted_with_unsupported_predictions"
      }
    }
    predictions[[j]] <- data.frame(Metabolite_raw = m, Year = grid, Estimate_loess = p)
    audit[[j]] <- data.frame(Metabolite_raw = m, n_valid_sets = nrow(d),
                             n_unique_times = length(unique(d$Year)),
                             n_finite_predictions = sum(is.finite(p)), fit_status = fit_status)
  }
  list(predictions = do.call(rbind, predictions), audit = do.call(rbind, audit))
}

# ------------------------- Descriptive clustering --------------------------
cluster_predictions <- function(pred, k = 3L, method = "ward.D") {
  ids <- sort(unique(pred$Metabolite_raw)); grid <- sort(unique(pred$Year))
  mat <- matrix(NA_real_, length(ids), length(grid), dimnames = list(ids, as.character(grid)))
  keys <- paste(pred$Metabolite_raw, pred$Year, sep = "|")
  if (anyDuplicated(keys)) stop("Duplicate metabolite/grid entries")
  mat[cbind(match(pred$Metabolite_raw, ids), match(pred$Year, grid))] <- pred$Estimate_loess
  complete <- apply(mat, 1L, function(x) all(is.finite(x)))
  excluded <- data.frame(Metabolite_raw = ids[!complete], reason = rep("missing_prediction_on_grid", sum(!complete)))
  mat <- mat[complete, , drop = FALSE]
  if (nrow(mat) < k) stop("Too few complete selected curves for k=", k,
                         ". Check selection, LOESS support, and exclusions; no fallback selection is used.")
  distance <- dist(mat, method = "euclidean")
  tree <- hclust(distance, method = method)
  raw <- cutree(tree, k = k)
  means <- tapply(rowMeans(mat), raw, mean)
  ord <- as.integer(names(means)[order(-means, as.integer(names(means)))])
  # Preserve raw labels; deterministic display labels only. For k=3: highest,
  # lowest, intermediate mean profile (without asserting a biological subtype).
  if (k == 3L) ord <- ord[c(1L, 3L, 2L)]
  display_map <- setNames(seq_len(k), as.character(ord))
  membership <- data.frame(Metabolite_raw = names(raw), cluster_raw = unname(raw),
                           cluster = unname(display_map[as.character(raw)]))
  joined <- merge(pred, membership, by = "Metabolite_raw", sort = FALSE)
  joined <- joined[order(joined$cluster, joined$Metabolite_raw, joined$Year), , drop = FALSE]
  summaries <- do.call(rbind, lapply(sort(unique(joined$cluster)), function(g) {
    d <- joined[joined$cluster == g, , drop = FALSE]
    data.frame(cluster = g, n_metabolites = length(unique(d$Metabolite_raw)),
               mean_profile = mean(d$Estimate_loess),
               far_12plus_years = safe_mean(d$Estimate_loess[d$Year <= -12]),
               middle_5_to_10_years = safe_mean(d$Estimate_loess[d$Year >= -10 & d$Year <= -5]),
               near_2_years = safe_mean(d$Estimate_loess[d$Year >= -2]),
               min_profile = min(d$Estimate_loess), max_profile = max(d$Estimate_loess))
  }))
  list(matrix = mat, distance = distance, hclust = tree, membership = membership,
       excluded = excluded, predictions = joined, summary = summaries)
}

# ------------------------- Plots --------------------------------------------
# Plot styles are illustrative; figure-source values remain untrimmed.
save_plot <- function(plot, name, dirs, cfg, width = 8, height = 5) {
  if (!cfg$make_plots) return(invisible(NULL))
  ggplot2::ggsave(file.path(dirs$figures, paste0(name, ".png")), plot = plot,
                 width = width, height = height, units = "in", dpi = cfg$plot_dpi)
  invisible(plot)
}

plot_quality <- function(matched, balance, z, completeness, dirs, cfg) {
  b <- balance[is.finite(balance$abs_smd), , drop = FALSE]
  b$variable <- factor(b$variable, levels = b$variable[order(b$abs_smd)])
  p <- ggplot2::ggplot(b, ggplot2::aes(x = abs_smd, y = variable)) +
    ggplot2::geom_point() + ggplot2::geom_vline(xintercept = cfg$balance_threshold, linetype = 2) +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "post-match balance"), x = "Absolute SMD", y = NULL) +
    ggplot2::theme_bw()
  save_plot(p, "figure1A_covariate_balance", dirs, cfg, height = 7)
  cases <- matched[matched$status == 1L, , drop = FALSE]
  write_table(cases[, c("participant_id", "subclass", "years_to_diagnosis"), drop = FALSE],
              "figure1B_interval_source", dirs)
  p <- ggplot2::ggplot(cases, ggplot2::aes(x = years_to_diagnosis)) +
    ggplot2::geom_histogram(bins = 25) +
    ggplot2::geom_vline(xintercept = median(cases$years_to_diagnosis), linetype = 2) +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "diagnosis intervals"), x = "Years from baseline to diagnosis", y = "Cases") +
    ggplot2::theme_bw()
  save_plot(p, "figure1B_diagnosis_intervals", dirs, cfg)
  zz <- z[is.finite(z$Estimate), , drop = FALSE]
  if (nrow(zz) > cfg$max_z_plot) zz <- zz[sample.int(nrow(zz), cfg$max_z_plot), , drop = FALSE]
  write_table(zz, "figure1C_z_score_plot_source", dirs)
  p <- ggplot2::ggplot(zz, ggplot2::aes(x = Estimate)) + ggplot2::geom_histogram(bins = 50) +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "matched-set Z-scores"), x = "Matched-set Z", y = "Observations") +
    ggplot2::theme_bw()
  save_plot(p, "figure1C_z_distribution", dirs, cfg)
  p <- ggplot2::ggplot(completeness, ggplot2::aes(x = Metabolite_raw, y = valid_fraction)) +
    ggplot2::geom_col() + ggplot2::coord_flip() +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "Z-score completeness"), x = NULL, y = "Fraction of matched sets with a valid Z") +
    ggplot2::theme_bw()
  save_plot(p, "figure1D_z_completeness", dirs, cfg, height = 6)
  p <- ggplot2::ggplot(z, ggplot2::aes(x = Metabolite_raw, y = control_sd)) +
    ggplot2::geom_boxplot(outlier.size = 0.5, na.rm = TRUE) + ggplot2::coord_flip() +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "denominator diagnostic"), x = NULL, y = "Within-set control residual SD") +
    ggplot2::theme_bw()
  save_plot(p, "diagnostic_control_sd", dirs, cfg, height = 6)
}

plot_cluster_results <- function(cl, dictionary, dirs, cfg, prefix = "figure2") {
  d <- cl$predictions
  d$Metabolite_name <- dictionary$Metabolite_name[match(d$Metabolite_raw, dictionary$Metabolite_raw)]
  d$cluster_label <- paste("Cluster", d$cluster)
  avg <- aggregate(Estimate_loess ~ cluster_label + Year, data = d, FUN = mean)
  write_table(d, paste0(prefix, "_curve_and_heatmap_source"), dirs)
  write_table(avg, paste0(prefix, "_cluster_mean_source"), dirs)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Year, y = Estimate_loess, group = Metabolite_raw)) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_line(alpha = 0.30) +
    ggplot2::geom_line(data = avg, ggplot2::aes(group = cluster_label), linewidth = 1) +
    ggplot2::facet_wrap(~cluster_label, ncol = 1) +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "descriptive gradients"), x = "Years relative to diagnosis (0 = diagnosis)", y = "Smoothed matched-set Z") +
    ggplot2::theme_bw()
  save_plot(p, paste0(prefix, "_cluster_curves"), dirs, cfg, height = 8)
  order_ids <- cl$membership$Metabolite_raw[order(cl$membership$cluster, cl$membership$Metabolite_raw)]
  order_names <- dictionary$Metabolite_name[match(order_ids, dictionary$Metabolite_raw)]
  d$Metabolite_name <- factor(d$Metabolite_name, levels = rev(order_names))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = Year, y = Metabolite_name, fill = Estimate_loess)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(midpoint = 0, limits = c(-cfg$heatmap_display_limit, cfg$heatmap_display_limit),
      oob = function(x, range, ...) pmin(pmax(x, range[1]), range[2]), name = "Smoothed Z") +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "selected-metabolite heatmap"), x = "Years relative to diagnosis", y = NULL) +
    ggplot2::theme_bw()
  save_plot(p, paste0(prefix, "_heatmap"), dirs, cfg, width = 10,
            height = max(5, min(14, 0.25 * length(order_ids) + 2)))
}

# ------------------------- Sensitivity analyses -----------------------------
run_sensitivities <- function(z, selected, primary, primary_cluster, grid, cfg, dirs) {
  base_z <- z[z$Metabolite_raw %in% selected, , drop = FALSE]
  scenarios <- list()
  for (s in cfg$sensitivity_spans) scenarios[[paste0("span_", format(s, nsmall = 2))]] <- list(data = base_z, span = s)
  scenarios$exclude_early <- list(data = base_z[base_z$Year <= -cfg$exclude_early_years, , drop = FALSE], span = cfg$loess_span)
  wz <- base_z
  for (m in selected) {
    ix <- which(wz$Metabolite_raw == m & is.finite(wz$Estimate))
    if (length(ix)) {
      q <- quantile(wz$Estimate[ix], cfg$winsor_probs, names = FALSE, type = 7)
      wz$Estimate[ix] <- pmin(pmax(wz$Estimate[ix], q[1]), q[2])
    }
  }
  scenarios$winsorized <- list(data = wz, span = cfg$loess_span)
  curves <- list(); agreement <- list(); outputs <- list()
  pm <- setNames(primary_cluster$membership$cluster_raw, primary_cluster$membership$Metabolite_raw)
  for (nm in names(scenarios)) {
    sc <- scenarios[[nm]]
    supported <- sc$data$Year[is.finite(sc$data$Year) & is.finite(sc$data$Estimate)]
    # Use only primary grid points also supported by the restricted scenario.
    sg <- if (length(supported)) grid[grid >= min(supported) & grid <= max(supported)] else numeric()
    if (length(sg) < 2L) {
      warning("Skipping sensitivity scenario with insufficient time support: ", nm)
      agreement[[nm]] <- data.frame(scenario = nm, ARI = NA_real_, n_common_metabolites = 0L, status = "insufficient_time_support")
      next
    }
    fit <- fit_gradients(sc$data, selected, sg, sc$span, cfg$loess_degree,
                         cfg$min_loess_sets, cfg$min_unique_times)
    write_table(fit$predictions, paste0("sensitivity_predictions_", nm), dirs)
    write_table(fit$audit, paste0("sensitivity_fit_audit_", nm), dirs)
    common <- merge(primary, fit$predictions, by = c("Metabolite_raw", "Year"), suffixes = c("_primary", "_sensitivity"))
    curves[[nm]] <- do.call(rbind, lapply(selected, function(m) {
      dd <- common[common$Metabolite_raw == m, , drop = FALSE]
      a <- dd$Estimate_loess_primary; b <- dd$Estimate_loess_sensitivity
      keep <- is.finite(a) & is.finite(b)
      data.frame(scenario = nm, Metabolite_raw = m, n_common_grid_points = sum(keep),
                 correlation = safe_cor(a, b), RMSE = if (any(keep)) sqrt(mean((a[keep] - b[keep])^2)) else NA_real_)
    }))
    cl <- tryCatch(cluster_predictions(fit$predictions, cfg$n_clusters, cfg$linkage), error = function(e) e)
    if (inherits(cl, "error")) {
      warning("Sensitivity clustering unavailable for ", nm, ": ", conditionMessage(cl))
      agreement[[nm]] <- data.frame(scenario = nm, ARI = NA_real_, n_common_metabolites = 0L, status = conditionMessage(cl))
    } else {
      cm <- setNames(cl$membership$cluster_raw, cl$membership$Metabolite_raw)
      ids <- intersect(names(pm), names(cm))
      agreement[[nm]] <- data.frame(scenario = nm,
        ARI = if (length(ids) >= 2L) mclust::adjustedRandIndex(pm[ids], cm[ids]) else NA_real_,
        n_common_metabolites = length(ids), status = "computed")
      write_table(cl$membership, paste0("sensitivity_membership_", nm), dirs)
    }
    outputs[[nm]] <- list(fit = fit, clustering = cl)
  }
  alt <- cluster_predictions(primary, cfg$n_clusters, cfg$alternative_linkage)
  cm <- setNames(alt$membership$cluster_raw, alt$membership$Metabolite_raw)
  ids <- intersect(names(pm), names(cm))
  agreement$alternative_linkage <- data.frame(scenario = cfg$alternative_linkage,
    ARI = mclust::adjustedRandIndex(pm[ids], cm[ids]), n_common_metabolites = length(ids), status = "computed")
  write_table(alt$membership, "sensitivity_alternative_linkage_membership", dirs)
  curve_table <- if (length(curves)) do.call(rbind, curves) else data.frame(
    scenario = character(), Metabolite_raw = character(), n_common_grid_points = integer(), correlation = numeric(), RMSE = numeric())
  agreement_table <- do.call(rbind, agreement)
  write_table(curve_table, "sensitivity_curve_similarity", dirs)
  write_table(agreement_table, "sensitivity_cluster_agreement", dirs)
  # Ward.D2 intentionally has no row in this curve-correlation table or plot.
  summaries <- lapply(unique(curve_table$scenario), function(nm) {
    x <- curve_table$correlation[curve_table$scenario == nm]
    x <- x[is.finite(x)]
    q <- if (length(x)) quantile(x, c(.25, .5, .75), names = FALSE) else rep(NA_real_, 3)
    data.frame(scenario = nm, q25 = q[1], median = q[2], q75 = q[3], n = length(x))
  })
  if (length(summaries)) {
    cs <- do.call(rbind, summaries)
    write_table(cs, "figure3_curve_similarity_source", dirs)
    p <- ggplot2::ggplot(cs, ggplot2::aes(x = scenario, y = median)) +
      ggplot2::geom_hline(yintercept = cfg$curve_reference, linetype = 2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = q25, ymax = q75), width = .15, na.rm = TRUE) +
      ggplot2::geom_point(na.rm = TRUE) + ggplot2::coord_flip() +
      ggplot2::labs(title = paste(cfg$plot_title_prefix, "curve sensitivity"), x = NULL, y = "Median Pearson correlation (IQR)") +
      ggplot2::theme_bw()
    save_plot(p, "figure3A_curve_sensitivity", dirs, cfg)
  }
  p <- ggplot2::ggplot(agreement_table, ggplot2::aes(x = scenario, y = ARI)) +
    ggplot2::geom_hline(yintercept = cfg$ari_reference, linetype = 2) + ggplot2::geom_point(na.rm = TRUE) +
    ggplot2::coord_flip() + ggplot2::labs(title = paste(cfg$plot_title_prefix, "cluster sensitivity"), x = NULL, y = "Adjusted Rand index") +
    ggplot2::theme_bw()
  save_plot(p, "figure3B_cluster_sensitivity", dirs, cfg)
  saveRDS(list(scenarios = outputs, curve_similarity = curve_table, cluster_agreement = agreement_table),
          file.path(dirs$data, "sensitivity_results.rds"))
  invisible(agreement_table)
}

cluster_count_diagnostics <- function(cl, cfg, dirs) {
  ks <- sort(unique(as.integer(cfg$silhouette_k)))
  ks <- ks[ks >= 2L & ks < nrow(cl$matrix)]
  if (!length(ks)) return(invisible(NULL))
  stats <- do.call(rbind, lapply(ks, function(k) {
    membership <- cutree(cl$hclust, k = k)
    silhouette <- cluster::silhouette(membership, cl$distance)
    data.frame(k = k, mean_silhouette_width = mean(silhouette[, "sil_width"]))
  }))
  write_table(stats, "cluster_count_silhouette", dirs)
  neighbour_k <- unique(c(cfg$n_clusters - 1L, cfg$n_clusters, cfg$n_clusters + 1L))
  neighbour_k <- neighbour_k[neighbour_k >= 2L & neighbour_k <= nrow(cl$matrix)]
  assignments <- data.frame(Metabolite_raw = rownames(cl$matrix))
  for (k in neighbour_k) assignments[[paste0("k", k)]] <- unname(cutree(cl$hclust, k = k)[assignments$Metabolite_raw])
  write_table(assignments, "cluster_count_membership_comparison", dirs)
  for (k in setdiff(neighbour_k, cfg$n_clusters)) {
    tab <- as.data.frame(table(primary = assignments[[paste0("k", cfg$n_clusters)]],
                               alternative = assignments[[paste0("k", k)]]))
    write_table(tab, paste0("cluster_merge_split_k", k), dirs)
  }
  p <- ggplot2::ggplot(stats, ggplot2::aes(x = k, y = mean_silhouette_width)) +
    ggplot2::geom_line() + ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = cfg$n_clusters, linetype = 2) +
    ggplot2::labs(title = paste(cfg$plot_title_prefix, "cluster-count diagnostic"), x = "Number of clusters", y = "Mean silhouette width") +
    ggplot2::theme_bw()
  save_plot(p, "diagnostic_cluster_count", dirs, cfg)
  invisible(stats)
}

# ------------------------- Cohort-agnostic analysis entry -------------------
# Pass already preprocessed data and an explicit metabolite dictionary.
# No UKB field numbers, private Excel files, or global 'dat' object are required.
run_gradient_pipeline <- function(cohort, metabolite_cols, dictionary,
                                  out_dir, config = protocol_config(), synthetic = TRUE) {
  check_packages(); check_config(config)
  if (length(synthetic) != 1L || is.na(synthetic)) stop("synthetic must be TRUE or FALSE")
  needed <- c("Metabolite_raw", "Metabolite_name")
  if (!all(needed %in% names(dictionary)) || anyDuplicated(dictionary$Metabolite_raw))
    stop("dictionary must have unique Metabolite_raw and Metabolite_name columns")
  ix <- match(metabolite_cols, dictionary$Metabolite_raw)
  if (anyNA(ix) || anyNA(dictionary$Metabolite_name[ix]) || any(!nzchar(dictionary$Metabolite_name[ix])) ||
      anyDuplicated(dictionary$Metabolite_name[ix])) stop("Every metabolite needs a unique, nonmissing display name")
  dirs <- make_directories(out_dir)
  cfg <- config
  cfg$plot_title_prefix <- if (isTRUE(synthetic)) "SYNTHETIC DEMONSTRATION:" else "Matched-gradient analysis:"
  cfg$data_type <- if (isTRUE(synthetic)) "FULLY_SYNTHETIC" else "USER_SUPPLIED"
  set.seed(cfg$seed)
  saveRDS(cfg, file.path(dirs$data, "analysis_configuration.rds"))
  dput(cfg, file = file.path(dirs$root, "analysis_configuration.R"))
  writeLines(capture.output(sessionInfo()), file.path(dirs$root, "sessionInfo.txt"))
  warnings_seen <- character(); completed <- FALSE
  on.exit({
    writeLines(if (length(warnings_seen)) unique(warnings_seen) else "No R warnings captured.", file.path(dirs$root, "warnings.log"))
    status_line <- if (!completed) {
      "NOT COMPLETED: inspect console error and partial outputs; do not claim successful execution."
    } else if (isTRUE(synthetic)) {
      "COMPLETED: synthetic execution only; inspect diagnostic warnings."
    } else {
      "COMPLETED: inspect diagnostic warnings and independently validate interpretation."
    }
    writeLines(status_line, file.path(dirs$root, "execution_status.txt"))
  }, add = TRUE)
  result <- withCallingHandlers({
    message("1/9 Constructing the eligible cohort")
    prep <- prepare_cohort(cohort, metabolite_cols, cfg); d <- prep$data
    write_table(prep$flow, "cohort_flow", dirs); write_table(prep$excluded, "cohort_exclusions", dirs)
    saveRDS(d, file.path(dirs$data, "analysis_cohort.rds"))
    message("2/9 Fitting Cox models (no external results file)")
    cox <- fit_cox_results(d, metabolite_cols, dictionary, cfg$cox_covariates)
    write_table(cox, if (isTRUE(synthetic)) "synthetic_cox_results" else "cox_results", dirs)
    message("3/9 Matching and assessing balance")
    match <- match_and_assess(d, cfg)
    saveRDS(match$matchit, file.path(dirs$data, "matchit_object.rds"))
    saveRDS(match$matched, file.path(dirs$data, "matched_dataset.rds"))
    write_table(match$integrity, "matched_set_integrity", dirs)
    write_table(match$balance, "covariate_balance", dirs)
    write_table(match$unmatched_cases, "unmatched_cases", dirs)
    saveRDS(list(set_qc = match$integrity, balance = match$balance), file.path(dirs$data, "matched_quality_control.rds"))
    message("4/9 Residualizing metabolites and calculating matched-set Z-scores")
    residuals <- residualize_metabolites(match$matched, metabolite_cols, cfg$adjustment_covariates)
    saveRDS(residuals, file.path(dirs$data, "metabolite_residuals.rds"))
    z <- calculate_z_scores(residuals, metabolite_cols, cfg$min_valid_controls, cfg$min_control_sd)
    write_table(z, "case_metabolite_z_scores_long", dirs)
    wide <- reshape(z[, c("subclass", "Year", "Metabolite_raw", "Estimate")],
                    idvar = c("subclass", "Year"), timevar = "Metabolite_raw", direction = "wide")
    names(wide) <- sub("^Estimate\\.", "", names(wide))
    write_table(wide, "case_metabolite_z_scores_wide", dirs)
    completeness <- do.call(rbind, lapply(metabolite_cols, function(m) {
      q <- z[z$Metabolite_raw == m, , drop = FALSE]
      data.frame(Metabolite_raw = m, n_sets = nrow(q), n_valid = sum(is.finite(q$Estimate)),
                 valid_fraction = mean(is.finite(q$Estimate)))
    }))
    write_table(completeness, "z_score_completeness", dirs)
    denominator <- do.call(rbind, lapply(metabolite_cols, function(m) {
      q <- z[z$Metabolite_raw == m, , drop = FALSE]
      s <- q$control_sd[is.finite(q$control_sd)]
      qs <- if (length(s)) quantile(s, c(cfg$denominator_diagnostic_quantile, .5), names = FALSE) else c(NA_real_, NA_real_)
      data.frame(Metabolite_raw = m, minimum_control_sd = if (length(s)) min(s) else NA_real_,
                 low_sd_quantile = qs[1], median_control_sd = qs[2],
                 invalid_z_fraction = mean(!is.finite(q$Estimate)),
                 n_sd_below_numeric_threshold = sum(is.finite(q$control_sd) & q$control_sd <= cfg$min_control_sd))
    }))
    write_table(denominator, "control_sd_diagnostics", dirs)
    message("5/9 Fitting all metabolite curves on an inward-rounded grid")
    grid <- make_prediction_grid(z, cfg$grid_interval)
    write_table(data.frame(Year = grid), "prediction_grid", dirs)
    fit <- fit_gradients(z, metabolite_cols, grid, cfg$loess_span, cfg$loess_degree,
                         cfg$min_loess_sets, cfg$min_unique_times)
    write_table(fit$predictions, "loess_predictions_all_metabolites", dirs)
    write_table(fit$audit, "loess_fit_audit", dirs)
    selected <- cox$Metabolite_raw[is.finite(cox$p_bonf) & cox$p_bonf < cfg$selection_alpha]
    write_table(cox[cox$Metabolite_raw %in% selected, , drop = FALSE], "selected_metabolites", dirs)
    if (length(selected) < cfg$n_clusters) stop("Only ", length(selected),
      " metabolites pass the prespecified Bonferroni threshold. No forced top-N fallback is used.")
    pred <- fit$predictions[fit$predictions$Metabolite_raw %in% selected, , drop = FALSE]
    write_table(pred, "loess_predictions_significant", dirs)
    message("6/9 Clustering selected curves")
    cl <- cluster_predictions(pred, cfg$n_clusters, cfg$linkage)
    write_table(cl$membership, "cluster_membership", dirs)
    write_table(cl$summary, "cluster_summary", dirs)
    write_table(cl$excluded, "clustering_excluded_metabolites", dirs)
    saveRDS(cl, file.path(dirs$data, "cluster_outputs.rds"))
    message("7/9 Writing diagnostic figures and figure-source tables")
    plot_quality(match$matched, match$balance, z, completeness, dirs, cfg)
    plot_cluster_results(cl, dictionary, dirs, cfg)
    message("8/9 Running prespecified sensitivity analyses")
    sensitivity <- run_sensitivities(z, selected, pred, cl, grid, cfg, dirs)
    message("9/9 Saving validation checks and optional cluster-count diagnostics")
    if (cfg$run_cluster_count_diagnostics) cluster_count_diagnostics(cl, cfg, dirs)
    finite <- is.finite(z$Estimate)
    formula_error <- max(abs(z$Estimate[finite] -
      (z$case_residual[finite] - z$control_mean[finite]) / z$control_sd[finite]))
    checks <- data.frame(
      check = c("unique_matched_participants", "one_case_per_set", "requested_controls_per_set",
                "exact_matching", "matched_z_formula", "primary_grid_no_extrapolation",
                "clustered_curves_complete", "cohort_and_dictionary_identifiers"),
      passed = c(!anyDuplicated(match$matched$participant_id), all(match$integrity$n_case == 1L),
                 all(match$integrity$n_control == cfg$match_ratio), all(match$integrity$exact_match_ok),
                 is.finite(formula_error) && formula_error < 1e-10,
                 min(grid) >= min(z$Year[finite]) && max(grid) <= max(z$Year[finite]),
                 all(is.finite(cl$matrix)), all(selected %in% dictionary$Metabolite_raw)))
    write_table(checks, "execution_checks", dirs)
    if (!all(checks$passed)) stop("An execution check failed; inspect execution_checks.csv")
    metrics <- data.frame(metric = c("eligible_participants", "incident_cases", "matched_cases", "matched_controls",
                                    "metabolites_tested", "metabolites_selected", "metabolites_clustered", "max_postmatch_absolute_smd"),
                         value = c(nrow(d), sum(d$status), sum(match$matched$status == 1L), sum(match$matched$status == 0L),
                                   length(metabolite_cols), length(selected), nrow(cl$matrix),
                                   max(match$balance$abs_smd, na.rm = TRUE)))
    write_table(metrics, "run_summary", dirs)
    list(summary = metrics, checks = checks, output_dir = normalizePath(out_dir),
         sensitivity = sensitivity)
  }, warning = function(w) {
    warnings_seen <<- c(warnings_seen, conditionMessage(w))
    message("DIAGNOSTIC WARNING: ", conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  completed <- TRUE
  message("Completed ", if (isTRUE(synthetic)) "synthetic demonstration" else "analysis", ". Results: ", normalizePath(out_dir))
  invisible(result)
}

# ------------------------- One-command synthetic driver ---------------------
run_synthetic_example <- function(out_dir = file.path(.protocol_demo_dir, "synthetic_demo_output"),
                                   config = protocol_config()) {
  check_packages(); check_config(config)
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    stop("Output directory is not empty. Use a new directory to avoid mixing runs.")
  generator <- file.path(.protocol_demo_dir, "generate_synthetic_data.R")
  if (!file.exists(generator)) stop("Place generate_synthetic_data.R beside this script")
  gen_env <- new.env(parent = globalenv())
  sys.source(generator, envir = gen_env)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "FULLY SYNTHETIC DATA AND RESULTS. No real participant records or fitted UKB parameters.",
    "This is a software execution demonstration, not external validation or a new biological analysis."
  ), file.path(out_dir, "SYNTHETIC_ONLY.txt"))
  generated <- gen_env$generate_synthetic_data(file.path(out_dir, "input"),
    n = config$n, n_metabolites = config$n_metabolites, seed = config$seed)
  # Data are already median-imputed, log-transformed, and standardized here.
  # Do NOT log-transform this processed cohort a second time.
  ans <- run_gradient_pipeline(generated$cohort,
    metabolite_cols = generated$dictionary$Metabolite_raw,
    dictionary = generated$dictionary, out_dir = file.path(out_dir, "analysis"), config = config)
  invisible(ans)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  cfg <- protocol_config()
  out <- if (length(args) >= 1L) args[1L] else file.path(.protocol_demo_dir, "synthetic_demo_output")
  if (length(args) >= 2L) cfg$n <- as.integer(args[2L])
  if (length(args) >= 3L) cfg$n_metabolites <- as.integer(args[3L])
  if (length(args) >= 4L) cfg$seed <- as.integer(args[4L])
  tryCatch(run_synthetic_example(out, cfg), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1L)
  })
}
