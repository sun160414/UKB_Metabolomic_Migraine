#!/usr/bin/env Rscript

## 02_trajectories_clusters.R
## Matched case-control time-to-diagnosis gradient analysis for migraine metabolomics
##
## Expected inputs:
##   1) data/mydata_base.csv
##      Individual-level analysis data with eid, migraine_status or incident_migraine,
##      migraine_years, covariates, and Meta_* metabolite columns.
##   2) data/metabolite_name_map.csv
##      Columns: Meta, Original_Metabolite.
##   3) results/cox_overall.csv
##      Output from 01_cox_linear.R. Should contain Metabolite_raw or raw Meta_* IDs
##      and p_bonf / p_fdr columns.
##
## Note: UK Biobank individual-level data are not included in this repository.

required_pkgs <- c(
  "data.table", "dplyr", "tidyr", "tibble", "purrr", "stringr",
  "MatchIt", "ggplot2", "patchwork", "ComplexHeatmap", "circlize", "grid"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Please install required packages: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(MatchIt)
  library(ggplot2)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
})

dir.create("results", showWarnings = FALSE, recursive = TRUE)
dir.create("figures", showWarnings = FALSE, recursive = TRUE)

##==========================
## User settings
##==========================
data_file <- "data/mydata_base.csv"
map_file  <- "data/metabolite_name_map.csv"
cox_file  <- "results/cox_overall.csv"

matching_ratio <- 10
loess_span <- 0.75
loess_step <- 0.2
n_clusters <- 3
significance_column <- "p_bonf"  # use "p_fdr" if FDR-significant metabolites are preferred

covars <- c(
  "age", "sex", "ethn", "Qualification", "bmi", "Socioeconomic",
  "Smoking_status", "Alcohol_consumption",
  "screen time (TV)", "screen time (computer)", "sleep duration",
  "diabetes_status", "CVD_status"
)

##==========================
## Helper functions
##==========================
read_or_get <- function(file, object_name) {
  if (file.exists(file)) {
    return(data.table::fread(file, data.table = FALSE, check.names = FALSE))
  }
  if (exists(object_name, envir = .GlobalEnv)) {
    return(get(object_name, envir = .GlobalEnv))
  }
  stop("Cannot find input file '", file, "' or object '", object_name, "'.")
}

quote_if_needed <- function(x) {
  ifelse(make.names(x) != x, paste0("`", x, "`"), x)
}

standardize_covariates <- function(df) {
  numeric_vars <- intersect(
    c("age", "bmi", "Socioeconomic", "Qualification",
      "screen time (TV)", "screen time (computer)", "sleep duration", "migraine_years"),
    names(df)
  )
  factor_vars <- intersect(
    c("sex", "ethn", "Smoking_status", "Alcohol_consumption",
      "diabetes_status", "CVD_status"),
    names(df)
  )

  df %>%
    mutate(across(all_of(numeric_vars), as.numeric)) %>%
    mutate(across(all_of(factor_vars), as.factor))
}

make_meta_map <- function(map_file) {
  if (file.exists(map_file)) {
    metabolite_name_map <- data.table::fread(map_file, data.table = FALSE, check.names = FALSE)
    if (!all(c("Meta", "Original_Metabolite") %in% names(metabolite_name_map))) {
      stop("metabolite_name_map must contain columns: Meta and Original_Metabolite.")
    }
    return(setNames(metabolite_name_map$Original_Metabolite, metabolite_name_map$Meta))
  }

  if (exists("meta_reverse_map", envir = .GlobalEnv)) {
    return(get("meta_reverse_map", envir = .GlobalEnv))
  }

  warning("No metabolite name map was found. Raw Meta_* IDs will be used as labels.")
  NULL
}

add_metabolite_names <- function(df, meta_reverse_map, raw_col = "Metabolite_raw") {
  if (is.null(meta_reverse_map)) {
    df$Metabolite_name <- df[[raw_col]]
  } else {
    mapped <- unname(meta_reverse_map[trimws(df[[raw_col]])])
    df$Metabolite_name <- ifelse(is.na(mapped), df[[raw_col]], mapped)
  }
  df
}

get_significant_metabolites <- function(cox_file, significance_column = "p_bonf") {
  cox_df <- read_or_get(cox_file, "cox_df_pcorrected")

  if (!significance_column %in% names(cox_df)) {
    stop("Column '", significance_column, "' was not found in the Cox result.")
  }

  if ("Metabolite_raw" %in% names(cox_df)) {
    raw_ids <- cox_df$Metabolite_raw
  } else if ("Metabolite" %in% names(cox_df) && all(grepl("^Meta_", cox_df$Metabolite))) {
    raw_ids <- cox_df$Metabolite
  } else {
    stop(
      "Cox result must contain raw metabolite IDs. ",
      "Please keep a 'Metabolite_raw' column in results/cox_overall.csv."
    )
  }

  raw_ids[!is.na(cox_df[[significance_column]]) & cox_df[[significance_column]] < 0.05]
}

prepare_matching_data <- function(df, covars) {
  df <- standardize_covariates(df)

  if ("incident_migraine" %in% names(df)) {
    df$migraine_status <- as.integer(df$incident_migraine == 1)
  } else if ("migraine_status" %in% names(df)) {
    df$migraine_status <- as.integer(df$migraine_status == 1)
  } else {
    stop("Input data must contain either incident_migraine or migraine_status.")
  }

  required_vars <- unique(c("eid", "migraine_status", "migraine_years", covars))
  missing_vars <- setdiff(required_vars, names(df))
  if (length(missing_vars) > 0) {
    stop("Missing variables in mydata_base: ", paste(missing_vars, collapse = ", "))
  }

  df %>%
    filter(!is.na(migraine_status), !is.na(migraine_years), migraine_years > 0)
}

run_matching <- function(df, covars, ratio = 10) {
  match_covars <- setdiff(covars, "sex")
  fml <- as.formula(
    paste("migraine_status ~", paste(quote_if_needed(match_covars), collapse = " + "))
  )

  m.out <- MatchIt::matchit(
    formula = fml,
    data = df,
    method = "nearest",
    distance = "mahalanobis",
    ratio = ratio,
    exact = ~ sex
  )

  matched <- MatchIt::match.data(m.out)

  case_time <- matched %>%
    group_by(subclass) %>%
    summarise(
      n_case = sum(migraine_status == 1, na.rm = TRUE),
      case_years = migraine_years[migraine_status == 1][1],
      .groups = "drop"
    )

  if (any(case_time$n_case != 1)) {
    warning("Some matched subclasses do not contain exactly one case.")
  }

  matched %>%
    left_join(case_time, by = "subclass") %>%
    mutate(time_scale = -case_years) %>%
    select(-n_case, -case_years)
}

residualize_metabolites <- function(df, meta_cols, covars, min_n = 10) {
  rhs <- paste(quote_if_needed(covars), collapse = " + ")
  out <- data.frame(eid = df$eid)

  for (m in meta_cols) {
    keep_cols <- unique(c("eid", covars, m))
    model_df <- as.data.frame(df[, keep_cols, drop = FALSE])
    model_df <- model_df[complete.cases(model_df), , drop = FALSE]

    if (nrow(model_df) < min_n) next

    fml <- as.formula(paste0(quote_if_needed(m), " ~ ", rhs))
    fit <- tryCatch(lm(fml, data = model_df), error = function(e) NULL)
    if (is.null(fit)) next

    res_df <- data.frame(
      eid = model_df$eid,
      residual = residuals(fit)
    )
    names(res_df)[2] <- paste0(m, "_res")
    out <- left_join(out, res_df, by = "eid")
  }

  info_cols <- intersect(
    c("eid", "migraine_status", "migraine_years", "subclass", "weights", "time_scale", covars),
    names(df)
  )

  left_join(out, as.data.frame(df[, info_cols, drop = FALSE]), by = "eid")
}

safe_z <- function(x, status) {
  case_val <- x[status == 1]
  ctrl_val <- x[status == 0]

  if (length(case_val) != 1) return(NA_real_)

  sd_ctrl <- sd(ctrl_val, na.rm = TRUE)
  if (is.na(sd_ctrl) || sd_ctrl == 0) return(NA_real_)

  (case_val - mean(ctrl_val, na.rm = TRUE)) / sd_ctrl
}

compute_case_control_z <- function(resid_df) {
  res_cols <- grep("^Meta_.*_res$", names(resid_df), value = TRUE)
  if (length(res_cols) == 0) stop("No residualized metabolite columns were found.")

  z_df <- resid_df %>%
    mutate(subclass = as.factor(subclass)) %>%
    group_by(subclass) %>%
    summarise(
      across(all_of(res_cols), ~ safe_z(.x, migraine_status)),
      Year = first(time_scale),
      .groups = "drop"
    ) %>%
    rename_with(~ str_remove(.x, "_res$"), all_of(res_cols))

  z_df
}

fit_loess_grid <- function(z_df, meta_cols, span = 0.75, step = 0.2, min_n = 10) {
  long_df <- z_df %>%
    pivot_longer(
      cols = all_of(meta_cols),
      names_to = "Metabolite_raw",
      values_to = "Estimate"
    ) %>%
    filter(!is.na(Year), !is.na(Estimate))

  year_grid <- seq(
    from = floor(min(long_df$Year, na.rm = TRUE)),
    to = ceiling(max(long_df$Year, na.rm = TRUE)),
    by = step
  )

  purrr::map_dfr(meta_cols, function(m) {
    df_m <- long_df %>% filter(Metabolite_raw == m)

    if (nrow(df_m) < min_n || length(unique(df_m$Year)) < 5) return(NULL)

    fit <- tryCatch(
      loess(Estimate ~ Year, data = df_m, span = span),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)

    pred_df <- data.frame(Year = year_grid)
    pred <- tryCatch(predict(fit, newdata = pred_df), error = function(e) rep(NA_real_, length(year_grid)))

    tibble(
      Year = year_grid,
      Metabolite_raw = m,
      Estimate_loess = as.numeric(pred)
    )
  })
}

cluster_loess_matrix <- function(loess_df, k = 3, method = "ward.D") {
  mat_df <- loess_df %>%
    mutate(Year = sprintf("%.1f", Year)) %>%
    select(Metabolite_raw, Year, Estimate_loess) %>%
    pivot_wider(names_from = Year, values_from = Estimate_loess)

  mat <- mat_df %>%
    column_to_rownames("Metabolite_raw") %>%
    as.matrix()

  mat <- mat[complete.cases(mat), , drop = FALSE]
  if (nrow(mat) < k) stop("Too few complete metabolite trajectories for clustering.")

  hc <- hclust(dist(mat, method = "euclidean"), method = method)

  cluster_df <- tibble(
    Metabolite_raw = rownames(mat),
    cluster = factor(cutree(hc, k = k), levels = seq_len(k))
  )

  list(matrix = mat, hclust = hc, clusters = cluster_df)
}

plot_heatmap <- function(mat, cluster_df, meta_reverse_map, file) {
  cluster_df <- add_metabolite_names(cluster_df, meta_reverse_map)

  row_order <- cluster_df %>%
    arrange(cluster, Metabolite_name) %>%
    pull(Metabolite_raw)

  mat_plot <- mat[row_order, , drop = FALSE]
  row_split <- cluster_df$cluster[match(row_order, cluster_df$Metabolite_raw)]

  mapped_names <- if (is.null(meta_reverse_map)) row_order else unname(meta_reverse_map[row_order])
  rownames(mat_plot) <- ifelse(is.na(mapped_names), row_order, mapped_names)

  col_fun <- circlize::colorRamp2(c(-1, 0, 1), c("#2b6cb0", "white", "#c53030"))

  ht <- ComplexHeatmap::Heatmap(
    mat_plot,
    name = "Z score",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = row_split,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 6),
    column_names_gp = grid::gpar(fontsize = 7),
    column_title = "Years before diagnosis",
    heatmap_legend_param = list(direction = "horizontal")
  )

  pdf(file, width = 7, height = max(8, 0.13 * nrow(mat_plot) + 3), useDingbats = FALSE)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "top")
  dev.off()
}

plot_cluster_trajectories <- function(loess_df, cluster_df, meta_reverse_map, file) {
  plot_df <- loess_df %>%
    inner_join(cluster_df, by = "Metabolite_raw") %>%
    add_metabolite_names(meta_reverse_map)

  cluster_colors <- c("#C0392B", "#E67E22", "#2980B9", "#16A085", "#8E44AD", "#7F8C8D")
  cluster_levels <- levels(plot_df$cluster)

  p_list <- lapply(cluster_levels, function(cl) {
    df_cl <- plot_df %>% filter(cluster == cl)
    n_meta <- n_distinct(df_cl$Metabolite_raw)
    col <- cluster_colors[(as.integer(cl) - 1) %% length(cluster_colors) + 1]

    ggplot(df_cl, aes(x = Year, y = Estimate_loess, group = Metabolite_raw)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_line(color = col, alpha = 0.30, linewidth = 0.3) +
      stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.2, color = col) +
      annotate(
        "text",
        x = min(df_cl$Year, na.rm = TRUE),
        y = max(df_cl$Estimate_loess, na.rm = TRUE),
        label = paste0("n = ", n_meta),
        hjust = 0,
        vjust = 1,
        size = 4
      ) +
      labs(
        title = paste0("Cluster ", cl),
        x = "Years before diagnosis",
        y = "Z score"
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid = element_blank())
  })

  combined <- patchwork::wrap_plots(p_list, ncol = 1)

  pdf(file, width = 6, height = 4.5 * length(p_list), useDingbats = FALSE)
  print(combined)
  dev.off()

  invisible(combined)
}

##==========================
## Main analysis
##==========================
mydata_base <- read_or_get(data_file, "mydata_base")
meta_reverse_map <- make_meta_map(map_file)

mydata_match <- prepare_matching_data(mydata_base, covars)
meta_cols <- grep("^Meta_", names(mydata_match), value = TRUE)

if (length(meta_cols) == 0) {
  stop("No Meta_* metabolite columns were found in mydata_base.")
}

message("Running nearest-neighbor Mahalanobis matching...")
matched_data <- run_matching(mydata_match, covars, ratio = matching_ratio)
data.table::fwrite(matched_data, "results/match_incident_migraine_metabolome.csv")

message("Residualizing metabolites...")
resid_df <- residualize_metabolites(matched_data, meta_cols, covars)
data.table::fwrite(resid_df, "results/resid_ukb_metabolome.csv")

message("Computing matched case-control Z scores...")
z_df <- compute_case_control_z(resid_df)
data.table::fwrite(z_df, "results/matched_case_control_z_scores.csv")

message("Selecting significant metabolites from Cox results...")
sig_meta <- get_significant_metabolites(cox_file, significance_column = significance_column)
sig_meta <- intersect(sig_meta, names(z_df))

if (length(sig_meta) == 0) {
  stop("No significant metabolites were available for trajectory analysis.")
}

message("Fitting LOESS curves for ", length(sig_meta), " significant metabolites...")
loess_sig <- fit_loess_grid(
  z_df = z_df,
  meta_cols = sig_meta,
  span = loess_span,
  step = loess_step
) %>%
  add_metabolite_names(meta_reverse_map)

data.table::fwrite(loess_sig, "results/loess_trajectories_significant_long.csv")

loess_wide <- loess_sig %>%
  select(Metabolite_raw, Metabolite_name, Year, Estimate_loess) %>%
  mutate(Year = sprintf("%.1f", Year)) %>%
  pivot_wider(names_from = Year, values_from = Estimate_loess)

data.table::fwrite(loess_wide, "results/loess_trajectories_significant_wide.csv")

message("Clustering LOESS-smoothed trajectories...")
clust <- cluster_loess_matrix(loess_sig, k = n_clusters, method = "ward.D")

cluster_df <- clust$clusters %>%
  add_metabolite_names(meta_reverse_map) %>%
  arrange(cluster, Metabolite_name)

data.table::fwrite(cluster_df, "results/cluster_metabolite_list_long.csv")

cluster_summary <- cluster_df %>%
  group_by(cluster) %>%
  summarise(
    n_metabolites = n(),
    metabolites = paste(Metabolite_name, collapse = ", "),
    .groups = "drop"
  )

data.table::fwrite(cluster_summary, "results/cluster_metabolite_collapsed.csv")

cluster_wide <- cluster_df %>%
  group_by(cluster) %>%
  mutate(idx = row_number()) %>%
  ungroup() %>%
  select(cluster, idx, Metabolite_name) %>%
  pivot_wider(names_from = cluster, values_from = Metabolite_name) %>%
  arrange(idx)

data.table::fwrite(cluster_wide, "results/cluster_metabolite_list_wide.csv")

message("Drawing figures...")
plot_heatmap(
  mat = clust$matrix,
  cluster_df = clust$clusters,
  meta_reverse_map = meta_reverse_map,
  file = "figures/Fig2A_heatmap_metabolome.pdf"
)

plot_cluster_trajectories(
  loess_df = loess_sig,
  cluster_df = clust$clusters,
  meta_reverse_map = meta_reverse_map,
  file = "figures/Fig2B_trajectories_metabolome.pdf"
)

message("Done.")
