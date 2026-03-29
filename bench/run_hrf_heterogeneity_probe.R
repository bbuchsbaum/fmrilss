#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}
source("bench/stglmnet_vs_fmrilss_harness.R")

out_dir <- file.path("bench", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
n_reps <- as.integer(Sys.getenv("HET_REPS", "12"))

regimes <- list(
  medium = list(
    noise_model = "ar1",
    ar_rho = 0.35,
    isi_model = "exponential",
    isi = 1.3,
    trial_amplitude_model = "lognormal",
    trial_amplitude_cv = 0.4
  ),
  hard = list(
    noise_model = "ar1",
    ar_rho = 0.6,
    isi_model = "uniform",
    isi = 0.9,
    isi_params = list(width = 0.5),
    trial_amplitude_model = "lognormal",
    trial_amplitude_cv = 0.8
  )
)

global_hrf <- list(
  nominal = list(true_hrf_lag_sec = 0, true_hrf_widen_sec = 0),
  mismatch = list(true_hrf_lag_sec = 2.0, true_hrf_widen_sec = 1.5)
)

hetero_levels <- list(
  none = list(
    true_hrf_lag_vox_sd = 0,
    true_hrf_widen_vox_sd = 0,
    true_hrf_lag_trial_sd = 0,
    true_hrf_widen_trial_sd = 0
  ),
  moderate = list(
    true_hrf_lag_vox_sd = 0.6,
    true_hrf_widen_vox_sd = 0.6,
    true_hrf_lag_trial_sd = 0.25,
    true_hrf_widen_trial_sd = 0.25
  ),
  high = list(
    true_hrf_lag_vox_sd = 1.0,
    true_hrf_widen_vox_sd = 1.0,
    true_hrf_lag_trial_sd = 0.45,
    true_hrf_widen_trial_sd = 0.45
  )
)

scenarios <- list()
for (r_nm in names(regimes)) {
  for (g_nm in names(global_hrf)) {
    for (h_nm in names(hetero_levels)) {
      nm <- paste(r_nm, g_nm, h_nm, sep = "__")
      scenarios[[nm]] <- c(regimes[[r_nm]], global_hrf[[g_nm]], hetero_levels[[h_nm]])
    }
  }
}

base_args <- list(
  n_reps = n_reps,
  n_tp = c(110, 110),
  n_runs = 2,
  n_trials = 36,
  n_vox = 12,
  cv_folds = 4,
  include_oasis = TRUE,
  oasis_grid = expand.grid(
    ridge_mode = c("absolute", "fractional"),
    ridge_x = c(0.02, 0.05),
    ridge_b = c(0.02, 0.05),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ),
  include_sbhm = TRUE,
  trace = FALSE
)

suite <- benchmark_stglmnet_vs_fmrilss_suite(
  scenarios = scenarios,
  base_args = base_args,
  base_seed = 20260415L,
  trace = TRUE
)

metrics <- suite$metrics_long
parts <- strsplit(as.character(metrics$scenario), "__", fixed = TRUE)
metrics$regime <- vapply(parts, `[`, character(1), 1L)
metrics$global_hrf <- vapply(parts, `[`, character(1), 2L)
metrics$hetero <- vapply(parts, `[`, character(1), 3L)
metrics$hetero <- factor(metrics$hetero, levels = c("none", "moderate", "high"))

key <- interaction(metrics$scenario, metrics$rep, drop = TRUE)
condensed <- do.call(rbind, lapply(split(metrics, key), function(d) {
  non <- d[!grepl("^fmrilss_oasis", d$method), , drop = FALSE]
  oo <- d[grepl("^fmrilss_oasis", d$method), , drop = FALSE]
  if (nrow(oo) > 0L) {
    oo <- oo[order(oo$cor, decreasing = TRUE), , drop = FALSE]
    best <- oo[1L, , drop = FALSE]
    best$oasis_variant <- as.character(best$method)
    best$method <- "fmrilss_oasis_best"
    non$oasis_variant <- NA_character_
    rbind(non, best)
  } else {
    non$oasis_variant <- NA_character_
    non
  }
}))
rownames(condensed) <- NULL
parts2 <- strsplit(as.character(condensed$scenario), "__", fixed = TRUE)
condensed$regime <- vapply(parts2, `[`, character(1), 1L)
condensed$global_hrf <- vapply(parts2, `[`, character(1), 2L)
condensed$hetero <- factor(vapply(parts2, `[`, character(1), 3L), levels = c("none", "moderate", "high"))

summary_rows <- do.call(rbind, lapply(split(condensed, interaction(condensed$scenario, condensed$method, drop = TRUE)), function(d) {
  data.frame(
    scenario = d$scenario[1],
    regime = d$regime[1],
    global_hrf = d$global_hrf[1],
    hetero = d$hetero[1],
    method = d$method[1],
    mean_cor = mean(d$cor, na.rm = TRUE),
    ci_cor_lo = .bmk_empirical_ci(d$cor)[1],
    ci_cor_hi = .bmk_empirical_ci(d$cor)[2],
    mean_rmse = mean(d$rmse, na.rm = TRUE),
    mean_sign_acc = mean(d$sign_acc, na.rm = TRUE),
    mean_abs_slope_bias = mean(abs(d$slope_bias), na.rm = TRUE),
    success_rate = mean(d$ok, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
rownames(summary_rows) <- NULL

winners <- do.call(rbind, lapply(split(summary_rows, summary_rows$scenario), function(d) {
  d <- d[order(d$mean_cor, decreasing = TRUE), , drop = FALSE]
  d[1L, c("scenario", "regime", "global_hrf", "hetero", "method", "mean_cor", "success_rate"), drop = FALSE]
}))
rownames(winners) <- NULL

pdf(file.path(out_dir, paste0("hrf-heterogeneity-probe-plots-", stamp, ".pdf")), width = 14, height = 10)

print(
  ggplot(summary_rows, aes(x = hetero, y = mean_cor, color = method, group = method)) +
    geom_line(linewidth = 0.9, alpha = 0.9) +
    geom_point(size = 2.2) +
    geom_errorbar(aes(ymin = ci_cor_lo, ymax = ci_cor_hi), width = 0.08, alpha = 0.5) +
    facet_grid(regime ~ global_hrf) +
    labs(
      title = "Beta Correlation Across HRF Heterogeneity Levels",
      x = "HRF heterogeneity level",
      y = "Mean beta correlation"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
)

print(
  ggplot(summary_rows, aes(x = hetero, y = mean_abs_slope_bias, color = method, group = method)) +
    geom_line(linewidth = 0.9, alpha = 0.9) +
    geom_point(size = 2.2) +
    facet_grid(regime ~ global_hrf) +
    labs(
      title = "Absolute Slope Bias Across HRF Heterogeneity Levels",
      x = "HRF heterogeneity level",
      y = "Mean |slope bias|"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
)

print(
  ggplot(winners, aes(x = hetero, y = regime, fill = method, label = method)) +
    geom_tile(color = "white") +
    geom_text(size = 3.2) +
    facet_wrap(~global_hrf) +
    labs(
      title = "Winner by Scenario (Highest Mean Correlation)",
      x = "HRF heterogeneity level",
      y = "Regime"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
)

dev.off()

write.csv(metrics, file.path(out_dir, paste0("hrf-heterogeneity-probe-metrics-long-", stamp, ".csv")), row.names = FALSE)
write.csv(condensed, file.path(out_dir, paste0("hrf-heterogeneity-probe-metrics-condensed-", stamp, ".csv")), row.names = FALSE)
write.csv(summary_rows, file.path(out_dir, paste0("hrf-heterogeneity-probe-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(winners, file.path(out_dir, paste0("hrf-heterogeneity-probe-winners-", stamp, ".csv")), row.names = FALSE)
if (is.data.frame(suite$true_hrf_diag_long) && nrow(suite$true_hrf_diag_long) > 0L) {
  write.csv(suite$true_hrf_diag_long, file.path(out_dir, paste0("hrf-heterogeneity-probe-true-hrf-diag-", stamp, ".csv")), row.names = FALSE)
}
if (is.data.frame(suite$sbhm_diag_long) && nrow(suite$sbhm_diag_long) > 0L) {
  write.csv(suite$sbhm_diag_long, file.path(out_dir, paste0("hrf-heterogeneity-probe-sbhm-diag-", stamp, ".csv")), row.names = FALSE)
}

cat("\nSaved heterogeneity probe outputs to:", normalizePath(out_dir), "\n")
