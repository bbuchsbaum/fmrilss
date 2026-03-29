#!/usr/bin/env Rscript

# Orientation benchmark: 3x3 scenario grid with all core methods.
# Rows: AR noise regimes
# Cols: ISI regimes
# Amplitude variability cycles across cells to sample low/medium/high regimes.
# We run this grid for two HRF regimes:
# - nominal canonical HRF
# - shifted+widened true HRF (model mismatch)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}
source("bench/benchmark_harness.R")
has_hrfals <- requireNamespace("hrfals", quietly = TRUE)
want_hrfals <- identical(Sys.getenv("INCLUDE_HRFALS", "0"), "1")
use_hrfals <- isTRUE(want_hrfals && has_hrfals)
if (!has_hrfals) {
  cat("Package 'hrfals' not available: running without hrfals methods.\n")
} else if (!use_hrfals) {
  cat("hrfals is installed but disabled by default (set INCLUDE_HRFALS=1 to include).\n")
}

ar_regimes <- list(
  list(label = "ar_iid", noise_model = "iid", ar_rho = 0),
  list(label = "ar_mod", noise_model = "ar1", ar_rho = 0.35),
  list(label = "ar_high", noise_model = "ar1", ar_rho = 0.60)
)

isi_regimes <- list(
  list(label = "isi_dense", isi_model = "uniform", isi = 1.0, isi_params = list(width = 0.4)),
  list(label = "isi_mid", isi_model = "exponential", isi = 3.0),
  list(label = "isi_sparse", isi_model = "lognormal", isi = 6.0, isi_params = list(sdlog = 0.7))
)

amp_regimes <- list(
  list(label = "amp_none", trial_amplitude_model = "none", trial_amplitude_cv = 0.0),
  list(label = "amp_mod", trial_amplitude_model = "lognormal", trial_amplitude_cv = 0.4),
  list(label = "amp_high", trial_amplitude_model = "lognormal", trial_amplitude_cv = 0.8)
)

hrf_regimes <- list(
  list(label = "hrf_nominal", true_hrf_lag_sec = 0, true_hrf_widen_sec = 0),
  list(label = "hrf_lag2_widen1p5", true_hrf_lag_sec = 2.0, true_hrf_widen_sec = 1.5)
)

scenarios <- list()
for (h in seq_along(hrf_regimes)) {
  hrf_cfg <- hrf_regimes[[h]]
  for (i in seq_along(ar_regimes)) {
    for (j in seq_along(isi_regimes)) {
      k <- ((i + j - 2) %% length(amp_regimes)) + 1L
      ar_cfg <- ar_regimes[[i]]
      isi_cfg <- isi_regimes[[j]]
      amp_cfg <- amp_regimes[[k]]

      scenario_name <- sprintf(
        "h%d_r%d_c%d__%s__%s__%s__%s",
        h, i, j,
        hrf_cfg$label,
        ar_cfg$label,
        isi_cfg$label,
        amp_cfg$label
      )

      scenarios[[scenario_name]] <- c(
        hrf_cfg[setdiff(names(hrf_cfg), "label")],
        ar_cfg[setdiff(names(ar_cfg), "label")],
        isi_cfg[setdiff(names(isi_cfg), "label")],
        amp_cfg[setdiff(names(amp_cfg), "label")]
      )
    }
  }
}

base_args <- list(
  n_reps = 6,
  n_tp = c(110, 110),
  n_runs = 2,
  n_trials = 36,
  n_vox = 12,
  n_conditions = 3,
  condition_effect = 1.0,
  cv_folds = 4,
  include_oasis = TRUE,
  oasis_grid = expand.grid(
    ridge_mode = c("absolute", "fractional"),
    ridge_x = c(0.01, 0.02, 0.05),
    ridge_b = c(0.01, 0.02, 0.05),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ),
  include_hrfals = use_hrfals,
  hrfals_modes = c("r", "cpp"),
  hrfals_fit_methods = c("ls_svd_only", "ls_svd_1als", "cf_als"),
  hrfals_basis = "spmg2",
  include_sbhm = TRUE,
  sbhm_rank = 6,
  sbhm_match = list(
    shrink = list(tau = 0.15, ref = NULL, snr = NULL),
    topK = 3,
    soft_blend = TRUE,
    blend_margin = 0.08,
    whiten = FALSE,
    sv_floor_rel = 0.05,
    whiten_power = 0.5,
    min_margin = 0.08,
    min_beta_norm = 1e-3,
    orient_ref = TRUE,
    alpha_source = "prepass"
  ),
  sbhm_oasis = list(
    ridge_mode = "fractional",
    ridge_x = 0.02,
    ridge_b = 0.02
  ),
  sbhm_amplitude = list(
    method = "oasis_voxel",
    ridge = list(mode = "fractional", lambda = 0.02),
    ridge_frac = list(x = 0.02, b = 0.02),
    cond_gate = NULL,
    adaptive = list(enable = FALSE, base = 0.02, k0 = 1000, max = 0.08),
    return_se = FALSE
  ),
  trace = FALSE
)

suite <- benchmark_lss_methods_suite(
  scenarios = scenarios,
  base_args = base_args,
  base_seed = 20260310L,
  trace = TRUE
)

sm <- suite$scenario_method_summary
winners <- do.call(rbind, lapply(split(sm, sm$scenario), function(d) {
  d <- d[order(d$mean_cor, decreasing = TRUE), , drop = FALSE]
  d[1L, c("scenario", "method", "mean_cor", "mean_rmse", "success_rate"), drop = FALSE]
}))
rownames(winners) <- NULL

# Condense OASIS variants into an OASIS-best row per scenario for easier method-level comparison.
condensed <- do.call(rbind, lapply(split(sm, sm$scenario), function(d) {
  non_oasis <- d[!grepl("^fmrilss_oasis", d$method), , drop = FALSE]
  oasis <- d[grepl("^fmrilss_oasis", d$method), , drop = FALSE]
  if (nrow(oasis) > 0L) {
    best <- oasis[order(oasis$mean_cor, decreasing = TRUE), , drop = FALSE][1L, , drop = FALSE]
    best$method <- "fmrilss_oasis_best"
    best$oasis_variant <- as.character(oasis$method[order(oasis$mean_cor, decreasing = TRUE)][1L])
    non_oasis$oasis_variant <- NA_character_
    rbind(non_oasis, best)
  } else {
    non_oasis$oasis_variant <- NA_character_
    non_oasis
  }
}))
rownames(condensed) <- NULL

condensed_winners <- do.call(rbind, lapply(split(condensed, condensed$scenario), function(d) {
  d <- d[order(d$mean_cor, decreasing = TRUE), , drop = FALSE]
  d[1L, c("scenario", "method", "mean_cor", "mean_rmse", "success_rate", "oasis_variant"), drop = FALSE]
}))
rownames(condensed_winners) <- NULL

win_counts <- as.data.frame(table(winners$method), stringsAsFactors = FALSE)
names(win_counts) <- c("method", "n_scenarios_won")
win_counts <- win_counts[order(win_counts$n_scenarios_won, decreasing = TRUE), , drop = FALSE]

overall <- aggregate(mean_cor ~ method, data = sm, FUN = function(x) mean(x, na.rm = TRUE))
overall <- overall[order(overall$mean_cor, decreasing = TRUE), , drop = FALSE]

condensed_win_counts <- as.data.frame(table(condensed_winners$method), stringsAsFactors = FALSE)
names(condensed_win_counts) <- c("method", "n_scenarios_won")
condensed_win_counts <- condensed_win_counts[order(condensed_win_counts$n_scenarios_won, decreasing = TRUE), , drop = FALSE]

overall_condensed <- aggregate(mean_cor ~ method, data = condensed, FUN = function(x) mean(x, na.rm = TRUE))
overall_condensed <- overall_condensed[order(overall_condensed$mean_cor, decreasing = TRUE), , drop = FALSE]

cat("\n=== 3x3 Orientation Grid x HRF Regimes: Top method by scenario (mean correlation) ===\n")
print(winners)

cat("\n=== Win counts across all scenarios (full method set) ===\n")
print(win_counts)

cat("\n=== Overall mean correlation (full method set) ===\n")
print(overall)

cat("\n=== Condensed winners (OASIS collapsed to best variant per scenario) ===\n")
print(condensed_winners)

cat("\n=== Condensed win counts ===\n")
print(condensed_win_counts)

cat("\n=== Overall mean correlation (condensed methods) ===\n")
print(overall_condensed)

cat("\n=== Pairwise stglmnet vs lss summary ===\n")
print(suite$pairwise_summary[order(suite$pairwise_summary$delta_cor_mean, decreasing = TRUE), ])

out_dir <- file.path("bench", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
write.csv(sm, file.path(out_dir, paste0("orientation-3x3-scenario-method-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(winners, file.path(out_dir, paste0("orientation-3x3-winners-", stamp, ".csv")), row.names = FALSE)
write.csv(win_counts, file.path(out_dir, paste0("orientation-3x3-win-counts-", stamp, ".csv")), row.names = FALSE)
write.csv(overall, file.path(out_dir, paste0("orientation-3x3-overall-", stamp, ".csv")), row.names = FALSE)
write.csv(condensed, file.path(out_dir, paste0("orientation-3x3-condensed-method-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(condensed_winners, file.path(out_dir, paste0("orientation-3x3-condensed-winners-", stamp, ".csv")), row.names = FALSE)
write.csv(condensed_win_counts, file.path(out_dir, paste0("orientation-3x3-condensed-win-counts-", stamp, ".csv")), row.names = FALSE)
write.csv(overall_condensed, file.path(out_dir, paste0("orientation-3x3-condensed-overall-", stamp, ".csv")), row.names = FALSE)
write.csv(suite$pairwise_summary, file.path(out_dir, paste0("orientation-3x3-pairwise-", stamp, ".csv")), row.names = FALSE)

cat("\nSaved results to:", normalizePath(out_dir), "\n")
