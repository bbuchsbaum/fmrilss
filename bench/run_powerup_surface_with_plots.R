#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", quiet = TRUE)
}
source("bench/stglmnet_vs_fmrilss_harness.R")
has_hrfals <- requireNamespace("hrfals", quietly = TRUE)
want_hrfals <- identical(Sys.getenv("INCLUDE_HRFALS", "0"), "1")
use_hrfals <- isTRUE(want_hrfals && has_hrfals)
if (!has_hrfals) {
  cat("Package 'hrfals' not available: running without hrfals methods.\n")
} else if (!use_hrfals) {
  cat("hrfals is installed but disabled by default (set INCLUDE_HRFALS=1 to include).\n")
}

`%||%` <- function(a, b) if (is.null(a)) b else a

out_dir <- file.path("bench", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
power_reps <- as.integer(Sys.getenv("POWER_REPS", "40"))
surface_reps <- as.integer(Sys.getenv("SURFACE_REPS", "20"))

# ---------- Shared helpers ----------------------------------------------------

build_orientation_18 <- function() {
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

        nm <- sprintf(
          "h%d_r%d_c%d__%s__%s__%s__%s",
          h, i, j,
          hrf_cfg$label,
          ar_cfg$label,
          isi_cfg$label,
          amp_cfg$label
        )

        scenarios[[nm]] <- c(
          hrf_cfg[setdiff(names(hrf_cfg), "label")],
          ar_cfg[setdiff(names(ar_cfg), "label")],
          isi_cfg[setdiff(names(isi_cfg), "label")],
          amp_cfg[setdiff(names(amp_cfg), "label")]
        )
      }
    }
  }

  scenarios
}

add_orientation_fields <- function(df) {
  df$hrf_regime <- ifelse(grepl("__hrf_nominal__", df$scenario), "nominal", "lag2_widen1p5")
  df$ar_regime <- sub(".*__(ar_[^_]+)__.*", "\\1", df$scenario)
  df$isi_regime <- sub(".*__(isi_[^_]+)__.*", "\\1", df$scenario)
  df$amp_regime <- sub(".*__(amp_[^_]+)$", "\\1", df$scenario)
  df$cell <- sub("^(h[0-9]+_[^_]+_[^_]+)__.*$", "\\1", df$scenario)
  df
}

condense_oasis_per_rep <- function(metrics_long) {
  key <- interaction(metrics_long$scenario, metrics_long$rep, drop = TRUE)
  spl <- split(metrics_long, key)
  out <- lapply(spl, function(d) {
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
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

pairwise_ci <- function(metrics_long, method_a, method_b) {
  da <- metrics_long[metrics_long$method == method_a, , drop = FALSE]
  db <- metrics_long[metrics_long$method == method_b, , drop = FALSE]
  ab <- merge(da, db, by = c("scenario", "rep"), suffixes = c("_a", "_b"), all = FALSE)
  if (nrow(ab) < 1L) return(data.frame())

  split(ab, ab$scenario) |> lapply(function(d) {
    d_cor <- d$cor_a - d$cor_b
    d_rmse <- d$rmse_a - d$rmse_b
    d_sign <- d$sign_acc_a - d$sign_acc_b
    d_abs_slope <- abs(d$slope_bias_a) - abs(d$slope_bias_b)
    d_abs_intercept <- abs(d$intercept_bias_a) - abs(d$intercept_bias_b)

    data.frame(
      scenario = d$scenario[1],
      method_a = method_a,
      method_b = method_b,
      n = nrow(d),
      delta_cor_mean = mean(d_cor, na.rm = TRUE),
      delta_cor_ci_lo = .bmk_empirical_ci(d_cor)[1],
      delta_cor_ci_hi = .bmk_empirical_ci(d_cor)[2],
      delta_rmse_mean = mean(d_rmse, na.rm = TRUE),
      delta_rmse_ci_lo = .bmk_empirical_ci(d_rmse)[1],
      delta_rmse_ci_hi = .bmk_empirical_ci(d_rmse)[2],
      delta_sign_acc_mean = mean(d_sign, na.rm = TRUE),
      delta_sign_acc_ci_lo = .bmk_empirical_ci(d_sign)[1],
      delta_sign_acc_ci_hi = .bmk_empirical_ci(d_sign)[2],
      delta_abs_slope_bias_mean = mean(d_abs_slope, na.rm = TRUE),
      delta_abs_slope_bias_ci_lo = .bmk_empirical_ci(d_abs_slope)[1],
      delta_abs_slope_bias_ci_hi = .bmk_empirical_ci(d_abs_slope)[2],
      delta_abs_intercept_bias_mean = mean(d_abs_intercept, na.rm = TRUE),
      delta_abs_intercept_bias_ci_lo = .bmk_empirical_ci(d_abs_intercept)[1],
      delta_abs_intercept_bias_ci_hi = .bmk_empirical_ci(d_abs_intercept)[2],
      stringsAsFactors = FALSE
    )
  }) |> do.call(what = rbind)
}

# ---------- 1) Power-up pass --------------------------------------------------

cat(sprintf("Running power-up pass (18 scenarios x n_reps=%d)...\n", power_reps))

power_scenarios <- build_orientation_18()
power_base_args <- list(
  n_reps = power_reps,
  n_tp = c(110, 110),
  n_runs = 2,
  n_trials = 36,
  n_vox = 12,
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
  sbhm_oasis = list(ridge_mode = "fractional", ridge_x = 0.02, ridge_b = 0.02),
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

power_suite <- benchmark_stglmnet_vs_fmrilss_suite(
  scenarios = power_scenarios,
  base_args = power_base_args,
  base_seed = 20260401L,
  trace = TRUE
)

power_long <- power_suite$metrics_long
power_long <- add_orientation_fields(power_long)
power_condensed_long <- condense_oasis_per_rep(power_long)
power_condensed_long <- add_orientation_fields(power_condensed_long)

# Summaries
power_summary <- .bmk_suite_method_summary(power_long)
power_condensed_summary <- .bmk_suite_method_summary(power_condensed_long)

power_pairwise <- do.call(rbind, list(
  pairwise_ci(power_condensed_long, "st_default", "fmrilss_lss"),
  pairwise_ci(power_condensed_long, "st_default", "fmrilss_oasis_best"),
  pairwise_ci(power_condensed_long, "st_default", "fmrilss_sbhm")
))
power_pairwise <- add_orientation_fields(power_pairwise)
power_pairwise$comparison <- paste0(power_pairwise$method_a, " - ", power_pairwise$method_b)

# Write power-up outputs
write.csv(power_long, file.path(out_dir, paste0("powerup-18x40-metrics-long-", stamp, ".csv")), row.names = FALSE)
write.csv(power_condensed_long, file.path(out_dir, paste0("powerup-18x40-metrics-condensed-long-", stamp, ".csv")), row.names = FALSE)
write.csv(power_summary, file.path(out_dir, paste0("powerup-18x40-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(power_condensed_summary, file.path(out_dir, paste0("powerup-18x40-condensed-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(power_pairwise, file.path(out_dir, paste0("powerup-18x40-pairwise-ci-", stamp, ".csv")), row.names = FALSE)

# ---------- 2) HRF mismatch surface ------------------------------------------

cat(sprintf("Running HRF mismatch surface (lag x widen x regime, n_reps=%d)...\n", surface_reps))

regimes <- list(
  easy = list(
    noise_model = "iid",
    ar_rho = 0,
    isi_model = "lognormal",
    isi = 6.0,
    isi_params = list(sdlog = 0.5),
    trial_amplitude_model = "none",
    trial_amplitude_cv = 0
  ),
  medium = list(
    noise_model = "ar1",
    ar_rho = 0.35,
    isi_model = "exponential",
    isi = 3.0,
    trial_amplitude_model = "lognormal",
    trial_amplitude_cv = 0.4
  ),
  hard = list(
    noise_model = "ar1",
    ar_rho = 0.6,
    isi_model = "uniform",
    isi = 1.0,
    isi_params = list(width = 0.4),
    trial_amplitude_model = "lognormal",
    trial_amplitude_cv = 0.8
  )
)

lags <- c(0, 1, 2, 3)
widens <- c(0, 1, 2)

surface_scenarios <- list()
for (rg in names(regimes)) {
  for (lag in lags) {
    for (w in widens) {
      nm <- sprintf("%s__lag%s__widen%s", rg, lag, w)
      surface_scenarios[[nm]] <- c(
        regimes[[rg]],
        list(true_hrf_lag_sec = lag, true_hrf_widen_sec = w)
      )
    }
  }
}

surface_base_args <- list(
  n_reps = surface_reps,
  n_tp = c(110, 110),
  n_runs = 2,
  n_trials = 36,
  n_vox = 12,
  cv_folds = 4,
  include_oasis = TRUE,
  oasis_grid = data.frame(
    ridge_mode = "fractional",
    ridge_x = 0.05,
    ridge_b = 0.05,
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
  sbhm_oasis = list(ridge_mode = "fractional", ridge_x = 0.02, ridge_b = 0.02),
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

surface_suite <- benchmark_stglmnet_vs_fmrilss_suite(
  scenarios = surface_scenarios,
  base_args = surface_base_args,
  base_seed = 20260501L,
  trace = TRUE
)

surface_long <- surface_suite$metrics_long
surface_long$regime <- sub("__.*$", "", surface_long$scenario)
surface_long$lag <- as.numeric(sub(".*__lag([0-9.]+)__widen.*$", "\\1", surface_long$scenario))
surface_long$widen <- as.numeric(sub(".*__widen([0-9.]+)$", "\\1", surface_long$scenario))

surface_summary <- .bmk_suite_method_summary(surface_long)
surface_summary$regime <- sub("__.*$", "", surface_summary$scenario)
surface_summary$lag <- as.numeric(sub(".*__lag([0-9.]+)__widen.*$", "\\1", surface_summary$scenario))
surface_summary$widen <- as.numeric(sub(".*__widen([0-9.]+)$", "\\1", surface_summary$scenario))

write.csv(surface_long, file.path(out_dir, paste0("hrf-surface-metrics-long-", stamp, ".csv")), row.names = FALSE)
write.csv(surface_summary, file.path(out_dir, paste0("hrf-surface-summary-", stamp, ".csv")), row.names = FALSE)

# ---------- 3) Plots ----------------------------------------------------------

pretty_methods <- c(
  st_default = "stglmnet",
  fmrilss_lss = "LSS",
  fmrilss_sbhm = "SBHM",
  fmrilss_oasis_best = "OASIS (best)",
  "hrfals_fastlss[r]" = "hrfals fastLSS (R)",
  "hrfals_fastlss[cpp]" = "hrfals fastLSS (C++)",
  "hrfals_fit[ls_svd_only]" = "hrfals fit (ls_svd_only)",
  "hrfals_fit[ls_svd_1als]" = "hrfals fit (ls_svd_1als)",
  "hrfals_fit[cf_als]" = "hrfals fit (cf_als)",
  `fmrilss_oasis[fractional|x=0.05|b=0.05]` = "OASIS (frac 0.05/0.05)"
)

# Power-up PDF
power_pdf <- file.path(out_dir, paste0("powerup-18x40-plots-", stamp, ".pdf"))
pdf(power_pdf, width = 13, height = 8)

pp <- power_condensed_summary
pp <- add_orientation_fields(pp)
pp$method_pretty <- pretty_methods[pp$method] %||% pp$method
pp$cell_pretty <- paste(pp$ar_regime, pp$isi_regime, pp$amp_regime, sep = " | ")

p1 <- ggplot(pp, aes(x = method_pretty, y = reorder(cell_pretty, cell_pretty), fill = mean_cor)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~hrf_regime, ncol = 1, scales = "free_y") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b") +
  labs(
    title = "Power-Up (18 scenarios, n_reps=40): Mean Correlation",
    x = NULL,
    y = "Scenario Cell (AR | ISI | Amplitude)",
    fill = "Mean Cor"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )
print(p1)

pw <- power_pairwise
pw$method_b_pretty <- pretty_methods[pw$method_b] %||% pw$method_b
pw$cell_pretty <- paste(pw$ar_regime, pw$isi_regime, pw$amp_regime, sep = " | ")

p2 <- ggplot(pw, aes(x = delta_cor_mean, y = reorder(cell_pretty, cell_pretty), color = method_b_pretty)) +
  geom_vline(xintercept = 0, linetype = 2, color = "gray40") +
  geom_errorbar(
    aes(xmin = delta_cor_ci_lo, xmax = delta_cor_ci_hi),
    width = 0.15,
    orientation = "y",
    alpha = 0.8
  ) +
  geom_point(size = 2.0) +
  facet_wrap(~hrf_regime, ncol = 1, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  labs(
    title = "Pairwise Delta Correlation with 95% Empirical CI",
    subtitle = "Delta = stglmnet - comparator",
    x = "Delta Correlation",
    y = "Scenario Cell",
    color = "Comparator"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"))
print(p2)

pm <- aggregate(cbind(mean_sign_acc, mean_abs_slope_bias, mean_abs_intercept_bias) ~ hrf_regime + method_pretty, data = pp, FUN = mean, na.rm = TRUE)

p3 <- ggplot(pm, aes(x = method_pretty, y = mean_sign_acc, fill = hrf_regime)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Diagnostic Metric: Sign Accuracy",
    x = NULL,
    y = "Mean Sign Accuracy",
    fill = "HRF Regime"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
print(p3)

p4 <- ggplot(pm, aes(x = method_pretty, y = mean_abs_slope_bias, fill = hrf_regime)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Diagnostic Metric: Absolute Slope Bias (lower is better)",
    x = NULL,
    y = "Mean |Slope Bias|",
    fill = "HRF Regime"
  ) +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
print(p4)

dev.off()

# HRF surface PDF
surface_pdf <- file.path(out_dir, paste0("hrf-surface-plots-", stamp, ".pdf"))
pdf(surface_pdf, width = 13, height = 8)

ss <- surface_summary
ss$method_pretty <- pretty_methods[ss$method] %||% ss$method

p5 <- ggplot(ss, aes(x = factor(lag), y = factor(widen), fill = mean_cor)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(method_pretty ~ regime) +
  scale_fill_gradient(low = "#fff5f0", high = "#67000d") +
  labs(
    title = "HRF Surface: Mean Correlation",
    x = "True HRF Lag (s)",
    y = "True HRF Widen (s)",
    fill = "Mean Cor"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid = element_blank())
print(p5)

p6 <- ggplot(ss, aes(x = factor(lag), y = factor(widen), fill = mean_sign_acc)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(method_pretty ~ regime) +
  scale_fill_gradient(low = "#f7fcf5", high = "#00441b") +
  labs(
    title = "HRF Surface: Mean Sign Accuracy",
    x = "True HRF Lag (s)",
    y = "True HRF Widen (s)",
    fill = "Sign Acc"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid = element_blank())
print(p6)

p7 <- ggplot(ss, aes(x = factor(lag), y = factor(widen), fill = mean_abs_slope_bias)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(method_pretty ~ regime) +
  scale_fill_gradient(low = "#f7f7f7", high = "#252525") +
  labs(
    title = "HRF Surface: Mean Absolute Slope Bias (lower is better)",
    x = "True HRF Lag (s)",
    y = "True HRF Widen (s)",
    fill = "|Slope Bias|"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid = element_blank())
print(p7)

p8 <- ggplot(ss, aes(x = factor(lag), y = factor(widen), fill = mean_rmse)) +
  geom_tile(color = "white", linewidth = 0.25) +
  facet_grid(method_pretty ~ regime) +
  scale_fill_gradient(low = "#ffffe5", high = "#662506") +
  labs(
    title = "HRF Surface: Mean RMSE",
    x = "True HRF Lag (s)",
    y = "True HRF Widen (s)",
    fill = "RMSE"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"), panel.grid = element_blank())
print(p8)

dev.off()

cat("\nCompleted. Files written:\n")
cat(" -", normalizePath(power_pdf), "\n")
cat(" -", normalizePath(surface_pdf), "\n")
cat(" -", normalizePath(file.path(out_dir, paste0("powerup-18x40-pairwise-ci-", stamp, ".csv"))), "\n")
cat(" -", normalizePath(file.path(out_dir, paste0("hrf-surface-summary-", stamp, ".csv"))), "\n")
