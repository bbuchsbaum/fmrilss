#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

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

`%||%` <- function(a, b) if (is.null(a)) b else a

lookup_pretty <- function(x, map) {
  out <- unname(map[as.character(x)])
  miss <- is.na(out) | !nzchar(out)
  out[miss] <- as.character(x)[miss]
  out
}

out_dir <- file.path("bench", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
power_reps <- as.integer(Sys.getenv("POWER_REPS", "40"))
surface_reps <- as.integer(Sys.getenv("SURFACE_REPS", "20"))
run_profile <- if (power_reps <= 2L && surface_reps <= 2L) "smoke" else "full"
run_profile_label <- if (run_profile == "smoke") {
  sprintf("Smoke run only (%d power reps; %d HRF-surface reps). Use for plumbing, not conclusions.", power_reps, surface_reps)
} else {
  sprintf("Full benchmark run (%d power reps; %d HRF-surface reps). Suitable for recommendation summaries.", power_reps, surface_reps)
}

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

add_orientation_labels <- function(df) {
  if (!is.data.frame(df) || nrow(df) < 1L) return(df)

  hrf_map <- c(
    nominal = "Well-specified HRF",
    lag2_widen1p5 = "Misspecified HRF (+2 s lag, +1.5 s wider)"
  )
  ar_map <- c(
    ar_iid = "Independent noise",
    ar_mod = "Moderate AR(1) noise",
    ar_high = "High AR(1) noise"
  )
  isi_map <- c(
    isi_dense = "Dense timing",
    isi_mid = "Intermediate timing",
    isi_sparse = "Sparse timing"
  )
  amp_map <- c(
    amp_none = "No amplitude variability",
    amp_mod = "Moderate amplitude variability",
    amp_high = "High amplitude variability"
  )

  df$hrf_label <- lookup_pretty(df$hrf_regime, hrf_map)
  df$ar_label <- lookup_pretty(df$ar_regime, ar_map)
  df$isi_label <- lookup_pretty(df$isi_regime, isi_map)
  df$amp_label <- lookup_pretty(df$amp_regime, amp_map)
  df$cell_label <- paste(df$ar_label, df$isi_label, df$amp_label, sep = " / ")
  df$scenario_label <- paste(df$cell_label, "|", df$hrf_label)
  df
}

add_surface_labels <- function(df) {
  if (!is.data.frame(df) || nrow(df) < 1L) return(df)

  regime_map <- c(
    easy = "Easy setting",
    medium = "Moderate setting",
    hard = "Hard setting"
  )

  df$regime_label <- lookup_pretty(df$regime, regime_map)
  df$widen_label <- paste0("HRF widen = ", df$widen, " s")
  df$lag_label <- paste0("HRF lag = ", df$lag, " s")
  df$mismatch_band <- ifelse(
    df$lag == 0 & df$widen == 0,
    "Well-specified",
    ifelse(df$lag + df$widen <= 2, "Moderate mismatch", "Strong mismatch")
  )
  df
}

top_recommendation_rows <- function(df, group_cols) {
  if (!is.data.frame(df) || nrow(df) < 1L) return(data.frame())
  key <- interaction(df[, group_cols, drop = FALSE], drop = TRUE, lex.order = TRUE)
  out <- lapply(split(df, key), function(d) {
    d <- d[order(-d$n, -d$mean_margin, -d$mean_score), , drop = FALSE]
    d[1L, , drop = FALSE]
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

write_benchmark_brief <- function(
  path,
  run_profile_label,
  power_recommend,
  surface_recommend,
  recommendation_summary,
  st_variant_summary = data.frame()
) {
  lines <- c(
    "# Benchmark Brief",
    "",
    run_profile_label,
    "",
    "## Terms",
    "- Well-specified HRF: the estimator uses the same HRF family that generated the signal.",
    "- Misspecified HRF (+2 s lag, +1.5 s wider): the true signal is slower and delayed relative to the estimation basis.",
    "- Dense/intermediate/sparse timing: events are packed closely together, moderately spaced, or well separated in time.",
    "- Balanced recommendation: winner from the benchmark's combined score over recovery, bias, stability, classification, and speed.",
    ""
  )

  power_bal <- power_recommend[power_recommend$profile == "balanced", , drop = FALSE]
  if (nrow(power_bal) > 0L) {
    power_counts <- aggregate(
      cbind(n = score_margin, mean_margin = score_margin, mean_score = overall_score) ~ hrf_label + recommended_pretty,
      data = transform(power_bal, n = 1),
      FUN = function(x) c(sum = length(x), mean = mean(x, na.rm = TRUE))
    )
    power_counts <- data.frame(
      hrf_label = power_counts$hrf_label,
      recommended_pretty = power_counts$recommended_pretty,
      n = power_counts$n[, "sum"],
      mean_margin = power_counts$mean_margin[, "mean"],
      mean_score = power_counts$mean_score[, "mean"],
      stringsAsFactors = FALSE
    )
    power_top <- top_recommendation_rows(power_counts, "hrf_label")

    lines <- c(lines, "## Power-Up Recommendations")
    for (i in seq_len(nrow(power_top))) {
      lines <- c(
        lines,
        sprintf(
          "- %s: %s is the most common balanced recommendation in %d scenarios, with mean winner margin %.3f.",
          power_top$hrf_label[i],
          power_top$recommended_pretty[i],
          power_top$n[i],
          power_top$mean_margin[i]
        )
      )
    }
    lines <- c(lines, "")
  }

  surface_bal <- surface_recommend[surface_recommend$profile == "balanced", , drop = FALSE]
  if (nrow(surface_bal) > 0L) {
    surface_counts <- aggregate(
      cbind(n = score_margin, mean_margin = score_margin, mean_score = overall_score) ~ regime_label + mismatch_band + recommended_pretty,
      data = transform(surface_bal, n = 1),
      FUN = function(x) c(sum = length(x), mean = mean(x, na.rm = TRUE))
    )
    surface_counts <- data.frame(
      regime_label = surface_counts$regime_label,
      mismatch_band = surface_counts$mismatch_band,
      recommended_pretty = surface_counts$recommended_pretty,
      n = surface_counts$n[, "sum"],
      mean_margin = surface_counts$mean_margin[, "mean"],
      mean_score = surface_counts$mean_score[, "mean"],
      stringsAsFactors = FALSE
    )
    surface_top <- top_recommendation_rows(surface_counts, c("regime_label", "mismatch_band"))

    lines <- c(lines, "## HRF-Mismatch Guidance")
    for (i in seq_len(nrow(surface_top))) {
      lines <- c(
        lines,
        sprintf(
          "- %s, %s: %s is the most common recommendation in %d cells, with mean winner margin %.3f.",
          surface_top$regime_label[i],
          surface_top$mismatch_band[i],
          surface_top$recommended_pretty[i],
          surface_top$n[i],
          surface_top$mean_margin[i]
        )
      )
    }
    lines <- c(lines, "")
  }

  if (is.data.frame(recommendation_summary) && nrow(recommendation_summary) > 0L) {
    rs <- recommendation_summary[recommendation_summary$profile == "balanced", , drop = FALSE]
    rs <- rs[order(-rs$n_scenarios, -rs$mean_score), , drop = FALSE]
    if (nrow(rs) > 0L) {
      lines <- c(lines, "## Overall Balanced Ranking")
      for (i in seq_len(min(5L, nrow(rs)))) {
        lines <- c(
          lines,
          sprintf(
            "- %s: selected in %d scenarios; mean score %.3f; mean elapsed %.3f s.",
            rs$method[i],
            rs$n_scenarios[i],
            rs$mean_score[i],
            rs$mean_elapsed_sec[i]
          )
        )
      }
      lines <- c(lines, "")
    }
  }

  if (is.data.frame(st_variant_summary) && nrow(st_variant_summary) > 0L) {
    lines <- c(lines, "## Internal stglmnet Variants")
    for (i in seq_len(min(5L, nrow(st_variant_summary)))) {
      lines <- c(
        lines,
        sprintf(
          "- %s: best internal stglmnet variant in %d scenarios; mean correlation %.3f.",
          st_variant_summary$method_pretty[i],
          st_variant_summary$n_scenarios[i],
          st_variant_summary$mean_cor[i]
        )
      )
    }
    lines <- c(lines, "")
  }

  writeLines(lines, con = path, useBytes = TRUE)
}

summarize_best_st_variant <- function(df, pretty_methods) {
  if (!is.data.frame(df) || nrow(df) < 1L) return(data.frame())
  st_df <- df[grepl("^st_", df$method), , drop = FALSE]
  if (nrow(st_df) < 1L) return(data.frame())

  split_by_scenario <- split(st_df, st_df$scenario)
  winners <- lapply(split_by_scenario, function(d) {
    d <- d[order(-d$mean_cor, d$mean_rmse, d$mean_elapsed_sec, na.last = TRUE), , drop = FALSE]
    d[1L, , drop = FALSE]
  })
  winners <- do.call(rbind, winners)
  rownames(winners) <- NULL
  winners$method_pretty <- lookup_pretty(winners$method, pretty_methods)

  agg <- aggregate(
    cbind(n_scenarios = mean_cor, mean_cor = mean_cor, mean_rmse = mean_rmse) ~ method_pretty,
    data = transform(winners, n_scenarios = 1),
    FUN = function(x) c(sum = length(x), mean = mean(x, na.rm = TRUE))
  )
  out <- data.frame(
    method_pretty = agg$method_pretty,
    n_scenarios = agg$n_scenarios[, "sum"],
    mean_cor = agg$mean_cor[, "mean"],
    mean_rmse = agg$mean_rmse[, "mean"],
    stringsAsFactors = FALSE
  )
  out[order(-out$n_scenarios, -out$mean_cor), , drop = FALSE]
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
  n_conditions = 3,
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

power_suite <- benchmark_lss_methods_suite(
  scenarios = power_scenarios,
  base_args = power_base_args,
  base_seed = 20260401L,
  trace = TRUE
)

power_long <- power_suite$metrics_long
power_long <- add_orientation_fields(power_long)
power_condensed_long <- condense_oasis_per_rep(power_long)
power_condensed_long <- add_orientation_fields(power_condensed_long)
primary_st_name <- unique(power_suite$pairwise_summary$primary_st)
primary_st_name <- primary_st_name[is.finite(nchar(primary_st_name))]
if (length(primary_st_name) < 1L) {
  primary_st_name <- unique(vapply(power_suite$results, function(x) x$primary_st, character(1)))
}
primary_st_name <- primary_st_name[1L]

# Summaries
power_summary <- .bmk_suite_method_summary(power_long)
power_condensed_summary <- .bmk_suite_method_summary(power_condensed_long)

power_pairwise <- .bmk_bind_rows(list(
  pairwise_ci(power_condensed_long, primary_st_name, "fmrilss_lss"),
  pairwise_ci(power_condensed_long, primary_st_name, "fmrilss_oasis_best"),
  pairwise_ci(power_condensed_long, primary_st_name, "fmrilss_sbhm")
))
if (is.null(power_pairwise)) {
  power_pairwise <- data.frame()
} else {
  power_pairwise <- add_orientation_fields(power_pairwise)
  power_pairwise$comparison <- paste0(power_pairwise$method_a, " - ", power_pairwise$method_b)
}
power_recommend <- power_suite$recommendations
power_recommend <- power_recommend[power_recommend$profile == "balanced", , drop = FALSE]
power_recommend <- add_orientation_fields(power_recommend)
power_recommend <- add_orientation_labels(power_recommend)

# Write power-up outputs
write.csv(power_long, file.path(out_dir, paste0("powerup-18x40-metrics-long-", stamp, ".csv")), row.names = FALSE)
write.csv(power_condensed_long, file.path(out_dir, paste0("powerup-18x40-metrics-condensed-long-", stamp, ".csv")), row.names = FALSE)
write.csv(power_summary, file.path(out_dir, paste0("powerup-18x40-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(power_condensed_summary, file.path(out_dir, paste0("powerup-18x40-condensed-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(power_pairwise, file.path(out_dir, paste0("powerup-18x40-pairwise-ci-", stamp, ".csv")), row.names = FALSE)
write.csv(power_suite$recommendations, file.path(out_dir, paste0("powerup-18x40-recommendations-", stamp, ".csv")), row.names = FALSE)
write.csv(power_suite$recommendation_summary, file.path(out_dir, paste0("powerup-18x40-recommendation-summary-", stamp, ".csv")), row.names = FALSE)

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
  n_conditions = 3,
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

surface_suite <- benchmark_lss_methods_suite(
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
surface_summary <- add_surface_labels(surface_summary)
surface_recommend <- surface_suite$recommendations
surface_recommend <- surface_recommend[surface_recommend$profile == "balanced", , drop = FALSE]
surface_recommend$regime <- sub("__.*$", "", surface_recommend$scenario)
surface_recommend$lag <- surface_recommend$true_hrf_lag_sec
surface_recommend$widen <- surface_recommend$true_hrf_widen_sec
surface_recommend <- add_surface_labels(surface_recommend)

write.csv(surface_long, file.path(out_dir, paste0("hrf-surface-metrics-long-", stamp, ".csv")), row.names = FALSE)
write.csv(surface_summary, file.path(out_dir, paste0("hrf-surface-summary-", stamp, ".csv")), row.names = FALSE)
write.csv(surface_suite$recommendations, file.path(out_dir, paste0("hrf-surface-recommendations-", stamp, ".csv")), row.names = FALSE)
write.csv(surface_suite$recommendation_summary, file.path(out_dir, paste0("hrf-surface-recommendation-summary-", stamp, ".csv")), row.names = FALSE)

# ---------- 3) Plots ----------------------------------------------------------

pretty_methods <- c(
  st_pool_compcv_wnever = "stglmnet: pooled CV (no overlap weighting)",
  st_pool_compcv = "stglmnet: pooled CV",
  st_pool_joint_mult = "stglmnet: joint multiplier",
  st_pool_joint_mult_1se = "stglmnet: joint multiplier (1SE)",
  st_pool_compcv_wB = "stglmnet: pooled CV + B weighting",
  st_pool_mean = "stglmnet: mean pooled",
  fmrilss_lss = "LSS",
  fmrilss_sbhm = "SBHM",
  fmrilss_oasis_best = "OASIS (best)",
  `fmrilss_oasis[K=1|fractional|x=0.05|b=0.05]` = "OASIS",
  "hrfals_fastlss[r]" = "hrfals fastLSS (R)",
  "hrfals_fastlss[cpp]" = "hrfals fastLSS (C++)",
  "hrfals_fit[ls_svd_only]" = "hrfals fit (ls_svd_only)",
  "hrfals_fit[ls_svd_1als]" = "hrfals fit (ls_svd_1als)",
  "hrfals_fit[cf_als]" = "hrfals fit (cf_als)",
  `fmrilss_oasis[fractional|x=0.05|b=0.05]` = "OASIS (frac 0.05/0.05)"
)
pretty_methods[primary_st_name] <- "stglmnet (default)"

power_display_methods <- intersect(
  c(primary_st_name, "fmrilss_lss", "fmrilss_oasis_best", "fmrilss_sbhm"),
  unique(power_condensed_summary$method)
)
surface_display_methods <- intersect(
  c(primary_st_name, "fmrilss_lss", "fmrilss_oasis[K=1|fractional|x=0.05|b=0.05]", "fmrilss_sbhm"),
  unique(surface_summary$method)
)

recommendation_brief <- file.path(out_dir, paste0("benchmark-brief-", stamp, ".md"))
power_st_variant_summary <- summarize_best_st_variant(power_condensed_summary, pretty_methods)
surface_st_variant_summary <- summarize_best_st_variant(surface_summary, pretty_methods)

# Power-up PDF
power_pdf <- file.path(out_dir, paste0("powerup-18x40-plots-", stamp, ".pdf"))
pdf(power_pdf, width = 13, height = 8)

pp <- power_condensed_summary[power_condensed_summary$method %in% power_display_methods, , drop = FALSE]
pp <- add_orientation_fields(pp)
pp <- add_orientation_labels(pp)
pp$method_pretty <- lookup_pretty(pp$method, pretty_methods)

power_recommend$recommended_pretty <- lookup_pretty(power_recommend$recommended_method, pretty_methods)
cell_levels <- unique(power_recommend$cell_label)
power_recommend$cell_label <- factor(power_recommend$cell_label, levels = rev(cell_levels))
power_recommend$hrf_label <- factor(power_recommend$hrf_label, levels = unique(power_recommend$hrf_label))

p1 <- ggplot(power_recommend, aes(x = score_margin, y = cell_label, color = recommended_pretty)) +
  geom_segment(aes(x = 0, xend = score_margin, y = cell_label, yend = cell_label), alpha = 0.45, linewidth = 0.8) +
  geom_point(size = 3.2) +
  facet_wrap(~ hrf_label, ncol = 1, scales = "free_y") +
  labs(
    title = "How decisive is the recommendation?",
    subtitle = run_profile_label,
    x = "Balanced-score margin over the runner-up",
    y = NULL,
    color = "Recommended method"
  ) +
  bmk_plot_theme(base_size = 12)
print(p1)

power_counts <- aggregate(
  list(n = power_recommend$scenario),
  by = list(hrf_label = power_recommend$hrf_label, recommended_pretty = power_recommend$recommended_pretty),
  FUN = length
)
p2 <- ggplot(power_counts, aes(x = n, y = recommended_pretty, fill = hrf_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  labs(
    title = "Which method gets picked most often?",
    subtitle = "Counts over the 18 power-up scenarios",
    x = "Number of scenarios",
    y = NULL,
    fill = NULL
  ) +
  bmk_plot_theme(base_size = 12)
print(p2)

if (nrow(power_st_variant_summary) > 0L) {
  p2b <- ggplot(power_st_variant_summary, aes(x = n_scenarios, y = reorder(method_pretty, n_scenarios), fill = method_pretty)) +
    geom_col(width = 0.65, show.legend = FALSE) +
    labs(
      title = "Which internal stglmnet variant wins most often?",
      subtitle = "Best internal variant by scenario, using mean correlation within the stglmnet family.",
      x = "Number of scenarios won",
      y = NULL
    ) +
    bmk_plot_theme(base_size = 12)
  print(p2b)
}

pw <- power_pairwise
if (nrow(pw) > 0L) {
  pw <- add_orientation_labels(pw)
  pw$method_b_pretty <- lookup_pretty(pw$method_b, pretty_methods)
  pw$cell_label <- factor(pw$cell_label, levels = rev(cell_levels))

  p3 <- bmk_plot_pairwise_delta(
    data = pw,
    x = "delta_cor_mean",
    y = "cell_label",
    color = "method_b_pretty",
    facet_rows = "hrf_label",
    xmin = "delta_cor_ci_lo",
    xmax = "delta_cor_ci_hi",
    title = "stglmnet versus the other contenders",
    subtitle = "Positive values favor stglmnet on mean correlation",
    x_lab = "Delta correlation",
    y_lab = NULL,
    base_size = 12
  )
  print(p3)
}

pm <- aggregate(
  cbind(mean_cor, mean_sign_acc, mean_elapsed_sec, success_rate) ~ hrf_label + method_pretty,
  data = pp,
  FUN = mean,
  na.rm = TRUE
)

p4 <- ggplot(pm, aes(x = mean_cor, y = method_pretty, color = method_pretty)) +
  geom_segment(aes(x = 0, xend = mean_cor, y = method_pretty, yend = method_pretty), alpha = 0.35, linewidth = 0.8) +
  geom_point(size = 3.2) +
  facet_wrap(~ hrf_label, ncol = 1) +
  labs(
    title = "Average recovery by method",
    subtitle = "Higher is better. This is the simplest read on which estimators recover trialwise betas most faithfully.",
    x = "Mean correlation with truth",
    y = NULL,
    color = "Method"
  ) +
  bmk_plot_theme(base_size = 12)
print(p4)

p5 <- ggplot(pm, aes(x = mean_elapsed_sec, y = mean_cor, color = method_pretty, label = method_pretty)) +
  geom_point(size = 3.2) +
  geom_text(check_overlap = TRUE, nudge_x = max(pm$mean_elapsed_sec, na.rm = TRUE) * 0.02, show.legend = FALSE) +
  facet_wrap(~ hrf_label) +
  labs(
    title = "Recovery versus compute time",
    subtitle = "Upper-left is better: more recovery for less runtime.",
    x = "Mean elapsed time (sec)",
    y = "Mean correlation",
    color = "Method"
  ) +
  bmk_plot_theme(base_size = 12)
print(p5)

p6 <- ggplot(pm, aes(x = method_pretty, y = mean_sign_acc, fill = hrf_label)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  coord_flip() +
  labs(
    title = "Sign accuracy check",
    subtitle = "Higher is better. This catches whether methods preserve effect direction reliably.",
    x = NULL,
    y = "Mean sign accuracy",
    fill = NULL
  ) +
  bmk_plot_theme(base_size = 12)
print(p6)

dev.off()

# HRF surface PDF
surface_pdf <- file.path(out_dir, paste0("hrf-surface-plots-", stamp, ".pdf"))
pdf(surface_pdf, width = 13, height = 8)

ss <- surface_summary[surface_summary$method %in% surface_display_methods, , drop = FALSE]
ss$method_pretty <- lookup_pretty(ss$method, pretty_methods)
ss$regime_label <- factor(ss$regime_label, levels = c("Easy setting", "Moderate setting", "Hard setting"))
ss$widen_label <- factor(ss$widen_label, levels = paste0("HRF widen = ", widens, " s"))

p7 <- ggplot(ss, aes(x = lag, y = mean_cor, color = method_pretty, group = method_pretty)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_grid(regime_label ~ widen_label) +
  labs(
    title = "Recovery as HRF lag mismatch grows",
    subtitle = run_profile_label,
    x = "True HRF lag shift (sec)",
    y = "Mean correlation",
    color = "Method"
  ) +
  bmk_plot_theme(base_size = 11)
print(p7)

p8 <- ggplot(ss, aes(x = lag, y = mean_sign_acc, color = method_pretty, group = method_pretty)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_grid(regime_label ~ widen_label) +
  labs(
    title = "Sign accuracy as HRF lag mismatch grows",
    subtitle = "Higher is better. This plot is often easier to read than a surface tile map.",
    x = "True HRF lag shift (sec)",
    y = "Mean sign accuracy",
    color = "Method"
  ) +
  bmk_plot_theme(base_size = 11)
print(p8)

p9 <- ggplot(ss, aes(x = lag, y = mean_rmse, color = method_pretty, group = method_pretty)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.4) +
  facet_grid(regime_label ~ widen_label) +
  labs(
    title = "RMSE as HRF lag mismatch grows",
    subtitle = "Lower is better.",
    x = "True HRF lag shift (sec)",
    y = "Mean RMSE",
    color = "Method"
  ) +
  bmk_plot_theme(base_size = 11)
print(p9)

surface_recommend$recommended_pretty <- lookup_pretty(surface_recommend$recommended_method, pretty_methods)
surface_recommend$regime_label <- factor(surface_recommend$regime_label, levels = c("Easy setting", "Moderate setting", "Hard setting"))
surface_recommend$mismatch_band <- factor(surface_recommend$mismatch_band, levels = c("Well-specified", "Moderate mismatch", "Strong mismatch"))
surface_counts <- aggregate(
  list(n = surface_recommend$scenario),
  by = list(
    regime_label = surface_recommend$regime_label,
    mismatch_band = surface_recommend$mismatch_band,
    recommended_pretty = surface_recommend$recommended_pretty
  ),
  FUN = length
)

p10 <- ggplot(surface_counts, aes(x = mismatch_band, y = n, fill = recommended_pretty)) +
  geom_col(position = "fill") +
  facet_wrap(~ regime_label) +
  labs(
    title = "How recommendations shift as HRF mismatch increases",
    subtitle = "Each bar shows the share of cells assigned to each method.",
    x = NULL,
    y = "Share of cells",
    fill = "Recommended method"
  ) +
  bmk_plot_theme(base_size = 11)
print(p10)

dev.off()

recommendation_summary_pretty <- power_suite$recommendation_summary
recommendation_summary_pretty$method <- lookup_pretty(recommendation_summary_pretty$method, pretty_methods)
write_benchmark_brief(
  path = recommendation_brief,
  run_profile_label = run_profile_label,
  power_recommend = power_recommend,
  surface_recommend = surface_recommend,
  recommendation_summary = recommendation_summary_pretty,
  st_variant_summary = power_st_variant_summary
)

cat("\nCompleted. Files written:\n")
cat(" -", normalizePath(power_pdf), "\n")
cat(" -", normalizePath(surface_pdf), "\n")
cat(" -", normalizePath(file.path(out_dir, paste0("powerup-18x40-pairwise-ci-", stamp, ".csv"))), "\n")
cat(" -", normalizePath(file.path(out_dir, paste0("powerup-18x40-recommendations-", stamp, ".csv"))), "\n")
cat(" -", normalizePath(file.path(out_dir, paste0("hrf-surface-summary-", stamp, ".csv"))), "\n")
cat(" -", normalizePath(file.path(out_dir, paste0("hrf-surface-recommendations-", stamp, ".csv"))), "\n")
cat(" -", normalizePath(recommendation_brief), "\n")
