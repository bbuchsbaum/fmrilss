if (!exists("benchmark_lss_methods", mode = "function")) {
  source(testthat::test_path("..", "..", "bench", "benchmark_harness.R"))
}

test_that("benchmark_lss_methods returns benchmark object", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3201)
  bmk <- benchmark_lss_methods(
    n_reps = 2,
    seed = 3201,
    n_tp = 70,
    n_trials = 10,
    n_vox = 3,
    isi = 1.4,
    cv_folds = 2
  )

  expect_s3_class(bmk, "lss_benchmark")
  expect_true(all(c("cor_st", "cor_lss", "rmse_st", "rmse_lss") %in% names(bmk$metrics)))
  expect_true("st_better_corr_rate" %in% names(bmk$summary))
  expect_true("metrics_long" %in% names(bmk))
  expect_true("method_summary" %in% names(bmk))
})

test_that("benchmark_lss_methods can include oasis variants", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3202)
  bmk <- benchmark_lss_methods(
    n_reps = 1,
    seed = 3202,
    n_tp = c(60, 60),
    n_runs = 2,
    n_trials = 12,
    n_vox = 3,
    cv_folds = 2,
    include_oasis = TRUE,
    oasis_grid = data.frame(
      ridge_mode = "absolute",
      ridge_x = 0,
      ridge_b = 0,
      stringsAsFactors = FALSE
    )
  )

  expect_s3_class(bmk, "lss_benchmark")
  expect_true(any(grepl("^fmrilss_oasis", bmk$metrics_long$method)))
  expect_true("fmrilss_lss" %in% bmk$metrics_long$method)
})

test_that("benchmark_lss_methods_suite aggregates scenarios", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sc <- default_benchmark_scenarios("quick")[1:2]
  suite <- benchmark_lss_methods_suite(
    scenarios = sc,
    base_args = list(
      n_reps = 1,
      seed = 3301,
      n_tp = 60,
      n_trials = 10,
      n_vox = 3,
      cv_folds = 2
    ),
    trace = FALSE
  )

  expect_s3_class(suite, "lss_benchmark_suite")
  expect_true(all(c("scenario", "method", "mean_cor") %in% names(suite$scenario_method_summary)))
  expect_equal(length(unique(suite$metrics_long$scenario)), 2)
})

test_that("suite carries flattened scenario metadata and recommendations", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sc <- list(
    baseline = list(name = "baseline", noise_model = "iid", true_hrf_lag_sec = 0),
    mismatch = list(name = "mismatch", noise_model = "ar1", ar_rho = 0.4, true_hrf_lag_sec = 2)
  )

  suite <- benchmark_lss_methods_suite(
    scenarios = sc,
    base_args = list(
      n_reps = 1,
      seed = 3310,
      n_tp = 60,
      n_trials = 8,
      n_vox = 2,
      cv_folds = 2
    ),
    trace = FALSE
  )

  expect_true(all(c("noise_model", "true_hrf_lag_sec") %in% names(suite$metrics_long)))
  expect_true(all(c("noise_model", "true_hrf_lag_sec") %in% names(suite$scenario_method_summary)))
  expect_true(all(c("recommended_method", "profile", "reason") %in% names(suite$recommendations)))
  expect_true("scenario_design" %in% names(suite))
  expect_true(nrow(suite$recommendation_summary) > 0)
})

test_that("benchmark_lss_methods can include sbhm method", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3302)
  bmk <- benchmark_lss_methods(
    n_reps = 1,
    seed = 3302,
    n_tp = 60,
    n_trials = 8,
    n_vox = 2,
    cv_folds = 2,
    include_sbhm = TRUE,
    sbhm_rank = 2
  )

  expect_s3_class(bmk, "lss_benchmark")
  expect_true("fmrilss_sbhm" %in% bmk$metrics_long$method)
})

test_that("benchmark_lss_methods can include hrfals fastlss methods", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("hrfals")

  set.seed(3304)
  bmk <- benchmark_lss_methods(
    n_reps = 1,
    seed = 3304,
    n_tp = 60,
    n_trials = 8,
    n_vox = 2,
    cv_folds = 2,
    include_hrfals = TRUE,
    hrfals_modes = "r",
    hrfals_fit_methods = "ls_svd_only",
    hrfals_basis = "spmg2"
  )

  expect_s3_class(bmk, "lss_benchmark")
  expect_true("hrfals_fastlss[r]" %in% bmk$metrics_long$method)
  expect_true("hrfals_fit[ls_svd_only]" %in% bmk$metrics_long$method)
})

test_that("benchmark_lss_methods supports HRF heterogeneity controls", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3303)
  bmk <- benchmark_lss_methods(
    n_reps = 1,
    seed = 3303,
    n_tp = 60,
    n_trials = 8,
    n_vox = 3,
    cv_folds = 2,
    true_hrf_lag_vox_sd = 0.5,
    true_hrf_widen_vox_sd = 0.3,
    true_hrf_lag_trial_sd = 0.2,
    true_hrf_widen_trial_sd = 0.2
  )

  expect_s3_class(bmk, "lss_benchmark")
  expect_true("true_hrf_diag_long" %in% names(bmk))
  expect_equal(nrow(bmk$true_hrf_diag_long), 1)
  expect_true(all(c("mean_abs_lag", "mean_widen") %in% names(bmk$true_hrf_diag_long)))
})

test_that("benchmark_lss_methods supports stress noise and amplitude knobs", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3305)
  bmk <- benchmark_lss_methods(
    n_reps = 1,
    seed = 3305,
    n_tp = 60,
    n_trials = 8,
    n_vox = 2,
    cv_folds = 2,
    noise_model = "arma11",
    ar_rho = 0.55,
    arma_theta = -0.35,
    noise_spike_prob = 0.01,
    noise_spike_scale = 6,
    trial_amplitude_model = "lognormal",
    trial_amplitude_cv = 1.2,
    trial_amplitude_outlier_prob = 0.08,
    trial_amplitude_outlier_scale = 5,
    true_hrf_lag_sec = 3,
    true_hrf_widen_sec = 2
  )

  expect_s3_class(bmk, "lss_benchmark")
  expect_equal(bmk$params$noise_model, "arma11")
  expect_equal(bmk$params$noise_spike_prob, 0.01)
  expect_equal(bmk$params$trial_amplitude_outlier_prob, 0.08)
  expect_true(nrow(bmk$method_summary) > 0)
})

test_that("benchmark results include elapsed_sec timing", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  bmk <- benchmark_lss_methods(
    n_reps = 1, seed = 4001, n_tp = 60, n_trials = 8, n_vox = 2, cv_folds = 2
  )
  expect_true("elapsed_sec" %in% names(bmk$metrics_long))
  expect_true(all(is.finite(bmk$metrics_long$elapsed_sec)))
  expect_true(all(bmk$metrics_long$elapsed_sec >= 0))
})

test_that("benchmark computes classification accuracy when n_conditions > 0", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  bmk <- benchmark_lss_methods(
    n_reps = 1, seed = 4002, n_tp = 60, n_trials = 12, n_vox = 4,
    cv_folds = 2, n_conditions = 3, condition_effect = 1.5
  )
  expect_true("classif_acc" %in% names(bmk$metrics_long))
  ok_rows <- bmk$metrics_long[bmk$metrics_long$ok == TRUE, ]
  expect_true(any(is.finite(ok_rows$classif_acc)))
  expect_true(all(ok_rows$classif_acc >= 0 & ok_rows$classif_acc <= 1, na.rm = TRUE))
})

test_that("fast_mode reduces benchmark parameters", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  bmk <- benchmark_lss_methods(
    seed = 4003, fast_mode = TRUE, cv_folds = 2
  )
  expect_equal(bmk$params$n_vox, 3L)
  expect_s3_class(bmk, "lss_benchmark")
})

test_that("score function returns per-voxel correlations", {
  # Test the internal scoring function directly
  set.seed(4004)
  beta_hat <- matrix(rnorm(30), 10, 3)
  beta_true <- beta_hat + matrix(rnorm(30, sd = 0.3), 10, 3)

  score <- .bmk_score_betas(beta_hat, beta_true)
  expect_true("per_voxel_cor" %in% names(score))
  expect_equal(length(score$per_voxel_cor), 3)
  expect_true(all(is.finite(score$per_voxel_cor)))
})

test_that("spatial_rho introduces cross-voxel correlation", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  bmk <- benchmark_lss_methods(
    n_reps = 1, seed = 4005, n_tp = 60, n_trials = 8, n_vox = 3,
    cv_folds = 2, spatial_rho = 0.5
  )
  expect_s3_class(bmk, "lss_benchmark")
  expect_equal(bmk$params$spatial_rho, 0.5)
})

test_that("suite pairwise summary includes Wilcoxon p-values", {
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sc <- default_benchmark_scenarios("quick")[1]
  suite <- benchmark_lss_methods_suite(
    scenarios = sc,
    base_args = list(n_reps = 3, seed = 4006, n_tp = 60, n_trials = 8, n_vox = 2, cv_folds = 2),
    trace = FALSE
  )
  expect_true("p_wilcox_cor" %in% names(suite$pairwise_summary))
  expect_true("p_wilcox_rmse" %in% names(suite$pairwise_summary))
})

test_that("benchmark plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")

  surf <- data.frame(
    method = c("stglmnet", "stglmnet"),
    regime = c("easy", "easy"),
    lag = c(0, 1),
    widen = c(0, 0),
    mean_cor = c(0.8, 0.7),
    mean_rmse = c(0.3, 0.4),
    stringsAsFactors = FALSE
  )
  surf$lag_f <- factor(surf$lag)
  surf$widen_f <- factor(surf$widen)
  curves <- data.frame(
    method = c("stglmnet", "lss"),
    hetero = c("low", "low"),
    mean_cor = c(0.8, 0.75),
    ci_lo = c(0.75, 0.7),
    ci_hi = c(0.85, 0.8),
    stringsAsFactors = FALSE
  )
  deltas <- data.frame(
    scenario = "baseline",
    comparison = "stglmnet - lss",
    delta_cor_mean = 0.03,
    delta_cor_ci_lo = 0.01,
    delta_cor_ci_hi = 0.05,
    stringsAsFactors = FALSE
  )
  recs <- data.frame(
    hrf_regime = "nominal",
    cell = "cell_1",
    recommended_method = "stglmnet",
    stringsAsFactors = FALSE
  )

  expect_s3_class(
    bmk_plot_surface_heatmap(surf, x = "lag_f", y = "widen_f", metric = "mean_cor", facet_rows = "method", facet_cols = "regime"),
    "ggplot"
  )
  expect_s3_class(
    bmk_plot_method_curves(curves, x = "hetero", y = "mean_cor", color = "method", group = "method", ymin = "ci_lo", ymax = "ci_hi"),
    "ggplot"
  )
  expect_s3_class(
    bmk_plot_pairwise_delta(deltas, x = "delta_cor_mean", y = "scenario", color = "comparison"),
    "ggplot"
  )
  expect_s3_class(
    bmk_plot_recommendation_tiles(recs, x = "hrf_regime", y = "cell"),
    "ggplot"
  )
})

test_that("wide scenario preset includes SNR and TR scenarios", {
  sc <- default_benchmark_scenarios("wide")
  names_sc <- vapply(sc, function(x) .bmk_or(x$name, ""), character(1))
  expect_true("snr_high" %in% names_sc)
  expect_true("snr_low" %in% names_sc)
  expect_true("tr_hcp" %in% names_sc)
  expect_true("ar1_prewhitened" %in% names_sc)
})
