if (!exists("benchmark_stglmnet_vs_fmrilss", mode = "function")) {
  source(testthat::test_path("..", "..", "bench", "stglmnet_vs_fmrilss_harness.R"))
}

test_that("benchmark_stglmnet_vs_fmrilss returns benchmark object", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3201)
  bmk <- benchmark_stglmnet_vs_fmrilss(
    n_reps = 2,
    seed = 3201,
    n_tp = 70,
    n_trials = 10,
    n_vox = 3,
    isi = 1.4,
    cv_folds = 2
  )

  expect_s3_class(bmk, "stglmnet_benchmark")
  expect_true(all(c("cor_st", "cor_lss", "rmse_st", "rmse_lss") %in% names(bmk$metrics)))
  expect_true("st_better_corr_rate" %in% names(bmk$summary))
  expect_true("metrics_long" %in% names(bmk))
  expect_true("method_summary" %in% names(bmk))
})

test_that("benchmark_stglmnet_vs_fmrilss can include oasis variants", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3202)
  bmk <- benchmark_stglmnet_vs_fmrilss(
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

  expect_s3_class(bmk, "stglmnet_benchmark")
  expect_true(any(grepl("^fmrilss_oasis", bmk$metrics_long$method)))
  expect_true("fmrilss_lss" %in% bmk$metrics_long$method)
})

test_that("benchmark_stglmnet_vs_fmrilss_suite aggregates scenarios", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  sc <- default_stglmnet_vs_fmrilss_scenarios("quick")[1:2]
  suite <- benchmark_stglmnet_vs_fmrilss_suite(
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

  expect_s3_class(suite, "stglmnet_benchmark_suite")
  expect_true(all(c("scenario", "method", "mean_cor") %in% names(suite$scenario_method_summary)))
  expect_equal(length(unique(suite$metrics_long$scenario)), 2)
})

test_that("benchmark_stglmnet_vs_fmrilss can include sbhm method", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3302)
  bmk <- benchmark_stglmnet_vs_fmrilss(
    n_reps = 1,
    seed = 3302,
    n_tp = 60,
    n_trials = 8,
    n_vox = 2,
    cv_folds = 2,
    include_sbhm = TRUE,
    sbhm_rank = 2
  )

  expect_s3_class(bmk, "stglmnet_benchmark")
  expect_true("fmrilss_sbhm" %in% bmk$metrics_long$method)
})

test_that("benchmark_stglmnet_vs_fmrilss can include hrfals fastlss methods", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("hrfals")

  set.seed(3304)
  bmk <- benchmark_stglmnet_vs_fmrilss(
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

  expect_s3_class(bmk, "stglmnet_benchmark")
  expect_true("hrfals_fastlss[r]" %in% bmk$metrics_long$method)
  expect_true("hrfals_fit[ls_svd_only]" %in% bmk$metrics_long$method)
})

test_that("benchmark_stglmnet_vs_fmrilss supports HRF heterogeneity controls", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3303)
  bmk <- benchmark_stglmnet_vs_fmrilss(
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

  expect_s3_class(bmk, "stglmnet_benchmark")
  expect_true("true_hrf_diag_long" %in% names(bmk))
  expect_equal(nrow(bmk$true_hrf_diag_long), 1)
  expect_true(all(c("mean_abs_lag", "mean_widen") %in% names(bmk$true_hrf_diag_long)))
})

test_that("benchmark_stglmnet_vs_fmrilss supports stress noise and amplitude knobs", {
  skip_if_not_installed("stglmnet")
  skip_if_not_installed("fmridesign")
  skip_if_not_installed("fmrihrf")

  set.seed(3305)
  bmk <- benchmark_stglmnet_vs_fmrilss(
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

  expect_s3_class(bmk, "stglmnet_benchmark")
  expect_equal(bmk$params$noise_model, "arma11")
  expect_equal(bmk$params$noise_spike_prob, 0.01)
  expect_equal(bmk$params$trial_amplitude_outlier_prob, 0.08)
  expect_true(nrow(bmk$method_summary) > 0)
})
