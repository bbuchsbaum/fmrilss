test_that("design_spec warns for likely run-relative onsets and injects run intercepts", {
  skip_if_not_installed("fmrihrf")

  set.seed(123)
  sframe <- fmrihrf::sampling_frame(blocklens = c(50, 50), TR = 2)
  # Onsets that look run-relative (<= one run duration)
  spec <- list(
    sframe = sframe,
    cond = list(onsets = c(10, 20, 30), duration = 0, span = 20, hrf = fmrihrf::HRF_SPMG1)
  )
  Y <- matrix(rnorm(100 * 6), 100, 6)

  # Expect a warning about run-relative onsets
  expect_warning(
    beta_auto <- lss(Y, X = NULL, method = "oasis", oasis = list(design_spec = spec)),
    "run-relative",
    ignore.case = TRUE
  )

  # Manual run-wise intercepts should match the auto-injected ones
  runs <- rep(1:2, c(50, 50))
  Z_run <- stats::model.matrix(~ 0 + factor(runs))
  beta_manual <- lss(Y, X = NULL, Z = Z_run, method = "oasis", oasis = list(design_spec = spec))

  expect_equal(dim(beta_auto), dim(beta_manual))
  expect_equal(beta_auto, beta_manual, tolerance = 1e-8)
})

