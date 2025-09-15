skip_on_cran()

test_that("fit_oasis_grid handles multi-basis HRFs without error", {
  skip_if_not_installed("fmrihrf")
  set.seed(7)
  n_time  <- 160
  V       <- 6
  onsets  <- sort(sample(10:(n_time-30), 12))
  sframe  <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)

  # Build synthetic Y with a multi-basis HRF aggregate
  times   <- fmrihrf::samples(sframe, global = TRUE)
  reg_all <- fmrihrf::regressor(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, duration = 0, span = 30)
  Xagg    <- fmrihrf::evaluate(reg_all, times)
  if (is.matrix(Xagg)) Xagg <- rowSums(Xagg)
  Xagg    <- matrix(Xagg, ncol = 1)
  B       <- matrix(rnorm(1*V, sd=0.6), 1, V)
  Y       <- Xagg %*% B + matrix(rnorm(n_time*V, sd=0.7), n_time, V)

  grid <- create_lwu_grid(n_tau = 3, n_sigma = 2, n_rho = 2)
  out  <- fit_oasis_grid(Y, onsets, sframe, grid,
                         ridge_x = 0.02, ridge_b = 0.02)

  expect_true(is.list(out))
  expect_true(all(c("best_idx","best_params","best_hrf","beta","scores","grid") %in% names(out)))
  expect_equal(nrow(out$beta), length(onsets))
  expect_equal(ncol(out$beta), V)
  expect_true(is.finite(out$scores[out$best_idx]))
})

