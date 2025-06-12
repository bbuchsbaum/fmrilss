library(testthat)
library(fmrilss)


test_that("estimate_voxel_hrf recovers known coefficients", {
  skip_if_not_installed("fmrihrf")

  set.seed(42)
  n_time <- 100
  events <- data.frame(
    onset = c(10, 40, 70),
    duration = c(1, 1, 1),
    condition = c("A", "A", "A")
  )
  basis <- fmrihrf::HRF_FIR(length = 10)
  X_basis <- fmrihrf::regressor_set(events, basis = basis, n = n_time)$X
  n_vox <- 3
  true_coef <- matrix(rnorm(ncol(X_basis) * n_vox), ncol(X_basis), n_vox)
  Y <- X_basis %*% true_coef
  est <- estimate_voxel_hrf(Y, events, basis)
  rmse <- sqrt(mean((est$coefficients - true_coef)^2))
  expect_lt(rmse, 1e-6)
})


test_that("estimate_voxel_hrf input validation", {
  skip_if_not_installed("fmrihrf")

  Y <- matrix(rnorm(20), 10, 2)
  events <- data.frame(onset = 1, duration = 1, condition = "A")
  basis <- fmrihrf::HRF_FIR(length = 5)

  expect_error(estimate_voxel_hrf("no", events, basis),
               "Y must be a numeric matrix")
  bad_events <- data.frame(time = 1)
  expect_error(estimate_voxel_hrf(Y, bad_events, basis),
               "events must be a data.frame")
  expect_error(estimate_voxel_hrf(Y, events, list()),
               "basis must be")
  bad_nuis <- matrix(1, 5, 1)
  expect_error(estimate_voxel_hrf(Y, events, basis, nuisance_regs = bad_nuis),
               "nuisance_regs")
})
