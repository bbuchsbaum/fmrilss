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

test_that("lss_with_hrf recovers trial betas", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")
  skip_if_not_installed("progress")

  set.seed(1)
  n_time <- 60
  events <- data.frame(
    onset = c(10, 30, 50),
    duration = c(1, 1, 1),
    condition = "A"
  )
  basis <- fmrihrf::HRF_FIR(length = 1)
  X_basis <- fmrihrf::regressor_set(events, basis = basis, n = n_time)$X
  n_vox <- 2
  hcoef <- runif(n_vox, 0.5, 1.5)
  true_beta <- matrix(rnorm(ncol(X_basis) * n_vox), ncol(X_basis), n_vox)
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    Y[, v] <- X_basis %*% true_beta[, v] * hcoef[v]
  }
  hrf_est <- list(coefficients = matrix(hcoef, 1, n_vox),
                  basis = basis,
                  conditions = "A")
  class(hrf_est) <- "VoxelHRF"

  res <- lss_with_hrf(Y, events, hrf_est, verbose = FALSE, chunk_size = n_vox)
  est <- as.matrix(res$betas[])
  expect_equal(est, true_beta, tolerance = 1e-6)
})

test_that("lss_with_hrf equivalent to lss when HRF identical", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")
  skip_if_not_installed("progress")

  set.seed(2)
  n_time <- 50
  events <- data.frame(onset = c(5, 25, 40), duration = 1, condition = "A")
  basis <- fmrihrf::HRF_FIR(length = 1)
  X_basis <- fmrihrf::regressor_set(events, basis = basis, n = n_time)$X
  n_vox <- 3
  coef_shared <- rep(1, n_vox)
  betas <- matrix(rnorm(ncol(X_basis) * n_vox), ncol(X_basis), n_vox)
  Y <- matrix(0, n_time, n_vox)
  for (v in 1:n_vox) {
    Y[, v] <- X_basis %*% betas[, v]
  }
  hrf_est <- list(coefficients = matrix(coef_shared, 1, n_vox),
                  basis = basis,
                  conditions = "A")
  class(hrf_est) <- "VoxelHRF"

  lss_res <- lss(Y, X_basis)
  res <- lss_with_hrf(Y, events, hrf_est, verbose = FALSE, chunk_size = n_vox)
  est <- as.matrix(res$betas[])
  expect_equal(est, lss_res, tolerance = 1e-6)
})

test_that("lss_with_hrf basis aliasing", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")
  skip_if_not_installed("progress")

  set.seed(3)
  n_time <- 80
  events <- data.frame(onset = c(10, 40), duration = 1, condition = "A")
  basis_true <- fmrihrf::HRF_BSPLINE(length = 8)
  basis_fit <- fmrihrf::HRF_FIR(length = 8)

  X_basis_true <- fmrihrf::regressor_set(events, basis = basis_true, n = n_time)$X
  hcoef <- rnorm(ncol(X_basis_true))
  hrf_kernel <- fmrihrf::hrf_from_coefficients(basis_true, hcoef)
  C <- fmrihrf::convolve_design(events$onset, events$duration, hrf_kernel, n_time)
  beta <- rnorm(ncol(C))
  Y <- C %*% beta
  hrf_est <- list(coefficients = matrix(hcoef, ncol(X_basis_true), 1),
                  basis = basis_fit,
                  conditions = "A")
  class(hrf_est) <- "VoxelHRF"
  res <- lss_with_hrf(matrix(Y, ncol = 1), events, hrf_est, verbose = FALSE, chunk_size = 1)
  est_kernel <- fmrihrf::hrf_from_coefficients(basis_fit, hrf_est$coefficients[,1])
  expect_equal(est_kernel, hrf_kernel, tolerance = 1e-6)
})
