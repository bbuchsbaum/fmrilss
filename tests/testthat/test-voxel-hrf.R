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
  basis <- fmrihrf::hrf_fir_generator(nbasis = 10)
  # Build regressor set with proper API
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  rset <- fmrihrf::regressor_set(onsets = events$onset, 
                                 fac = factor(1:nrow(events)),
                                 hrf = basis,
                                 duration = events$duration,
                                 span = 30)
  X_basis <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
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
  basis <- fmrihrf::hrf_fir_generator(nbasis = 5)

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
  
  set.seed(1)
  n_time <- 60
  n_trials <- 3
  n_vox <- 2
  
  # Simple setup: create event onsets
  events <- data.frame(
    onset = c(10, 30, 50),
    duration = c(1, 1, 1),
    condition = "A"
  )
  
  # Create a simple HRF kernel directly (bypass fmrihrf complexity)
  # This is a simple gamma-like shape
  hrf_kernel <- c(0, 0.1, 0.3, 0.6, 0.9, 1.0, 0.8, 0.5, 0.3, 0.1, 0.05, 0)
  L <- length(hrf_kernel)
  
  # Single basis (K=1) 
  hrf_basis_kernels <- matrix(hrf_kernel, nrow = L, ncol = 1)
  
  # Voxel-specific HRF weights (all 1 for simplicity)
  coefficients <- matrix(1, nrow = 1, ncol = n_vox)
  
  # Build event matrix manually
  onset_idx <- as.integer(events$onset)
  durations <- as.integer(events$duration)
  
  # Create true betas
  true_beta <- matrix(rnorm(n_trials * n_vox, mean = 1, sd = 0.5), n_trials, n_vox)
  
  # Generate Y by convolving events with HRF and scaling by betas
  # This mimics what would happen in real fMRI data
  Y <- matrix(0, n_time, n_vox)
  for (trial in 1:n_trials) {
    # Create impulse for this trial
    impulse <- rep(0, n_time)
    impulse[onset_idx[trial]:(onset_idx[trial] + durations[trial])] <- 1
    
    # Convolve with HRF
    conv_signal <- stats::convolve(impulse, rev(hrf_kernel), type = "open")[1:n_time]
    
    # Add to Y scaled by trial beta
    for (v in 1:n_vox) {
      Y[, v] <- Y[, v] + conv_signal * true_beta[trial, v]
    }
  }
  
  # Add small noise
  Y <- Y + matrix(rnorm(n_time * n_vox, sd = 0.01), n_time, n_vox)
  
  # Create HRF estimates object
  hrf_est <- list(
    coefficients = coefficients,
    basis = structure(list(), class = "HRF"),  # Dummy HRF object
    conditions = "A"
  )
  class(hrf_est) <- "VoxelHRF"
  
  # Use pure R implementation directly
  res <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = durations,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    Z = NULL,
    Nuisance = NULL,
    verbose = FALSE,
    method = "r"
  )
  
  # Check that we recover the betas reasonably well
  # With some tolerance for numerical differences
  expect_equal(dim(res), c(n_trials, n_vox))
  
  # The recovered betas should be close to true betas
  # Allow higher tolerance due to LSS estimation
  correlation <- cor(as.vector(res), as.vector(true_beta))
  expect_gt(correlation, 0.95)  # High correlation expected
})

test_that("lss_with_hrf equivalent to lss when HRF identical", {
  set.seed(2)
  n_time <- 50
  n_trials <- 3
  n_vox <- 3
  
  # Simple event setup
  events <- data.frame(
    onset = c(5, 25, 40),
    duration = c(1, 1, 1),
    condition = "A"
  )
  
  # Create a simple HRF kernel
  hrf_kernel <- c(0, 0.2, 0.5, 0.8, 1.0, 0.7, 0.4, 0.2, 0.1, 0)
  L <- length(hrf_kernel)
  
  # Single basis, same coefficients for all voxels (identity HRF)
  hrf_basis_kernels <- matrix(hrf_kernel, nrow = L, ncol = 1)
  coefficients <- matrix(1, nrow = 1, ncol = n_vox)  # All voxels have same HRF
  
  # Build design matrix manually (what standard LSS would use)
  X <- matrix(0, n_time, n_trials)
  onset_idx <- as.integer(events$onset)
  durations <- as.integer(events$duration)
  
  for (trial in 1:n_trials) {
    # Create impulse for this trial
    impulse <- rep(0, n_time)
    if (onset_idx[trial] <= n_time) {
      end_idx <- min(n_time, onset_idx[trial] + durations[trial])
      impulse[onset_idx[trial]:end_idx] <- 1
    }
    # Convolve with HRF
    X[, trial] <- stats::convolve(impulse, rev(hrf_kernel), type = "open")[1:n_time]
  }
  
  # Generate data
  true_betas <- matrix(rnorm(n_trials * n_vox, mean = 1, sd = 0.5), n_trials, n_vox)
  Y <- X %*% true_betas + matrix(rnorm(n_time * n_vox, sd = 0.1), n_time, n_vox)
  
  # Run standard LSS
  lss_res <- lss(Y, X, method = "r_optimized")
  
  # Run lss_with_hrf with identical HRFs for all voxels
  hrf_res <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = durations,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    Z = NULL,
    Nuisance = NULL,
    verbose = FALSE,
    method = "r"
  )
  
  # They should be very similar (not identical due to numerical differences)
  expect_equal(dim(hrf_res), dim(lss_res))
  
  # Check correlation is very high
  correlation <- cor(as.vector(hrf_res), as.vector(lss_res))
  expect_gt(correlation, 0.99)  # Should be nearly identical
  
  # Check mean squared difference is small
  mse <- mean((hrf_res - lss_res)^2)
  expect_lt(mse, 0.01)
})

test_that("lss_with_hrf basis aliasing", {
  skip("Test needs redesign - convolve_design not exported from fmrihrf")
})
