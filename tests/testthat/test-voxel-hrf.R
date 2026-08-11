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
                                 fac = factor(rep("all events", nrow(events))),
                                 hrf = basis,
                                 duration = events$duration,
                                 span = 30)
  X_basis <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  n_vox <- 3
  true_coef <- matrix(rnorm(ncol(X_basis) * n_vox), ncol(X_basis), n_vox)
  Y <- X_basis %*% true_coef
  est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
  recovered_raw <- sweep(est$coefficients, 2, est$amplitude_scale, "*")
  rmse <- sqrt(mean((recovered_raw - true_coef)^2))
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

test_that("voxel-HRF public inputs fail closed with semantic errors", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
  events <- data.frame(
    onset = c(5, 25, 45), duration = 0, condition = "A"
  )
  basis <- fmrihrf::HRF_SPMG1
  design <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)$X
  Y <- design %*% matrix(c(1, 1, 1, 2, 2, 2), 3, 2)
  colnames(Y) <- c("voxel_a", "voxel_b")
  estimate <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)

  duplicated_y <- Y
  colnames(duplicated_y) <- c("voxel_a", "voxel_a")
  expect_error(
    estimate_voxel_hrf(duplicated_y, events, basis, sframe = sframe),
    "voxel names must be complete and unique"
  )
  nonfinite_y <- Y
  nonfinite_y[1, 1] <- NA_real_
  expect_error(
    estimate_voxel_hrf(nonfinite_y, events, basis, sframe = sframe),
    "finite values"
  )

  empty <- events[0, , drop = FALSE]
  missing_onset <- events
  missing_onset$onset[1] <- NA_real_
  late_onset <- events
  late_onset$onset[1] <- 1000
  negative_duration <- events
  negative_duration$duration[1] <- -1
  for (bad in list(empty, missing_onset, late_onset, negative_duration)) {
    expect_error(estimate_voxel_hrf(Y, bad, basis, sframe = sframe))
    expect_error(lss_with_hrf(Y, bad, estimate, sframe = sframe, verbose = FALSE))
  }

  character_coefficients <- estimate
  character_coefficients$coefficients <- matrix("x", 1, 2)
  expect_error(
    lss_with_hrf(Y, events, character_coefficients, verbose = FALSE),
    "numeric coefficient matrix"
  )
  missing_coefficients <- estimate
  missing_coefficients$coefficients[1, 1] <- NA_real_
  expect_error(
    lss_with_hrf(Y, events, missing_coefficients, verbose = FALSE),
    "coefficients must be finite"
  )
})

test_that("voxel-HRF condition pooling and engine metadata are explicit", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")

  sframe <- fmrihrf::sampling_frame(blocklens = 70, TR = 1)
  events <- data.frame(
    onset = c(5, 25, 45), duration = 0, condition = "A"
  )
  basis <- fmrihrf::HRF_SPMG1
  X <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)$X
  Y <- X %*% matrix(rep(c(2, 3, 4), each = nrow(events)), nrow(events), 3)
  colnames(Y) <- c("voxel_a", "voxel_b", "voxel_c")

  one_condition <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
  relabeled <- events
  relabeled$condition <- c("A", "B", "B")
  two_conditions <- estimate_voxel_hrf(Y, relabeled, basis, sframe = sframe)
  expect_identical(one_condition$condition_pooling, "all-events")
  expect_identical(two_conditions$condition_pooling, "all-events")
  expect_equal(one_condition$coefficients, two_conditions$coefficients)

  fit_r <- lss_with_hrf(Y, events, one_condition, engine = "R", verbose = FALSE)
  fit_cpp <- lss_with_hrf(
    Y, events, one_condition, engine = "C++", chunk_size = 1,
    verbose = FALSE
  )
  dense_cpp <- as.matrix(fit_cpp)
  expect_lte(max(abs(fit_r - dense_cpp)), 1e-10)
  expect_identical(attr(fit_r, "engine_requested"), "R")
  expect_identical(attr(fit_r, "engine_used"), "r")
  expect_identical(fit_cpp$engine_requested, "C++")
  expect_true(fit_cpp$engine_used %in% c("cpp_arma", "cpp", "r"))
  expect_identical(fit_cpp$chunk_size, 1L)
  expect_identical(attr(dense_cpp, "engine_used"), fit_cpp$engine_used)
  expect_identical(attr(dense_cpp, "chunk_size"), 1L)
  expect_identical(attr(dense_cpp, "event_amplitude"), rep(1, nrow(events)))
  expect_identical(attr(dense_cpp, "event_duration"), rep(0, nrow(events)))
})

test_that("non-unit event amplitudes have explicit per-unit coefficient semantics", {
  skip_if_not_installed("fmrihrf")
  skip_if_not_installed("bigmemory")

  sframe <- fmrihrf::sampling_frame(blocklens = 80, TR = 1)
  events <- data.frame(
    onset = c(5, 25, 50), duration = 0, condition = "A",
    amplitude = c(2, 3, 4)
  )
  basis <- fmrihrf::HRF_SPMG1
  built <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)
  Y <- built$X %*% rep(3, nrow(events))
  colnames(Y) <- "voxel_a"
  estimate <- structure(
    list(
      coefficients = matrix(1, 1, 1, dimnames = list(NULL, "voxel_a")),
      basis = basis, sframe = sframe, normalization = "positive-peak"
    ),
    class = "VoxelHRF"
  )

  fit <- lss_with_hrf(Y, events, estimate, engine = "R", verbose = FALSE)
  fit_cpp <- as.matrix(lss_with_hrf(
    Y, events, estimate, engine = "C++", chunk_size = 1, verbose = FALSE
  ))
  doubled_events <- events
  doubled_events$amplitude <- 2 * doubled_events$amplitude
  doubled_fit <- lss_with_hrf(
    Y, doubled_events, estimate, engine = "R", verbose = FALSE
  )

  expect_equal(unname(fit), matrix(3, nrow(events), 1), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(fit_cpp, fit, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(unname(doubled_fit), matrix(1.5, nrow(events), 1),
               tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(2 * unname(doubled_fit), unname(fit), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_identical(
    attr(fit, "units"),
    "coefficient on supplied event design with unit-peak HRF shape"
  )
  expect_identical(attr(fit, "event_amplitude"), events$amplitude)
  expect_identical(attr(fit_cpp, "event_amplitude"), events$amplitude)
  expect_identical(attr(fit, "event_duration"), events$duration)
  expect_identical(attr(fit_cpp, "event_duration"), events$duration)

  zero_events <- events
  zero_events$amplitude[1] <- 0
  expect_error(
    lss_with_hrf(Y, zero_events, estimate, engine = "R", verbose = FALSE),
    "must be nonzero"
  )
})

test_that("nonzero durations return event-design coefficients, not peak amplitudes", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1)
  events <- data.frame(
    onset = c(5, 25, 50, 75), duration = c(0, 2, 4, 8), condition = "A"
  )
  basis <- fmrihrf::HRF_SPMG1
  built <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)
  Y <- built$X %*% rep(3, nrow(events))
  colnames(Y) <- "voxel_a"
  estimate <- structure(
    list(
      coefficients = matrix(1, 1, 1, dimnames = list(NULL, "voxel_a")),
      basis = basis, sframe = sframe, normalization = "positive-peak"
    ),
    class = "VoxelHRF"
  )

  fit <- lss_with_hrf(Y, events, estimate, engine = "R", verbose = FALSE)
  modeled_peaks <- 3 * apply(built$X, 2L, max)

  expect_equal(unname(fit), matrix(3, nrow(events), 1), tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_gt(max(abs(modeled_peaks - 3)), 0.1)
  expect_identical(
    attr(fit, "units"),
    "coefficient on supplied event design with unit-peak HRF shape"
  )
  expect_identical(attr(fit, "event_duration"), events$duration)
})

test_that("voxel-HRF estimation fails closed when the common span identifies the HRF", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 120, TR = 0.8)
  events <- data.frame(
    onset = seq(8, 72, length.out = 8), duration = 0, condition = "A"
  )
  basis <- fmrihrf::HRF_SPMG3
  built <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)
  pooled_basis <- do.call(cbind, lapply(built$basis_convolved, rowSums))
  raw_coefficients <- matrix(c(1, 0.15, -0.08, 0.8, -0.1, 0.05), 3, 2)
  Y <- pooled_basis %*% raw_coefficients

  expect_error(
    estimate_voxel_hrf(
      Y, events, basis, sframe = sframe, nuisance_regs = pooled_basis
    ),
    "not identifiable after projecting the common design"
  )
})

test_that("voxel-HRF estimation is invariant to redundant rescaled nuisance bases", {
  skip_if_not_installed("fmrihrf")

  sframe <- fmrihrf::sampling_frame(blocklens = 140, TR = 1)
  events <- data.frame(
    onset = seq(8, 108, length.out = 9), duration = 0, condition = "A"
  )
  basis <- fmrihrf::HRF_SPMG3
  built <- fmrilss:::.voxhrf_trial_basis(events, basis, sframe)
  pooled_basis <- do.call(cbind, lapply(built$basis_convolved, rowSums))
  raw_coefficients <- matrix(c(1, 0.12, -0.05, 0.7, -0.08, 0.04), 3, 2)
  time <- seq_len(nrow(pooled_basis))
  nuisance <- cbind(sin(time / 13), cos(time / 17))
  Y <- pooled_basis %*% raw_coefficients + nuisance %*% matrix(c(2, -1, -3, 0.5), 2, 2)

  reference <- estimate_voxel_hrf(
    Y, events, basis, sframe = sframe, nuisance_regs = nuisance
  )
  equivalent <- cbind(
    nuisance[, 1] * 1e-8,
    nuisance[, 1] * 2e8,
    nuisance[, 2] * 1e8,
    nuisance[, 2] * -3e-8
  )
  mutated <- estimate_voxel_hrf(
    Y, events, basis, sframe = sframe, nuisance_regs = equivalent
  )
  fit_reference <- lss_with_hrf(
    Y, events, reference, sframe = sframe,
    nuisance_regs = nuisance, engine = "R", verbose = FALSE
  )
  fit_mutated <- lss_with_hrf(
    Y, events, reference, sframe = sframe,
    nuisance_regs = equivalent, engine = "R", verbose = FALSE
  )
  recovered_raw <- sweep(
    reference$coefficients, 2L, reference$amplitude_scale, "*"
  )

  expect_equal(mutated$coefficients, reference$coefficients, tolerance = 1e-10)
  expect_equal(mutated$amplitude_scale, reference$amplitude_scale, tolerance = 1e-10)
  expect_equal(recovered_raw, raw_coefficients, tolerance = 1e-10)
  expect_equal(fit_mutated, fit_reference, tolerance = 1e-10)
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

test_that("lss_with_hrf handles multi-basis HRF correctly", {
  skip_if_not_installed("fmrihrf")
  set.seed(42)

  # Setup parameters
  n_time <- 100
  n_trials <- 5
  n_vox <- 10
  TR <- 1.0

  # Create sampling frame
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = TR)
  times <- fmrihrf::samples(sframe, global = TRUE)

  # Define events
  events <- data.frame(
    onset = c(10, 25, 40, 60, 80),
    duration = rep(1, n_trials),
    condition = rep("A", n_trials)
  )

  # Create multi-basis HRF (SPMG3 has 3 basis functions)
  hrf_multi <- fmrihrf::HRF_SPMG3
  K <- 3  # Number of basis functions

  # Build design matrix using fmrihrf
  rset <- fmrihrf::regressor_set(
    onsets = events$onset,
    fac = factor(1:n_trials),
    hrf = hrf_multi,
    duration = events$duration,
    span = 30,
    summate = TRUE
  )

  # Evaluate to get convolved design matrix
  X_multi <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  if (inherits(X_multi, "Matrix")) X_multi <- as.matrix(X_multi)

  # X_multi should have dimensions: n_time x (n_trials * K)
  expect_equal(dim(X_multi), c(n_time, n_trials * K))

  # Generate synthetic data with multi-basis structure
  # Each trial has K coefficients (one per basis)
  true_coefs <- matrix(rnorm(n_trials * K * n_vox, mean = 0, sd = 0.5),
                       nrow = n_trials * K, ncol = n_vox)
  # Add stronger signal to first basis to make it dominant
  true_coefs[seq(1, n_trials * K, by = K), ] <-
    true_coefs[seq(1, n_trials * K, by = K), ] + 1

  Y <- X_multi %*% true_coefs + matrix(rnorm(n_time * n_vox, sd = 0.2), n_time, n_vox)

  # Extract basis kernels from HRF object
  hrf_eval_times <- seq(0, 30, by = TR)
  hrf_basis_kernels <- fmrihrf::evaluate(hrf_multi, hrf_eval_times)
  if (inherits(hrf_basis_kernels, "Matrix")) {
    hrf_basis_kernels <- as.matrix(hrf_basis_kernels)
  }

  # Create coefficients matrix (K x n_vox)
  # For testing, use identity - each voxel uses all basis functions equally
  coefficients <- matrix(1/K, nrow = K, ncol = n_vox)

  # Run lss_with_hrf with multi-basis HRF
  onset_idx <- as.integer(events$onset)
  durations <- as.integer(events$duration)

  result <- fmrilss:::lss_with_hrf_pure_r(
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

  # Check dimensions
  expect_equal(dim(result), c(n_trials, n_vox))

  # Results should be finite
  expect_true(all(is.finite(result)))

  # Check that we can recover some signal
  # Since first basis has stronger signal, trial estimates should correlate with it
  trial_means <- rowMeans(result)
  expect_true(var(trial_means) > 0)  # Should have variation across trials

  # Alternative test using OASIS with multi-basis
  if (is.null(colnames(X_multi))) colnames(X_multi) <- paste0("x", seq_len(ncol(X_multi)))
  oasis_result <- lss(
    Y = Y,
    X = X_multi,
    method = "oasis",
    oasis = list(
      K = K, ntrials = n_trials,
      trial_basis_map = data.frame(
        column = colnames(X_multi),
        trial = rep(seq_len(n_trials), each = K),
        basis = rep(seq_len(K), times = n_trials)
      )
    )
  )

  # OASIS should return K*n_trials rows for multi-basis
  expect_equal(nrow(oasis_result), n_trials * K)
  expect_equal(ncol(oasis_result), n_vox)

  # Compare aggregated results (sum across basis functions per trial)
  oasis_aggregated <- matrix(0, n_trials, n_vox)
  for (i in 1:n_trials) {
    idx <- ((i-1)*K + 1):(i*K)
    oasis_aggregated[i, ] <- colSums(oasis_result[idx, , drop = FALSE])
  }

  # Results should be somewhat correlated
  correlation <- cor(as.vector(result), as.vector(oasis_aggregated))
  expect_gt(correlation, 0.5)  # Should have positive correlation
})
