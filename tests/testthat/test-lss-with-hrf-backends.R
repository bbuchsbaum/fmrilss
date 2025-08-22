library(testthat)
library(fmrilss)

test_that("all backends produce equivalent results", {
  skip_on_cran()
  
  set.seed(2025)
  n_time <- 60L
  n_trials <- 5L
  n_vox <- 10L
  
  # Generate test data
  onset_idx <- as.integer(seq(5, n_time - 10, length.out = n_trials))
  durations <- rep(2L, n_trials)
  
  # Simple HRF basis (K=2)
  hrf_basis_kernels <- cbind(
    c(0, 0.2, 0.5, 0.8, 1.0, 0.7, 0.4, 0.2, 0.1, 0),
    c(0, 0.1, 0.3, 0.4, 0.3, 0.2, 0.1, 0.05, 0, 0)
  )
  
  # Random voxel weights
  coefficients <- matrix(runif(2 * n_vox, 0.5, 1.5), nrow = 2L, ncol = n_vox)
  
  # Build event matrix
  Xev <- matrix(0, n_time, n_trials)
  for (i in seq_len(n_trials)) {
    i1 <- onset_idx[i]
    i2 <- min(n_time, i1 + durations[i])
    Xev[i1:i2, i] <- 1
  }
  
  # Generate Y data
  conv_open_trim <- function(x, k) {
    as.numeric(stats::convolve(x, rev(as.numeric(k)), type = "open"))[seq_len(length(x))]
  }
  
  # Create average HRF for simulation
  avg_hrf <- drop(hrf_basis_kernels %*% rowMeans(coefficients))
  X_ref <- vapply(seq_len(n_trials), function(j) conv_open_trim(Xev[, j], avg_hrf),
                  numeric(n_time))
  
  # Add experimental regressors and nuisance
  Z <- cbind(1, scale(seq_len(n_time)))
  Nuisance <- cbind(scale(sin(seq_len(n_time) / 7)), scale(cos(seq_len(n_time) / 8)))
  
  # Generate data
  true_betas <- matrix(rnorm(n_trials * n_vox, mean = 1, sd = 0.5), n_trials, n_vox)
  Y <- X_ref %*% true_betas + 
       Z %*% matrix(rnorm(ncol(Z) * n_vox, sd = 0.3), ncol(Z), n_vox) +
       Nuisance %*% matrix(rnorm(ncol(Nuisance) * n_vox, sd = 0.2), ncol(Nuisance), n_vox) +
       matrix(rnorm(n_time * n_vox, sd = 0.4), n_time, n_vox)
  
  # Test all available backends
  available_methods <- c("r")
  results <- list()
  
  # Always test R backend
  results$r <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = durations,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    Z = Z,
    Nuisance = Nuisance,
    verbose = FALSE,
    method = "r"
  )
  
  # Test cpp if available
  cpp_available <- FALSE
  try({
    get("lss_engine_vox_hrf_cpp", envir = asNamespace("fmrilss"))
    cpp_available <- TRUE
  }, silent = TRUE)
  
  if (cpp_available) {
    available_methods <- c(available_methods, "cpp")
    results$cpp <- fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = onset_idx,
      durations = durations,
      hrf_basis_kernels = hrf_basis_kernels,
      coefficients = coefficients,
      Z = Z,
      Nuisance = Nuisance,
      verbose = FALSE,
      method = "cpp"
    )
  }
  
  # Test cpp_arma if available
  arma_available <- FALSE
  try({
    get("lss_engine_vox_hrf_arma", envir = asNamespace("fmrilss"))
    arma_available <- TRUE
  }, silent = TRUE)
  
  if (arma_available) {
    available_methods <- c(available_methods, "cpp_arma")
    results$cpp_arma <- fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = onset_idx,
      durations = durations,
      hrf_basis_kernels = hrf_basis_kernels,
      coefficients = coefficients,
      Z = Z,
      Nuisance = Nuisance,
      verbose = FALSE,
      method = "cpp_arma"
    )
  }
  
  # Test cpp_omp if available
  omp_available <- FALSE
  try({
    get("lss_engine_vox_hrf_omp", envir = asNamespace("fmrilss"))
    omp_available <- TRUE
  }, silent = TRUE)
  
  if (omp_available) {
    available_methods <- c(available_methods, "cpp_omp")
    results$cpp_omp <- fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = onset_idx,
      durations = durations,
      hrf_basis_kernels = hrf_basis_kernels,
      coefficients = coefficients,
      Z = Z,
      Nuisance = Nuisance,
      verbose = FALSE,
      method = "cpp_omp"
    )
  }
  
  # Compare all results to R backend
  if (length(available_methods) > 1) {
    for (method in available_methods[-1]) {  # Skip "r" since it's the reference
      expect_equal(results[[method]], results$r, tolerance = 1e-8,
                   info = paste("Backend", method, "should match R backend"))
    }
  }
  
  # Verify at least R backend works
  expect_equal(dim(results$r), c(n_trials, n_vox))
  expect_false(any(is.na(results$r)))
})

test_that("fallback chain works correctly", {
  skip_on_cran()
  
  set.seed(123)
  n_time <- 30L
  n_trials <- 2L
  n_vox <- 3L
  
  # Minimal test data
  onset_idx <- c(5L, 20L)
  durations <- c(1L, 1L)
  hrf_basis_kernels <- matrix(c(0, 0.5, 1, 0.5, 0), ncol = 1)
  coefficients <- matrix(1, nrow = 1L, ncol = n_vox)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  
  # Test that requesting unavailable backends falls back gracefully
  # This should always work (falls back to R)
  result <- fmrilss:::lss_with_hrf_pure_r(
    Y = Y,
    onset_idx = onset_idx,
    durations = durations,
    hrf_basis_kernels = hrf_basis_kernels,
    coefficients = coefficients,
    method = "cpp_omp",  # May not be available
    verbose = FALSE
  )
  
  expect_equal(dim(result), c(n_trials, n_vox))
  expect_false(any(is.na(result)))
})

test_that("backends handle edge cases correctly", {
  skip_on_cran()
  
  set.seed(456)
  n_time <- 50L
  n_vox <- 5L
  
  # Single trial case
  onset_idx <- 25L
  durations <- 3L
  hrf_basis_kernels <- matrix(c(0, 0.5, 1, 0.5, 0), ncol = 1)
  coefficients <- matrix(runif(n_vox, 0.5, 1.5), nrow = 1L, ncol = n_vox)
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  
  # Test each available backend
  for (method in c("r", "cpp", "cpp_arma", "cpp_omp")) {
    result <- try(fmrilss:::lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = onset_idx,
      durations = durations,
      hrf_basis_kernels = hrf_basis_kernels,
      coefficients = coefficients,
      method = method,
      verbose = FALSE
    ), silent = TRUE)
    
    if (!inherits(result, "try-error")) {
      expect_equal(dim(result), c(1L, n_vox))
      expect_false(any(is.na(result)))
    }
  }
})