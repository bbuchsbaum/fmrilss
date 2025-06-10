test_that("All LSS implementations produce equivalent results", {
  skip_if_not_installed("testthat")
  
  # Test parameters
  set.seed(123)
  n_timepoints <- 120
  n_trials <- 20
  n_voxels <- 50
  
  # Create trial design matrix: each trial has a 5-timepoint boxcar
  X <- matrix(0, n_timepoints, n_trials)
  trial_onsets <- seq(10, n_timepoints - 15, length.out = n_trials)
  for (i in 1:n_trials) {
    onset <- round(trial_onsets[i])
    X[onset:(onset + 4), i] <- 1
  }
  
  # Fixed effects: intercept + linear trend + quadratic trend
  Z <- cbind(
    intercept = rep(1, n_timepoints),
    linear = scale(1:n_timepoints)[,1],
    quadratic = scale((1:n_timepoints)^2)[,1]
  )
  
  # Generate synthetic data with known signal
  true_betas <- matrix(rnorm(n_trials * n_voxels, mean = 0, sd = 2), n_trials, n_voxels)
  base_effects <- matrix(rnorm(3 * n_voxels, mean = 1, sd = 0.5), 3, n_voxels)
  noise <- matrix(rnorm(n_timepoints * n_voxels, mean = 0, sd = 1), n_timepoints, n_voxels)
  
  # Construct data: Y = Z * base_effects + X * true_betas + noise
  Y <- Z %*% base_effects + X %*% true_betas + noise
  
  # Test all methods with new interface
  beta_r_optimized <- lss(Y, X, Z, method = "r_optimized")
  beta_cpp_optimized <- lss(Y, X, Z, method = "cpp_optimized") 
  beta_r_vectorized <- lss(Y, X, Z, method = "r_vectorized")
  beta_cpp <- lss(Y, X, Z, method = "cpp")
  beta_naive <- lss(Y, X, Z, method = "naive")
  
  # Check dimensions
  expect_equal(dim(beta_r_optimized), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp_optimized), c(n_trials, n_voxels))
  expect_equal(dim(beta_r_vectorized), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp), c(n_trials, n_voxels))
  expect_equal(dim(beta_naive), c(n_trials, n_voxels))
  
  # Check row/column names
  expect_equal(rownames(beta_r_optimized), paste0("Trial_", 1:n_trials))
  expect_equal(colnames(beta_r_optimized), paste0("Voxel_", 1:n_voxels))
  
  # Check numerical equivalence (allowing for small numerical differences)
  expect_equal(beta_r_vectorized, beta_cpp, tolerance = 1e-8, 
               info = "R vectorized and C++ implementations should be numerically equivalent")
  
  expect_equal(beta_r_vectorized, beta_naive, tolerance = 1e-6,
               info = "R vectorized and naive implementations should be numerically equivalent")
  
  expect_equal(beta_cpp, beta_naive, tolerance = 1e-6,
               info = "C++ and naive implementations should be numerically equivalent")
  
  expect_equal(beta_r_optimized, beta_r_vectorized, tolerance = 1e-8,
               info = "Optimized R should be equivalent to vectorized R implementation")
  
  expect_equal(beta_cpp_optimized, beta_cpp, tolerance = 1e-6,
               info = "Optimized C++ should be equivalent to standard C++ implementation")
})

test_that("LSS implementations handle edge cases correctly", {
  set.seed(456)
  
  # Test 1: Single trial case
  n_timepoints <- 50
  n_voxels <- 10
  
  Z <- cbind(rep(1, n_timepoints), 1:n_timepoints)
  X <- matrix(0, n_timepoints, 1)
  X[10:15, 1] <- 1
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  beta_r_optimized <- lss(Y, X, Z, method = "r_optimized")
  beta_cpp_optimized <- lss(Y, X, Z, method = "cpp_optimized")
  beta_r_vectorized <- lss(Y, X, Z, method = "r_vectorized")
  beta_cpp <- lss(Y, X, Z, method = "cpp")
  beta_naive <- lss(Y, X, Z, method = "naive")
  
  expect_equal(beta_r_vectorized, beta_cpp, tolerance = 1e-8)
  expect_equal(beta_r_vectorized, beta_naive, tolerance = 1e-6)
  expect_equal(beta_r_optimized, beta_r_vectorized, tolerance = 1e-8)
  expect_equal(beta_cpp_optimized, beta_cpp, tolerance = 1e-6)
  
  # Test 2: Default Z matrix (intercept only)
  beta_default_z <- lss(Y, X, method = "r_optimized")
  Z_intercept <- matrix(1, n_timepoints, 1)
  beta_explicit_z <- lss(Y, X, Z_intercept, method = "r_optimized")
  
  expect_equal(beta_default_z, beta_explicit_z, tolerance = 1e-12,
               info = "Default Z (intercept) should match explicit intercept matrix")
})

test_that("Nuisance regression works correctly", {
  set.seed(789)
  
  n_timepoints <- 100
  n_trials <- 15
  n_voxels <- 30
  
  # Create trial design matrix
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i-1) * 6 + 5
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset+3), i] <- 1
    }
  }
  
  # Fixed effects
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  
  # Nuisance regressors (e.g., motion parameters)
  Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)
  
  # Create data with nuisance signal
  true_betas <- matrix(rnorm(n_trials * n_voxels, 0, 1), n_trials, n_voxels)
  nuisance_effects <- matrix(rnorm(6 * n_voxels, 0, 2), 6, n_voxels)
  Y <- X %*% true_betas + Z %*% matrix(rnorm(2 * n_voxels), 2, n_voxels) + 
       Nuisance %*% nuisance_effects + matrix(rnorm(n_timepoints * n_voxels, 0, 0.5), n_timepoints, n_voxels)
  
  # Test with and without nuisance regression
  beta_no_nuisance <- lss(Y, X, Z, method = "r_optimized")
  beta_with_nuisance <- lss(Y, X, Z, Nuisance = Nuisance, method = "r_optimized")
  
  # Results should be different (nuisance regression should change estimates)
  expect_false(isTRUE(all.equal(beta_no_nuisance, beta_with_nuisance, tolerance = 1e-3)),
               info = "Nuisance regression should change beta estimates")
  
  # Test that all methods give same result with nuisance regression
  beta_cpp_nuisance <- lss(Y, X, Z, Nuisance = Nuisance, method = "cpp_optimized")
  beta_naive_nuisance <- lss(Y, X, Z, Nuisance = Nuisance, method = "naive")
  
  expect_equal(beta_with_nuisance, beta_cpp_nuisance, tolerance = 1e-6)
  expect_equal(beta_with_nuisance, beta_naive_nuisance, tolerance = 1e-6)
})

test_that("Input validation works correctly", {
  n_timepoints <- 50
  n_trials <- 5
  n_voxels <- 10
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  X <- matrix(rnorm(n_timepoints * n_trials), n_timepoints, n_trials)
  Z <- matrix(1, n_timepoints, 1)
  
  # Test dimension mismatches
  Y_wrong <- matrix(rnorm(40 * n_voxels), 40, n_voxels)  # Wrong timepoints
  X_wrong <- matrix(rnorm(40 * n_trials), 40, n_trials)  # Wrong timepoints  
  Z_wrong <- matrix(1, 40, 1)  # Wrong timepoints
  Nuisance_wrong <- matrix(rnorm(40 * 3), 40, 3)  # Wrong timepoints
  
  expect_error(lss(Y_wrong, X), "Y and X must have the same number of rows")
  expect_error(lss(Y, X_wrong), "Y and X must have the same number of rows")
  expect_error(lss(Y, X, Z_wrong), "Z must be a numeric matrix with the same number of rows as Y")
  expect_error(lss(Y, X, Nuisance = Nuisance_wrong), "Nuisance must be a numeric matrix with the same number of rows as Y")
  
  # Test non-matrix inputs
  expect_error(lss("not_matrix", X), "Y must be a numeric matrix")
  expect_error(lss(Y, "not_matrix"), "X must be a numeric matrix")
  expect_error(lss(Y, X, "not_matrix"), "Z must be a numeric matrix")
  expect_error(lss(Y, X, Nuisance = "not_matrix"), "Nuisance must be a numeric matrix")
  
  # Test invalid method
  expect_error(lss(Y, X, method = "invalid_method"), "should be one of")
})

test_that("Column and row naming works correctly", {
  set.seed(999)
  
  n_timepoints <- 60
  n_trials <- 8
  n_voxels <- 15
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  X <- matrix(0, n_timepoints, n_trials)
  
  # Add some signal
  for (i in 1:n_trials) {
    onset <- i * 7
    if (onset + 2 <= n_timepoints) {
      X[onset:(onset+2), i] <- 1
    }
  }
  
  # Test with named matrices
  colnames(Y) <- paste0("Region_", 1:n_voxels)
  colnames(X) <- paste0("Condition_", 1:n_trials)
  
  beta <- lss(Y, X)
  
  expect_equal(rownames(beta), paste0("Condition_", 1:n_trials))
  expect_equal(colnames(beta), paste0("Region_", 1:n_voxels))
  
  # Test without names (should get defaults)
  colnames(Y) <- NULL
  colnames(X) <- NULL
  
  beta_no_names <- lss(Y, X)
  
  expect_equal(rownames(beta_no_names), paste0("Trial_", 1:n_trials))
  expect_equal(colnames(beta_no_names), paste0("Voxel_", 1:n_voxels))
})

test_that("LSS handles numerical edge cases", {
  set.seed(131415)
  
  n_timepoints <- 60
  n_voxels <- 20
  n_trials <- 8
  
  # Create design with some nearly collinear trials
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  X <- matrix(0, n_timepoints, n_trials)
  
  # First few trials with some overlap
  X[10:15, 1] <- 1
  X[12:17, 2] <- 0.8  # Overlapping with trial 1
  X[25:30, 3] <- 1
  X[40:45, 4] <- 1
  
  # Remaining trials non-overlapping
  for (i in 5:n_trials) {
    onset <- 50 + (i-5)*2
    if (onset + 2 <= n_timepoints) {
      X[onset:(onset+2), i] <- 1
    }
  }
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # All methods should handle this without error
  expect_no_error(beta_r <- lss(Y, X, Z, method = "r_optimized"))
  expect_no_error(beta_cpp <- lss(Y, X, Z, method = "cpp_optimized"))
  expect_no_error(beta_naive <- lss(Y, X, Z, method = "naive"))
  
  # Results should still be similar (allowing for larger tolerance due to numerical issues)
  expect_equal(beta_r, beta_cpp, tolerance = 1e-6)
  expect_equal(beta_r, beta_naive, tolerance = 1e-4)
})

test_that("Backward compatibility with old interface works", {
  # Test that old bdes interface still works for existing code
  set.seed(246)
  
  n_timepoints <- 80
  n_trials <- 10
  n_voxels <- 25
  
  # Old interface
  dmat_base <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i-1) * 7 + 5
    if (onset + 3 <= n_timepoints) {
      dmat_ran[onset:(onset+3), i] <- 1
    }
  }
  
  bdes <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    fixed_ind = NULL
  )
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Test old functions still work
  beta_old_r <- lss_fast(dset = NULL, bdes = bdes, Y = Y, use_cpp = FALSE)
  beta_old_cpp <- lss_fast(dset = NULL, bdes = bdes, Y = Y, use_cpp = TRUE)
  beta_naive_old <- lss_naive(Y = Y, bdes = bdes)
  beta_optimized_old <- lss_optimized(Y = Y, bdes = bdes)
  
  # Test new interface gives same results
  beta_new_r <- lss(Y, dmat_ran, dmat_base, method = "r_vectorized")
  beta_new_cpp <- lss(Y, dmat_ran, dmat_base, method = "cpp")
  beta_new_naive <- lss(Y, dmat_ran, dmat_base, method = "naive")
  beta_new_optimized <- lss(Y, dmat_ran, dmat_base, method = "r_optimized")
  
  # Should be equivalent (ignoring row/column names)
  expect_equal(unname(as.matrix(beta_old_r)), unname(as.matrix(beta_new_r)), tolerance = 1e-12)
  expect_equal(unname(as.matrix(beta_old_cpp)), unname(as.matrix(beta_new_cpp)), tolerance = 1e-12)
  expect_equal(unname(as.matrix(beta_naive_old)), unname(as.matrix(beta_new_naive)), tolerance = 1e-12)
  expect_equal(unname(as.matrix(beta_optimized_old)), unname(as.matrix(beta_new_optimized)), tolerance = 1e-12)
}) 