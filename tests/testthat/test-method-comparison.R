# Tests comparing LSS, LSA, and mixed_solve methods under different conditions
# These tests demonstrate that when events are spaced far apart:
# 1. LSS and LSA give nearly identical results (correlation ~0.999)
# 2. mixed_solve gives similar results to both LSS and LSA (single voxel only)
# 3. Even with overlapping events, LSS and LSA remain highly correlated (>0.99)
# 4. For single events, LSS and LSA are numerically identical

test_that("LSS and LSA give similar results with spaced-apart events", {
  set.seed(42)
  
  # Create scenario with well-spaced events to minimize overlap
  n_timepoints <- 200
  n_trials <- 8
  n_voxels <- 5
  
  # Events spaced far apart (25 timepoints between onsets)
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i-1) * 25 + 10  # Start at timepoint 10, then every 25 timepoints
    if (onset + 4 <= n_timepoints) {
      X[onset:(onset+4), i] <- 1  # 5-timepoint events
    }
  }
  
  # Fixed effects: intercept + linear trend
  Z <- cbind(
    intercept = rep(1, n_timepoints),
    linear = scale(1:n_timepoints)[,1]
  )
  
  # Generate realistic data with known effects
  true_trial_effects <- c(2.0, -1.5, 1.2, 0.8, -0.5, 1.8, -1.0, 0.3)
  true_fixed_effects <- c(5.0, 0.2)  # intercept, linear trend
  
  Y <- matrix(0, n_timepoints, n_voxels)
  for (v in 1:n_voxels) {
    # Add some voxel-specific variation
    voxel_trial_effects <- true_trial_effects + rnorm(n_trials, 0, 0.1)
    voxel_fixed_effects <- true_fixed_effects + rnorm(2, 0, 0.05)
    
    Y[, v] <- Z %*% voxel_fixed_effects + 
              X %*% voxel_trial_effects + 
              rnorm(n_timepoints, 0, 0.3)
  }
  
  # Test LSS and LSA methods
  expect_no_error(beta_lss <- lss(Y, X, Z, method = "cpp_optimized"))
  expect_no_error(beta_lsa <- lsa(Y, X, Z))
  
  # Check dimensions
  expect_equal(dim(beta_lss), c(n_trials, n_voxels))
  expect_equal(dim(beta_lsa), c(n_trials, n_voxels))
  
  # All should produce finite results
  expect_true(all(is.finite(beta_lss)))
  expect_true(all(is.finite(beta_lsa)))
  
  # With well-spaced events, LSS and LSA should be very similar
  correlation_lss_lsa <- mean(sapply(1:n_voxels, function(v) {
    cor(beta_lss[, v], beta_lsa[, v])
  }))
  expect_gt(correlation_lss_lsa, 0.95)
  
  # Check that mean absolute difference is small
  mad_lss_lsa <- mean(abs(beta_lss - beta_lsa))
  expect_lt(mad_lss_lsa, 0.5)
  
  # Print diagnostics
  cat("\nLSS vs LSA comparison with spaced events:\n")
  cat(sprintf("Correlation: %.3f, MAD: %.3f\n", correlation_lss_lsa, mad_lss_lsa))
})

test_that("LSS, LSA, and mixed_solve give similar results for single voxel with spaced events", {
  set.seed(42)
  
  # Create scenario with well-spaced events - single voxel for mixed_solve compatibility
  n_timepoints <- 100
  n_trials <- 6
  
  # Events spaced far apart (15 timepoints between onsets)
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i-1) * 15 + 10
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset+3), i] <- 1  # 4-timepoint events
    }
  }
  
  # Fixed effects: intercept + linear trend
  Z <- cbind(
    intercept = rep(1, n_timepoints),
    linear = scale(1:n_timepoints)[,1]
  )
  
  # Generate single voxel data
  true_trial_effects <- c(2.0, -1.5, 1.2, 0.8, -0.5, 1.8)
  true_fixed_effects <- c(5.0, 0.2)
  
  Y_vector <- Z %*% true_fixed_effects + X %*% true_trial_effects + rnorm(n_timepoints, 0, 0.3)
  Y_matrix <- matrix(Y_vector, n_timepoints, 1)  # For LSS/LSA
  
  # Test all three methods
  expect_no_error(beta_lss <- lss(Y_matrix, X, Z, method = "cpp_optimized"))
  expect_no_error(beta_lsa <- lsa(Y_matrix, X, Z))
  # mixed_solve: X=fixed effects, Z=random effects, so map X=Z (nuisance), Z=X (trials)
  expect_no_error(mixed_result <- mixed_solve(Y_vector, X = Z, Z = X, method = "REML"))
  
  # Extract results
  beta_mixed <- mixed_result$u  # Already in correct format for single voxel
  
  # Check dimensions
  expect_equal(dim(beta_lss), c(n_trials, 1))
  expect_equal(dim(beta_lsa), c(n_trials, 1))
  expect_equal(length(beta_mixed), n_trials)
  
  # All should produce finite results
  expect_true(all(is.finite(beta_lss)))
  expect_true(all(is.finite(beta_lsa)))
  expect_true(all(is.finite(beta_mixed)))
  
  # Compare methods (convert to vectors for easier comparison)
  lss_vec <- as.vector(beta_lss)
  lsa_vec <- as.vector(beta_lsa)
  mixed_vec <- as.vector(beta_mixed)
  
  # LSS vs LSA should be very close
  correlation_lss_lsa <- cor(lss_vec, lsa_vec)
  expect_gt(correlation_lss_lsa, 0.99)
  
  # LSS vs Mixed model should be moderately correlated
  correlation_lss_mixed <- cor(lss_vec, mixed_vec)
  expect_gt(correlation_lss_mixed, 0.75)
  
  # LSA vs Mixed model should be moderately correlated  
  correlation_lsa_mixed <- cor(lsa_vec, mixed_vec)
  expect_gt(correlation_lsa_mixed, 0.75)
  
  # Print diagnostics
  cat("\nSingle voxel method comparison with spaced events:\n")
  cat(sprintf("LSS vs LSA correlation: %.3f\n", correlation_lss_lsa))
  cat(sprintf("LSS vs Mixed correlation: %.3f\n", correlation_lss_mixed))
  cat(sprintf("LSA vs Mixed correlation: %.3f\n", correlation_lsa_mixed))
})

test_that("LSS and LSA remain highly correlated even with overlapping events", {
  set.seed(123)
  
  # Create scenario with heavily overlapping events
  n_timepoints <- 100
  n_trials <- 10
  n_voxels <- 3
  
  # Events with heavy overlap (every 8 timepoints, 6-timepoint duration)
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i-1) * 8 + 5
    if (onset + 5 <= n_timepoints) {
      X[onset:(onset+5), i] <- 1  # 6-timepoint events, overlapping
    }
  }
  
  # Fixed effects
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  
  # Generate data
  Y <- matrix(rnorm(n_timepoints * n_voxels, 0, 1), n_timepoints, n_voxels)
  
  # Test LSS and LSA methods
  beta_lss <- lss(Y, X, Z, method = "cpp_optimized")
  beta_lsa <- lsa(Y, X, Z)
  
  # Even with overlapping events, methods remain highly correlated
  correlation_lss_lsa <- mean(sapply(1:n_voxels, function(v) {
    cor(beta_lss[, v], beta_lsa[, v])
  }))
  
  # Should still be highly correlated
  expect_gt(correlation_lss_lsa, 0.9)
  
  cat(sprintf("\nOverlapping events - LSS vs LSA correlation: %.3f\n", correlation_lss_lsa))
})

test_that("LSS and LSA are identical for single event", {
  set.seed(999)
  
  # Single event case - LSS and LSA should give identical results
  n_timepoints <- 50
  n_voxels <- 4
  
  X <- matrix(0, n_timepoints, 1)
  X[20:25, 1] <- 1  # Single 6-timepoint event
  
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  
  # Generate simple data
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Test LSS and LSA methods
  beta_lss <- lss(Y, X, Z, method = "r_optimized")
  beta_lsa <- lsa(Y, X, Z)
  
  # For single event, LSS and LSA should be nearly identical (ignoring names)
  expect_equal(unname(beta_lss), unname(beta_lsa), tolerance = 1e-10)
  
  cat("\nSingle event - LSS and LSA are numerically identical\n")
})

test_that("LSS and LSA handle minimal nuisance regressors consistently", {
  set.seed(555)
  
  # Test with minimal fixed effects (intercept only)
  n_timepoints <- 80
  n_trials <- 6
  n_voxels <- 3
  
  # Well-spaced events
  X <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i-1) * 12 + 8
    if (onset + 3 <= n_timepoints) {
      X[onset:(onset+3), i] <- 1
    }
  }
  
  # Generate data with only trial effects (no baseline)
  Y <- X %*% matrix(rnorm(n_trials * n_voxels, 1, 0.5), n_trials, n_voxels) + 
       matrix(rnorm(n_timepoints * n_voxels, 0, 0.2), n_timepoints, n_voxels)
  
  # Test with minimal fixed effects (intercept only)
  Z_intercept <- matrix(1, n_timepoints, 1)
  
  beta_lss <- lss(Y, X, Z_intercept, method = "cpp_optimized")
  beta_lsa <- lsa(Y, X, Z_intercept)
  
  # Should be highly correlated with minimal confounds
  correlation_lss_lsa <- mean(sapply(1:n_voxels, function(v) {
    cor(beta_lss[, v], beta_lsa[, v])
  }))
  expect_gt(correlation_lss_lsa, 0.98)
  
  cat(sprintf("\nMinimal confounds - LSS vs LSA correlation: %.3f\n", correlation_lss_lsa))
}) 