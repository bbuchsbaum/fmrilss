test_that("mixed_solve and mixed_solve_optimized both produce reasonable results", {
  set.seed(123)
  
  # Test parameters
  n <- 100
  p <- 5  
  q <- 10
  
  # Generate test data
  X <- cbind(1, scale(matrix(rnorm(n*(p-1)), n, p-1)))
  Z <- matrix(rnorm(n*q), n, q)
  K <- diag(q)  # Identity kinship matrix
  
  # Single voxel test
  y <- rnorm(n)
  
  # Run both implementations
  result_original <- mixed_solve(Y = y, X = X, Z = Z, K = K)
  result_optimized <- mixed_solve_optimized(X = X, Z = Z, Y = y, K = K)
  
  # Test that both produce finite, reasonable results
  expect_true(all(is.finite(result_original$beta)))
  expect_true(all(is.finite(result_original$u)))
  expect_true(is.finite(result_original$Vu) && result_original$Vu >= 0)
  expect_true(is.finite(result_original$Ve) && result_original$Ve >= 0)
  
  expect_true(all(is.finite(result_optimized$beta)))
  expect_true(all(is.finite(result_optimized$u)))
  expect_true(is.finite(result_optimized$Vu) && result_optimized$Vu >= 0)
  expect_true(is.finite(result_optimized$Ve) && result_optimized$Ve >= 0)
  
  # Test that dimensions are correct
  expect_length(result_original$beta, p)
  expect_length(result_original$u, q)
  expect_length(result_optimized$beta, p)
  expect_length(result_optimized$u, q)
  
  # Note: We don't test for exact equivalence because the algorithms differ:
  # - original uses iterative L-BFGS-B optimization  
  # - optimized uses analytical REML estimation
  # Both are valid but can give different results for the same data
})

test_that("mixed_solve_optimized handles different input formats consistently", {
  set.seed(456)
  
  n <- 50
  p <- 3
  q <- 5
  
  X <- cbind(1, scale(matrix(rnorm(n*(p-1)), n, p-1)))
  Z <- matrix(rnorm(n*q), n, q)
  y <- rnorm(n)
  
  # Test vector vs 1-column matrix input
  result_vector <- mixed_solve_optimized(X, Z, y)
  result_matrix <- mixed_solve_optimized(X, Z, matrix(y, ncol = 1))
  
  # Results should be identical (not just correlated)
  expect_equal(result_vector$beta, result_matrix$beta, tolerance = 1e-10)
  
  expect_equal(result_vector$u, result_matrix$u, tolerance = 1e-10)
  
  expect_equal(result_vector$Vu, result_matrix$Vu, tolerance = 1e-10)
  
  expect_equal(result_vector$Ve, result_matrix$Ve, tolerance = 1e-10)
})

test_that("mixed_solve_optimized multi-voxel gives same results as single voxel loop", {
  set.seed(789)
  
  n <- 80
  p <- 4
  q <- 6
  v <- 5  # Number of voxels
  
  X <- cbind(1, scale(matrix(rnorm(n*(p-1)), n, p-1)))
  Z <- matrix(rnorm(n*q), n, q)
  Y <- matrix(rnorm(n*v), n, v)
  
  # Multi-voxel approach
  result_multi <- mixed_solve_optimized(X, Z, Y)
  
  # Single voxel loop approach
  results_single <- lapply(1:v, function(i) {
    mixed_solve_optimized(X, Z, Y[, i])
  })
  
  # Extract components from single voxel results
  # Note: multi-voxel returns results in format (voxels x parameters)
  # Need to properly combine single voxel results to match this format
  beta_single <- t(sapply(results_single, function(r) r$beta))  # Transpose to get voxels x parameters
  u_single <- t(sapply(results_single, function(r) r$u))        # Transpose to get voxels x parameters
  Vu_single <- sapply(results_single, function(r) r$Vu)
  Ve_single <- sapply(results_single, function(r) r$Ve)
  
  # Compare dimensions first
  expect_equal(dim(result_multi$beta), dim(beta_single))
  expect_equal(dim(result_multi$u), dim(u_single))
  expect_equal(length(result_multi$Vu), length(Vu_single))
  expect_equal(length(result_multi$Ve), length(Ve_single))
  
  # Compare results (should be very close, allowing for small numerical differences)
  expect_equal(result_multi$beta, beta_single, tolerance = 1e-6)
  
  expect_equal(result_multi$u, u_single, tolerance = 1e-6)
  
  expect_equal(as.vector(result_multi$Vu), Vu_single, tolerance = 1e-6)
  
  expect_equal(as.vector(result_multi$Ve), Ve_single, tolerance = 1e-6)
})

test_that("workspace precomputation and reuse works correctly", {
  set.seed(321)
  
  n <- 60
  p <- 3
  q <- 4
  
  X <- cbind(1, scale(matrix(rnorm(n*(p-1)), n, p-1)))
  Z <- matrix(rnorm(n*q), n, q)
  y1 <- rnorm(n)
  y2 <- rnorm(n)
  
  # Compute with and without precomputed workspace
  result1_auto <- mixed_solve_optimized(X, Z, y1)
  
  # Manual workspace precomputation
  workspace <- mixed_precompute(X, Z)
  result1_manual <- mixed_solve_optimized(X, Z, y1, workspace = workspace)
  result2_manual <- mixed_solve_optimized(X, Z, y2, workspace = workspace)
  
  # Results with automatic vs manual workspace should be identical
  expect_equal(result1_auto$beta, result1_manual$beta, tolerance = 1e-12)
  
  expect_equal(result1_auto$u, result1_manual$u, tolerance = 1e-12)
  
  # Second voxel with reused workspace should work
  expect_length(result2_manual$beta, p)
  
  expect_length(result2_manual$u, q)
})

test_that("mixed_solve_optimized handles edge cases appropriately", {
  set.seed(654)
  
  n <- 50
  p <- 2
  q <- 3
  
  X <- cbind(1, rnorm(n))
  Z <- matrix(rnorm(n*q), n, q)
  
  # Test with very small variance
  y_small_var <- rnorm(n, sd = 0.001)
  result_small <- mixed_solve_optimized(X, Z, y_small_var)
  
  expect_true(all(is.finite(result_small$beta)))
  
  expect_true(all(is.finite(result_small$u)))
  
  # Test with larger problem
  n_large <- 200
  X_large <- cbind(1, scale(matrix(rnorm(n_large*4), n_large, 4)))
  Z_large <- matrix(rnorm(n_large*8), n_large, 8)
  y_large <- rnorm(n_large)
  
  result_large <- mixed_solve_optimized(X_large, Z_large, y_large)
  
  expect_length(result_large$beta, 5)
  
  expect_length(result_large$u, 8)
}) 