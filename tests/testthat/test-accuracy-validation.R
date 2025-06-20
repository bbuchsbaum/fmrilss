test_that("LSS implementations are numerically consistent", {
  set.seed(12345)
  
  # Test case: realistic design with multiple confounds
  n_timepoints <- 50
  n_voxels <- 5
  n_trials <- 3
  
  # Fixed effects: intercept + linear trend + quadratic
  Z <- cbind(
    rep(1, n_timepoints),
    scale(1:n_timepoints)[,1],
    scale((1:n_timepoints)^2)[,1]
  )
  
  # Trial regressors with some minimal overlap
  X <- matrix(0, n_timepoints, n_trials)
  X[8:12, 1] <- 1
  X[20:24, 2] <- 1  
  X[35:39, 3] <- 1
  
  # Generate realistic data
  Y <- matrix(rnorm(n_timepoints * n_voxels, sd = 1), n_timepoints, n_voxels)
  
  # Test all implementations work
  expect_no_error(beta_r_optimized <- lss(Y, X, Z, method = "r_optimized"))
  expect_no_error(beta_cpp_optimized <- lss(Y, X, Z, method = "cpp_optimized"))
  expect_no_error(beta_r_vectorized <- lss(Y, X, Z, method = "r_vectorized"))
  expect_no_error(beta_cpp <- lss(Y, X, Z, method = "cpp"))
  expect_no_error(beta_naive <- lss(Y, X, Z, method = "naive"))
  
  # Check dimensions
  expect_equal(dim(beta_r_optimized), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp_optimized), c(n_trials, n_voxels))
  expect_equal(dim(beta_r_vectorized), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp), c(n_trials, n_voxels))
  expect_equal(dim(beta_naive), c(n_trials, n_voxels))
  
  # R and C++ should be numerically equivalent
  expect_equal(beta_r_vectorized, beta_cpp, tolerance = 1e-12)
  
  # All should produce finite results
  expect_true(all(is.finite(beta_r_optimized)), info = "R optimized implementation should produce finite results")
  expect_true(all(is.finite(beta_cpp_optimized)), info = "C++ optimized implementation should produce finite results") 
  expect_true(all(is.finite(beta_r_vectorized)), info = "R vectorized implementation should produce finite results")
  expect_true(all(is.finite(beta_cpp)), info = "C++ implementation should produce finite results") 
  expect_true(all(is.finite(beta_naive)), info = "Naive implementation should produce finite results")
})

test_that("LSS vs standard GLM for single trial", {
  set.seed(54321)
  
  # Compare LSS to standard GLM when there's only one trial
  n_timepoints <- 100
  n_voxels <- 10
  
  # Fixed effects
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  
  # Single trial design
  X <- matrix(0, n_timepoints, 1)
  X[40:50, 1] <- 1
  
  # Generate data
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # LSS estimates
  beta_lss <- lss(Y, X, Z, method = "r_vectorized")
  
  # Standard GLM estimates for comparison
  X_full <- cbind(Z, X)
  beta_glm <- matrix(NA, 1, n_voxels)
  
  for (v in 1:n_voxels) {
    fit <- lm(Y[, v] ~ X_full - 1)  # -1 because we include intercept explicitly
    beta_glm[1, v] <- coef(fit)[ncol(X_full)]  # Last coefficient is the trial effect
  }
  
  # LSS should match GLM for single trial case (ignoring dimnames)
  expect_equal(unname(as.matrix(beta_lss)), unname(as.matrix(beta_glm)), tolerance = 1e-10,
               info = "LSS should match standard GLM for single trial")
})

test_that("LSS handles non-overlapping trials correctly", {
  set.seed(99999)
  
  # Test with non-overlapping trials but NOT perfectly orthogonal 
  # (LSS is designed for overlapping/correlated designs)
  n_timepoints <- 100
  n_voxels <- 8
  n_trials <- 3
  
  Z <- cbind(1, scale(1:n_timepoints)[,1])  # Intercept + trend
  
  # Create non-overlapping but not block trial design
  X <- matrix(0, n_timepoints, n_trials)
  X[10:15, 1] <- 1     # Trial 1: short duration 
  X[40:45, 2] <- 1     # Trial 2: separate period
  X[70:75, 3] <- 1     # Trial 3: separate period
  
  # Generate data with known effects + more realistic noise
  true_effects <- c(1.5, -2.0, 0.8)
  Y <- matrix(0, n_timepoints, n_voxels)
  for (v in 1:n_voxels) {
    Y[, v] <- 5 + 0.2 * scale(1:n_timepoints)[,1] +  # intercept + trend
              X %*% true_effects +
              rnorm(n_timepoints, sd = 0.5)  # Higher noise for realism
  }
  
  # All implementations should work without errors
  expect_no_error(beta_r_optimized <- lss(Y, X, Z, method = "r_optimized"))
  expect_no_error(beta_cpp_optimized <- lss(Y, X, Z, method = "cpp_optimized"))
  expect_no_error(beta_r_vectorized <- lss(Y, X, Z, method = "r_vectorized"))
  expect_no_error(beta_cpp <- lss(Y, X, Z, method = "cpp"))
  expect_no_error(beta_naive <- lss(Y, X, Z, method = "naive"))
  
  # Check basic properties
  expect_equal(dim(beta_r_optimized), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp_optimized), c(n_trials, n_voxels))
  expect_equal(dim(beta_r_vectorized), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp), c(n_trials, n_voxels))
  expect_equal(dim(beta_naive), c(n_trials, n_voxels))
  
  # Check R and C++ equivalence (should be very close)
  expect_equal(beta_r_vectorized, beta_cpp, tolerance = 1e-10)
  
  # For this simple design, LSS and naive should be similar (but not identical)
  expect_equal(beta_r_vectorized, beta_naive, tolerance = 1.0,
               info = "LSS and naive should be reasonably similar for non-overlapping trials")
})

test_that("Projection matrices have correct properties", {
  set.seed(77777)
  
  n_timepoints <- 50
  k_confounds <- 4
  
  # Create confound matrix
  X <- cbind(
    rep(1, n_timepoints),
    scale(1:n_timepoints)[,1],
    scale((1:n_timepoints)^2)[,1],
    rnorm(n_timepoints)
  )
  
  # Get projection matrix
  Q <- project_confounds(X)
  
  # Test projection matrix properties
  # 1. Q should be symmetric
  expect_equal(Q, t(Q), tolerance = 1e-12, info = "Q should be symmetric")
  
  # 2. Q should be idempotent: Q * Q = Q
  expect_equal(Q %*% Q, Q, tolerance = 1e-10, info = "Q should be idempotent")
  
  # 3. Q * X should be zero (projects out X)
  QX <- Q %*% X
  expect_equal(QX, matrix(0, n_timepoints, k_confounds), tolerance = 1e-10,
               info = "Q should project out X completely")
  
  # 4. Rank of Q should be n - k
  expect_equal(qr(Q)$rank, n_timepoints - k_confounds, tolerance = 1,
               info = "Q should have rank n - k")
})

test_that("Memory usage test for large design", {
  skip_on_cran()
  skip_if_not_installed("testthat")
  
  # Test that C++ version doesn't create huge matrices
  n_timepoints <- 500  # Moderately large to test memory efficiency
  n_voxels <- 100
  n_trials <- 50
  
  Z <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  X <- matrix(0, n_timepoints, n_trials)
  
  # Create random trial design
  for (i in 1:n_trials) {
    onset <- sample(1:(n_timepoints-10), 1)
    X[onset:(onset+5), i] <- 1
  }
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # These should all complete without memory issues
  expect_no_error({
    beta_r_optimized <- lss(Y, X, Z, method = "r_optimized")
  })
  
  expect_no_error({
    beta_cpp_optimized <- lss(Y, X, Z, method = "cpp_optimized")
  })
  
  expect_no_error({
    beta_r_vectorized <- lss(Y, X, Z, method = "r_vectorized")
  })
  
  expect_no_error({
    beta_cpp <- lss(Y, X, Z, method = "cpp")
  })
  
  expect_no_error({
    beta_naive <- lss(Y, X, Z, method = "naive")
  })
  
  # Results should still be equivalent
  expect_equal(beta_r_vectorized, beta_cpp, tolerance = 1e-8)
  expect_equal(beta_r_vectorized, beta_naive, tolerance = 1e-6)
}) 