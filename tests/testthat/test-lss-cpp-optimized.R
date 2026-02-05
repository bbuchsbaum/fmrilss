test_that("lss_cpp_optimized produces correct results", {
  set.seed(123)
  n_timepoints <- 80
  n_trials <- 12
  n_voxels <- 25

  # Create bdes structure (legacy interface)
  dmat_base <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 6 + 5
    if (onset + 3 <= n_timepoints) {
      dmat_ran[onset:(onset + 3), i] <- 1
    }
  }

  bdes <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    dmat_fixed = NULL
  )

  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  beta <- lss_cpp_optimized(Y, bdes)

  # Check dimensions
  expect_equal(dim(beta), c(n_trials, n_voxels))
})

test_that("lss_cpp_optimized matches other implementations", {
  set.seed(234)
  n_timepoints <- 60
  n_trials <- 10
  n_voxels <- 20

  # Create design matrices
  dmat_base <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 5 + 5
    if (onset + 3 <= n_timepoints) {
      dmat_ran[onset:(onset + 3), i] <- 1
    }
  }

  bdes <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    dmat_fixed = NULL
  )

  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Legacy interface
  beta_cpp_opt <- lss_cpp_optimized(Y, bdes)

  # Modern interface
  beta_modern <- lss(Y, dmat_ran, dmat_base, method = "cpp_optimized")

  expect_equal(unname(as.matrix(beta_cpp_opt)), unname(as.matrix(beta_modern)),
               tolerance = 1e-8)
})

test_that("lss_cpp_optimized handles dmat_fixed", {
  set.seed(345)
  n_timepoints <- 70
  n_trials <- 10
  n_voxels <- 15

  # Create design matrices
  dmat_base <- cbind(rep(1, n_timepoints))
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 6 + 5
    if (onset + 3 <= n_timepoints) {
      dmat_ran[onset:(onset + 3), i] <- 1
    }
  }

  # Additional fixed effects (e.g., linear trend)
  dmat_fixed <- cbind(scale(1:n_timepoints)[, 1])

  bdes_with_fixed <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    dmat_fixed = dmat_fixed
  )

  bdes_combined <- list(
    dmat_base = cbind(dmat_base, dmat_fixed),
    dmat_ran = dmat_ran,
    dmat_fixed = NULL
  )

  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  beta_with_fixed <- lss_cpp_optimized(Y, bdes_with_fixed)
  beta_combined <- lss_cpp_optimized(Y, bdes_combined)

  # Should produce same results
  expect_equal(beta_with_fixed, beta_combined, tolerance = 1e-8)
})

test_that("lss_cpp_optimized handles single trial", {
  set.seed(456)
  n_timepoints <- 50
  n_voxels <- 10

  dmat_base <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[, 1])
  dmat_ran <- matrix(0, n_timepoints, 1)
  dmat_ran[10:15, 1] <- 1

  bdes <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    dmat_fixed = NULL
  )

  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  beta <- lss_cpp_optimized(Y, bdes)

  expect_equal(nrow(beta), 1)
  expect_equal(ncol(beta), n_voxels)
})

test_that("project_confounds_cpp and lss_beta_cpp work correctly", {
  set.seed(567)
  n_timepoints <- 60
  n_trials <- 8
  n_voxels <- 15

  # Create confounds
  X_confounds <- cbind(
    intercept = rep(1, n_timepoints),
    linear = scale(1:n_timepoints)[, 1]
  )

  # Create trial design
  C_trials <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- (i - 1) * 7 + 5
    if (onset + 3 <= n_timepoints) {
      C_trials[onset:(onset + 3), i] <- 1
    }
  }

  Y_data <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

  # Project confounds
  result <- project_confounds_cpp(X_confounds, Y_data, C_trials)

  expect_true(is.list(result))
  expect_true("residual_data" %in% names(result))
  expect_true("Q_dmat_ran" %in% names(result))
  expect_equal(dim(result$residual_data), c(n_timepoints, n_voxels))
  expect_equal(dim(result$Q_dmat_ran), c(n_timepoints, n_trials))

  # Compute betas
  betas <- lss_beta_cpp(result$Q_dmat_ran, result$residual_data)

  expect_equal(dim(betas), c(n_trials, n_voxels))
})

test_that("project_confounds_cpp removes confound influence", {
  set.seed(678)
  n_timepoints <- 100
  n_voxels <- 10

  # Create confounds with strong signal
  X_confounds <- cbind(
    intercept = rep(1, n_timepoints),
    linear = 1:n_timepoints
  )

  # Create data with confound signal
  true_intercept <- rep(5, n_voxels)
  true_slope <- seq(0.1, 1, length.out = n_voxels)

  Y_data <- X_confounds %*% rbind(true_intercept, true_slope) +
    matrix(rnorm(n_timepoints * n_voxels, sd = 0.1), n_timepoints, n_voxels)

  C_trials <- matrix(rnorm(n_timepoints * 5), n_timepoints, 5)

  result <- project_confounds_cpp(X_confounds, Y_data, C_trials)

  # Residuals should have mean close to zero and no linear trend
  residual_means <- colMeans(result$residual_data)
  expect_true(all(abs(residual_means) < 0.5))
})
