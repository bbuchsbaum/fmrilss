test_that("All LSS implementations produce equivalent results", {
  skip_if_not_installed("testthat")
  
  # Test parameters
  set.seed(123)
  n_timepoints <- 120
  n_trials <- 20
  n_voxels <- 50
  
  # Create synthetic design matrices
  # Base design: intercept + linear trend + quadratic trend
  dmat_base <- cbind(
    intercept = rep(1, n_timepoints),
    linear = scale(1:n_timepoints)[,1],
    quadratic = scale((1:n_timepoints)^2)[,1]
  )
  
  # Trial design matrix: each trial has a 5-timepoint boxcar
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  trial_onsets <- seq(10, n_timepoints - 15, length.out = n_trials)
  for (i in 1:n_trials) {
    onset <- round(trial_onsets[i])
    dmat_ran[onset:(onset + 4), i] <- 1
  }
  
  # Create design list
  bdes <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    fixed_ind = NULL,
    dmat_fixed = NULL
  )
  
  # Generate synthetic data with known signal
  true_betas <- matrix(rnorm(n_trials * n_voxels, mean = 0, sd = 2), n_trials, n_voxels)
  base_effects <- matrix(rnorm(3 * n_voxels, mean = 1, sd = 0.5), 3, n_voxels)
  noise <- matrix(rnorm(n_timepoints * n_voxels, mean = 0, sd = 1), n_timepoints, n_voxels)
  
  # Construct data: Y = X_base * base_effects + dmat_ran * true_betas + noise
  Y <- dmat_base %*% base_effects + dmat_ran %*% true_betas + noise
  
  # Test 1: Main LSS functions equivalence
  beta_r <- lss(Y = Y, bdes = bdes, use_cpp = FALSE)
  beta_cpp <- lss(Y = Y, bdes = bdes, use_cpp = TRUE)
  beta_naive <- lss_naive(Y = Y, bdes = bdes)
  
  # Check dimensions
  expect_equal(dim(beta_r), c(n_trials, n_voxels))
  expect_equal(dim(beta_cpp), c(n_trials, n_voxels))
  expect_equal(dim(beta_naive), c(n_trials, n_voxels))
  
  # Check numerical equivalence (allowing for small numerical differences)
  expect_equal(beta_r, beta_cpp, tolerance = 1e-8, 
               info = "R and C++ implementations should be numerically equivalent")
  
  expect_equal(beta_r, beta_naive, tolerance = 1e-6,
               info = "R and naive implementations should be numerically equivalent")
  
  expect_equal(beta_cpp, beta_naive, tolerance = 1e-6,
               info = "C++ and naive implementations should be numerically equivalent")
})

test_that("LSS implementations handle edge cases correctly", {
  set.seed(456)
  
  # Test 2: Single trial case
  n_timepoints <- 50
  n_voxels <- 10
  
  dmat_base <- cbind(rep(1, n_timepoints), 1:n_timepoints)
  dmat_ran <- matrix(0, n_timepoints, 1)
  dmat_ran[10:15, 1] <- 1
  
  bdes_single <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    fixed_ind = NULL
  )
  
  Y_single <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  beta_r_single <- lss(Y = Y_single, bdes = bdes_single, use_cpp = FALSE)
  beta_cpp_single <- lss(Y = Y_single, bdes = bdes_single, use_cpp = TRUE)
  beta_naive_single <- lss_naive(Y = Y_single, bdes = bdes_single)
  
  expect_equal(beta_r_single, beta_cpp_single, tolerance = 1e-8)
  expect_equal(beta_r_single, beta_naive_single, tolerance = 1e-6)
  
  # Test 3: Case with fixed effects
  dmat_fixed <- matrix(rnorm(n_timepoints * 2), n_timepoints, 2)
  bdes_fixed <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    dmat_fixed = dmat_fixed,
    fixed_ind = 1:2
  )
  
  beta_r_fixed <- lss(Y = Y_single, bdes = bdes_fixed, use_cpp = FALSE)
  beta_cpp_fixed <- lss(Y = Y_single, bdes = bdes_fixed, use_cpp = TRUE)
  beta_naive_fixed <- lss_naive(Y = Y_single, bdes = bdes_fixed)
  
  expect_equal(beta_r_fixed, beta_cpp_fixed, tolerance = 1e-8)
  expect_equal(beta_r_fixed, beta_naive_fixed, tolerance = 1e-6)
})

test_that("Projection functions are equivalent", {
  set.seed(789)
  
  n_timepoints <- 100
  n_voxels <- 30
  n_trials <- 15
  
  # Create test matrices
  X_confounds <- cbind(
    rep(1, n_timepoints),
    scale(1:n_timepoints)[,1],
    rnorm(n_timepoints)
  )
  Y_data <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  C_trials <- matrix(rnorm(n_timepoints * n_trials), n_timepoints, n_trials)
  
  # Test R vs C++ projection
  Q_r <- project_confounds(X_confounds)
  Y_proj_r <- Q_r %*% Y_data
  C_proj_r <- Q_r %*% C_trials
  
  result_cpp <- project_confounds_cpp(X_confounds, Y_data, C_trials)
  Y_proj_cpp <- result_cpp$residual_data
  C_proj_cpp <- result_cpp$Q_dmat_ran
  
  expect_equal(Y_proj_r, Y_proj_cpp, tolerance = 1e-10,
               info = "R and C++ projection of Y should be equivalent")
  expect_equal(C_proj_r, C_proj_cpp, tolerance = 1e-10,
               info = "R and C++ projection of C should be equivalent")
})

test_that("C++ functions match internal implementations", {
  set.seed(101112)
  
  n_timepoints <- 80
  n_voxels <- 25
  n_trials <- 12
  
  # Create test data
  dmat_base <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  for (i in 1:n_trials) {
    onset <- sample(1:(n_timepoints-8), 1)
    dmat_ran[onset:(onset+3), i] <- 1
  }
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Test using internal C++ functions directly
  result_cpp <- project_confounds_cpp(dmat_base, Y, dmat_ran)
  beta_cpp_direct <- lss_beta_cpp(result_cpp$Q_dmat_ran, result_cpp$residual_data)
  
  # Test using main lss function with C++
  bdes <- list(dmat_base = dmat_base, dmat_ran = dmat_ran, fixed_ind = NULL)
  beta_cpp_main <- lss(Y = Y, bdes = bdes, use_cpp = TRUE)
  
  expect_equal(beta_cpp_direct, beta_cpp_main, tolerance = 1e-12,
               info = "Direct C++ functions should match main lss() C++ path")
})

test_that("LSS handles numerical edge cases", {
  set.seed(131415)
  
  n_timepoints <- 60
  n_voxels <- 20
  n_trials <- 8
  
  # Create design with some nearly collinear trials
  dmat_base <- cbind(rep(1, n_timepoints), scale(1:n_timepoints)[,1])
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  
  # First few trials with some overlap
  dmat_ran[10:15, 1] <- 1
  dmat_ran[12:17, 2] <- 0.8  # Overlapping with trial 1
  dmat_ran[25:30, 3] <- 1
  dmat_ran[40:45, 4] <- 1
  
  # Remaining trials non-overlapping
  for (i in 5:n_trials) {
    onset <- 50 + (i-5)*2
    if (onset + 2 <= n_timepoints) {
      dmat_ran[onset:(onset+2), i] <- 1
    }
  }
  
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  bdes <- list(dmat_base = dmat_base, dmat_ran = dmat_ran, fixed_ind = NULL)
  
  # All methods should handle this without error
  expect_no_error(beta_r <- lss(Y = Y, bdes = bdes, use_cpp = FALSE))
  expect_no_error(beta_cpp <- lss(Y = Y, bdes = bdes, use_cpp = TRUE))
  expect_no_error(beta_naive <- lss_naive(Y = Y, bdes = bdes))
  
  # Results should still be similar (allowing for larger tolerance due to numerical issues)
  expect_equal(beta_r, beta_cpp, tolerance = 1e-6)
  expect_equal(beta_r, beta_naive, tolerance = 1e-4)
})

test_that("Data validation works correctly", {
  n_timepoints <- 50
  n_trials <- 5
  
  dmat_base <- cbind(rep(1, n_timepoints))
  dmat_ran <- matrix(rnorm(n_timepoints * n_trials), n_timepoints, n_trials)
  bdes <- list(dmat_base = dmat_base, dmat_ran = dmat_ran, fixed_ind = NULL)
  
  # Test dimension mismatch
  Y_wrong <- matrix(rnorm(40 * 10), 40, 10)  # Wrong number of timepoints
  
  expect_error(lss(Y = Y_wrong, bdes = bdes), 
               "Y has .* timepoints, design has .*")
  expect_error(lss_naive(Y = Y_wrong, bdes = bdes),
               "Y must have .* rows to match design matrix dimensions")
  
  # Test data.frame conversion
  Y_df <- data.frame(matrix(rnorm(n_timepoints * 10), n_timepoints, 10))
  expect_warning(beta_df <- lss(Y = Y_df, bdes = bdes),
                 "Converting Y from data.frame to matrix")
  expect_equal(dim(beta_df), c(n_trials, 10))
}) 