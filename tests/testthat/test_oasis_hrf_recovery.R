# Test OASIS HRF Recovery
# Tests whether OASIS can recover non-canonical HRF shapes better than standard models

library(testthat)
library(fmrilss)

test_that("OASIS recovers LWU HRF parameters better than SPMG1 in rapid designs", {
  skip_on_cran()
  set.seed(42)
  
  # Source helper functions
  source(system.file("R", "oasis_hrf_recovery.R", package = "fmrilss"))
  
  # Test parameters
  true_tau <- 7      # Non-canonical peak time
  true_sigma <- 3    # Wider than typical
  true_rho <- 0.45   # Stronger undershoot
  
  # Generate rapid event-related design
  onsets <- generate_rapid_design(
    n_events = 20,
    total_time = 200,
    min_isi = 2,
    max_isi = 4,
    seed = 123
  )
  
  # Generate synthetic data with low noise
  data <- generate_lwu_data(
    onsets = onsets,
    tau = true_tau,
    sigma = true_sigma,
    rho = true_rho,
    TR = 1.0,
    total_time = 200,
    n_voxels = 5,
    noise_sd = 0.1,  # Low noise for initial test
    seed = 456
  )
  
  # Create HRF grid for OASIS
  hrf_grid <- create_lwu_grid(
    tau_range = c(5, 9),
    sigma_range = c(2, 4),
    rho_range = c(0.2, 0.6),
    n_tau = 5,
    n_sigma = 3,
    n_rho = 3
  )
  
  # Compare recovery methods
  results <- compare_hrf_recovery(data, hrf_grid)
  
  # Calculate metrics
  metrics <- calculate_recovery_metrics(results, data$true_hrf)
  
  # Tests
  # 1. OASIS should have lowest MSE
  oasis_mse <- metrics$mse[metrics$method == "OASIS"]
  spmg1_mse <- metrics$mse[metrics$method == "SPMG1"]
  expect_lt(oasis_mse, spmg1_mse, 
            label = "OASIS MSE should be lower than SPMG1")
  
  # 2. OASIS should have highest correlation
  oasis_cor <- metrics$correlation[metrics$method == "OASIS"]
  spmg1_cor <- metrics$correlation[metrics$method == "SPMG1"]
  expect_gt(oasis_cor, spmg1_cor,
            label = "OASIS correlation should be higher than SPMG1")
  
  # 3. OASIS should recover tau within 1.5 seconds
  oasis_tau_error <- metrics$peak_time_error[metrics$method == "OASIS"]
  expect_lt(oasis_tau_error, 1.5,
            label = "OASIS should recover tau within 1.5 seconds")
  
  # 4. OASIS parameters should be close to ground truth
  expect_lt(abs(results$oasis$best_params$tau - true_tau), 2,
            label = "Recovered tau should be within 2s of truth")
  expect_lte(abs(results$oasis$best_params$sigma - true_sigma), 1,
             label = "Recovered sigma should be within 1s of truth (<=)")
  expect_lt(abs(results$oasis$best_params$rho - true_rho), 0.2,
            label = "Recovered rho should be within 0.2 of truth")
})

test_that("OASIS performance degrades gracefully with noise", {
  skip_on_cran()
  set.seed(789)
  
  source(system.file("R", "oasis_hrf_recovery.R", package = "fmrilss"))
  
  # Test parameters
  true_tau <- 6.5
  true_sigma <- 2.8
  true_rho <- 0.4
  
  # Test with different noise levels
  noise_levels <- c(0.05, 0.2, 0.5, 1.0)
  correlations <- numeric(length(noise_levels))
  
  for (i in seq_along(noise_levels)) {
    # Generate design
    onsets <- generate_rapid_design(
      n_events = 25,
      total_time = 250,
      min_isi = 2,
      max_isi = 3.5,
      seed = 100 + i
    )
    
    # Generate data
    data <- generate_lwu_data(
      onsets = onsets,
      tau = true_tau,
      sigma = true_sigma,
      rho = true_rho,
      TR = 1.0,
      total_time = 250,
      n_voxels = 10,
      noise_sd = noise_levels[i],
      seed = 200 + i
    )
    
    # Fit OASIS
    hrf_grid <- create_lwu_grid(
      tau_range = c(5, 8),
      sigma_range = c(2, 3.5),
      rho_range = c(0.25, 0.55),
      n_tau = 4,
      n_sigma = 3,
      n_rho = 3
    )
    
    results <- compare_hrf_recovery(data, hrf_grid)
    metrics <- calculate_recovery_metrics(results, data$true_hrf)
    
    correlations[i] <- metrics$correlation[metrics$method == "OASIS"]
  }
  
  # Test that correlation decreases with noise
  expect_true(all(diff(correlations) <= 0.1),  # Allow small increases due to randomness
              label = "Correlation should generally decrease with noise")
  
  # Even with high noise, should maintain reasonable correlation
  expect_gt(correlations[length(correlations)], 0.5,
            label = "Should maintain correlation > 0.5 even with high noise")
})

test_that("OASIS outperforms SPMG models for overlapping HRFs", {
  skip_on_cran()
  set.seed(321)
  
  source(system.file("R", "oasis_hrf_recovery.R", package = "fmrilss"))
  
  # Very rapid design to ensure overlapping HRFs
  onsets <- generate_rapid_design(
    n_events = 30,
    total_time = 180,
    min_isi = 1.5,  # Very short ISI
    max_isi = 2.5,
    seed = 999
  )
  
  # Non-canonical HRF shape
  data <- generate_lwu_data(
    onsets = onsets,
    tau = 5,      # Earlier peak
    sigma = 3.5,  # Wider
    rho = 0.5,    # Strong undershoot
    TR = 1.0,
    total_time = 180,
    n_voxels = 8,
    noise_sd = 0.15,
    seed = 888
  )
  
  # Fit models
  hrf_grid <- create_lwu_grid(
    tau_range = c(4, 7),
    sigma_range = c(2.5, 4),
    rho_range = c(0.3, 0.6),
    n_tau = 4,
    n_sigma = 3,
    n_rho = 3
  )
  
  results <- compare_hrf_recovery(data, hrf_grid)
  metrics <- calculate_recovery_metrics(results, data$true_hrf)
  
  # OASIS should significantly outperform both SPMG models
  oasis_mse <- metrics$mse[metrics$method == "OASIS"]
  spmg1_mse <- metrics$mse[metrics$method == "SPMG1"]
  spmg3_mse <- metrics$mse[metrics$method == "SPMG3"]
  
  # OASIS should have at least 30% lower MSE than SPMG1
  improvement <- (spmg1_mse - oasis_mse) / spmg1_mse
  expect_gt(improvement, 0.3,
            label = "OASIS should show >30% MSE improvement over SPMG1")
  
  # OASIS should also beat SPMG3
  expect_lt(oasis_mse, spmg3_mse,
            label = "OASIS should have lower MSE than SPMG3")
  
  # Beta recovery should be better with OASIS
  if ("beta_correlation" %in% names(metrics)) {
    oasis_beta_cor <- metrics$beta_correlation[metrics$method == "OASIS"]
    spmg1_beta_cor <- metrics$beta_correlation[metrics$method == "SPMG1"]
    expect_gt(oasis_beta_cor, spmg1_beta_cor,
              label = "OASIS should have better beta recovery")
  }
})

test_that("OASIS grid search selects appropriate HRF parameters", {
  skip_on_cran()
  set.seed(555)
  
  source(system.file("R", "oasis_hrf_recovery.R", package = "fmrilss"))
  
  # Test with canonical-like parameters
  true_tau <- 6
  true_sigma <- 2.5
  true_rho <- 0.35
  
  onsets <- generate_rapid_design(
    n_events = 20,
    total_time = 200,
    min_isi = 2.5,
    max_isi = 4,
    seed = 777
  )
  
  data <- generate_lwu_data(
    onsets = onsets,
    tau = true_tau,
    sigma = true_sigma,
    rho = true_rho,
    TR = 1.0,
    total_time = 200,
    n_voxels = 10,
    noise_sd = 0.1,
    seed = 666
  )
  
  # Create grid that includes true parameters
  hrf_grid <- create_lwu_grid(
    tau_range = c(5, 7),
    sigma_range = c(2, 3),
    rho_range = c(0.25, 0.45),
    n_tau = 5,
    n_sigma = 5,
    n_rho = 5
  )
  
  # Fit OASIS
  oasis_result <- fit_oasis_grid(
    Y = data$Y,
    onsets = onsets,
    sframe = data$sframe,
    hrf_grid = hrf_grid,
    ridge_x = 0.01,
    ridge_b = 0.01
  )
  
  # Check that selected parameters are close to truth
  expect_lt(abs(oasis_result$best_params$tau - true_tau), 1,
            label = "Selected tau should be within 1s of truth")
  expect_lte(abs(oasis_result$best_params$sigma - true_sigma), 0.5,
             label = "Selected sigma should be within 0.5s of truth (<=)")
  expect_lt(abs(oasis_result$best_params$rho - true_rho), 0.15,
            label = "Selected rho should be within 0.15 of truth")
  
  # Check that best model has high score
  expect_gt(max(oasis_result$scores), 0.7,
            label = "Best model should have R-squared > 0.7")
})

# Run comprehensive simulation test
test_that("OASIS shows consistent improvement across parameter space", {
  skip_on_cran()
  skip_if_not(interactive(), "Skipping comprehensive test in non-interactive mode")
  
  set.seed(1234)
  source(system.file("R", "oasis_hrf_recovery.R", package = "fmrilss"))
  
  # Test multiple parameter combinations
  test_params <- expand.grid(
    tau = c(5, 6, 7),
    sigma = c(2, 2.5, 3),
    rho = c(0.3, 0.4, 0.5)
  )
  
  n_tests <- nrow(test_params)
  improvements <- numeric(n_tests)
  
  for (i in seq_len(min(n_tests, 5))) {  # Limit to 5 for speed
    params <- test_params[i, ]
    
    # Generate data
    onsets <- generate_rapid_design(
      n_events = 15,
      total_time = 150,
      min_isi = 2,
      max_isi = 3.5,
      seed = 1000 + i
    )
    
    data <- generate_lwu_data(
      onsets = onsets,
      tau = params$tau,
      sigma = params$sigma,
      rho = params$rho,
      TR = 1.0,
      total_time = 150,
      n_voxels = 5,
      noise_sd = 0.2,
      seed = 2000 + i
    )
    
    # Quick grid search
    hrf_grid <- create_lwu_grid(
      tau_range = c(params$tau - 2, params$tau + 2),
      sigma_range = c(params$sigma - 1, params$sigma + 1),
      rho_range = c(max(0.1, params$rho - 0.2), min(0.7, params$rho + 0.2)),
      n_tau = 3,
      n_sigma = 3,
      n_rho = 3
    )
    
    results <- compare_hrf_recovery(data, hrf_grid)
    metrics <- calculate_recovery_metrics(results, data$true_hrf)
    
    # Calculate improvement
    oasis_mse <- metrics$mse[metrics$method == "OASIS"]
    spmg1_mse <- metrics$mse[metrics$method == "SPMG1"]
    improvements[i] <- (spmg1_mse - oasis_mse) / spmg1_mse
  }
  
  # OASIS should show improvement in most cases
  expect_gt(mean(improvements[1:min(n_tests, 5)]), 0.2,
            label = "Average improvement should be > 20%")
  expect_gt(sum(improvements[1:min(n_tests, 5)] > 0), min(n_tests, 5) * 0.8,
            label = "OASIS should improve in >80% of cases")
})
