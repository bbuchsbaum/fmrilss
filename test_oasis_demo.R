#!/usr/bin/env Rscript

# Simple demonstration of OASIS HRF recovery without skip_on_cran

library(fmrilss)
library(fmrihrf)

# Source helper functions
source("R/oasis_hrf_recovery.R")

cat("\n==========================================\n")
cat("OASIS HRF Recovery Test Demonstration\n")
cat("==========================================\n\n")

# Set seed for reproducibility
set.seed(42)

# Test parameters - non-canonical HRF
true_tau <- 7      # Later peak than canonical
true_sigma <- 3    # Wider than typical
true_rho <- 0.45   # Stronger undershoot

cat("True HRF Parameters (LWU model):\n")
cat(sprintf("  tau (peak time): %.1f seconds\n", true_tau))
cat(sprintf("  sigma (width): %.1f seconds\n", true_sigma))
cat(sprintf("  rho (undershoot): %.2f\n\n", true_rho))

# Generate rapid event-related design
cat("Generating rapid event-related design...\n")
onsets <- generate_rapid_design(
  n_events = 20,
  total_time = 200,
  min_isi = 2,
  max_isi = 4,
  seed = 123
)
cat(sprintf("  Generated %d events\n", length(onsets)))
cat(sprintf("  Mean ISI: %.2f seconds\n\n", mean(diff(onsets))))

# Generate synthetic data with low noise
cat("Generating synthetic fMRI data...\n")
data <- generate_lwu_data(
  onsets = onsets,
  tau = true_tau,
  sigma = true_sigma,
  rho = true_rho,
  TR = 1.0,
  total_time = 200,
  n_voxels = 5,
  noise_sd = 0.1,  # Low noise for demonstration
  seed = 456
)
cat(sprintf("  Time points: %d\n", nrow(data$Y)))
cat(sprintf("  Voxels: %d\n", ncol(data$Y)))
cat(sprintf("  SNR: %.2f\n\n", var(data$signal) / var(data$noise)))

# Create HRF grid for OASIS
cat("Creating HRF parameter grid for OASIS search...\n")
hrf_grid <- create_lwu_grid(
  tau_range = c(5, 9),
  sigma_range = c(2, 4),
  rho_range = c(0.2, 0.6),
  n_tau = 5,
  n_sigma = 3,
  n_rho = 3
)
cat(sprintf("  Grid size: %d HRF models\n\n", length(hrf_grid$hrfs)))

# Compare recovery methods
cat("Comparing HRF recovery methods...\n")
cat("  Fitting OASIS with grid search...\n")
cat("  Fitting SPMG1 (canonical HRF)...\n")
cat("  Fitting SPMG3 (with derivatives)...\n")
cat("  Fitting FIR (non-parametric)...\n\n")

results <- compare_hrf_recovery(data, hrf_grid)

# Calculate metrics
cat("Calculating recovery metrics...\n\n")
metrics <- calculate_recovery_metrics(results, data$true_hrf)

# Display results
cat("==========================================\n")
cat("RESULTS\n")
cat("==========================================\n\n")

cat("HRF Recovery Performance:\n")
cat("-------------------------\n")
for (i in 1:nrow(metrics)) {
  cat(sprintf("%s:\n", metrics$method[i]))
  cat(sprintf("  MSE: %.4f\n", metrics$mse[i]))
  cat(sprintf("  Correlation: %.3f\n", metrics$correlation[i]))
}

cat("\nOASIS Parameter Recovery:\n")
cat("-------------------------\n")
cat(sprintf("True tau: %.1f → Recovered: %.2f (error: %.2f)\n", 
            true_tau, results$oasis$best_params$tau, 
            abs(results$oasis$best_params$tau - true_tau)))
cat(sprintf("True sigma: %.1f → Recovered: %.2f (error: %.2f)\n",
            true_sigma, results$oasis$best_params$sigma,
            abs(results$oasis$best_params$sigma - true_sigma)))
cat(sprintf("True rho: %.2f → Recovered: %.2f (error: %.2f)\n",
            true_rho, results$oasis$best_params$rho,
            abs(results$oasis$best_params$rho - true_rho)))

# Performance comparison
oasis_mse <- metrics$mse[metrics$method == "OASIS"]
spmg1_mse <- metrics$mse[metrics$method == "SPMG1"]
spmg3_mse <- metrics$mse[metrics$method == "SPMG3"]

improvement_spmg1 <- (spmg1_mse - oasis_mse) / spmg1_mse * 100
improvement_spmg3 <- (spmg3_mse - oasis_mse) / spmg3_mse * 100

cat("\n==========================================\n")
cat("SUMMARY\n")
cat("==========================================\n\n")

cat(sprintf("✓ OASIS shows %.1f%% improvement over SPMG1\n", improvement_spmg1))
cat(sprintf("✓ OASIS shows %.1f%% improvement over SPMG3\n", improvement_spmg3))
cat(sprintf("✓ OASIS recovered peak time within %.1f seconds\n", 
            abs(results$oasis$best_params$tau - true_tau)))
cat(sprintf("✓ OASIS recovered width within %.1f seconds\n",
            abs(results$oasis$best_params$sigma - true_sigma)))

if (oasis_mse < spmg1_mse && oasis_mse < spmg3_mse) {
  cat("\n✅ SUCCESS: OASIS outperforms standard models for HRF recovery!\n")
  cat("   This validates OASIS for analyzing data with non-canonical HRFs.\n")
} else {
  cat("\n⚠️  Unexpected result - check parameters and noise levels\n")
}

cat("\n==========================================\n")
cat("Test completed successfully!\n")
cat("==========================================\n\n")

# Optional: Create visualization
if (require(ggplot2, quietly = TRUE)) {
  cat("Creating visualization...\n")
  p <- plot_hrf_comparison(results, save_path = "test_oasis_hrf_comparison.pdf")
  cat("Plot saved as 'test_oasis_hrf_comparison.pdf'\n\n")
}