## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----simulate-----------------------------------------------------------------
library(fmrilss)
set.seed(123)

# Create simple event design
n_time <- 100
n_trials <- 4
n_vox <- 3

# Event onsets
events <- data.frame(
  onset = c(10, 30, 50, 70),
  duration = 1,
  condition = "A"
)

# Create a simple HRF kernel (single basis for simplicity)
hrf_kernel <- c(0, 0.1, 0.3, 0.6, 0.9, 1.0, 0.8, 0.5, 0.3, 0.1, 0.05, 0)

# Generate data with different HRF scaling per voxel
Y <- matrix(0, n_time, n_vox)
true_betas <- matrix(rnorm(n_trials * n_vox, mean = 1, sd = 0.5), n_trials, n_vox)

for (trial in 1:n_trials) {
  # Create impulse for this trial
  impulse <- rep(0, n_time)
  onset_idx <- events$onset[trial]
  impulse[onset_idx:(onset_idx + events$duration[trial])] <- 1
  
  # Convolve with HRF
  conv_signal <- stats::convolve(impulse, rev(hrf_kernel), type = "open")[1:n_time]
  
  # Add to Y with trial-specific beta and voxel-specific scaling
  for (v in 1:n_vox) {
    # Simulate voxel-specific HRF scaling
    voxel_scale <- 0.5 + v * 0.3  # Different HRF amplitude per voxel
    Y[, v] <- Y[, v] + conv_signal * true_betas[trial, v] * voxel_scale
  }
}

# Add noise
Y <- Y + matrix(rnorm(n_time * n_vox, sd = 0.2), n_time, n_vox)

## ----standard-lss-------------------------------------------------------------
# Create design matrix with convolved events
X <- matrix(0, n_time, n_trials)
for (trial in 1:n_trials) {
  impulse <- rep(0, n_time)
  onset_idx <- events$onset[trial]
  impulse[onset_idx:(onset_idx + events$duration[trial])] <- 1
  X[, trial] <- stats::convolve(impulse, rev(hrf_kernel), type = "open")[1:n_time]
}

# Run standard LSS
standard_betas <- lss(Y, X, method = "r_optimized")
print(round(standard_betas, 2))

## ----voxel-hrf-lss------------------------------------------------------------
# Simulate estimated HRF coefficients (in real analysis, these would come from estimate_voxel_hrf)
# Here we use the known voxel scalings with some estimation error
coefficients <- matrix(c(0.8, 1.1, 1.4) + rnorm(3, sd = 0.1), nrow = 1, ncol = n_vox)

# Create HRF estimates object
hrf_est <- list(
  coefficients = coefficients,
  basis = structure(list(), class = "HRF"),  # Dummy HRF object
  conditions = "A"
)
class(hrf_est) <- "VoxelHRF"

# Convert events to onset indices
onset_idx <- as.integer(events$onset)
durations <- as.integer(events$duration)

# Use the pure R implementation with voxel-wise HRFs
hrf_basis_kernels <- matrix(hrf_kernel, ncol = 1)

voxel_betas <- fmrilss:::lss_with_hrf_pure_r(
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

print(round(voxel_betas, 2))

## ----compare------------------------------------------------------------------
# Calculate correlation between estimated and true betas
cor_standard <- cor(as.vector(standard_betas), as.vector(true_betas))
cor_voxel <- cor(as.vector(voxel_betas), as.vector(true_betas))

cat("Correlation with true betas:\n")
cat("  Standard LSS:", round(cor_standard, 3), "\n")
cat("  Voxel-wise HRF LSS:", round(cor_voxel, 3), "\n")

# Plot comparison
par(mfrow = c(1, 2))
plot(true_betas, standard_betas, 
     xlab = "True Betas", ylab = "Estimated Betas",
     main = paste("Standard LSS (r =", round(cor_standard, 2), ")"),
     pch = 19, col = rep(1:n_vox, each = n_trials))
abline(0, 1, lty = 2)

plot(true_betas, voxel_betas,
     xlab = "True Betas", ylab = "Estimated Betas", 
     main = paste("Voxel-wise HRF LSS (r =", round(cor_voxel, 2), ")"),
     pch = 19, col = rep(1:n_vox, each = n_trials))
abline(0, 1, lty = 2)

