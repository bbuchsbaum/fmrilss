# Examples of using the new prewhitening features in fmrilss with fmriAR

library(fmrilss)

# Generate example data
set.seed(123)
n_time <- 200
n_voxels <- 50
n_trials <- 12

# Create AR(1) data with temporal autocorrelation
phi <- 0.5
Y <- matrix(0, n_time, n_voxels)
for (v in 1:n_voxels) {
  e <- rnorm(n_time)
  Y[1, v] <- e[1]
  for (t in 2:n_time) {
    Y[t, v] <- phi * Y[t-1, v] + e[t]
  }
}

# Create trial design
X <- matrix(0, n_time, n_trials)
for (i in 1:n_trials) {
  start <- (i-1) * 15 + 5
  if (start + 5 <= n_time) {
    X[start:(start+5), i] <- 1
  }
}

# Add some signal
true_betas <- matrix(rnorm(n_trials * n_voxels, 2, 1), n_trials, n_voxels)
Y <- Y + X %*% true_betas

# ============================================================================
# Basic Prewhitening Examples
# ============================================================================

# Example 1: Simple AR(1) prewhitening
cat("Example 1: AR(1) prewhitening\n")
beta_ar1 <- lss(Y, X, prewhiten = list(method = "ar", p = 1))
cat("  Dimensions:", dim(beta_ar1), "\n\n")

# Example 2: Automatic AR order selection
cat("Example 2: Auto AR order selection\n")
beta_auto <- lss(Y, X, prewhiten = list(
  method = "ar",
  p = "auto",
  p_max = 6  # Maximum order to consider
))
cat("  Dimensions:", dim(beta_auto), "\n\n")

# Example 3: No prewhitening (for comparison)
cat("Example 3: No prewhitening\n")
beta_none <- lss(Y, X, prewhiten = list(method = "none"))
# Or equivalently: beta_none <- lss(Y, X)
cat("  Dimensions:", dim(beta_none), "\n\n")

# ============================================================================
# Advanced Pooling Strategies
# ============================================================================

# Example 4: Voxel-specific AR parameters
cat("Example 4: Voxel-specific AR parameters\n")
beta_voxel <- lss(Y, X, prewhiten = list(
  method = "ar",
  p = 1,
  pooling = "voxel"  # Separate AR(1) for each voxel
))
cat("  Dimensions:", dim(beta_voxel), "\n\n")

# Example 5: Run-aware estimation for multi-run data
cat("Example 5: Run-aware AR estimation\n")
runs <- rep(1:2, each = n_time/2)
beta_run <- lss(Y, X, prewhiten = list(
  method = "ar",
  p = "auto",
  pooling = "run",
  runs = runs
))
cat("  Dimensions:", dim(beta_run), "\n\n")

# Example 6: Parcel-based pooling
cat("Example 6: Parcel-based pooling\n")
# Create example parcels (e.g., from clustering)
parcels <- sample(1:5, n_voxels, replace = TRUE)
beta_parcel <- lss(Y, X, prewhiten = list(
  method = "ar",
  p = 1,
  pooling = "parcel",
  parcels = parcels
))
cat("  Dimensions:", dim(beta_parcel), "\n\n")

# ============================================================================
# ARMA Models for Complex Noise
# ============================================================================

# Example 7: ARMA(2,1) model
cat("Example 7: ARMA(2,1) model\n")
beta_arma <- lss(Y, X, prewhiten = list(
  method = "arma",
  p = 2,  # AR order
  q = 1   # MA order
))
cat("  Dimensions:", dim(beta_arma), "\n\n")

# ============================================================================
# Integration with Different LSS Methods
# ============================================================================

# Example 8: Prewhitening with OASIS method
cat("Example 8: OASIS with prewhitening\n")
beta_oasis <- lss(Y, X, method = "oasis",
                  prewhiten = list(method = "ar", p = 1))
cat("  Dimensions:", dim(beta_oasis), "\n\n")

# Example 9: Prewhitening with C++ optimized method
cat("Example 9: C++ optimized with prewhitening\n")
beta_cpp <- lss(Y, X, method = "cpp_optimized",
                prewhiten = list(method = "ar", p = "auto"))
cat("  Dimensions:", dim(beta_cpp), "\n\n")

# ============================================================================
# With Additional Regressors
# ============================================================================

# Example 10: Prewhitening with experimental and nuisance regressors
cat("Example 10: Full model with prewhitening\n")
Z <- cbind(1, scale(1:n_time))  # Intercept + linear trend
Nuisance <- matrix(rnorm(n_time * 6), n_time, 6)  # Motion parameters

beta_full <- lss(Y, X, Z = Z, Nuisance = Nuisance,
                prewhiten = list(
                  method = "ar",
                  p = 1,
                  compute_residuals = TRUE  # Compute residuals after nuisance regression
                ))
cat("  Dimensions:", dim(beta_full), "\n\n")

# ============================================================================
# Backward Compatibility
# ============================================================================

# Example 11: Old syntax still works (with deprecation warning)
cat("Example 11: Backward compatibility (old syntax)\n")
suppressMessages({
  beta_old <- lss(Y, X, method = "oasis",
                 oasis = list(whiten = "ar1"))
})
cat("  Dimensions:", dim(beta_old), "\n")
cat("  Note: This syntax is deprecated. Use prewhiten parameter instead.\n\n")

# ============================================================================
# Performance Comparison
# ============================================================================

cat("Performance comparison:\n")
cat("  Correlation between no-whitening and AR(1):\n")
cat("    Mean:", mean(diag(cor(t(beta_none), t(beta_ar1)))), "\n")
cat("  Correlation between AR(1) and auto-selected:\n")
cat("    Mean:", mean(diag(cor(t(beta_ar1), t(beta_auto)))), "\n")
cat("  Correlation between global and voxel-specific:\n")
cat("    Mean:", mean(diag(cor(t(beta_ar1), t(beta_voxel)))), "\n")