## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(fmrilss)

set.seed(42)
n_timepoints <- 150
n_trials <- 12
n_voxels <- 25

## ----design-------------------------------------------------------------------
# Trial design matrix (X)
X <- matrix(0, n_timepoints, n_trials)
onsets <- seq(from = 10, to = n_timepoints - 12, length.out = n_trials)
for(i in 1:n_trials) {
  X[onsets[i]:(onsets[i] + 5), i] <- 1
}

# Experimental regressors (Z) - intercept and condition effects
# These are regressors we want to model and get estimates for, but not trial-wise
Z <- cbind(Intercept = 1, LinearTrend = scale(1:n_timepoints, center = TRUE, scale = FALSE))

# Nuisance regressors - e.g., 6 motion parameters
Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)

## ----data---------------------------------------------------------------------
# Simulate effects for each component
true_trial_betas <- matrix(rnorm(n_trials * n_voxels, 0, 1.2), n_trials, n_voxels)
true_fixed_effects <- matrix(rnorm(2 * n_voxels, c(5, -0.1), 0.5), 2, n_voxels)
true_nuisance_effects <- matrix(rnorm(6 * n_voxels, 0, 2), 6, n_voxels)

# Combine signals and add noise
Y <- (Z %*% true_fixed_effects) + 
     (X %*% true_trial_betas) +
     (Nuisance %*% true_nuisance_effects) +
     matrix(rnorm(n_timepoints * n_voxels, 0, 1), n_timepoints, n_voxels)

## ----lss_basic----------------------------------------------------------------
beta_basic <- lss(Y, X)
# The result is a trials-by-voxels matrix
dim(beta_basic)

## ----lss_advanced-------------------------------------------------------------
beta_full <- lss(Y, X, Z = Z, Nuisance = Nuisance)
# The output dimensions remain the same
dim(beta_full)

## ----lss_cpp, eval=require("Rcpp") && require("RcppArmadillo")----------------
# Run the same analysis with the high-performance C++ engine
beta_fast <- lss(Y, X, Z = Z, Nuisance = Nuisance, method = "cpp_optimized")

# The results are numerically identical to the R version
all.equal(beta_full, beta_fast, tolerance = 1e-8)

