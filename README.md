# fmrilss

[![R-CMD-check](https://github.com/bbuchsbaum/fmrilss/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bbuchsbaum/fmrilss/actions/workflows/R-CMD-check.yaml)

Least Squares Separate (LSS) Analysis for fMRI Data.

## Overview

The `fmrilss` package provides an efficient and flexible implementation of the Least Squares Separate (LSS) method for fMRI analysis, as proposed by Mumford et al. (2012). This approach models each trial with a separate GLM, making it a powerful technique for multivariate pattern analysis (MVPA) and connectivity studies where trial-specific estimates are needed.

The package offers multiple backends, from a simple reference implementation to a highly optimized, parallel C++ engine, all accessible through a clean, unified interface.

## Features

- **Five Implementations**: Includes a highly optimized C++ version, an optimized R version, standard vectorized R and C++ versions, and a simple R loop for testing and validation.
- **Parallel Processing**: The optimized C++ version uses OpenMP for multi-threaded execution to maximize performance on modern hardware.
- **Flexible & Modern Interface**: A clean `lss(Y, X, Z, Nuisance)` signature that is powerful and intuitive.
- **Nuisance Regression**: Built-in support for projecting out nuisance regressors (e.g., motion parameters, physiological noise) before the LSS analysis.
- **CRAN-Compliant**: Built with portable configurations suitable for CRAN submission.

## Installation

You can install the development version of `fmrilss` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("bbuchsbaum/fmrilss")
```

## Usage

The primary function is `lss()`, which takes your data `Y`, trial design `X`, experimental regressors `Z`, and optional nuisance regressors.

```r
library(fmrilss)

# 1. Generate synthetic data
set.seed(123)
n_timepoints <- 120
n_trials <- 15
n_voxels <- 50

# Trial design matrix (X): one column per trial
X <- matrix(0, n_timepoints, n_trials)
onsets <- seq(from = 5, to = n_timepoints - 10, length.out = n_trials)
for(i in 1:n_trials) {
  X[onsets[i]:(onsets[i] + 4), i] <- 1
}
colnames(X) <- paste0("Trial_", 1:n_trials)

# Experimental regressors (Z): intercept and condition-specific effects
# These are experimental regressors we want to model and get beta estimates for,
# but not trial-wise (e.g., condition differences, block effects)
Z <- cbind(Intercept = 1, LinearTrend = scale(1:n_timepoints, center = TRUE, scale = FALSE))

# Nuisance regressors: e.g., 6 motion parameters
Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)

# Data (Y): timepoints x voxels
# (Simulate some effects for demonstration)
true_betas <- matrix(rnorm(n_trials * n_voxels, 0, 1.5), n_trials, n_voxels)
Y <- Z %*% matrix(c(10, -0.2), 2, n_voxels) + 
     X %*% true_betas +
     Nuisance %*% matrix(rnorm(6 * n_voxels, 0, 2), 6, n_voxels) +
     matrix(rnorm(n_timepoints * n_voxels, 0, 0.8), n_timepoints, n_voxels)
colnames(Y) <- paste0("Voxel_", 1:n_voxels)


# 2. Run LSS analysis

# Example 1: Basic LSS with default intercept
# If Z is NULL, an intercept is automatically added.
beta_estimates <- lss(Y, X)

# Example 2: LSS with experimental regressors (intercept + condition effects)
beta_fixed <- lss(Y, X, Z = Z)

# Example 3: LSS with experimental regressors and nuisance regression
beta_clean <- lss(Y, X, Z = Z, Nuisance = Nuisance)

# Example 4: Use the super-fast, parallelized C++ implementation
beta_fast <- lss(Y, X, Z = Z, Nuisance = Nuisance, method = "cpp_optimized")

# The result is a (trials x voxels) matrix of beta estimates
print(dim(beta_fast))
#> [1] 15 50
print(beta_fast[1:5, 1:4])
```

## References

Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012). Deconvolving BOLD activation in event-related designs for multivoxel pattern classification analyses. *NeuroImage*, 59(3), 2636-2643.

## License

GPL-3 