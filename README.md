# fmrilss

Least Squares Separate (LSS) Analysis for fMRI Data

## Overview

The `fmrilss` package implements efficient least squares separate (LSS) analysis for functional magnetic resonance imaging (fMRI) data. LSS is used to estimate trial-by-trial activation patterns in event-related fMRI designs.

## Features

- Fast C++ implementation using RcppArmadillo
- Fallback R implementation for compatibility  
- Support for various design matrix configurations
- Robust numerical handling for edge cases

## Installation

```r
# Install development version from GitHub
# devtools::install_github("your-username/fmrilss")

# Or install from local source
# install.packages(".", repos = NULL, type = "source")
```

## Usage

```r
library(fmrilss)

# Create example design matrices
n_timepoints <- 100
n_trials <- 20
n_voxels <- 1000

# Base design (intercept + linear trend)
dmat_base <- cbind(1, 1:n_timepoints)

# Trial design matrix (one column per trial)
dmat_ran <- matrix(0, n_timepoints, n_trials)
for(i in 1:n_trials) {
  trial_onset <- sample(1:(n_timepoints-10), 1)
  dmat_ran[trial_onset:(trial_onset+5), i] <- 1
}

# Create design list
bdes <- list(
  dmat_base = dmat_base,
  dmat_ran = dmat_ran,
  fixed_ind = NULL
)

# Simulate data
Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

# Run LSS analysis
beta_estimates <- lss(Y = Y, bdes = bdes)
```

## References

Mumford, J.A., et al. (2012). Deconvolving BOLD activation in event-related designs for multivoxel pattern classification analyses. NeuroImage, 59(3), 2636-2643.

## License

GPL-3 