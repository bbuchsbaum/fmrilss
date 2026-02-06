# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

`fmrilss` is an R package implementing efficient Least Squares Separate (LSS) analysis for functional magnetic resonance imaging (fMRI) data. LSS is used to estimate trial-by-trial activation patterns in event-related fMRI designs, critical for multivariate pattern analysis (MVPA) and connectivity studies. The package provides multiple backends, from R implementations to highly optimized C++ with OpenMP parallelization.

## Core Architecture

### LSS Implementations

The package provides multiple implementations of the LSS algorithm, accessible through a unified interface:
- **r_optimized**: Optimized R implementation (default, recommended)
- **cpp_optimized**: Parallelized C++ using OpenMP (fastest for large datasets)
- **r_vectorized**: Standard vectorized R
- **cpp**: Standard C++ implementation
- **naive**: Simple loop-based R (reference implementation for testing)
- **oasis**: Mathematically equivalent reformulation with ridge regularization support

Key design pattern: All methods use the same `lss(Y, X, Z, Nuisance)` signature where:
- `Y`: data matrix (timepoints x voxels)
- `X`: trial design matrix (one column per trial)
- `Z`: experimental regressors (intercept, trends, blocks)
- `Nuisance`: regressors to project out (motion, physiology)

### Advanced Features

**OASIS Method**: Reformulates LSS as a single matrix operation, eliminating per-trial GLM redundancy. Supports:
- Multi-basis HRF models (K > 1)
- Ridge regularization (absolute or fractional modes)
- Automatic design construction from event onsets via `fmrihrf`
- Standard error computation
- AR(1) prewhitening

**Voxel-wise HRF Estimation**: `estimate_voxel_hrf()` fits per-voxel HRF basis coefficients, enabling data-driven HRF modeling before LSS.

**Shared-Basis HRF Models (SBHM)**: Learn low-rank shared time bases from parameterized HRF libraries via SVD, enabling efficient multi-voxel HRF estimation with reduced parameters.

**Prewhitening**: Integration with `fmriAR` package for AR/ARMA noise modeling with flexible pooling strategies (global, voxel-wise, run-aware, parcel-based).

## Essential Commands

```r
# Build and install package
devtools::build()
devtools::install()

# Run tests
devtools::test()                    # Run all tests
devtools::test_active_file()       # Test current file
testthat::test_file("tests/testthat/test-lss-equivalence.R")  # Single test file

# Check package
devtools::check()                  # Full R CMD check
devtools::check_examples()         # Check examples only

# Documentation
devtools::document()               # Update roxygen documentation
pkgdown::build_site()             # Build package website

# Load for development
devtools::load_all()              # Load package without installing
```

## C++ Compilation

The package uses Rcpp with Armadillo for matrix operations. OpenMP is configured via `src/Makevars`:
- Automatic OpenMP detection and configuration
- Falls back gracefully if OpenMP unavailable
- Links: Rcpp, RcppArmadillo, roptim, bigmemory, BH

## Testing Strategy

Tests verify:
1. Mathematical equivalence across all implementations
2. Correctness against ground truth
3. Edge cases (missing data, singular matrices)
4. HRF convolution integration
5. OASIS deconvolution methods

Test files follow pattern `test-{feature}.R` in `tests/testthat/`.

## Key Dependencies

- **fmrihrf**: HRF modeling and convolution (GitHub: bbuchsbaum/fmrihrf)
- **fmriAR**: AR/ARMA prewhitening (GitHub: bbuchsbaum/fmriAR)
- **RcppArmadillo**: Matrix operations in C++
- **roptim**: Optimization routines for HRF parameter estimation
- **MASS**: Statistical functions (used in mixed models)
- **bigmemory**: Memory-efficient matrix operations

## Vignettes

Five main vignettes demonstrate package usage:
- `fmrilss.Rmd`: Basic LSS concepts and usage (pkgdown "Get started" entry)
- `oasis_method.Rmd`: OASIS method with ridge regularization
- `oasis_theory.Rmd`: Mathematical foundations of OASIS
- `voxel-wise-hrf.Rmd`: Voxel-wise HRF estimation and integration with LSS
- `sbhm.Rmd`: Shared-Basis HRF Matching for efficient voxel-specific HRF estimation