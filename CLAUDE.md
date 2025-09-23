# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

`fmrilss` is an R package implementing efficient Least Squares Separate (LSS) analysis for functional magnetic resonance imaging (fMRI) data. LSS is used to estimate trial-by-trial activation patterns in event-related fMRI designs, critical for multivariate pattern analysis (MVPA) and connectivity studies. The package provides multiple backends, from R implementations to highly optimized C++ with OpenMP parallelization.

## Core Architecture

The package provides five implementations of the LSS algorithm, accessible through a unified interface:
- **cpp_optimized**: Parallelized C++ using OpenMP (fastest)
- **optimized**: Optimized R version
- **vectorized**: Standard vectorized R
- **cpp_vectorized**: Standard C++
- **loop**: Simple R loop (reference implementation)

Key design pattern: All methods use the same `lss(Y, X, Z, Nuisance)` signature where:
- `Y`: data matrix (timepoints x voxels)
- `X`: trial design matrix (one column per trial)
- `Z`: experimental regressors (intercept, trends, blocks)
- `Nuisance`: regressors to project out (motion, physiology)

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

- **fmrihrf**: Custom HRF modeling (GitHub: bbuchsbaum/fmrihrf)
- **RcppArmadillo**: Matrix operations in C++
- **roptim**: Optimization routines
- **MASS**: Statistical functions

## Vignettes

Three main vignettes demonstrate package usage:
- `getting_started.Rmd`: Basic LSS concepts and usage
- `oasis_method.Rmd`: OASIS deconvolution integration
- `performance_optimization.Rmd`: Backend comparison and optimization