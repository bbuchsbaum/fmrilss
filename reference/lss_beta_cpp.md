# Vectorized LSS Beta Computation Using C++

Fast C++ implementation of least squares separate (LSS) beta estimation
using vectorized matrix operations. Computes all trial betas in a single
pass without loops.

## Usage

``` r
lss_beta_cpp(C_projected, Y_projected)
```

## Arguments

- C_projected:

  Projected trial regressors (n x T) from project_confounds_cpp

- Y_projected:

  Projected data (n x V) from project_confounds_cpp

## Value

Beta matrix (T x V) with LSS estimates for each trial and voxel

## Details

This vectorized implementation computes all LSS betas simultaneously
using matrix algebra. It's significantly faster than per-trial loops and
automatically benefits from BLAS multithreading. The algorithm handles
numerical edge cases by applying finite lower bounds to degenerate
denominators.

For best performance on large datasets, ensure your R installation uses
optimized BLAS (like OpenBLAS or Intel MKL).

## Examples

``` r
# \donttest{
n <- 40; V <- 3; T <- 4
X_confounds <- cbind(1, seq_len(n))
Y_data <- matrix(rnorm(n * V), n, V)
C_trials <- matrix(rnorm(n * T), n, T)
result <- project_confounds_cpp(X_confounds, Y_data, C_trials)
betas <- lss_beta_cpp(result$Q_dmat_ran, result$residual_data)
# }
```
