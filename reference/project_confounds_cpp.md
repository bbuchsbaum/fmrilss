# Project Out Confounds Using C++

Fast C++ implementation for projecting out confound variables from data
and trial design matrices. This uses Cholesky decomposition for
numerical stability and avoids creating large projection matrices.

## Usage

``` r
project_confounds_cpp(X_confounds, Y_data, C_trials)
```

## Arguments

- X_confounds:

  Confound design matrix (n x k)

- Y_data:

  Data matrix (n x V) where V is number of voxels

- C_trials:

  Trial design matrix (n x T) where T is number of trials

## Value

List with projected data (residual_data) and projected trials
(Q_dmat_ran)

## Details

This function computes residuals Y - X(X'X)^(-1)X'Y and C -
X(X'X)^(-1)X'C without explicitly forming the projection matrix Q = I -
X(X'X)^(-1)X'. This approach uses ~100x less memory for large n and is
numerically more stable.

## Examples

``` r
# \donttest{
n <- 100; V <- 50; T <- 10
X_confounds <- cbind(1, seq_len(n), matrix(rnorm(n * 3), n, 3))
Y_data <- matrix(rnorm(n*V), n, V)
C_trials <- matrix(rnorm(n*T), n, T)

result <- project_confounds_cpp(X_confounds, Y_data, C_trials)
# }
```
