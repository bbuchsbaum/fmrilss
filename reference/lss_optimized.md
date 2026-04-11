# Optimized LSS Analysis (Pure R)

An optimized version of the LSS analysis that avoids creating large
intermediate matrices, providing a significant speedup and lower memory
usage for the pure R implementation.

## Usage

``` r
lss_optimized(Y = NULL, bdes, dset = NULL, use_cpp = TRUE)
```

## Arguments

- Y:

  A numeric matrix where rows are timepoints and columns are
  voxels/features.

- bdes:

  A list containing the design matrices.

- dset:

  Optional dataset object.

- use_cpp:

  Logical. If TRUE (default), uses the C++ implementation. If FALSE,
  uses the new optimized R implementation.

## Value

A numeric matrix of LSS beta estimates.

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(16), 8, 2)
X_trials <- matrix(0, 8, 2)
X_trials[2:3, 1] <- 1
X_trials[5:6, 2] <- 1
bdes <- list(
  dmat_base = matrix(1, 8, 1),
  dmat_ran = X_trials,
  dmat_fixed = NULL,
  fixed_ind = NULL
)
lss_optimized(Y, bdes, use_cpp = FALSE)
#>            [,1]        [,2]
#> [1,] -0.8746378  0.09179092
#> [2,] -0.7941255 -1.92937571
```
