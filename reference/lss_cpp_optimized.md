# A wrapper for the optimized C++ LSS implementation

A wrapper for the optimized C++ LSS implementation

## Usage

``` r
lss_cpp_optimized(Y, bdes)
```

## Arguments

- Y:

  the voxel by time data matrix

- bdes:

  the block design list created by `block_design`

## Value

a matrix of beta estimates

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(16), 8, 2)
X_trials <- matrix(0, 8, 2)
X_trials[2:3, 1] <- 1
X_trials[5:6, 2] <- 1
bdes <- list(
  dmat_base = matrix(1, 8, 1),
  dmat_fixed = NULL,
  dmat_ran = X_trials
)
lss_cpp_optimized(Y, bdes)
#>            [,1]        [,2]
#> [1,] -0.8746378  0.09179092
#> [2,] -0.7941255 -1.92937571
```
