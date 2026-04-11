# Optimized LSS with modern signature

This is a convenience wrapper around
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) that
selects one of the optimized implementations.

## Usage

``` r
lss_optimized_fit(
  Y,
  X,
  Z = NULL,
  Nuisance = NULL,
  engine = c("cpp", "r"),
  block_size = 96,
  prewhiten = NULL
)
```

## Arguments

- Y:

  Numeric matrix (timepoints x voxels).

- X:

  Trial design matrix (timepoints x trials).

- Z:

  Optional experimental regressors.

- Nuisance:

  Optional nuisance regressors to project out.

- engine:

  `"cpp"` (default) for `method="cpp_optimized"` or `"r"` for
  `method="r_optimized"`.

- block_size:

  Block size used by the C++ optimized path.

- prewhiten:

  Optional prewhitening options list (see
  [`prewhiten_options()`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)).

## Value

A numeric matrix (trials x voxels) of beta estimates.

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(16), 8, 2)
X <- matrix(0, 8, 2)
X[2:3, 1] <- 1
X[5:6, 2] <- 1
lss_optimized_fit(Y, X, engine = "r")
#>            Voxel_1     Voxel_2
#> Trial_1 -0.8746378  0.09179092
#> Trial_2 -0.7941255 -1.92937571
```
