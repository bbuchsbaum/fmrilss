# Naive LSS with modern signature

Equivalent to calling `lss(..., method = "naive")`.

## Usage

``` r
lss_naive_fit(Y, X, Z = NULL, Nuisance = NULL, prewhiten = NULL)
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
lss_naive_fit(Y, X)
#>            Voxel_1     Voxel_2
#> Trial_1 -0.8746378  0.09179092
#> Trial_2 -0.7941255 -1.92937571
```
