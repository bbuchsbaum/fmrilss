# Optimized Mixed Model Solver

An optimized implementation of mixed model estimation that precomputes
expensive matrix operations and can be reused across multiple voxels for
significant performance improvements.

## Usage

``` r
mixed_solve_optimized(
  X,
  Z,
  Y,
  K = NULL,
  workspace = NULL,
  compute_se = FALSE,
  n_threads = 0
)
```

## Arguments

- X:

  Fixed effects design matrix (n × p)

- Z:

  Random effects design matrix (n × q)

- Y:

  Response data - can be a vector (single voxel) or matrix (n × V for
  multiple voxels)

- K:

  Kinship/covariance matrix for random effects (q × q). Defaults to
  identity.

- workspace:

  Precomputed workspace (optional, will compute if NULL)

- compute_se:

  Whether to compute standard errors (default: FALSE)

- n_threads:

  Number of OpenMP threads for multi-voxel (0 = auto)

## Value

List with estimated parameters and variance components

## Examples

``` r
set.seed(1)
X <- matrix(1, 6, 1)
Z <- diag(6)
Y <- matrix(rnorm(12), 6, 2)
ws <- mixed_precompute(X, Z)
#> Workspace precomputed successfully:
#>   - n=6, p=1, q=6
#>   - Effective rank: 5
#>   - Using identity K: yes
fit <- mixed_solve_optimized(X, Z, Y, workspace = ws)
#> Using 4 OpenMP threads
names(fit)
#> [1] "beta"   "u"      "Vu"     "Ve"     "lambda"
```
