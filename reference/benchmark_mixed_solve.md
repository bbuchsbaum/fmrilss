# Benchmark Mixed Model Implementations

Compare performance between the standard `mixed_solve` implementation
and the optimized `mixed_solve_optimized` version.

## Usage

``` r
benchmark_mixed_solve(X, Z, K = NULL, Y, n_reps = 5)
```

## Arguments

- X:

  Fixed effects design matrix

- Z:

  Random effects design matrix

- K:

  Kinship matrix (optional, defaults to identity)

- Y:

  Response matrix (n x V)

- n_reps:

  Number of repetitions for benchmarking

## Value

Data frame with timing results

## Examples

``` r
if (FALSE) { # \dontrun{
X <- matrix(rnorm(100 * 2), 100, 2)
Z <- matrix(rnorm(100 * 3), 100, 3)
Y <- matrix(rnorm(100 * 5), 100, 5)
benchmark_mixed_solve(X, Z, Y = Y, n_reps = 2)
} # }
```
