# Mixed Model Solver using C++

C++ implementation of the mixed model solver. This function is typically
called through the main `mixed_solve` function rather than directly.

## Usage

``` r
.mixed_solve_cpp(
  Y,
  X = NULL,
  Z = NULL,
  K = NULL,
  method = "REML",
  bounds = c(1e-09, 1e+09),
  SE = FALSE,
  return_Hinv = FALSE
)
```

## Arguments

- Y:

  Response vector.

- X:

  Design matrix for fixed effects (default: intercept only).

- Z:

  Design matrix for random effects (default: identity matrix).

- K:

  Kinship matrix (default: identity matrix).

- method:

  Optimization method, either "REML" or "ML".

- bounds:

  Bounds for the optimizer.

- SE:

  Logical, whether to return standard errors.

- return_Hinv:

  Logical, whether to return the inverse of H.

## Value

A list with mixed model results.
