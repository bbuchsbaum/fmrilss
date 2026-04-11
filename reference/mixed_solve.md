# Mixed Model Solver

Solves mixed models with random effects using REML or ML estimation.
This function provides a unified interface to mixed model estimation,
similar to the lss/lsa functions in this package.

## Usage

``` r
mixed_solve(
  Y,
  X = NULL,
  Z = NULL,
  K = NULL,
  Nuisance = NULL,
  method = c("REML", "ML"),
  bounds = c(1e-09, 1e+09),
  SE = FALSE,
  return_Hinv = FALSE
)

mixed_solve_cpp(
  Y,
  X = NULL,
  Z = NULL,
  K = NULL,
  Nuisance = NULL,
  method = c("REML", "ML"),
  bounds = c(1e-09, 1e+09),
  SE = FALSE,
  return_Hinv = FALSE
)
```

## Arguments

- Y:

  Response vector or matrix. If a matrix, each column is treated as a
  separate response variable.

- X:

  Design matrix for fixed effects. If NULL, defaults to intercept only.

- Z:

  Design matrix for random effects. If NULL, defaults to identity
  matrix.

- K:

  Kinship matrix for random effects. If NULL, defaults to identity
  matrix.

- Nuisance:

  An alias for X, provided for consistency with lss/lsa interface. If
  both X and Nuisance are provided, X takes precedence.

- method:

  Character string specifying the estimation method:

  - "REML" - Restricted Maximum Likelihood (default)

  - "ML" - Maximum Likelihood

- bounds:

  Numeric vector of length 2 specifying bounds for variance component
  optimization. Defaults to c(1e-9, 1e9).

- SE:

  Logical, whether to compute and return standard errors. Defaults to
  FALSE.

- return_Hinv:

  Logical, whether to return the inverse of the H matrix. Defaults to
  FALSE.

## Value

A list containing:

- Vu:

  Estimated variance component for random effects.

- Ve:

  Estimated variance component for residuals.

- beta:

  Estimated fixed effects coefficients.

- u:

  Estimated random effects coefficients.

- LL:

  Log-likelihood of the model.

- beta.SE:

  Standard errors of fixed effects coefficients (if SE = TRUE).

- u.SE:

  Standard errors of random effects coefficients (if SE = TRUE).

- Hinv:

  Inverse of H matrix (if return_Hinv = TRUE).

## Details

This function fits the mixed model: Y = X*beta + Z*u + error, where u ~
N(0, Vu*K) and error ~ N(0, Ve*I). The variance components Vu and Ve are
estimated using REML or ML.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(123)
n <- 100
Y <- rnorm(n)
Z <- matrix(rnorm(n * 5), n, 5)
K <- diag(5)
X <- matrix(1, n, 1)

result <- mixed_solve(Y, X, Z, K)
} # }
```
