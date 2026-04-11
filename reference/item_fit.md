# Fit ITEM decoder weights

Fit ITEM weights with a ridge-stabilized generalized least-squares
solve.

## Usage

``` r
item_fit(
  Gamma_train,
  T_train,
  U_train,
  ridge = 0,
  method = c("chol", "svd", "pinv"),
  tol = sqrt(.Machine$double.eps)
)
```

## Arguments

- Gamma_train:

  Numeric matrix (`n_train x n_features`).

- T_train:

  Numeric target matrix (`n_train x p`).

- U_train:

  Trial covariance (`n_train x n_train`) or run-block list.

- ridge:

  Non-negative ridge term added to the left-hand system.

- method:

  Preferred solver path (`"chol"`, `"svd"`, `"pinv"`).

- tol:

  Numerical tolerance for rank/solver fallbacks.

## Value

Numeric weight matrix `W_hat` (`n_features x p`) with `item_diagnostics`
attribute.

## Examples

``` r
Gamma_train <- matrix(
  c(1, 0,
    0.9, 0.1,
    0.1, 0.9),
  ncol = 2,
  byrow = TRUE
)
T_train <- rbind(c(1, 0), c(1, 0), c(0, 1))
W_hat <- item_fit(Gamma_train, T_train, diag(3))
item_predict(Gamma_train, W_hat)
#>             [,1]        [,2]
#> [1,] 1.054794521 -0.05479452
#> [2,] 0.938356164  0.06164384
#> [3,] 0.006849315  0.99315068
```
