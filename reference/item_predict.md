# Predict targets from ITEM weights

Compute `T_hat = Gamma_test %*% W_hat`.

## Usage

``` r
item_predict(Gamma_test, W_hat)
```

## Arguments

- Gamma_test:

  Numeric matrix (`n_test x n_features`).

- W_hat:

  Numeric matrix (`n_features x p`).

## Value

Numeric matrix of predictions (`n_test x p`).

## Examples

``` r
Gamma_test <- matrix(c(1, 0, 0, 1), ncol = 2, byrow = TRUE)
W_hat <- diag(2)
item_predict(Gamma_test, W_hat)
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
```
