# Compute ITEM trial covariance matrix

Compute the trial covariance term as the inverse of a ridge-stabilized
normal-equation system, using stable solve paths with fallbacks.

## Usage

``` r
item_compute_u(
  X_t,
  V = NULL,
  v_type = c("cov", "precision"),
  ridge = 0,
  method = c("chol", "svd", "pinv"),
  run_id = NULL,
  output = c("matrix", "by_run"),
  tol = sqrt(.Machine$double.eps)
)
```

## Arguments

- X_t:

  Numeric trial-wise design matrix (`n_time x n_trials`).

- V:

  Optional temporal covariance/precision object. Accepted forms:

  - `NULL` (identity),

  - dense/sparse matrix (`n_time x n_time`),

  - run-block list of square matrices whose block sizes sum to `n_time`.

- v_type:

  Whether `V` is covariance (`"cov"`) or precision (`"precision"`).

- ridge:

  Non-negative ridge term added to the trial-wise system matrix before
  inversion.

- method:

  Preferred solver path (`"chol"`, `"svd"`, or `"pinv"`).

- run_id:

  Optional trial-level run ids (length `n_trials`). Required only when
  `output = "by_run"`.

- output:

  Return full `U` (`"matrix"`) or split run blocks (`"by_run"`).

- tol:

  Numerical tolerance for rank/solver fallbacks.

## Value

Numeric matrix `U` by default, or a named list `U_by_run` when
`output = "by_run"`. Returned object includes `item_diagnostics`
attribute with rank/condition/solver details.

## Examples

``` r
X_t <- diag(4)
item_compute_u(X_t)
#>      [,1] [,2] [,3] [,4]
#> [1,]    1    0    0    0
#> [2,]    0    1    0    0
#> [3,]    0    0    1    0
#> [4,]    0    0    0    1
#> attr(,"item_diagnostics")
#> attr(,"item_diagnostics")$rank
#> [1] 4
#> 
#> attr(,"item_diagnostics")$condition_number
#> [1] 1
#> 
#> attr(,"item_diagnostics")$solver_path
#> [1] "chol"
#> 
#> attr(,"item_diagnostics")$ridge
#> [1] 0
#> 
#> attr(,"item_diagnostics")$warnings
#> character(0)
#> 
#> attr(,"item_diagnostics")$v_type
#> [1] "cov"
#> 
#> attr(,"item_diagnostics")$n_time
#> [1] 4
#> 
#> attr(,"item_diagnostics")$n_trials
#> [1] 4
#> 
item_compute_u(X_t, run_id = c(1, 1, 2, 2), output = "by_run")
#> $`1`
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> $`2`
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> attr(,"item_diagnostics")
#> attr(,"item_diagnostics")$rank
#> [1] 4
#> 
#> attr(,"item_diagnostics")$condition_number
#> [1] 1
#> 
#> attr(,"item_diagnostics")$solver_path
#> [1] "chol"
#> 
#> attr(,"item_diagnostics")$ridge
#> [1] 0
#> 
#> attr(,"item_diagnostics")$warnings
#> character(0)
#> 
#> attr(,"item_diagnostics")$v_type
#> [1] "cov"
#> 
#> attr(,"item_diagnostics")$n_time
#> [1] 4
#> 
#> attr(,"item_diagnostics")$n_trials
#> [1] 4
#> 
```
