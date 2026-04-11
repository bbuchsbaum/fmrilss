# Crossvalidated ITEM decoding

Run deterministic leave-one-run-out (LOSO) crossvalidation for
classification or regression using ITEM decoder weights.

## Usage

``` r
item_cv(
  Gamma,
  T_target = NULL,
  U = NULL,
  run_id = NULL,
  mode = c("classification", "regression"),
  metric = NULL,
  ridge = 0,
  method = c("chol", "svd", "pinv"),
  class_levels = NULL,
  trial_id = NULL,
  trial_hash = NULL,
  check_hash = FALSE
)
```

## Arguments

- Gamma:

  Trial-wise beta matrix (`n_trials x n_features`) or an `item_bundle`
  object.

- T_target:

  Supervised targets (`n_trials x p`) or target vector. Ignored when
  `Gamma` is an `item_bundle`.

- U:

  Trial covariance as full matrix (`n_trials x n_trials`) or run-block
  list. Ignored when `Gamma` is an `item_bundle`.

- run_id:

  Run/session id vector (`n_trials`). Ignored when `Gamma` is an
  `item_bundle`.

- mode:

  Decoding mode: `"classification"` or `"regression"`.

- metric:

  Optional metric name. Classification: `"accuracy"` (default),
  `"balanced_accuracy"`. Regression: `"correlation"` (default),
  `"rmse"`.

- ridge:

  Ridge passed to
  [`item_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/item_fit.md).

- method:

  Solver preference for
  [`item_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/item_fit.md).

- class_levels:

  Optional fixed class order for classification.

- trial_id:

  Optional trial id vector (used if `Gamma` is not a bundle).

- trial_hash:

  Optional trial hash (used if `Gamma` is not a bundle).

- check_hash:

  Logical; validate stored trial hash before CV.

## Value

Object of class `item_cv_result` with per-fold metrics, aggregate metric
summary, and trial-level predictions.

## Examples

``` r
Gamma <- matrix(
  c(1, 0,
    0.9, 0.1,
    0.1, 0.9,
    0, 1),
  ncol = 2,
  byrow = TRUE
)
item_cv(
  Gamma = Gamma,
  T_target = factor(c("A", "A", "B", "B")),
  U = diag(4),
  run_id = c(1, 1, 2, 2)
)
#> $mode
#> [1] "classification"
#> 
#> $metric
#> [1] "accuracy"
#> 
#> $folds
#>   fold test_run n_train n_test metric
#> 1    1        1       2      2      0
#> 2    2        2       2      2      0
#> 
#> $aggregate
#> $aggregate$metric
#> [1] "accuracy"
#> 
#> $aggregate$mean
#> [1] 0
#> 
#> $aggregate$sd
#> [1] 0
#> 
#> $aggregate$n_folds
#> [1] 2
#> 
#> 
#> $predictions
#> $predictions$T_hat
#>         A B
#> Trial_1 0 1
#> Trial_2 0 1
#> Trial_3 1 0
#> Trial_4 1 0
#> 
#> $predictions$T_true
#>   A B
#> 1 1 0
#> 2 1 0
#> 3 0 1
#> 4 0 1
#> attr(,"assign")
#> [1] 1 1
#> attr(,"contrasts")
#> attr(,"contrasts")$fac
#> [1] "contr.treatment"
#> 
#> 
#> $predictions$predicted_class
#> [1] "B" "B" "A" "A"
#> 
#> $predictions$true_class
#> [1] "A" "A" "B" "B"
#> 
#> $predictions$class_levels
#> [1] "A" "B"
#> 
#> 
#> $diagnostics
#> $diagnostics$fold_order
#> [1] "1" "2"
#> 
#> $diagnostics$solver
#> [1] "chol"
#> 
#> $diagnostics$ridge
#> [1] 0
#> 
#> 
#> attr(,"class")
#> [1] "item_cv_result"
```
