# Slice an ITEM bundle into train/test fold objects

Create deterministic leave-one-run-out slices for `Gamma`, `T_target`,
and trial covariance (`U` or `U_by_run`).

## Usage

``` r
item_slice_fold(bundle, test_run, check_hash = FALSE)
```

## Arguments

- bundle:

  Object of class `item_bundle`.

- test_run:

  Run/session id to hold out for testing.

- check_hash:

  Logical; if `TRUE`, validate stored `trial_hash`.

## Value

A list with train/test slices: `Gamma_train`, `Gamma_test`, `T_train`,
`T_test`, `U_train`, `U_test`, `train_idx`, `test_idx`, `train_runs`,
`test_run`.

## Examples

``` r
bundle <- item_build_design(
  X_t = diag(4),
  T_target = factor(c("A", "B", "A", "B")),
  run_id = c(1, 1, 2, 2)
)
bundle$Gamma <- matrix(c(1, 0, 0.8, 0.2, 0.2, 0.8, 0, 1), 4, 2, byrow = TRUE)
bundle$U <- diag(4)
fold <- item_slice_fold(bundle, test_run = 2)
fold$test_idx
#> [1] 3 4
```
