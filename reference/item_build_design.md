# Build ITEM design metadata

Build and validate trial-wise design objects used by the ITEM helper
layer. The returned object has class `item_bundle` and carries
trial-level metadata needed by
[`item_cv()`](https://bbuchsbaum.github.io/fmrilss/reference/item_cv.md)
and
[`item_slice_fold()`](https://bbuchsbaum.github.io/fmrilss/reference/item_slice_fold.md).

## Usage

``` r
item_build_design(
  X_t,
  T_target = NULL,
  run_id = NULL,
  C_transform = NULL,
  trial_id = NULL,
  trial_hash = NULL,
  meta = list(),
  diagnostics = list(),
  validate = TRUE
)
```

## Arguments

- X_t:

  Numeric trial-wise design matrix with shape `n_time x n_trials`.

- T_target:

  Optional supervised targets with `n_trials` rows. Accepts: numeric
  matrix/data.frame, numeric vector (regression), or
  factor/character/logical vector (classification labels).

- run_id:

  Optional run/session identifier of length `n_trials`. If `NULL`, all
  trials are assigned to a single run.

- C_transform:

  Optional transformation matrix used to map `X_t` to the working design
  matrix `X`. Must have `n_trials` rows when provided.

- trial_id:

  Optional trial identifier vector of length `n_trials`. Defaults to
  `colnames(X_t)` when available, else `Trial_1..Trial_n`.

- trial_hash:

  Optional hash used by alignment guards. If supplied,
  `item_cv(..., check_hash = TRUE)` validates it.

- meta:

  Optional metadata list.

- diagnostics:

  Optional diagnostics list to attach to the bundle.

- validate:

  Logical; when `TRUE` run strict structure checks.

## Value

An object of class `item_bundle` with fields: `Gamma`, `X_t`,
`C_transform`, `T_target`, `U`, `U_by_run`, `run_id`, `trial_id`,
`trial_hash`, `trial_info`, `meta`, and `diagnostics`.

## Examples

``` r
bundle <- item_build_design(
  X_t = diag(4),
  T_target = factor(c("A", "B", "A", "B")),
  run_id = c(1, 1, 2, 2)
)
bundle$trial_info
#>   trial_index trial_id run_id
#> 1           1  Trial_1      1
#> 2           2  Trial_2      1
#> 3           3  Trial_3      2
#> 4           4  Trial_4      2
```
