# Build an ITEM bundle from LS-A estimates

Convenience wrapper that runs
[`lsa()`](https://bbuchsbaum.github.io/fmrilss/reference/lsa.md) to
estimate trial-wise amplitudes, computes `U`, and returns an
`item_bundle` ready for crossvalidation.

## Usage

``` r
item_from_lsa(
  Y,
  X_t,
  T_target,
  run_id,
  Z = NULL,
  Nuisance = NULL,
  V = NULL,
  v_type = c("cov", "precision"),
  ridge = 0,
  lsa_method = c("r", "cpp"),
  solver = c("chol", "svd", "pinv"),
  u_output = c("matrix", "by_run"),
  C_transform = NULL,
  trial_id = NULL,
  trial_hash = NULL,
  meta = list(),
  validate = TRUE
)
```

## Arguments

- Y:

  Numeric data matrix (`n_time x n_features`).

- X_t:

  Numeric trial-wise design matrix (`n_time x n_trials`).

- T_target:

  Supervised targets with `n_trials` rows.

- run_id:

  Trial-level run/session identifier (`n_trials`).

- Z:

  Optional nuisance matrix passed to
  [`lsa()`](https://bbuchsbaum.github.io/fmrilss/reference/lsa.md).

- Nuisance:

  Alias for `Z` in
  [`lsa()`](https://bbuchsbaum.github.io/fmrilss/reference/lsa.md).

- V:

  Optional covariance/precision for
  [`item_compute_u()`](https://bbuchsbaum.github.io/fmrilss/reference/item_compute_u.md).

- v_type:

  Whether `V` is covariance or precision.

- ridge:

  Ridge passed to
  [`item_compute_u()`](https://bbuchsbaum.github.io/fmrilss/reference/item_compute_u.md).

- lsa_method:

  LS-A backend (`"r"` or `"cpp"`).

- solver:

  Solver preference for
  [`item_compute_u()`](https://bbuchsbaum.github.io/fmrilss/reference/item_compute_u.md).

- u_output:

  Return full `U` matrix or `U_by_run` blocks.

- C_transform:

  Optional transform matrix used to map `X_t` to the working design
  matrix `X`.

- trial_id:

  Optional trial id vector.

- trial_hash:

  Optional trial hash.

- meta:

  Optional metadata list.

- validate:

  Logical; enforce strict checks before returning.

## Value

`item_bundle` with fields `Gamma`, `X_t`, `T_target`, `U`/`U_by_run`,
`run_id`, `meta`, and `diagnostics`.

## Examples

``` r
set.seed(1)
X_t <- diag(4)
Y <- X_t %*% matrix(
  c(1, 0,
    0.8, 0.2,
    0.2, 0.8,
    0, 1),
  ncol = 2,
  byrow = TRUE
) + matrix(rnorm(8, sd = 0.01), 4, 2)
bundle <- item_from_lsa(
  Y = Y,
  X_t = X_t,
  T_target = factor(c("A", "B", "A", "B")),
  run_id = c(1, 1, 2, 2),
  lsa_method = "r"
)
names(bundle)
#>  [1] "Gamma"       "X_t"         "C_transform" "T_target"    "U"          
#>  [6] "run_id"      "trial_id"    "trial_hash"  "trial_info"  "meta"       
#> [11] "diagnostics"
```
