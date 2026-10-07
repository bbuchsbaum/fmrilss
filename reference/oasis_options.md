# Construct OASIS options

Convenience constructor for the `oasis=` list accepted by
`lss(method="oasis")`. Unknown fields are rejected so misspelled
scientific controls cannot be silently ignored.

## Usage

``` r
oasis_options(
  design_spec = NULL,
  K = NULL,
  ntrials = NULL,
  trial_basis_map = NULL,
  ridge_mode = c("fractional", "absolute"),
  ridge_x = 0.05,
  ridge_b = 0.05,
  block_cols = 4096L,
  return_se = FALSE,
  return_diag = FALSE,
  add_intercept = TRUE,
  hrf_mode = NULL,
  ...
)
```

## Arguments

- design_spec:

  Optional design spec list used to build `X` via `fmrihrf`.

- K:

  Optional basis dimension override.

- ntrials:

  Number of trials for a raw multi-basis design.

- trial_basis_map:

  For a raw multi-basis design, a data frame with one row per `X` column
  and fields `column`, `trial`, and `basis`.

- ridge_mode:

  `"fractional"` (default) or `"absolute"`.

- ridge_x, ridge_b:

  Non-negative ridge penalties (defaults 0.05 each).

- block_cols:

  Voxel block size for blocked products.

- return_se:

  Logical; return standard errors.

- return_diag:

  Logical; return diagnostics.

- add_intercept:

  Logical; add intercept when `Z` is NULL.

- hrf_mode:

  Optional mode (e.g. `"voxhrf"`); advanced use.

- ...:

  Certified advanced options: `infer_K_from_X`, `lambda_shape`,
  `mu_rough`, `ref_hrf`, `shrink_global`, `orient_ref`, and the
  deprecated `whiten` compatibility field. Unknown names are rejected.

## Value

A list with class `"fmrilss_oasis_options"`.

## Examples

``` r
oasis_options(ridge_mode = "fractional", ridge_x = 0.1, ridge_b = 0.1)
#> $design_spec
#> NULL
#> 
#> $K
#> NULL
#> 
#> $ntrials
#> NULL
#> 
#> $trial_basis_map
#> NULL
#> 
#> $ridge_mode
#> [1] "fractional"
#> 
#> $ridge_x
#> [1] 0.1
#> 
#> $ridge_b
#> [1] 0.1
#> 
#> $block_cols
#> [1] 4096
#> 
#> $return_se
#> [1] FALSE
#> 
#> $return_diag
#> [1] FALSE
#> 
#> $add_intercept
#> [1] TRUE
#> 
#> $hrf_mode
#> NULL
#> 
#> attr(,"class")
#> [1] "fmrilss_oasis_options" "list"                 
```
