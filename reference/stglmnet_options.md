# Construct stglmnet backend options

Convenience constructor for the `stglmnet=` list accepted by
`lss(method = "stglmnet")`. Unknown fields are allowed via `...` for
forward compatibility.

## Usage

``` r
stglmnet_options(
  mode = c("cv", "fixed"),
  alpha = 0.2,
  lambda = NULL,
  overlap_strategy = c("none", "multiplicative", "additive", "hybrid", "threshold"),
  pool_to_mean = FALSE,
  pool_strength = 1,
  pool_mean_penalty = 0,
  whiten = c("inherit", "auto", "never", "always"),
  cv_folds = 5L,
  cv_type.measure = c("auto", "mse", "correlation", "reliability", "composite"),
  cv_select = c("optimal", "1se"),
  return_fit = FALSE,
  ...
)
```

## Arguments

- mode:

  `"cv"` (default) selects lambda by internal cross-validation, while
  `"fixed"` uses the supplied lambda sequence or the smallest fitted
  lambda when no scalar is provided.

- alpha:

  Elastic-net mixing parameter passed to `glmnet`.

- lambda:

  Optional lambda sequence (or scalar in fixed mode).

- overlap_strategy:

  Trial-overlap penalty mapping. One of `"none"`, `"multiplicative"`,
  `"additive"`, `"hybrid"`, or `"threshold"`.

- pool_to_mean:

  Logical; reparameterize trial effects into a pooled mean plus
  orthogonal contrasts.

- pool_strength:

  Penalty multiplier applied to pooled contrasts.

- pool_mean_penalty:

  Penalty applied to the pooled mean coefficient.

- whiten:

  One of `"inherit"` (default), `"auto"`, `"never"`, or `"always"`.
  `"inherit"` uses the top-level `prewhiten=` argument only.

- cv_folds:

  Number of folds used when `mode = "cv"`.

- cv_type.measure:

  Cross-validation objective.

- cv_select:

  Lambda selection rule in CV mode. `"optimal"` uses the best-scoring
  lambda, `"1se"` applies the one-standard-error rule.

- return_fit:

  Logical; when `TRUE`, `lss(method="stglmnet")` returns a list
  containing `beta`, fit metadata, and the selected lambda.

- ...:

  Additional backend options.

## Value

A list with class `"fmrilss_stglmnet_options"`.

## Examples

``` r
stglmnet_options(mode = "fixed", lambda = 0.05, alpha = 0.5)
#> $mode
#> [1] "fixed"
#> 
#> $alpha
#> [1] 0.5
#> 
#> $lambda
#> [1] 0.05
#> 
#> $overlap_strategy
#> [1] "none"
#> 
#> $pool_to_mean
#> [1] FALSE
#> 
#> $pool_strength
#> [1] 1
#> 
#> $pool_mean_penalty
#> [1] 0
#> 
#> $whiten
#> [1] "inherit"
#> 
#> $cv_folds
#> [1] 5
#> 
#> $cv_type.measure
#> [1] "auto"
#> 
#> $cv_select
#> [1] "optimal"
#> 
#> $return_fit
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "fmrilss_stglmnet_options" "list"                    
```
