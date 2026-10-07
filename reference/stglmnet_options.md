# Construct stglmnet backend options

Convenience constructor for the `stglmnet=` list accepted by
`lss(method = "stglmnet")`. Advanced fields accepted through `...` are
validated; unknown fields are rejected so misspelled scientific controls
cannot be silently ignored.

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

  Certified advanced backend fields such as graph-pooling,
  overlap-strength, fold, and whitening controls. Unknown names are
  rejected.

## Value

A list with class `"fmrilss_stglmnet_options"`.

## Examples

``` r
stglmnet_options(mode = "fixed", lambda = 0.05, alpha = 0.5)
#> $mode
#> [1] "fixed"
#> 
#> $run_id
#> NULL
#> 
#> $alpha
#> [1] 0.5
#> 
#> $lambda
#> [1] 0.05
#> 
#> $family
#> NULL
#> 
#> $standardize
#> [1] FALSE
#> 
#> $intercept
#> [1] FALSE
#> 
#> $overlap_strategy
#> [1] "none"
#> 
#> $overlap_strength
#> [1] 1
#> 
#> $overlap_mix
#> [1] 0.5
#> 
#> $overlap_threshold
#> [1] 0.75
#> 
#> $overlap_exponent
#> [1] 1
#> 
#> $graph_pool
#> [1] FALSE
#> 
#> $graph_strength
#> [1] 1
#> 
#> $graph_exponent
#> [1] 1
#> 
#> $graph_mean_penalty
#> [1] 0.2
#> 
#> $graph_metric
#> [1] "corr"
#> 
#> $graph_scale_by_overlap
#> [1] TRUE
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
#> $pool_scale_by_overlap
#> [1] TRUE
#> 
#> $nuisance_penalty
#> [1] 0
#> 
#> $whiten
#> [1] "inherit"
#> 
#> $whiten_threshold
#> [1] 0.15
#> 
#> $prewhiten_args
#> $prewhiten_args$method
#> [1] "ar"
#> 
#> $prewhiten_args$p
#> [1] 1
#> 
#> $prewhiten_args$pooling
#> [1] "global"
#> 
#> 
#> $cv_folds
#> [1] 5
#> 
#> $cv_type.measure
#> [1] "auto"
#> 
#> $cv_fold_scheme
#> [1] "run"
#> 
#> $cv_select
#> [1] "optimal"
#> 
#> $overlap_low_threshold
#> [1] 0.12
#> 
#> $composite_weights
#>         mse correlation reliability 
#>         0.4         0.3         0.3 
#> 
#> $return_fit
#> [1] FALSE
#> 
#> attr(,"class")
#> [1] "fmrilss_stglmnet_options" "list"                    
```
