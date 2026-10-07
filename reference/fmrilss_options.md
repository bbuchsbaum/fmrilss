# Option constructors for nested interfaces

These helpers create validated option lists for
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) and
friends.

## Value

No value itself. This topic groups the documented constructors
[`stglmnet_options()`](https://bbuchsbaum.github.io/fmrilss/reference/stglmnet_options.md),
[`oasis_options()`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md),
and
[`prewhiten_options()`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md).

## Examples

``` r
stglmnet_options(mode = "fixed", lambda = 0.1)
#> $mode
#> [1] "fixed"
#> 
#> $run_id
#> NULL
#> 
#> $alpha
#> [1] 0.2
#> 
#> $lambda
#> [1] 0.1
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
oasis_options(ridge_x = 0.1, ridge_b = 0.1)
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
prewhiten_options(method = "ar", p = 1)
#> $method
#> [1] "ar"
#> 
#> $p
#> [1] 1
#> 
#> $q
#> [1] 0
#> 
#> $p_max
#> [1] 6
#> 
#> $pooling
#> [1] "global"
#> 
#> $runs
#> NULL
#> 
#> $parcels
#> NULL
#> 
#> $exact_first
#> [1] "ar1"
#> 
#> $compute_residuals
#> [1] TRUE
#> 
#> $design
#> NULL
#> 
#> $acvf_correction
#> NULL
#> 
#> $correction_max_lag
#> [1] 25
#> 
#> attr(,"class")
#> [1] "fmrilss_prewhiten_options" "list"                     
```
