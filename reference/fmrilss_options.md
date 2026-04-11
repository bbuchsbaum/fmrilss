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
#> $alpha
#> [1] 0.2
#> 
#> $lambda
#> [1] 0.1
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
oasis_options(ridge_x = 0.1, ridge_b = 0.1)
#> $design_spec
#> NULL
#> 
#> $K
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
#> attr(,"class")
#> [1] "fmrilss_prewhiten_options" "list"                     
```
