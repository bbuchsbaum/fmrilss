# Plot HRF Recovery Comparison

Creates visualization comparing true vs recovered HRFs

## Usage

``` r
plot_hrf_comparison(results, save_path = NULL)
```

## Arguments

- results:

  Output from compare_hrf_recovery

- save_path:

  Optional path to save plot

## Value

A `ggplot2` plot object. When `save_path` is supplied, the same plot is
also written to disk.

## Examples

``` r
# \donttest{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
res <- compare_hrf_recovery(sim, hrf_grid = grid)
#> Fitting OASIS with HRF grid search...
#> Fitting SPMG1...
#> Fitting SPMG3...
#> Fitting FIR...
plot_hrf_comparison(res)
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the fmrilss package.
#>   Please report the issue at <https://github.com/bbuchsbaum/fmrilss/issues>.

# }
```
