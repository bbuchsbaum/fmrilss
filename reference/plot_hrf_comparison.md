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
if (FALSE) { # \dontrun{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
res <- compare_hrf_recovery(sim, hrf_grid = grid)
plot_hrf_comparison(res)
} # }
```
