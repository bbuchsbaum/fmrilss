# Calculate HRF Recovery Metrics

Evaluates how well each method recovered the true HRF

## Usage

``` r
calculate_recovery_metrics(results, true_hrf)
```

## Arguments

- results:

  Output from compare_hrf_recovery

- true_hrf:

  Ground truth HRF

## Value

Data frame with recovery metrics

## Examples

``` r
if (FALSE) { # \dontrun{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
res <- compare_hrf_recovery(sim, hrf_grid = grid)
calculate_recovery_metrics(res, sim$true_hrf)
} # }
```
