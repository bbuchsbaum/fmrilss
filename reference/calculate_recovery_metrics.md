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

  Ground truth HRF sampled on `results$hrf_times`. Shape comparisons
  normalize each curve to unit maximum absolute height.

## Value

Data frame with recovery metrics

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
calculate_recovery_metrics(res, sim$true_hrf)
#>   method       mse correlation peak_time_error width_error beta_correlation
#> 1  OASIS 0.1250698   0.6063464               2           1       -0.8164999
#> 2  SPMG1 0.0180997   0.9411653               1          NA       -0.3952613
#> 3  SPMG3 0.0180997   0.9411653               1          NA               NA
# }
```
