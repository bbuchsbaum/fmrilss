# Compare HRF Recovery Methods

Compares OASIS, SPMG1, SPMG3, and FIR for HRF recovery

## Usage

``` r
compare_hrf_recovery(data, hrf_grid = NULL)
```

## Arguments

- data:

  Synthetic data from generate_lwu_data

- hrf_grid:

  Optional pre-computed HRF grid for OASIS

## Value

List with results from all methods

## Examples

``` r
if (FALSE) { # \dontrun{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
res <- compare_hrf_recovery(sim, hrf_grid = grid)
names(res)
} # }
```
