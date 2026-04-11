# Fit OASIS with HRF Grid Search

Fits OASIS models with different HRF parameters and selects best

## Usage

``` r
fit_oasis_grid(Y, onsets, sframe, hrf_grid, ridge_x = 0.01, ridge_b = 0.01)
```

## Arguments

- Y:

  Data matrix (time x voxels)

- onsets:

  Event onset times

- sframe:

  Sampling frame

- hrf_grid:

  List of HRF models to test

- ridge_x:

  Ridge parameter for design matrix

- ridge_b:

  Ridge parameter for aggregator

## Value

List with best HRF index, parameters, and beta estimates

## Examples

``` r
if (FALSE) { # \dontrun{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
fit <- fit_oasis_grid(sim$Y, sim$onsets, sim$sframe, grid)
fit$best_params
} # }
```
