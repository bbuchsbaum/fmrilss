# Fit OASIS with HRF Grid Search

Selects an LWU HRF by joint-model profile fit, then estimates OASIS
betas.

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

## Details

Each candidate is scored by pooled R-squared from unpenalized joint
least squares with a voxel-specific intercept and all trial columns. The
denominator sums squared deviations from each voxel's own mean. Ridge
parameters affect only the final LSS estimates, not HRF selection.
Simulation, scoring, and fitting use the same event construction at
0.1-second precision. This is an in-sample selection criterion;
overlapping events can leave HRF parameters weakly identifiable.

## Examples

``` r
# \donttest{
onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
fit <- fit_oasis_grid(sim$Y, sim$onsets, sim$sframe, grid)
fit$best_params
#>   tau sigma rho
#> 4   8   3.5 0.1
# }
```
