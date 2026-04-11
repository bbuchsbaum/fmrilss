# Create LWU HRF Grid for OASIS Search

Generates a grid of LWU HRF models with varying parameters

## Usage

``` r
create_lwu_grid(
  tau_range = c(4, 8),
  sigma_range = c(1.5, 3.5),
  rho_range = c(0.1, 0.6),
  n_tau = 5,
  n_sigma = 3,
  n_rho = 3
)
```

## Arguments

- tau_range:

  Range of tau values to test

- sigma_range:

  Range of sigma values to test

- rho_range:

  Range of rho values to test

- n_tau:

  Number of tau values in grid

- n_sigma:

  Number of sigma values in grid

- n_rho:

  Number of rho values in grid

## Value

List of HRF models and their parameters

## Examples

``` r
grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
nrow(grid$parameters)
#> [1] 8
```
