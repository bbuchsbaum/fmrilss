# Project Trial-wise SBHM Coefficients to Scalar Amplitudes

Given trial-wise coefficients in the shared basis (rxntrialsxV) and the
voxel-specific matched library coordinates `alpha_hat` (rxV), compute
scalar amplitudes per trial and voxel via least-squares projection:
`a = (alpha' beta) / (alpha' alpha)`.

## Usage

``` r
sbhm_project(beta_rt, alpha_hat)
```

## Arguments

- beta_rt:

  3D array of shape r x ntrials x V containing per-trial coefficients in
  the SBHM basis (as returned by OASIS with K=r, reshaped).

- alpha_hat:

  Numeric matrix rxV of matched library coordinates per voxel (e.g.,
  `sbhm_match()$alpha_hat`). These should be in the same coordinate
  system as `beta_rt` (unwhitened, not L2-normalized) for interpretable
  amplitudes.

## Value

Numeric matrix ntrials x V of scalar amplitudes.

## Examples

``` r
# \donttest{
  set.seed(1)
  r <- 2; ntrials <- 3; nvox <- 2
  alpha_hat <- matrix(rnorm(r * nvox), r, nvox)
  beta_rt <- array(rnorm(r * ntrials * nvox), c(r, ntrials, nvox))
  amps <- sbhm_project(beta_rt, alpha_hat)
  dim(amps)
#> [1] 3 2
# }
```
