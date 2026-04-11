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
if (FALSE) { # \dontrun{
  r <- nrow(alpha_hat)
  ntrials <- nrow(beta_mat) / r
  beta_rt <- array(beta_mat, dim = c(r, ntrials, ncol(beta_mat)))
  amps <- sbhm_project(beta_rt, alpha_hat)
} # }
```
