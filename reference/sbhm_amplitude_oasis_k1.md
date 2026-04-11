# OASIS K=1 amplitudes per voxel with matched HRF columns. Returns list(beta, se)

OASIS K=1 amplitudes per voxel with matched HRF columns. Returns
list(beta, se)

## Usage

``` r
sbhm_amplitude_oasis_k1(
  Y,
  sbhm,
  design_spec,
  alpha_hat,
  Nuisance = NULL,
  ridge_frac = list(x = 0.02, b = 0.02),
  prewhiten = NULL,
  return_se = TRUE
)
```
