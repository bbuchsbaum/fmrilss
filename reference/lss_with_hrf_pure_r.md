# Least Squares Separate with voxel-wise HRF (basis-weighted)

Compute LSS trial-wise betas when each voxel has its own HRF formed as a
linear combination of K basis kernels sampled on the TR grid.

## Usage

``` r
lss_with_hrf_pure_r(
  Y,
  onset_idx,
  durations = NULL,
  hrf_basis_kernels,
  coefficients,
  Z = NULL,
  Nuisance = NULL,
  verbose = FALSE,
  method = c("r", "cpp", "cpp_arma", "cpp_omp")
)
```

## Arguments

- Y:

  numeric matrix (n_time x n_vox)

- onset_idx:

  integer vector (length n_trials), 1-based TR indices

- durations:

  numeric vector (length n_trials), in TRs; 0 means an impulse. Uses
  inclusive indexing: if duration = d, samples `o:(o+d)` are 1.

- hrf_basis_kernels:

  numeric matrix (L x K), K basis kernels on TR grid

- coefficients:

  numeric matrix (K x n_vox), voxel-wise HRF weights

- Z:

  optional numeric matrix (n_time x F) of experimental regressors; if
  NULL, an intercept (column of 1s) is used.

- Nuisance:

  optional numeric matrix (n_time x q) of confounds to project out

- verbose:

  logical; print progress every 1000 voxels

- method:

  character: "r" (default, pure R), "cpp" (C++ backend), "cpp_arma"
  (Armadillo backend), or "cpp_omp" (OpenMP parallel backend). Falls
  back automatically: cpp_omp -\> cpp_arma -\> cpp -\> r.

## Value

numeric matrix (n_trials x n_vox) of trial-wise beta estimates

## Details

**Design & nuisance handling match
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)**:

- The trial-of-interest (Xi) and the sum of all other trials (Xother)
  are included in each per-trial GLM.

- If `Nuisance` is supplied, it is projected out of **Y** and the trial
  regressors before LSS (standard residualization). Experimental
  regressors `Z` are *not* residualized, matching
  [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)
  documentation.

- If `Z` is `NULL`, an intercept-only design is used.

## Examples

``` r
if (FALSE) { # \dontrun{
betas <- lss_with_hrf_pure_r(Y, onset_idx, durations, basis, coeffs, Z = NULL, Nuisance = NULL)
betas <- lss_with_hrf_pure_r(Y, onset_idx, durations, basis, coeffs, method = "cpp")
} # }
```
