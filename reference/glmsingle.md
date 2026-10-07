# GLMsingle single-trial response estimation

Estimates single-trial response amplitudes with the GLMsingle procedure
(Prince et al., 2022): a library-based HRF per voxel, GLMdenoise noise
regressors chosen by cross-validation, and voxel-wise fractional ridge
regression chosen by cross-validation over repeated conditions. The
implementation reorganises the computation (per-run factorisation, a
single pass of data products per stage, compiled cross-validation and
sufficient-statistic R^2) so it is much faster than the reference
implementation while computing the same estimator.

## Usage

``` r
glmsingle(
  Y,
  design,
  tr,
  stimdur,
  runs = NULL,
  hrf_library = NULL,
  want_library = TRUE,
  want_glmdenoise = TRUE,
  want_fracridge = TRUE,
  fracs = seq(1, 0.05, by = -0.05),
  n_pcs = 10L,
  pcstop = 1.05,
  xval_scheme = NULL,
  session_indicator = NULL,
  extra_regressors = NULL,
  max_poly_deg = NULL,
  brain_thresh = c(99, 0.1),
  brain_r2 = NULL,
  brain_exclude = NULL,
  pc_r2_cutoff = NULL,
  pc_r2_cutoff_mask = NULL,
  want_percent_bold = TRUE,
  want_autoscale = TRUE,
  extras_in_denoise = c("always", "with_pcs"),
  zero_sd_cv = c("zero", "python"),
  singular = c("error", "pinv"),
  frac_alpha = c("fracridge", "exact"),
  full_glmbadness = FALSE,
  chunk_size = 50000L,
  n_threads = 1L,
  verbose = TRUE
)
```

## Arguments

- Y:

  Data: a list of time x voxel matrices (one per run), or one time x
  voxel matrix with `runs` giving the run of each row.

- design:

  Either a list of time x condition 0/1 matrices, one per run
  (GLMsingle's format; a 1 marks a trial onset), or a data frame with
  columns `run`, `onset` (seconds; must be on the TR grid) and
  `condition`. Repeated conditions drive the cross-validation. For
  matrix `Y`, event run IDs are matched to `runs` in data order. For
  list `Y`, ascending event run IDs correspond to the list order.

- tr:

  Repetition time in seconds.

- stimdur:

  Trial duration in seconds.

- runs:

  Run identifier per row of `Y` when `Y` is a single matrix.

- hrf_library:

  Optional time x HRF matrix sampled at the TR (first row at onset).
  Default:
  [`glmsingle_hrf_library()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle_hrf_library.md).

- want_library:

  Fit the HRF library per voxel (otherwise the canonical HRF is used for
  every voxel).

- want_glmdenoise:

  Derive and select GLMdenoise noise regressors.

- want_fracridge:

  Fit the fractional ridge model (type D).

- fracs:

  Ridge fractions in (0, 1\] to evaluate. A single value skips
  cross-validation and uses that fraction.

- n_pcs:

  Maximum number of noise PCs to evaluate. Capped, with a warning, at
  the available noise-pool rank across runs. Empty or rank-zero pools
  use zero PCs and skip PC-count cross-validation.

- pcstop:

  Stopping factor for choosing the number of PCs. A value `<= 0` uses
  `-pcstop` PCs without cross-validation.

- xval_scheme:

  List of integer vectors of runs held out together in each
  cross-validation fold. Default: leave one run out.

- session_indicator:

  Integer session per run, used to z-score betas within session for
  cross-validation. Default: one session.

- extra_regressors:

  Optional list (one per run) of time x regressor nuisance matrices
  (e.g. motion), or `NULL`.

- max_poly_deg:

  Polynomial drift degree, scalar or one per run. Default:
  `round(run_seconds / 120)` as in GLMsingle.

- brain_thresh:

  Two numbers: a percentile of the mean volume and a fraction of it;
  voxels brighter than their product may enter the noise pool.

- brain_r2:

  ON-OFF R^2 (percent) below which bright voxels enter the noise pool.
  Default: estimated tail threshold. With fewer than two distinct finite
  ON-OFF R^2 values, uses the common value (or zero if none are finite).
  Thresholds are only estimated when denoising is used.

- brain_exclude:

  Optional logical vector of voxels to exclude from the noise pool
  (`FALSE` excludes).

- pc_r2_cutoff:

  ON-OFF R^2 above which voxels summarise the PC-count cross-validation.
  Default: estimated tail threshold.

- pc_r2_cutoff_mask:

  Optional logical vector restricting those voxels.

- want_percent_bold:

  Express betas as percent signal change of the mean volume.

- want_autoscale:

  Rescale ridge betas to best match the unregularised betas (type D).

- extras_in_denoise:

  How user extra regressors enter the GLMdenoise and ridge fits:
  `"always"` (default, recommended) or `"with_pcs"` (only when at least
  one PC is used; GLMsingle's behaviour).

- zero_sd_cv:

  Treatment of voxels whose reference betas have zero variance within a
  session during cross-validation: `"zero"` (default) or `"python"`
  (reproduce GLMsingle's Python behaviour).

- singular:

  What to do when a run's trial regressors are linearly dependent:
  `"error"` (default, as GLMsingle) or `"pinv"` (minimum-norm solution
  with a warning).

- frac_alpha:

  How fractions map to ridge penalties: `"fracridge"` (default;
  GLMsingle's grid interpolation) or `"exact"` (solve the fraction
  equation exactly; experimental).

- full_glmbadness:

  Compute the PC-count cross-validation for every voxel instead of only
  the voxels that decide the PC count. Estimates are identical; only the
  `glmbadness` diagnostic is filled for all voxels.

- chunk_size:

  Number of voxels processed at a time (GLMsingle's `chunklen`). Lower
  it to reduce peak memory.

- n_threads:

  Threads for the per-voxel C++ loops (`0` = OpenMP default). Results do
  not depend on the thread count. Matrix products use the BLAS library's
  own threads; use `n_threads > 1` only with a single-threaded BLAS,
  because the two thread pools compete for cores.

- verbose:

  Print progress messages.

## Value

An object of class `glmsingle_fit`: a list with elements `typea`,
`typeb`, `typec`, `typed` (each a list using GLMsingle's field names;
betas are trials x voxels), `meanvol`, `hrf_library`, `design` (trial
bookkeeping), `settings` and `timing`.

## Details

Four models are returned, matching GLMsingle:

- typea:

  ON-OFF model: one canonical-HRF regressor for all trials.

- typeb:

  Single-trial OLS with the best library HRF per voxel.

- typec:

  Type B plus GLMdenoise noise regressors.

- typed:

  Type C plus voxel-wise fractional ridge regression.

## Argument names

Arguments use snake_case; the GLMsingle names are `wantlibrary`
(`want_library`), `wantglmdenoise` (`want_glmdenoise`), `wantfracridge`
(`want_fracridge`), `xvalscheme` (`xval_scheme`), `sessionindicator`
(`session_indicator`), `maxpolydeg` (`max_poly_deg`), `brainthresh`
(`brain_thresh`), `brainR2` (`brain_r2`), `brainexclude`
(`brain_exclude`), `pcR2cutoff` (`pc_r2_cutoff`), `pcR2cutoffmask`
(`pc_r2_cutoff_mask`), `wantpercentbold` (`want_percent_bold`),
`wantautoscale` (`want_autoscale`) and `chunklen` (`chunk_size`). Output
fields keep GLMsingle's names.

## Differences from GLMsingle

Defaults reproduce GLMsingle (Python, commit 1ab54a6) except where
GLMsingle is internally inconsistent:

- `extras_in_denoise = "always"` keeps user extra regressors in every
  GLMdenoise and ridge fit; GLMsingle (Python) drops them when zero
  noise PCs are used, including from the cross-validation reference.
  `"with_pcs"` reproduces GLMsingle.

- `zero_sd_cv = "zero"` gives voxels with zero beta variance no weight
  in cross-validation; GLMsingle's in-place division helper treats them
  inconsistently across candidates. `"python"` reproduces it.

- The ON-OFF R^2 threshold uses a deterministic two-Gaussian mixture fit
  with a small variance floor (as in GLMsingle's MATLAB code);
  GLMsingle's Python fit is unseeded.

- `want_percent_bold` applies to all model types (GLMsingle always
  scales types C and D).

- Computation is in double precision (GLMsingle uses single).

- HRF indices are 1-based; betas are trials x voxels.

The FIR diagnostic model, `hrfmodel = "optimize"`, `wantlss`, bootstrap
modes, figures and file outputs are not implemented.

## References

Prince, J. S., Charest, I., Kurzawski, J. W., Pyles, J. A., Tarr, M. J.,
& Kay, K. N. (2022). Improving the accuracy of single-trial fMRI
response estimates using GLMsingle. eLife, 11, e77599.

## Examples

``` r
set.seed(1)
tr <- 1; n_time <- 120; n_vox <- 40
design <- lapply(1:3, function(r) {
  D <- matrix(0, n_time, 6)
  onsets <- seq(5, 95, by = 8)
  D[cbind(onsets, rep_len(sample(6), length(onsets)))] <- 1
  D
})
Y <- lapply(design, function(D) {
  X <- apply(D, 2, function(s) stats::filter(s, glmsingle_hrf(3, tr),
                                              sides = 1, circular = FALSE))
  X[is.na(X)] <- 0
  100 + X %*% matrix(rnorm(6 * n_vox, 2), 6) + matrix(rnorm(n_time * n_vox), n_time)
})
fit <- glmsingle(Y, design, tr = tr, stimdur = 3, n_pcs = 2, verbose = FALSE)
dim(fit$typed$betasmd)
#> [1] 36 40
```
