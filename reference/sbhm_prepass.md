# SBHM Prepass: Aggregate Fit in a Shared Basis

Fit the trial-aggregated SBHM basis to every voxel after projecting
run-specific intercepts, supplied nuisance columns, and complete modeled
other-condition spans. The same run-safe design builder is used by the
full SBHM pipeline.

## Usage

``` r
sbhm_prepass(
  Y,
  sbhm,
  design_spec,
  Nuisance = NULL,
  prewhiten = NULL,
  ridge = list(mode = "fractional", lambda = 0.01, alpha_ref = NULL),
  data_fac = NULL
)
```

## Arguments

- Y:

  Finite numeric T by V response matrix.

- sbhm:

  Object returned by
  [`sbhm_build()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_build.md).

- design_spec:

  Event specification with `sframe` and `cond`; multi-run inputs use
  run-relative onsets and an explicit `cond$run` vector.

- Nuisance:

  Optional T by P nuisance matrix.

- prewhiten:

  Optional prewhitening options passed to the package prewhitening
  layer. Run labels are inferred from `design_spec$sframe` when omitted;
  supplied labels must match those sampling-frame boundaries.

- ridge:

  List with `mode`, nonnegative `lambda`, and optional rank-length
  `alpha_ref`.

- data_fac:

  Optional exact or approximate factorization with `scores` T by q and
  `loadings` q by V. Active prewhitening is unsupported with this
  shortcut, and the original full `Y` remains required by
  [`lss_sbhm()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md).
  If factor axes are named, scores and loadings must provide the same
  complete unique names and are aligned before multiplication. When `Y`
  is named, loadings must provide the same complete unique voxel names
  and are reordered to `Y`.

## Value

A list containing named `beta_bar` (rank by voxel), the residualized
aggregate design `A_agg`, its Gram matrix `G`, and design diagnostics
and identity maps. Active whitening is recorded in the `whiten_plan`
attribute.

## Examples

``` r
times <- seq(0, 30, by = 0.5)
H <- cbind(stats::dgamma(times, 5, 1), stats::dgamma(times, 7, 1))
basis <- sbhm_build(library_H = H, tgrid = times, span = 30, r = 2)
spec <- list(sframe = fmrihrf::sampling_frame(80L, TR = 1),
             cond = list(onsets = c(5, 20, 35, 50), duration = 0))
set.seed(1)
pre <- sbhm_prepass(matrix(rnorm(80 * 3), 80, 3), basis, spec)
dim(pre$beta_bar)
#> [1] 2 3
```
