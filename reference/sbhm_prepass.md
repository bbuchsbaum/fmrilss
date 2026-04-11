# SBHM Prepass: Aggregate Fit in Shared Basis

Compute per-voxel coefficients in the shared SBHM basis by fitting a
single aggregate GLM with one regressor per basis column (trials
summed), optionally residualizing by nuisances and prewhitening. This
produces `beta_bar` (r×V) that you can feed to
[`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md).

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

  Numeric matrix T×V of fMRI time series.

- sbhm:

  SBHM object from
  [`sbhm_build()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_build.md)
  (must contain B, S, A, tgrid, span).

- design_spec:

  List describing events (same shape as `oasis$design_spec`). Must
  contain `sframe` and `cond` with `onsets` (and optional `duration`,
  `amplitude`, `span`). `cond$hrf` is ignored and replaced with
  `sbhm_hrf`. Optional `others` (list of other conditions) will be
  aggregated as nuisances.

- Nuisance:

  Optional T×P nuisance matrix (motion, drift, etc.).

- prewhiten:

  Optional fmriAR prewhitening options (see
  [`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)). If
  provided, Y and design are prewhitened together.

- ridge:

  Optional list for targeted ridge shrinkage in the prepass solve:

  - `mode`: "fractional" (default) or "absolute". Fractional scales by
    mean(diag(G)).

  - `lambda`: nonnegative scalar (default 0.01 in fractional mode).

  - `alpha_ref`: r-vector to shrink towards (default zero vector).

- data_fac:

  Optional list for external factorization: `scores` (T×q), `loadings`
  (q×V). If provided, computes X'Y via (X'Scores) × Loadings. In this
  PR2 version, prewhitening is not applied when `data_fac` is used.

## Value

List with:

- `beta_bar` r×V aggregate coefficients

- `A_agg` T×r aggregated per-basis design (after any
  residualization/whitening)

- `G` r×r crossprod of A_agg

- `diag` list with K=r, ntrials, times, used_prewhiten

## Details

Notes:

- Aggregated per-basis regressors can be highly collinear, making G = A'
  A ill-conditioned. A small ridge is recommended for stability. The
  default uses fractional mode with `lambda = 0.01` (scaled by
  mean(diag(G))).

- When `data_fac` is provided (factorized data path), prewhitening is
  skipped in this version; both dense and factorized paths perform
  nuisance residualization consistently.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(fmrihrf)
  set.seed(1)
  Tlen <- 120; V <- 5; r <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(exp(-seq(0, 30, length.out = Tlen)/4),
             exp(-seq(0, 30, length.out = Tlen)/6))
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  onsets <- seq(5, 95, by = 10)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
  X <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  betas_true <- matrix(rnorm(r), r)
  Y <- matrix(rnorm(Tlen*V, sd = 0.5), Tlen, V)
  Y[,1] <- Y[,1] + X %*% betas_true
  pre <- sbhm_prepass(Y, sbhm, design_spec)
  str(pre)
} # }

```
