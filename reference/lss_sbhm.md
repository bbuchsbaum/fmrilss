# End-to-End LSS with Shared-Basis HRF Matching (SBHM)

Orchestrates the SBHM pipeline: (1) prepass aggregate fit in the learned
shared basis, (2) cosine matching to a library of HRFs represented in
the same basis, (3) Least Squares Separate (OASIS) with the SBHM basis
to obtain trial-wise r-dimensional coefficients, and (4) projection of
those coefficients onto the matched coordinates to produce scalar
amplitudes.

## Usage

``` r
lss_sbhm(
  Y,
  sbhm,
  design_spec,
  Nuisance = NULL,
  prewhiten = NULL,
  prepass = list(),
  match = list(shrink = list(tau = 0, ref = NULL, snr = NULL), topK = 3, soft_blend =
    TRUE, blend_margin = 0.08, whiten = FALSE, sv_floor_rel = 0.05, whiten_power = 0.5,
    min_margin = NULL, min_beta_norm = NULL, fallback_ref = NULL, orient_ref = TRUE,
    alpha_source = "prepass", rank1_min = 0),
  oasis = list(),
  amplitude = list(method = "lss1", ridge = list(mode = "fractional", lambda = 0.02),
    ridge_frac = list(x = 0.02, b = 0.02), cond_gate = NULL, adaptive = list(enable =
    FALSE, base = 0.02, k0 = 1000, max = 0.08), return_se = FALSE),
  return = c("amplitude", "coefficients", "both")
)
```

## Arguments

- Y:

  Numeric matrix T×V of fMRI time series.

- sbhm:

  SBHM object from
  [`sbhm_build()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_build.md).

- design_spec:

  List for design construction (same as `oasis$design_spec`):
  `list(sframe=..., cond=list(onsets=..., duration=0, span=...), others=list(...))`.
  The HRF in `cond$hrf` is ignored and replaced with the SBHM basis HRF.

- Nuisance:

  Optional T×P nuisance regressors.

- prewhiten:

  Optional prewhitening options (see
  [`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)).

- prepass:

  Optional list forwarded to
  [`sbhm_prepass()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_prepass.md)
  (e.g., ridge, data_fac).

- match:

  Optional list forwarded to
  [`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md)
  (e.g., shrink, topK, whiten, orient_ref). Additional fields handled
  here:

  - `alpha_source`: one of `"prepass"` (default), `"trial_projection"`,
    or `"oasis_rank1"`. `"trial_projection"` estimates voxel shape from
    per-trial projection coefficients in the shared basis.

  - `rank1_min`: optional minimum rank-1 variance fraction in `[0,1]`
    when `alpha_source="oasis_rank1"`. Voxels below threshold fall back
    to prepass.

  - `soft_blend` logical (default TRUE): when `topK > 1`, blend the
    top-K library coordinates per voxel using softmax weights returned
    by
    [`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md).
    If `blend_margin` is provided, blending is only applied to voxels
    with `margin < blend_margin`; others use the hard top-1 assignment.

  - `blend_margin` optional numeric threshold on the matching margin for
    conditional blending.

  - `whiten_power` numeric in `[0,1]` for partial singular-value
    whitening (`1`=full, `0.5`=partial).

  - `min_margin` optional minimum matching margin. Voxels below
    threshold fall back to `fallback_ref`.

  - `min_beta_norm` optional minimum norm of the shape summary used for
    matching. Voxels below threshold fall back to `fallback_ref`.

  - `fallback_ref` optional r-vector fallback coordinate (default
    `sbhm$ref$alpha_ref`).

- oasis:

  Optional list forwarded to `lss(..., method="oasis")`. `K` is set to
  `ncol(sbhm$B)` if not provided, and `design_spec` is injected
  automatically.

- amplitude:

  List controlling the scalar amplitude stage. Fields:

  - `method`: one of "lss1" (default), "global_ls", "oasis_voxel".

  - `ridge`: for `global_ls`, either numeric (absolute) or list(mode,
    lambda).

  - `ridge_frac`: for `lss1`/`oasis_voxel`, list(x, b) fractional ridge.

  - `cond_gate`: optional auto-fallback rule, e.g., list(metric="rho",
    thr=0.999, fallback="lss1").

- return:

  One of `"amplitude"`, `"coefficients"`, or `"both"` (default
  `"amplitude"`).

## Value

A list with components:

- `amplitude` ntrials×V matrix (when requested)

- `coeffs_r` r×ntrials×V array of trial-wise coefficients (when
  requested)

- `matched_idx` length-V integer indices into the library

- `margin` length-V confidence margins (top1 - top2 cosine)

- `alpha_coords` r×V matched coordinates per voxel

- `diag` list with `r`, `ntrials`, and `times`

## Details

Most users should treat the `prepass`, `match`, `oasis`, and `amplitude`
inputs as optional *override lists*: you can provide only the fields you
want to change, and rely on defaults for everything else.

If you already use `fmridesign`, prefer
[`lss_sbhm_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm_design.md)
to avoid manually assembling an OASIS `design_spec`.

## See also

[`lss_sbhm_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm_design.md),
[`sbhm_prepass()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_prepass.md),
[`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  library(fmrihrf)
  set.seed(3)
  Tlen <- 180; V <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(exp(-seq(0, 30, length.out = Tlen)/5),
             exp(-seq(0, 30, length.out = Tlen)/7))
  sbhm <- sbhm_build(library_H = H, r = 4, sframe = sframe, normalize = TRUE)
  onsets <- seq(8, 140, by = 12)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
  Xr <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  alpha_true <- rnorm(ncol(sbhm$B))
  Y <- matrix(rnorm(Tlen*V, sd = .6), Tlen, V)
  Y[,1] <- Y[,1] + Xr %*% alpha_true
  out <- lss_sbhm(Y, sbhm, design_spec)
  out2 <- lss_sbhm(Y, sbhm, design_spec,
                  match = list(topK = 3, soft_blend = TRUE),
                  return = "amplitude")
  names(out)
} # }
```
