# End-to-End LSS with Shared-Basis HRF Matching (SBHM)

Orchestrates the SBHM pipeline: (1) run-safe trial design and OASIS fit
in the shared basis, (2) a voxel shape summary from the requested
source, (3) hard library matching or explicit blending, and (4) a
separate scalar-coefficient refit using the selected voxel shape.

## Usage

``` r
lss_sbhm(
  Y,
  sbhm,
  design_spec,
  Nuisance = NULL,
  prewhiten = NULL,
  prepass = list(),
  match = list(),
  oasis = list(),
  amplitude = list(),
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

  List for design construction:
  `list(sframe=..., cond=list(onsets=..., duration=0, span=...), others=list(...))`.
  The target HRF is replaced by the SBHM basis. Multi-run inputs use
  run-relative onsets and an explicit `cond$run` vector; other
  conditions use the same event fields and may supply their own HRF.

- Nuisance:

  Optional T×P nuisance regressors.

- prewhiten:

  Optional prewhitening options (see
  [`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)). Run
  labels are inferred from `design_spec$sframe` when omitted; supplied
  labels must match those sampling-frame boundaries.

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

  Optional list forwarded to `lss(..., method="oasis")`. The basis
  dimension, trial count, run-specific intercept span, and trial/basis
  map are supplied from the run-safe SBHM design.

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

- `matched_name` and `matched_idx`: named top-scoring library identities

- `margin` named length-V score differences (top1 - top2 cosine); these
  are not calibrated confidence measures

- `alpha_coords` r×V matched coordinates per voxel

- `shape_mode` and `fallback_low_conf`: named vectors recording the
  shape policy actually used

- `prepass_fallback`: named logical vector identifying voxels whose
  requested shape source fell back to the aggregate prepass

- `trial_basis_map`: complete trial/basis output identity

- `event_amplitude`, `event_duration`, and `event_run`: named
  event-design metadata

- `diag` list with `r`, `ntrials`, and `times`

Scalar amplitudes are coefficients on the supplied event design using
the selected rank-truncated library shape. They are not automatically
peak BOLD responses. SBHM standard errors are not calibrated and
`return_se=TRUE` fails explicitly for every amplitude method.

## Details

Most users should treat the `prepass`, `match`, `oasis`, and `amplitude`
inputs as optional nested *override lists*: you can provide only the
fields you want to change. Unspecified nested defaults are preserved.

If you already use `fmridesign`, prefer
[`lss_sbhm_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm_design.md)
to avoid manually assembling an OASIS `design_spec`.

## See also

[`lss_sbhm_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm_design.md),
[`sbhm_prepass()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_prepass.md),
[`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md)

## Examples

``` r
# \donttest{
  library(fmrihrf)
  set.seed(3)
  Tlen <- 180; V <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(exp(-seq(0, 30, length.out = Tlen)/5),
             exp(-seq(0, 30, length.out = Tlen)/7))
  sbhm <- sbhm_build(library_H = H, r = 2, sframe = sframe, normalize = TRUE)
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
                  match = list(topK = 2, soft_blend = TRUE),
                  return = "amplitude")
  names(out)
#>  [1] "matched_idx"       "matched_name"      "margin"           
#>  [4] "alpha_coords"      "shape_mode"        "prepass_fallback" 
#>  [7] "fallback_low_conf" "trial_basis_map"   "event_amplitude"  
#> [10] "event_duration"    "event_run"         "diag"             
#> [13] "topK_idx"          "weights"           "alpha_mode"       
#> [16] "amplitude"        
# }
```
