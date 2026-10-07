# LSS Analysis with fmridesign Objects

Perform Least Squares Separate (LSS) analysis using event_model and
baseline_model objects from the fmridesign package. This provides a
streamlined interface for complex designs with multi-condition,
parametric modulators, and structured nuisance handling.

## Usage

``` r
lss_design(
  Y,
  event_model,
  baseline_model = NULL,
  method = "oasis",
  oasis = list(),
  prewhiten = NULL,
  blockids = NULL,
  validate = TRUE,
  ...,
  diagnostics = FALSE
)
```

## Arguments

- Y:

  Numeric matrix of fMRI data (timepoints × voxels).

- event_model:

  An event_model object from
  [`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html).
  It must contain exactly one
  [`trialwise()`](https://bbuchsbaum.github.io/fmridesign/reference/trialwise.html)
  target term. Any additional condition-level event terms are included
  as common fixed regressors.

- baseline_model:

  Optional baseline_model object from
  [`fmridesign::baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html).
  Defines drift correction, block intercepts, and nuisance regressors.
  If NULL, basic baseline intercepts are auto-injected: per-run
  intercepts derived from `blockids` (or the sampling frame) are used to
  ensure proper baseline modeling.

- method:

  LSS method to use. All methods accepted by
  [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) are
  supported for a one-basis trialwise event model. Multi-basis event
  models require `method = "oasis"`.

- oasis:

  List of OASIS-specific options: ridge regularization (`ridge_x`,
  `ridge_b`, `ridge_mode`), standard errors (`return_se`), etc. See
  [`oasis_options`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
  and the Details section of
  [`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) for the
  full list. Note: `design_spec` must not be supplied when providing
  event_model and is rejected at the adapter boundary; `oasis$whiten` is
  deprecated — use `prewhiten` instead.

- prewhiten:

  Optional prewhitening specification as a list (or `NULL` for no
  whitening). Controls temporal autocorrelation correction via the
  fmriAR package. Key fields: `method` (`"ar"`, `"arma"`, `"none"`), `p`
  (AR order or `"auto"`), `pooling` (`"global"`, `"voxel"`, `"run"`,
  `"parcel"`), and `runs`/`parcels`. See
  [`prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
  and [`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) for
  full details and examples.

- blockids:

  Optional complete exact-integer block/run identifiers, one per scan.
  If NULL, run intercepts are derived from the sampling frame.

- validate:

  Logical. If TRUE (default), performs validation checks on design
  compatibility, collinearity, and temporal alignment.

- ...:

  Additional arguments passed to the underlying LSS method.

- diagnostics:

  Logical. If TRUE, attach a `diagnostics` attribute with trial/basis
  and voxel identities, retained/excluded counts, non-finite output
  locations, and rank/conditioning of the assembled input design. The
  design diagnostics precede whitening and describe the full design, not
  each trial-specific LSS model. No trials or voxels are silently
  dropped; invalid inputs raise errors. This option does not change the
  return type.

## Value

Normally a trial-by-voxel beta matrix, or a (trial × basis)-by-voxel
matrix for multi-basis HRFs. When OASIS `return_diag = TRUE` or
`return_se = TRUE`, returns `list(beta, diag?, se?)` with the
matrix/list shapes documented in
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md).
Multi-basis rows are always trial-major with basis varying within trial.
Adapter metadata are attached to the returned object; multi-basis beta
and SE matrices retain the canonical `trial_basis_map`. With active
estimated prewhitening, `attr(result, "whiten_plan")` records the fitted
`fmriAR_plan`.

## Details

**Design Specification:**

The `event_model` should typically use
[`trialwise()`](https://bbuchsbaum.github.io/fmridesign/reference/trialwise.html)
for LSS:


      emod <- event_model(onset ~ trialwise(basis = "spmg1"),
                          data = events,
                          block = ~run,
                          sampling_frame = sframe)

Non-trial event terms are supported as common fixed regressors, for
example a parametric modulator alongside the unique trialwise target
term:


      emod <- event_model(onset ~ trialwise(basis = "spmg1") + hrf(RT),
                          data = events,
                          block = ~run,
                          sampling_frame = sframe)

**Baseline Model:**

If provided, baseline_model components are mapped as follows:

- `drift` and `block` terms → Z parameter (fixed effects)

- `nuisance` term → Nuisance parameter (confounds)

**Multi-Run Handling:**

Both event_model and baseline_model must use the same `sampling_frame`.
Event onsets should be run-relative (resetting to 0 each run) as per
fmridesign convention; conversion to global time and run boundaries in
the design are handled automatically. A joint call still uses one OASIS
other-trial aggregate across all runs, so its coefficients are cross-run
pooled. Fit each run separately when run-local LSS coefficients are
required.

**Prewhitening:**

Use the `prewhiten` parameter (not the `oasis` list) for temporal
autocorrelation correction. For active multi-run whitening, this adapter
infers the scan-level run segmentation from the sampling frame. An
explicit `prewhiten$runs` vector is allowed only when it encodes those
same boundaries. The event model's `blockids` usually has one value per
event and is not the required scan-level vector. Residual-autocovariance
bias correction remains an explicit opt-in: `lss_design()` does not
silently populate `prewhiten$design`. If supplied, it must be the
assembled residual-forming design, including trialwise, fixed event,
baseline, nuisance, and intercept columns as applicable. This preserves
fmriAR's requirement that the correction design be the one that actually
produced the residuals. See
[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) and
[`prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
for full details.

**Validation:**

When `validate = TRUE`, the function checks:

- Temporal alignment: nrow(Y) matches total scans in sampling_frame

- Collinearity: emits a warning when the full assembled design has a
  condition number above 30 and all effective ridge penalties are zero.
  The condition number is scale-dependent and is not returned.

- Compatibility: event_model and baseline_model use same sampling_frame

## See also

[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) for the
traditional matrix-based interface,
[`fmridesign::event_model`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
for event model creation,
[`fmridesign::baseline_model`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html)
for baseline model creation

## Examples

``` r
# \donttest{
library(fmridesign)
#> 
#> Attaching package: ‘fmridesign’
#> The following objects are masked from ‘package:stats’:
#> 
#>     contrasts, convolve
library(fmrihrf)
#> 
#> Attaching package: ‘fmrihrf’
#> The following object is masked from ‘package:stats’:
#> 
#>     deriv

sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)

trials <- data.frame(
  onset = c(10, 30, 50, 70, 90, 110,
            10, 30, 50, 70, 90, 110),
  run = rep(1:2, each = 6)
)

emod <- event_model(
  onset ~ trialwise(basis = "spmg1"),
  data = trials,
  block = ~run,
  sampling_frame = sframe
)

motion <- list(
  matrix(rnorm(150 * 6), 150, 6),
  matrix(rnorm(150 * 6), 150, 6)
)
bmodel <- baseline_model(
  basis = "bs",
  degree = 5,
  sframe = sframe,
  nuisance_list = motion
)

Y <- matrix(rnorm(300 * 1000), 300, 1000)
beta <- lss_design(Y, emod, bmodel, method = "oasis")

dim(beta)
#> [1]   12 1000
# }
```
