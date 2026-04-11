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
  ...
)
```

## Arguments

- Y:

  Numeric matrix of fMRI data (timepoints × voxels).

- event_model:

  An event_model object from
  [`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html).
  This defines the trial-wise or condition-wise task design. For LSS,
  typically created with
  [`trialwise()`](https://bbuchsbaum.github.io/fmridesign/reference/trialwise.html)
  to generate one regressor per trial.

- baseline_model:

  Optional baseline_model object from
  [`fmridesign::baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html).
  Defines drift correction, block intercepts, and nuisance regressors.
  If NULL, basic baseline intercepts are auto-injected: per-run
  intercepts derived from `blockids` (or the sampling frame) are used to
  ensure proper baseline modeling.

- method:

  LSS method to use. Currently only "oasis" is supported for event_model
  integration.

- oasis:

  List of OASIS-specific options: ridge regularization (`ridge_x`,
  `ridge_b`, `ridge_mode`), standard errors (`return_se`), etc. See
  [`oasis_options`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
  and the Details section of
  [`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) for the
  full list. Note: `design_spec` is not used when providing event_model,
  and `oasis$whiten` is deprecated — use `prewhiten` instead.

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

  Optional block/run identifiers for event_model. If NULL, extracted
  from event_model\$blockids.

- validate:

  Logical. If TRUE (default), performs validation checks on design
  compatibility, collinearity, and temporal alignment.

- ...:

  Additional arguments passed to the underlying LSS method.

## Value

Matrix of trial-wise beta estimates (trials × voxels), or (trials ×
basis_functions) × voxels for multi-basis HRFs.

## Details

**Design Specification:**

The `event_model` should typically use
[`trialwise()`](https://bbuchsbaum.github.io/fmridesign/reference/trialwise.html)
for LSS:

      emod <- event_model(onset ~ trialwise(basis = "spmg1"),
                          data = events,
                          block = ~run,
                          sampling_frame = sframe)

For factorial designs (e.g., estimating condition-level betas
separately):

      emod <- event_model(onset ~ hrf(condition),
                          data = events,
                          block = ~run,
                          sampling_frame = sframe)

**Baseline Model:**

If provided, baseline_model components are mapped as follows:

- `drift` and `block` terms → Z parameter (fixed effects)

- `nuisance` term → Nuisance parameter (confounds)

**Multi-Run Handling:**

Both event_model and baseline_model must use the same `sampling_frame`.
Run structure is automatically respected. Event onsets should be
run-relative (resetting to 0 each run) as per fmridesign convention -
conversion to global time is handled automatically.

**Prewhitening:**

Use the `prewhiten` parameter (not the `oasis` list) for temporal
autocorrelation correction. For multi-run data, pass
`prewhiten = list(method = "ar", p = 1, pooling = "run", runs = blockids)`
so that whitening respects run boundaries. See
[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) and
[`prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
for full details.

**Validation:**

When `validate = TRUE`, the function checks:

- Temporal alignment: nrow(Y) matches total scans in sampling_frame

- Collinearity: Design matrix condition number \< 30 (suppressed when
  ridge is already configured via `oasis$ridge_x` or `oasis$ridge_b`)

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
if (FALSE) { # \dontrun{
library(fmridesign)
library(fmrihrf)

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
} # }
```
