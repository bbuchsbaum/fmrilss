# SBHM Pipeline with fmridesign Models

Run the SBHM end-to-end pipeline using fmridesign's `event_model` and
optional `baseline_model`, mirroring the convenience of
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
but producing SBHM shape coefficients and (optionally) scalar
event-design trial coefficients in the component named `amplitude`.

## Usage

``` r
lss_sbhm_design(
  Y,
  sbhm,
  event_model,
  baseline_model = NULL,
  prewhiten = NULL,
  prepass = list(),
  match = list(),
  oasis = list(),
  amplitude = list(),
  return = c("amplitude", "coefficients", "both"),
  validate = TRUE,
  ...
)
```

## Arguments

- Y:

  Numeric matrix T×V of fMRI time series (timepoints × voxels).

- sbhm:

  SBHM object as returned by
  [`sbhm_build()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_build.md).

- event_model:

  An `event_model` from
  [`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
  defining the trial structure (typically created with
  [`trialwise()`](https://bbuchsbaum.github.io/fmridesign/reference/trialwise.html)).

- baseline_model:

  Optional `baseline_model` from
  [`fmridesign::baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html).
  Its drift, block, and nuisance terms are projected out as confounds.

- prewhiten:

  Optional prewhitening options (see
  [`?lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md)).

- prepass:

  Optional list forwarded to
  [`sbhm_prepass()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_prepass.md).

- match:

  Optional list forwarded to
  [`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md).

- oasis:

  Optional SBHM OASIS override list. Supported fields are `ridge_mode`,
  `ridge_x`, `ridge_b`, and `block_cols`; false-valued `return_se` and
  `return_diag` are accepted, while true values fail because SBHM
  uncertainty/design diagnostics are not exposed. The basis rank, trial
  count, intercept policy, and trial/basis map are fixed internally from
  `sbhm` and the event model.

- amplitude:

  Amplitude options (see
  [`?lss_sbhm`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md)).

- return:

  One of `"amplitude"`, `"coefficients"`, or `"both"`.

- validate:

  Logical; when TRUE, performs basic checks (sampling frame
  compatibility, temporal alignment) analogous to
  [`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md).

- ...:

  Reserved for future use.

## Value

Same return contract as
[`lss_sbhm()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md).

## Details

This function wraps
[`lss_sbhm()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md)
by locating the unique trialwise event term, rebuilding it with the SBHM
basis, and retaining run identities and heterogeneous durations.
Non-trial event terms and `baseline_model` columns enter the fixed
nuisance span. Multi-run regressors are constructed within runs so HRF
tails cannot cross run boundaries.

## Examples

``` r
# \donttest{
  library(fmridesign)
  sframe <- fmrihrf::sampling_frame(blocklens = 60, TR = 1)
  trials <- data.frame(onset = c(5, 20, 35), run = 1)
  emod <- event_model(onset ~ trialwise(basis = "spmg1"), data = trials,
                      block = ~run, sampling_frame = sframe)
  times <- fmrihrf::samples(sframe, global = TRUE)
  H <- cbind(exp(-times / 5), exp(-times / 7))
  sbhm <- sbhm_build(library_H = H, r = 2, sframe = sframe,
                     normalize = TRUE, baseline = NULL)
  Y <- matrix(rnorm(60 * 4), 60, 4)
  out <- lss_sbhm_design(Y, sbhm, emod)
# }
```
