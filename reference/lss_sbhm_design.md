# SBHM Pipeline with fmridesign Models

Run the SBHM end-to-end pipeline using fmridesign's \`event_model\` and
optional \`baseline_model\`, mirroring the convenience of
\`lss_design()\` but producing SBHM coefficients and (optionally) scalar
amplitudes.

## Usage

``` r
lss_sbhm_design(
  Y,
  sbhm,
  event_model,
  baseline_model = NULL,
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
  return = c("amplitude", "coefficients", "both"),
  validate = TRUE,
  ...
)
```

## Arguments

- Y:

  Numeric matrix T×V of fMRI time series (timepoints × voxels).

- sbhm:

  SBHM object as returned by \`sbhm_build()\`.

- event_model:

  An \`event_model\` from \`fmridesign::event_model()\` defining the
  trial structure (typically created with \`trialwise()\`).

- baseline_model:

  Optional \`baseline_model\` from \`fmridesign::baseline_model()\`. Its
  drift, block, and nuisance terms are projected out as confounds.

- prewhiten:

  Optional prewhitening options (see \`?lss\`).

- prepass:

  Optional list forwarded to \`sbhm_prepass()\`.

- match:

  Optional list forwarded to \`sbhm_match()\`.

- oasis:

  Optional list forwarded to \`lss(..., method = "oasis")\`. \`K\`
  defaults to \`ncol(sbhm\$B)\`.

- amplitude:

  Amplitude options (see \`?lss_sbhm\`).

- return:

  One of \`"amplitude"\`, \`"coefficients"\`, or \`"both"\`.

- validate:

  Logical; when TRUE, performs basic checks (sampling frame
  compatibility, temporal alignment) analogous to \`lss_design()\`.

- ...:

  Reserved for future use.

## Value

Same return contract as \`lss_sbhm()\`.

## Details

This function wraps \`lss_sbhm()\` by converting the \`event_model\`
into an OASIS \`design_spec\` that uses the SBHM basis HRF, and by
mapping \`baseline_model\` terms to nuisance regressors for projection.

## Examples

``` r
if (FALSE) { # \dontrun{
  library(fmridesign)
  sframe <- fmrihrf::sampling_frame(blocklens = c(150, 150), TR = 2)
  trials <- data.frame(onset = c(10,30,50, 10,30,50), run = rep(1:2, each=3))
  emod <- event_model(onset ~ trialwise(basis = "spmg1"), data = trials,
                      block = ~run, sampling_frame = sframe)
  Y <- matrix(rnorm(300*100), 300, 100)
  out <- lss_sbhm_design(Y, sbhm, emod)
} # }
```
