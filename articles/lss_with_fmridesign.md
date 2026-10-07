# Run-aware LSS with fmridesign

[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
is the run-aware adapter between `fmridesign` event and baseline models
and the OASIS estimator. Read this after
[`vignette("fmrilss")`](https://bbuchsbaum.github.io/fmrilss/articles/fmrilss.md)
and
[`vignette("oasis_method")`](https://bbuchsbaum.github.io/fmrilss/articles/oasis_method.md).
The adapter is useful when event tables have run-relative onsets, the
baseline is structured by run, or non-trial event terms must remain in
the common model.

This article requires the suggested `fmridesign` package. Load the three
packages used below before copying the workflow into a fresh session.

``` r

suppressPackageStartupMessages({
  library(fmrilss)
  library(fmridesign)
  library(fmrihrf)
})
```

## Know the adapter contract

| Question | Answer |
|:---|:---|
| Which columns become LSS targets? | Exactly one trialwise event term |
| What happens to other event terms? | They are fixed common regressors, not extra trials |
| How are baseline terms mapped? | drift + block go to Z; nuisance is projected as Nuisance |
| What does multi-basis output mean? | K coefficient rows per trial, in trial-major order |

The fmridesign-to-LSS mapping used by lss_design(). {.table}

The function requires one—and only one—trialwise event term. A
parametric or condition-level event term can accompany it, but that term
is part of every trial’s common design. It does not create another set
of trial targets.

## Start from a two-run event table

Onsets are relative to the start of each run. Here, the same five onset
times occur in both runs; the `run` column supplies their distinct
identities.

``` r

set.seed(20260821)

run_lengths <- c(110L, 130L)
TR <- 1
sframe <- sampling_frame(blocklens = run_lengths, TR = TR)

events <- data.frame(
  event_id = seq_len(10L),
  onset = rep(c(10, 30, 50, 70, 90), times = 2),
  run = rep(1:2, each = 5),
  RT = seq(0.4, 0.9, length.out = 10)
)
events$RT_c <- events$RT - mean(events$RT)
events
#>    event_id onset run        RT        RT_c
#> 1         1    10   1 0.4000000 -0.25000000
#> 2         2    30   1 0.4555556 -0.19444444
#> 3         3    50   1 0.5111111 -0.13888889
#> 4         4    70   1 0.5666667 -0.08333333
#> 5         5    90   1 0.6222222 -0.02777778
#> 6         6    10   2 0.6777778  0.02777778
#> 7         7    30   2 0.7333333  0.08333333
#> 8         8    50   2 0.7888889  0.13888889
#> 9         9    70   2 0.8444444  0.19444444
#> 10       10    90   2 0.9000000  0.25000000
```

Build one trialwise target term and one reaction-time event term. The
latter is a common regressor: it adjusts every LSS model for the
specified amplitude modulation but is not reported as a trial beta.

``` r

emod <- event_model(
  onset ~ trialwise(basis = "spmg1") + hrf(RT_c),
  data = events,
  block = ~run,
  sampling_frame = sframe
)

event_dm <- design_matrix(emod)
event_meta <- attr(event_dm, "col_metadata")
table(event_meta$term_tag)
#> 
#>  RT_c trial 
#>     1    10
```

The ten `trial` columns are the LSS targets; the one `RT_c` column is
fixed.
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
uses the design metadata and stable column identities rather than
assuming that columns arrived in a particular order.

## Verify the run boundary

The target matrix is block diagonal by run: trials 1–5 have no energy in
run 2, and trials 6–10 have no energy in run 1.

| Check                     | Maximum_absolute_design_value |
|:--------------------------|------------------------------:|
| run-1 trials inside run 2 |                             0 |
| run-2 trials inside run 1 |                             0 |

Cross-run leakage in the trialwise design. {.table}

![Heatmap with activity for trials 1 through 5 only before scan 110 and
trials 6 through 10 only after scan 110; no trial regressor crosses the
run
boundary.](lss_with_fmridesign_files/figure-html/run-boundary-plot-1.png)

The trialwise design is block diagonal by run; the vertical line marks
the boundary after scan 110.

## Add a structured baseline

[`baseline_model()`](https://bbuchsbaum.github.io/fmridesign/reference/baseline_model.html)
keeps run intercepts, drift, and nuisance inputs attached to the same
sampling frame. The example uses two synthetic motion columns per run.

``` r

motion <- list(
  matrix(rnorm(run_lengths[1] * 2, sd = 0.15), run_lengths[1], 2),
  matrix(rnorm(run_lengths[2] * 2, sd = 0.15), run_lengths[2], 2)
)

bmodel <- baseline_model(
  basis = "poly",
  degree = 1,
  sframe = sframe,
  intercept = "runwise",
  nuisance_list = motion
)

baseline_terms <- term_matrices(bmodel)
vapply(baseline_terms, ncol, integer(1))
#>    drift    block nuisance 
#>        2        2        4
```

The adapter sends `drift` and `block` to `Z`, sends `nuisance` to
`Nuisance`, and appends the fixed `RT_c` event column to `Z`. If no
baseline model is supplied, it inserts run-wise intercepts instead.

## Fit and check every trial against a full GLM

The simulation below is correctly specified for every LSS model: target
trials share a coefficient within each voxel, while drift, run
intercepts, motion, and the RT modulator all have their own common
coefficients. This is an adapter court, not evidence about recovery in
arbitrary designs.

Use explicit zero ridge when the goal is ordinary LSS. The default OASIS
configuration is penalized.

``` r

fit <- lss_design(
  Y,
  emod,
  bmodel,
  method = "oasis",
  oasis = oasis_options(
    ridge_mode = "absolute",
    ridge_x = 0,
    ridge_b = 0
  ),
  validate = FALSE
)
dim(fit)
#> [1] 10  4
fit[1:4, ]
#>                                           Voxel_1    Voxel_2    Voxel_3
#> trial_.trial_factor.length.onsets...01  0.8588198  1.7142706  0.2521402
#> trial_.trial_factor.length.onsets...02 -1.0983198  1.2285833 -3.0476933
#> trial_.trial_factor.length.onsets...03 -4.7123311  1.7288773  5.1002103
#> trial_.trial_factor.length.onsets...04  3.7497051 -0.9466139 -4.1816692
#>                                          Voxel_4
#> trial_.trial_factor.length.onsets...01  4.100424
#> trial_.trial_factor.length.onsets...02 -3.187782
#> trial_.trial_factor.length.onsets...03  0.736539
#> trial_.trial_factor.length.onsets...04  5.679483
```

`validate = FALSE` above skips the adapter’s optional full-design
condition-number warning; it does not disable the mandatory trial/basis
identity mapping. The full trialwise matrix contains `RT_c` as a
weighted sum of trial columns, so that screening matrix is rank
deficient even though each target-specific LSS model below is full rank.
The independent court checks the actual per-target designs and every
returned cell.

| Maximum_absolute_error | Minimum_target_model_rank | Target_model_columns |
|-----------------------:|--------------------------:|---------------------:|
|                      0 |                        11 |                   11 |

lss_design() versus independently assembled full GLMs. {.table}

The zero discrepancy also proves that the RT event term and baseline
nuisance span entered the model without becoming extra trial rows.

## Decide whether other-trial effects are pooled across runs

Run-aware onset placement does not make the default coefficient
estimator run-local. A single
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
call passes all target columns to one OASIS fit, whose other-trial
regressor aggregates trials from both runs. A run-1 target can therefore
change when the run-2 response changes. This is a pooled multi-run
estimand, not leakage in the event convolution.

The next mutation adds signal only in run 2. It changes the joint fit’s
run-1 coefficients, while two explicitly separate run fits leave run 1
unchanged.

| Estimator              | Maximum_change_in_run1_beta |
|:-----------------------|----------------------------:|
| one pooled two-run fit |                    13.83231 |
| two separate run fits  |                     0.00000 |

Effect of a run-2-only response mutation on run-1 coefficients. {.table}

Use one joint call only when that cross-run pooling is the intended LSS
model. For run-local coefficients, fit each run with its own event and
baseline models, then restore global event identifiers in both row names
and the trial/basis map before combining rows, as the helper above does.
Run-specific prewhitening does not by itself change the pooled
other-trial estimand.

## Preserve output identity

The result carries the event model, baseline model, sampling frame, and
an explicit trial/basis map. Use the map or row names; do not
reconstruct trial identity from a presumed source-column order.

``` r

attributes_kept <- c(
  event_model = !is.null(attr(fit, "event_model")),
  baseline_model = !is.null(attr(fit, "baseline_model")),
  sampling_frame = !is.null(attr(fit, "sampling_frame")),
  trial_basis_map = !is.null(attr(fit, "trial_basis_map"))
)
attributes_kept
#>     event_model  baseline_model  sampling_frame trial_basis_map 
#>            TRUE            TRUE            TRUE            TRUE
identity_map <- attr(fit, "trial_basis_map")
stopifnot(identical(identity_map$trial, events$event_id))
identity_map$input_event_id <- events$event_id
identity_map[1:3, c("input_event_id", "trial", "basis", "output_row")]
#>   input_event_id trial basis output_row
#> 1              1     1     1          1
#> 2              2     2     1          2
#> 3              3     3     1          3
```

## Multi-basis designs return coefficients, not amplitudes

With SPMG3,
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
detects three basis columns per trial and returns rows in trial-major,
basis-minor order. The next simulation is correctly specified in that
three-dimensional basis and is checked against independent full GLMs.

``` r

fit_spmg3 <- lss_design(
  Y_spmg3,
  emod_spmg3,
  method = "oasis",
  oasis = oasis_options(
    ridge_mode = "absolute",
    ridge_x = 0,
    ridge_b = 0
  ),
  validate = FALSE
)
c(rows = nrow(fit_spmg3), trials = n_trials, basis_dimension = K)
#>            rows          trials basis_dimension 
#>              30              10               3
rownames(fit_spmg3)[1:6]
#> [1] "trial_.trial_factor.length.onsets...01:basis_1"
#> [2] "trial_.trial_factor.length.onsets...01:basis_2"
#> [3] "trial_.trial_factor.length.onsets...01:basis_3"
#> [4] "trial_.trial_factor.length.onsets...02:basis_1"
#> [5] "trial_.trial_factor.length.onsets...02:basis_2"
#> [6] "trial_.trial_factor.length.onsets...02:basis_3"
attr(fit_spmg3, "trial_basis_map")[
  1:6, c("trial", "basis", "source_column", "output_row")
]
#>   trial basis source_column output_row
#> 1     1     1             1          1
#> 2     1     2             2          2
#> 3     1     3             3          3
#> 4     2     1             4          4
#> 5     2     2             5          5
#> 6     2     3             6          6
```

| Maximum_absolute_error |
|-----------------------:|
|                      0 |

Multi-basis lss_design() versus independent full GLMs. {.table}

Each trial now has canonical, temporal-derivative, and
dispersion-derivative coefficients. The canonical coefficient alone is
not a normalized response amplitude. Keep all three rows unless a
separately defined shape and amplitude estimand justifies a reduction.

## Ridge, standard errors, and whitening retain OASIS semantics

The adapter does not change the inference contract:

- [`oasis_options()`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
  uses fractional ridge by default.
- `return_se = TRUE` requires zero ridge, a fixed full-rank design,
  positive residual degrees of freedom, and no estimated prewhitening.
- A non-trial event term remains common whether or not ridge is used.

For multiple runs, whitening needs scan-level segmentation.
[`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
infers it from the sampling frame. The explicit vector below is
equivalent and shows the underlying contract; if supplied, it must
encode the same boundaries. The event model’s `blockids` describe
events, not scans, so do not pass that shorter vector as
`prewhiten$runs`.

``` r

scan_run <- rep(seq_along(run_lengths), times = run_lengths)

fit_whitened <- lss_design(
  Y,
  emod,
  bmodel,
  method = "oasis",
  oasis = oasis_options(),
  prewhiten = prewhiten_options(
    method = "ar",
    p = 1,
    pooling = "run",
    runs = scan_run
  )
)
```

This is a recipe, not a recommendation that AR(1) fits every dataset.
Choose the noise model from residual diagnostics. Shared-design OASIS
calls reject voxel- or parcel-specific whitening operators because those
require correspondingly voxel- or parcel-specific filtered designs.

## Validation and failure modes

With `validate = TRUE`, the adapter checks time dimensions and
sampling-frame agreement. It computes a scale-dependent condition number
for the full assembled design and emits a warning only when that number
exceeds 30 and the effective ridge penalties are all zero; it does not
return the number as a diagnostic. Treat the warning as screening, not
as a characterization of every residualized LSS target model. Mandatory
semantic checks remain active regardless of the flag: the event design
must have metadata, exactly one trialwise term, a complete trial/basis
rectangle, and stable identities.

Common failures are direct:

- `nrow(Y)` must equal `sum(blocklens(sframe))`.
- The event and baseline models must use the same sampling frame.
- Multi-basis metadata must agree with the HRF basis dimension.
- Parametric and condition-level terms are fixed regressors; they are
  not additional trialwise targets.

## Next steps

- [`vignette("voxel-wise-hrf")`](https://bbuchsbaum.github.io/fmrilss/articles/voxel-wise-hrf.md)
  — normalized voxel-specific HRF shapes
- [`vignette("sbhm")`](https://bbuchsbaum.github.io/fmrilss/articles/sbhm.md)
  — library-constrained voxel-specific HRFs
