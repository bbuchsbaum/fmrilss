# fmridesign front end for glmsingle()

Fits
[`glmsingle()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle.md)
using an
[`fmridesign::event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html)
to describe the trials. The event model must contain one event term with
a single factor whose levels are the conditions (repeated conditions
drive cross-validation). Onsets must lie on the TR grid.

## Usage

``` r
glmsingle_design(Y, event_model, baseline_model = NULL, stimdur = NULL, ...)
```

## Arguments

- Y:

  Time x voxel data matrix covering all runs of the sampling frame.

- event_model:

  An `event_model` from fmridesign.

- baseline_model:

  Optional `baseline_model` from fmridesign.

- stimdur:

  Trial duration in seconds. Default: the event durations, which must
  then be a single positive value.

- ...:

  Further arguments passed to
  [`glmsingle()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle.md).

## Value

A `glmsingle_fit`; see
[`glmsingle()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle.md).

## Details

A `baseline_model` contributes its nuisance term (e.g. motion) as
`extra_regressors`. Its drift and block terms are not used: GLMsingle
models drift with its own per-run polynomials (`max_poly_deg`).

## Examples

``` r
set.seed(1)
sf <- fmrihrf::sampling_frame(c(80L, 80L), TR = 1)
events <- data.frame(onset = rep(c(5, 20, 35, 50), 2),
                     run = rep(1:2, each = 4),
                     stimulus = factor(rep(c("A", "B"), 4)))
em <- fmridesign::event_model(
  onset ~ fmridesign::hrf(stimulus), data = events, block = ~run,
  sampling_frame = sf, durations = rep(2, nrow(events)))
Y <- matrix(100 + rnorm(160 * 4), 160, 4)
fit <- glmsingle_design(Y, em, want_glmdenoise = FALSE,
                        want_fracridge = FALSE, verbose = FALSE)
dim(coef(fit, type = "b"))
#> [1] 8 4
```
