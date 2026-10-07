# Estimate Voxel-wise HRF Basis Coefficients

Fits a common-amplitude GLM to estimate one pooled HRF shape for every
voxel. The HRF basis and response are residualized against the complete
fixed and nuisance span before fitting. Estimation fails when that span
makes the HRF basis unidentifiable or rank deficient.

## Usage

``` r
estimate_voxel_hrf(
  Y,
  events,
  basis,
  nuisance_regs = NULL,
  sframe = NULL,
  fixed_regs = NULL
)
```

## Arguments

- Y:

  Numeric matrix of BOLD data (time x voxels).

- events:

  Data frame with `onset`, `duration` and `condition` columns. For a
  multi-run sampling frame, a `run` column is required and onsets are
  physical times relative to that run.

- basis:

  HRF object from the `fmrihrf` package.

- nuisance_regs:

  Optional finite numeric matrix of nuisance regressors.

- sframe:

  Explicit `fmrihrf` sampling frame defining scan times and TR.
  Required; onset and duration values use its physical-time units.

- fixed_regs:

  Optional finite numeric matrix of fixed/common regressors. An
  intercept is added when it is not already in their span.

## Value

A [VoxelHRF](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
object containing at least:

- coefficients:

  Matrix of positive-peak-normalized HRF shape weights with one row per
  basis function and one column per voxel.

- amplitude_scale:

  The signed scale removed from each raw pooled-fit coefficient column.

- basis:

  The HRF basis object used.

- conditions:

  Observed event labels. Labels are metadata: all events are pooled into
  one shape per voxel.

- condition_pooling:

  The literal string "all-events".

## Examples

``` r
# \donttest{
set.seed(1)
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1,
                     condition = "A")
basis <- fmrihrf::HRF_SPMG1
sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
times <- fmrihrf::samples(sframe, global = TRUE)
rset <- fmrihrf::regressor_set(onsets = events$onset,
                               fac = factor(rep("all events", nrow(events))),
                               hrf = basis, duration = events$duration,
                               span = 30)
X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
Y <- X %*% coef + Y * 0.1
est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
str(est)
#> List of 8
#>  $ coefficients     : num [1, 1:2] 5.7 5.7
#>  $ amplitude_scale  : num [1:2] -0.1018 0.0121
#>  $ basis            :function (t)  
#>   ..- attr(*, "class")= chr [1:2] "HRF" "function"
#>   ..- attr(*, "name")= chr "SPMG1"
#>   ..- attr(*, "nbasis")= int 1
#>   ..- attr(*, "span")= num 24
#>   ..- attr(*, "param_names")= chr [1:3] "P1" "P2" "A1"
#>   ..- attr(*, "params")=List of 3
#>   .. ..$ P1: num 5
#>   .. ..$ P2: num 15
#>   .. ..$ A1: num 0.00833
#>  $ conditions       : chr "A"
#>  $ sframe           :List of 4
#>   ..$ blocklens : int 50
#>   ..$ TR        : num 1
#>   ..$ start_time: num 0.5
#>   ..$ precision : num 0.1
#>   ..- attr(*, "class")= chr "sampling_frame"
#>  $ condition_pooling: chr "all-events"
#>  $ normalization    : chr "positive-peak"
#>  $ coefficient_units: chr "unit-peak HRF shape weights"
#>  - attr(*, "class")= chr "VoxelHRF"
# }
```
