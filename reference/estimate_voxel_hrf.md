# Estimate Voxel-wise HRF Basis Coefficients

Fits a GLM to estimate HRF basis coefficients for every voxel.

## Usage

``` r
estimate_voxel_hrf(Y, events, basis, nuisance_regs = NULL)
```

## Arguments

- Y:

  Numeric matrix of BOLD data (time x voxels).

- events:

  Data frame with `onset`, `duration` and `condition` columns.

- basis:

  HRF object from the `fmrihrf` package.

- nuisance_regs:

  Optional numeric matrix of nuisance regressors.

## Value

A [VoxelHRF](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
object containing at least:

- coefficients:

  Matrix of HRF basis coefficients.

- basis:

  The HRF basis object used.

- conditions:

  Character vector of modeled conditions.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1,
                     condition = "A")
basis <- fmrihrf::hrf_gamma()
sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
times <- fmrihrf::samples(sframe, global = TRUE)
rset <- fmrihrf::regressor_set(onsets = events$onset,
                               fac = factor(1:nrow(events)),
                               hrf = basis, duration = events$duration,
                               span = 30)
X <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
Y <- X %*% coef + Y * 0.1
est <- estimate_voxel_hrf(Y, events, basis)
str(est)
} # }
```
