# LSSBeta object

Simple list-based S3 class returned by `lss_with_hrf` containing
trial-wise beta estimates.

## Value

No value itself. This topic documents the object returned by
`lss_with_hrf(..., engine = "C++")`.

## Examples

``` r
if (FALSE) { # \dontrun{
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1, condition = "A")
basis <- fmrihrf::HRF_SPMG1
est <- estimate_voxel_hrf(Y, events, basis)
fit <- lss_with_hrf(Y, events, est, engine = "C++", verbose = FALSE)
class(fit)
} # }
```
