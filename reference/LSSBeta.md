# LSSBeta object

Simple list-based S3 class returned by `lss_with_hrf` containing
trial-wise beta estimates.
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) materializes
C++-backed estimates as a dense matrix while preserving voxel names,
sampling-frame, normalization, and unit metadata.

## Value

No value itself. This topic documents the object returned by
`lss_with_hrf(..., engine = "C++")`.

## Stored fields

`betas` is the file-backed trial-by-voxel matrix; `dimnames` preserves
trial and voxel identity; `sframe`, `normalization`, `units`,
`degenerate`, and `event_amplitude` and `event_duration` define timing
and coefficient interpretation; and `engine_requested`, `engine_used`,
and `chunk_size` record execution. Use
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) for the documented
dense extraction path; these metadata are retained as matrix attributes.

## Examples

``` r
# \donttest{
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1, condition = "A")
basis <- fmrihrf::HRF_SPMG1
sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
fit <- lss_with_hrf(Y, events, est, engine = "C++", verbose = FALSE)
class(fit)
#> [1] "LSSBeta"
# }
```
