# VoxelHRF object

Simple list-based S3 class returned by `estimate_voxel_hrf` containing
voxel-wise HRF basis coefficients and related metadata.

## Value

No value itself. This topic documents the structure returned by
[`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md).

## Examples

``` r
if (FALSE) { # \dontrun{
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1, condition = "A")
basis <- fmrihrf::HRF_SPMG1
est <- estimate_voxel_hrf(Y, events, basis)
class(est)
} # }
```
