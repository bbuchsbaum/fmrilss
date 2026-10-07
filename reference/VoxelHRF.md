# VoxelHRF object

Simple list-based S3 class returned by `estimate_voxel_hrf` containing
voxel-wise HRF basis coefficients and related metadata.

## Value

No value itself. This topic documents the structure returned by
[`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md).

## Stored fields

`coefficients` contains one normalized HRF-shape column per voxel;
`amplitude_scale` records the removed positive-peak scale; `basis`
stores the HRF basis; `conditions` records observed labels while
`condition_pooling` states that all events estimate one pooled shape;
`sframe` preserves physical scan timing; and `normalization` plus
`coefficient_units` state the coefficient scale. Voxel names, when
supplied on `Y`, identify coefficient columns and amplitude-scale
entries.

## Examples

``` r
# \donttest{
Y <- matrix(rnorm(100), 50, 2)
events <- data.frame(onset = c(5, 25), duration = 1, condition = "A")
basis <- fmrihrf::HRF_SPMG1
sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
est <- estimate_voxel_hrf(Y, events, basis, sframe = sframe)
class(est)
#> [1] "VoxelHRF"
# }
```
