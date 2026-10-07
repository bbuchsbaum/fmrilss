# Construct prewhitening options

Convenience constructor for the `prewhiten=` list accepted by
[`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md).

## Usage

``` r
prewhiten_options(
  method = c("none", "ar", "arma"),
  p = "auto",
  q = 0L,
  p_max = 6L,
  pooling = c("global", "voxel", "run", "parcel"),
  runs = NULL,
  parcels = NULL,
  exact_first = c("ar1", "none"),
  compute_residuals = TRUE,
  design = NULL,
  acvf_correction = NULL,
  correction_max_lag = 25L
)
```

## Arguments

- method:

  `"none"`, `"ar"`, or `"arma"`.

- p:

  AR order or `"auto"`.

- q:

  MA order for ARMA.

- p_max:

  Maximum AR order when `p="auto"`.

- pooling:

  `"global"`, `"voxel"`, `"run"`, or `"parcel"`.

- runs:

  Optional complete exact-integer run identifiers. Required for
  `pooling = "run"`; execution requires one value per timepoint.

- parcels:

  Optional complete exact-integer parcel identifiers. Required for
  `pooling = "parcel"`; execution requires one value per voxel.

- exact_first:

  `"ar1"` or `"none"`.

- compute_residuals:

  Logical.

- design:

  Optional numeric design matrix whose projection produced the residuals
  used to estimate the noise model. Supplying it opts in to fmriAR's
  residual-autocovariance bias correction. It must represent the same
  column space as the design used for residualization.

- acvf_correction:

  Optional bias-correction matrix, or list of matrices, produced by
  [`fmriAR::acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.html).
  This is a cached alternative to `design`; the two options are mutually
  exclusive.

- correction_max_lag:

  Positive integer lag budget used when `design` is supplied. See
  [`fmriAR::fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.html)
  for the computational and filtering requirements of the correction.

## Value

A list with class `"fmrilss_prewhiten_options"`.

## Examples

``` r
prewhiten_options(method = "ar", p = 1, pooling = "run",
                  runs = rep(1:2, each = 50))
#> $method
#> [1] "ar"
#> 
#> $p
#> [1] 1
#> 
#> $q
#> [1] 0
#> 
#> $p_max
#> [1] 6
#> 
#> $pooling
#> [1] "run"
#> 
#> $runs
#>   [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#>  [38] 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#>  [75] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
#> 
#> $parcels
#> NULL
#> 
#> $exact_first
#> [1] "ar1"
#> 
#> $compute_residuals
#> [1] TRUE
#> 
#> $design
#> NULL
#> 
#> $acvf_correction
#> NULL
#> 
#> $correction_max_lag
#> [1] 25
#> 
#> attr(,"class")
#> [1] "fmrilss_prewhiten_options" "list"                     
```
