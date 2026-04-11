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
  compute_residuals = TRUE
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

  Optional run identifiers.

- parcels:

  Optional parcel ids.

- exact_first:

  `"ar1"` or `"none"`.

- compute_residuals:

  Logical.

## Value

A list with class `"fmrilss_prewhiten_options"`.

## Examples

``` r
prewhiten_options(method = "ar", p = 1, pooling = "run")
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
#> NULL
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
#> attr(,"class")
#> [1] "fmrilss_prewhiten_options" "list"                     
```
