# Wrap a Learned Basis as an HRF (SBHM HRF)

Convert a shared time basis matrix `B` (T×r) into an
[`fmrihrf::HRF`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF-class.html)
object so it can be used directly in OASIS design construction. The HRF
returns the r basis columns evaluated at arbitrary times by
piecewise-linear interpolation on the provided time grid.

## Usage

``` r
sbhm_hrf(B, tgrid, span)
```

## Arguments

- B:

  A numeric matrix (T×r) with orthonormal columns (shared basis).

- tgrid:

  Numeric vector of length T giving the global times (seconds)
  corresponding to the rows of `B`.

- span:

  Numeric HRF span passed to
  [`fmrihrf::HRF()`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF-class.html)
  metadata.

## Value

An
[`fmrihrf::HRF`](https://bbuchsbaum.github.io/fmrihrf/reference/HRF-class.html)
object with `nbasis = ncol(B)`.

## Examples

``` r
if (FALSE) { # \dontrun{
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
} # }
```
