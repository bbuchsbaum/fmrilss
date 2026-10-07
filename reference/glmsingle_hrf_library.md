# GLMsingle canonical HRF and HRF library

Reproduces GLMsingle's canonical HRF (`getcanonicalhrf`) and its library
of 20 HRFs (`getcanonicalhrflibrary`) for a stimulus of duration
`stimdur` sampled at `tr`. Each HRF is peak-normalised to 1, and the
first sample is coincident with stimulus onset.

## Usage

``` r
glmsingle_hrf_library(stimdur, tr)

glmsingle_hrf(stimdur, tr)
```

## Arguments

- stimdur:

  Stimulus duration in seconds (rounded to 0.1 s).

- tr:

  Repetition time in seconds.

## Value

`glmsingle_hrf()` returns a numeric vector; `glmsingle_hrf_library()`
returns a time x 20 matrix.

## References

Prince, J. S., et al. (2022). Improving the accuracy of single-trial
fMRI response estimates using GLMsingle. eLife, 11, e77599.

## Examples

``` r
lib <- glmsingle_hrf_library(stimdur = 3, tr = 1)
dim(lib)
#> [1] 53 20
```
