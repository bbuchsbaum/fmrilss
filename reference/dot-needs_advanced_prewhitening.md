# Check if advanced prewhitening is needed

Determines whether to use fmriAR (advanced) or keep simple AR(1) based
on the requested features.

## Usage

``` r
.needs_advanced_prewhitening(prewhiten)
```

## Arguments

- prewhiten:

  List of prewhitening options

## Value

Logical: TRUE if fmriAR is needed, FALSE for simple AR(1)
