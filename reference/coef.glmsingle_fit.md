# Single-trial betas from a glmsingle fit

Single-trial betas from a glmsingle fit

## Usage

``` r
# S3 method for class 'glmsingle_fit'
coef(object, type = c("d", "c", "b", "a"), ...)
```

## Arguments

- object:

  A `glmsingle_fit`.

- type:

  Model type: `"d"` (fractional ridge, default), `"c"` (GLMdenoise),
  `"b"` (HRF library) or `"a"` (ON-OFF, one value per voxel).

- ...:

  Unused.

## Value

Trials x voxels matrix (a vector for type `"a"`).
