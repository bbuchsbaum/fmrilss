# Convenience wrappers for modern `lss()` usage

These functions provide a modern-signature entry point for methods that
historically required a block-design (`bdes`) object.

## Value

No value itself. This topic groups
[`lss_naive_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_naive_fit.md)
and
[`lss_optimized_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_optimized_fit.md).

## Examples

``` r
set.seed(1)
Y <- matrix(rnorm(16), 8, 2)
X <- matrix(0, 8, 2)
X[2:3, 1] <- 1
X[5:6, 2] <- 1
lss_naive_fit(Y, X)
#>            Voxel_1     Voxel_2
#> Trial_1 -0.8746378  0.09179092
#> Trial_2 -0.7941255 -1.92937571
```
