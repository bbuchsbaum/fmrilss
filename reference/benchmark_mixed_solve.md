# Benchmark Mixed Model Implementations

Compare performance between the standard `mixed_solve` implementation
and the optimized `mixed_solve_optimized` version.

## Usage

``` r
benchmark_mixed_solve(X, Z, K = NULL, Y, n_reps = 5)
```

## Arguments

- X:

  Fixed effects design matrix

- Z:

  Random effects design matrix

- K:

  Kinship matrix (optional, defaults to identity)

- Y:

  Response matrix (n x V)

- n_reps:

  Number of repetitions for benchmarking

## Value

Data frame with timing results

## Examples

``` r
# \donttest{
X <- matrix(rnorm(100 * 2), 100, 2)
Z <- matrix(rnorm(100 * 3), 100, 3)
Y <- matrix(rnorm(100 * 5), 100, 5)
benchmark_mixed_solve(X, Z, Y = Y, n_reps = 2)
#> Benchmarking mixed model implementations...
#> Data: n = 100 , p = 2 , q = 3 , voxels = 5 
#> Testing standard mixed_solve...
#> Testing optimized mixed_solve...
#> 
#> Results:
#>      method   mean_time median_time    min_time    max_time      sd_time
#> 1  standard 0.002838492 0.002838492 0.002372265 0.003304720 0.0006593453
#> 2 optimized 0.006560564 0.006560564 0.003158092 0.009963036 0.0048118214
#> 
#> Per-voxel timing:
#> Standard: 0.0006 sec/voxel
#> Optimized: 0.0013 sec/voxel
#> Speedup: 0.43x 
#>      method   mean_time median_time    min_time    max_time      sd_time
#> 1  standard 0.002838492 0.002838492 0.002372265 0.003304720 0.0006593453
#> 2 optimized 0.006560564 0.006560564 0.003158092 0.009963036 0.0048118214
# }
```
