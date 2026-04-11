# Extract Data Matrix from Dataset

Helper function to extract data matrix from various dataset formats.
This is a placeholder that should be customized based on your data
format.

## Usage

``` r
get_data_matrix(dset)
```

## Arguments

- dset:

  Dataset object (format depends on your specific use case)

## Value

A numeric matrix where rows are timepoints and columns are voxels

## Examples

``` r
get_data_matrix(matrix(1:6, nrow = 3))
#>      [,1] [,2]
#> [1,]    1    4
#> [2,]    2    5
#> [3,]    3    6
get_data_matrix(data.frame(a = 1:3, b = 4:6))
#>      a b
#> [1,] 1 4
#> [2,] 2 5
#> [3,] 3 6
```
