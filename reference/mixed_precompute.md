# Precompute Workspace for Optimized Mixed Model

Performs expensive matrix computations that don't depend on the response
vector, allowing for efficient reuse across multiple voxels.

## Usage

``` r
mixed_precompute(X, Z, K = NULL)
```

## Arguments

- X:

  Fixed effects design matrix (n × p)

- Z:

  Random effects design matrix (n × q)

- K:

  Kinship/covariance matrix for random effects (q × q)

## Value

Workspace object for use with mixed_solve_optimized

## Examples

``` r
X <- matrix(1, 6, 1)
Z <- diag(6)
ws <- mixed_precompute(X, Z)
#> Workspace precomputed successfully:
#>   - n=6, p=1, q=6
#>   - Effective rank: 5
#>   - Using identity K: yes
names(ws)
#>  [1] "Q"              "U"              "phi"            "theta"         
#>  [5] "X_orig"         "Z_orig"         "XtXinv"         "n"             
#>  [9] "p"              "q"              "use_identity_K"
```
