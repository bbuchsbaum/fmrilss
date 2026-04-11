# Project Out Confound Variables

Computes the orthogonal projection matrix Q = I - X(X'X)^(-1)X' that
projects out the space spanned by confound regressors X. This is useful
for advanced users who want to cache and reuse projection matrices.

## Usage

``` r
project_confounds(X)
```

## Arguments

- X:

  Confound design matrix (n x p) where n is number of timepoints and p
  is number of confound regressors

## Value

Projection matrix Q (n x n) that projects out the column space of X

## Details

This function uses QR decomposition for numerical stability instead of
computing the Moore-Penrose pseudoinverse directly. The resulting matrix
Q can be applied to data to remove the influence of confound regressors.

## Examples

``` r
if (FALSE) { # \dontrun{
n <- 100
X_confounds <- cbind(1, 1:n)

Q <- project_confounds(X_confounds)

Y_clean <- Q %*% Y_raw
} # }
```
