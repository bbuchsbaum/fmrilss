# Naive Least Squares Separate (LSS) Analysis

Performs LSS analysis using the naive approach where each trial model is
fit separately. This is the conceptually simplest implementation but
less efficient than the optimized
[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) function.

## Usage

``` r
lss_naive(Y = NULL, bdes, dset = NULL)
```

## Arguments

- Y:

  A numeric matrix where rows are timepoints and columns are
  voxels/features. If NULL, the function will attempt to extract data
  from `dset`.

- bdes:

  A list containing design matrices with components:

  - `dmat_base`: Base design matrix (e.g., intercept, drift terms)

  - `dmat_fixed`: Fixed effects design matrix (optional)

  - `dmat_ran`: Random/trial design matrix for LSS analysis

  - `fixed_ind`: Indices for fixed effects (optional)

- dset:

  Optional dataset object. If provided and Y is NULL, data will be
  extracted using `get_data_matrix`.

## Value

A numeric matrix with dimensions (n_events x n_voxels) containing the
LSS beta estimates for each trial and voxel.

## Details

This function implements the LSS approach by fitting a separate GLM for
each trial. Following the method described by Mumford et al. (2012), the
model for each trial includes:

- The regressor for the trial of interest.

- A single regressor representing all other trials (the sum of their
  individual regressors).

- All base regressors (e.g., intercept, drift terms).

- All fixed effects regressors (if any).

While less efficient than the optimized
[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) function,
this implementation is conceptually clear and serves as a reference for
validation.

## See also

[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) for the
optimized implementation

## Examples

``` r
# \donttest{
beta_estimates_naive <- lss_naive(Y = Y, bdes = bdes)
#> Error: object 'Y' not found

beta_estimates_fast <- lss(Y = Y, bdes = bdes)
#> Error in lss(Y = Y, bdes = bdes): unused argument (bdes = bdes)
max(abs(beta_estimates_naive - beta_estimates_fast))
#> Error: object 'beta_estimates_naive' not found
# }
```
