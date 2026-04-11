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

This function implements the naive LSS approach where for each trial, a
separate GLM is fitted that includes:

- All base regressors (intercept, drift, etc.)

- All fixed effects regressors (if any)

- Only the current trial's regressor from the trial design matrix

While less efficient than the optimized
[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) function,
this implementation is conceptually simpler and can serve as a reference
or for validation purposes.

## See also

[`lss`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) for the
optimized implementation

## Examples

``` r
if (FALSE) { # \dontrun{
# Using same setup as lss() examples
beta_estimates_naive <- lss_naive(Y = Y, bdes = bdes)

# Compare with optimized version
beta_estimates_fast <- lss(Y = Y, bdes = bdes)
max(abs(beta_estimates_naive - beta_estimates_fast))
} # }
```
