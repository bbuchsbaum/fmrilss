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
set.seed(1)
n_timepoints <- 20
n_trials <- 4
n_voxels <- 3

X <- matrix(0, n_timepoints, n_trials)
for (i in seq_len(n_trials)) {
  onset <- 1 + (i - 1) * 4
  X[onset:(onset + 2), i] <- 1
}

Z <- cbind(Intercept = 1, Linear = seq_len(n_timepoints))
Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

for (i in seq_len(n_trials)) {
  Y <- Y + X[, i] %*% matrix(rnorm(n_voxels, sd = 0.2), 1, n_voxels)
}

bdes <- list(
  dmat_base = Z,
  dmat_ran = X,
  dmat_fixed = NULL,
  fixed_ind = NULL
)

beta_estimates_naive <- lss_naive(Y = Y, bdes = bdes)

beta_estimates_fast <- lss(Y = Y, X = X, Z = Z)
max(abs(beta_estimates_naive - beta_estimates_fast))
#> [1] 8.881784e-16
# }
```
