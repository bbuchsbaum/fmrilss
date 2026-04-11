# Least Squares Separate (LSS) Analysis

Computes trial-wise beta estimates using the Least Squares Separate
approach of Mumford et al. (2012). This method fits a separate GLM for
each trial, with the trial of interest and all other trials as separate
regressors.

## Usage

``` r
lss(
  Y,
  X,
  Z = NULL,
  Nuisance = NULL,
  method = c("r_optimized", "cpp_optimized", "r_vectorized", "cpp", "naive", "oasis",
    "stglmnet"),
  block_size = 96,
  oasis = list(),
  stglmnet = list(),
  prewhiten = NULL
)
```

## Arguments

- Y:

  A numeric matrix of size n × V where n is the number of timepoints and
  V is the number of voxels/variables

- X:

  A numeric matrix of size n × T where T is the number of trials. Each
  column represents the design for one trial

- Z:

  A numeric matrix of size n × F representing experimental regressors to
  include in all trial-wise models. These are regressors we want to
  model and get beta estimates for, but not trial-wise (e.g., intercept,
  condition effects, block effects). If NULL, an intercept-only design
  is used. Defaults to NULL

- Nuisance:

  A numeric matrix of size n × N representing nuisance regressors to be
  projected out before LSS analysis (e.g., motion parameters,
  physiological noise). If NULL, no nuisance projection is performed.
  Defaults to NULL

- method:

  Character string specifying which implementation to use. Options are:

  - "r_optimized" - Optimized R implementation (recommended, default)

  - "cpp_optimized" - Optimized C++ implementation with parallel support

  - "r_vectorized" - Standard R vectorized implementation

  - "cpp" - Standard C++ implementation

  - "naive" - Simple loop-based R implementation (for testing)

  - "oasis" - OASIS method with HRF support and ridge regularization

  - "stglmnet" - overlap-aware elastic-net backend using `glmnet`

- block_size:

  An integer specifying the voxel block size for parallel processing,
  only applicable when `method = "cpp_optimized"`. Defaults to 96.

- oasis:

  A list of options for the OASIS method (ridge, SE, design
  construction, etc.). See Details and
  [`oasis_options`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
  for the full list. **Note:** `oasis$whiten` is deprecated and ignored.
  Use the `prewhiten` parameter instead for all temporal whitening.

- stglmnet:

  A list of options for the `method = "stglmnet"` backend. See Details
  and
  [`stglmnet_options`](https://bbuchsbaum.github.io/fmrilss/reference/stglmnet_options.md)
  for the common fields.

- prewhiten:

  A list of prewhitening options using the fmriAR package, or `NULL` (no
  whitening, the default). See Details and
  [`prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
  for the full list.

## Value

A numeric matrix of size T × V containing the trial-wise beta estimates.
Note: Currently only returns estimates for the trial regressors (X).
Beta estimates for the experimental regressors (Z) are computed but not
returned.

## Details

The LSS approach fits a separate GLM for each trial, where each model
includes:

- The trial of interest (from column i of X)

- All other trials combined (sum of all other columns of X)

- Experimental regressors (Z matrix) - these are modeled to get beta
  estimates but not trial-wise

If Nuisance regressors are provided, they are first projected out from
both Y and X using standard linear regression residualization.

When using method="oasis", the following options are available in the
oasis list (see also
[`oasis_options`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
for a validated constructor):

- `design_spec`: A list for building trial-wise designs from event
  onsets using fmrihrf. Must contain: `sframe` (sampling frame), `cond`
  (list with `onsets`, `hrf`, and optionally `span`), and optionally
  `others` (list of other conditions to be modeled as nuisances). When
  provided, X can be NULL and will be constructed automatically.

- `K`: Explicit basis dimension for multi-basis HRF models (e.g., 3 for
  SPMG3). If not provided, it's auto-detected from X dimensions or
  defaults to 1 for single-basis HRFs.

- `ridge_mode`: Either "fractional" (default) or "absolute". In absolute
  mode, ridge_x and ridge_b are used directly as regularization
  parameters. In fractional mode, they represent fractions of the mean
  design energy for adaptive regularization.

- `ridge_x`: Ridge parameter for trial-specific regressors (default
  0.05). Controls regularization strength for individual trial
  estimates.

- `ridge_b`: Ridge parameter for the aggregator regressor (default
  0.05). Controls regularization strength for the sum of all other
  trials.

- `return_se`: Logical, whether to return standard errors (default
  FALSE). When TRUE, returns a list with `beta` (trial estimates) and
  `se` (standard errors) components.

- `return_diag`: Logical, whether to return design diagnostics (default
  FALSE). When TRUE, includes diagnostic information about the design
  matrix structure.

- `block_cols`: Integer, voxel block size for memory-efficient
  processing (default 4096). Larger values use more memory but may be
  faster for systems with sufficient RAM.

- `ntrials`: Explicit number of trials (used when K \> 1 to determine
  output dimensions). If not provided, calculated as ncol(X) / K.

- `hrf_grid`: Vector of HRF indices for grid-based HRF selection
  (advanced use). Allows testing multiple HRF shapes simultaneously.

**Prewhitening (temporal autocorrelation correction):**

Use the top-level `prewhiten` parameter for all temporal whitening. This
replaces the old `oasis$whiten = "ar1"` syntax, which is now deprecated
and ignored. Do *not* put AR options inside the `oasis` list; they
belong in `prewhiten`.

When using `method = "stglmnet"`, the backend accepts an additional
nested `stglmnet=` list for lambda selection, overlap-adaptive
penalties, and optional pooled trial parameterizations. The common
pattern is `stglmnet = stglmnet_options(mode = "cv")` to select lambda
by cross-validation, or
`stglmnet = stglmnet_options(mode = "fixed", lambda = 0.01)` for a fixed
elastic-net fit. The backend reuses fmrilss prewhitening and
nuisance-projection utilities rather than maintaining a separate
whitening path.

The `prewhiten` list accepts the following fields (see also
[`prewhiten_options`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
for a validated constructor):

- `method`: Character, `"ar"` (default when the list is non-NULL),
  `"arma"`, or `"none"`. `"ar"` fits a pure autoregressive model;
  `"arma"` adds a moving-average component (requires `q > 0`).

- `p`: AR order. An integer, or `"auto"` (default) to select via AIC/BIC
  up to `p_max`. Use `p = 1` for a simple AR(1) model (the most common
  choice for fMRI); higher orders are rarely needed but may help with
  short TRs or multi-band sequences.

- `q`: Integer MA order for ARMA models (default 0). Only relevant when
  `method = "arma"`.

- `p_max`: Integer, maximum AR order when `p = "auto"` (default 6).

- `pooling`: How AR coefficients are estimated across voxels. One of:

  - `"global"`:

    (default) A single set of AR coefficients is estimated from the
    median autocorrelation across all voxels. Fast and usually adequate.

  - `"voxel"`:

    Fit a separate AR model per voxel. Most accurate but slow; consider
    `"parcel"` instead.

  - `"run"`:

    Fit one AR model per run (requires `runs`). Useful when noise
    structure differs between runs.

  - `"parcel"`:

    Fit one AR model per parcel (requires `parcels`). Good compromise
    between `"global"` and `"voxel"`.

- `runs`: Integer vector of length `nrow(Y)` giving run/block labels.
  Required for `pooling = "run"` and recommended whenever data span
  multiple runs so that whitening respects run boundaries.

- `parcels`: Integer vector of length `ncol(Y)` giving parcel labels.
  Required for `pooling = "parcel"`.

- `exact_first`: Character, `"ar1"` (default) or `"none"`. When `"ar1"`,
  the first observation of each segment is scaled by \\\sqrt{1 -
  \phi_1^2}\\ for the exact likelihood; `"none"` drops the first
  observation instead.

- `compute_residuals`: Logical (default TRUE). When TRUE, OLS residuals
  from the full design are computed before fitting the noise model. Set
  to FALSE only if Y is already residualized.

**Typical prewhiten recipes:**

      # Simple AR(1) — good default for most fMRI data
      prewhiten = list(method = "ar", p = 1)

      # Auto-select AR order (AIC), global pooling
      prewhiten = list(method = "ar", p = "auto")

      # Per-run AR(1) for multi-run data
      prewhiten = list(method = "ar", p = 1, pooling = "run",
                       runs = blockids)

      # Parcel-based AR with atlas labels
      prewhiten = list(method = "ar", p = 1, pooling = "parcel",
                       parcels = atlas_labels)

      # Or use the validated constructor:
      prewhiten = prewhiten_options(method = "ar", p = 1, pooling = "run",
                                    runs = blockids)

Prewhitening is applied before the LSS analysis to account for temporal
autocorrelation in the fMRI time series. Both Y and all design matrices
(X, Z, Nuisance) are filtered through the same whitening operator so
that OLS on the whitened system is equivalent to GLS on the original
data.

The OASIS method provides a mathematically equivalent but
computationally optimized version of standard LSS. It reformulates the
per-trial GLM fitting as a single matrix operation, eliminating
redundant computations. This is particularly beneficial for designs with
many trials or when processing large datasets. When K \> 1 (multi-basis
HRFs), the output will have K\*ntrials rows, with basis functions for
each trial arranged sequentially.

## References

Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012).
Deconvolving BOLD activation in event-related designs for multivoxel
pattern classification analyses. NeuroImage, 59(3), 2636-2643.

## Examples

``` r
n_timepoints <- 100
n_trials <- 10
n_voxels <- 50

X <- matrix(0, n_timepoints, n_trials)
for(i in 1:n_trials) {
  start <- (i-1) * 8 + 1
  if(start + 5 <= n_timepoints) {
    X[start:(start+5), i] <- 1
  }
}

Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
true_betas <- matrix(rnorm(n_trials * n_voxels, 0, 0.5), n_trials, n_voxels)
for(i in 1:n_trials) {
  Y <- Y + X[, i] %*% matrix(true_betas[i, ], 1, n_voxels)
}

beta_estimates <- lss(Y, X)

Z <- cbind(1, scale(1:n_timepoints))
beta_estimates_with_regressors <- lss(Y, X, Z = Z)

Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)
beta_estimates_clean <- lss(Y, X, Z = Z, Nuisance = Nuisance)

if (FALSE) { # \dontrun{
beta_oasis <- lss(Y, X, method = "oasis",
                  oasis = list(ridge_x = 0.1, ridge_b = 0.1,
                              ridge_mode = "fractional"))

result_with_se <- lss(Y, X, method = "oasis",
                     oasis = list(return_se = TRUE))
beta_estimates <- result_with_se$beta
standard_errors <- result_with_se$se

  sframe <- sampling_frame(blocklens = 200, TR = 1.0)

  beta_auto <- lss(Y, X = NULL, method = "oasis",
                   oasis = list(
                     design_spec = list(
                       sframe = sframe,
                       cond = list(
                         onsets = c(10, 30, 50, 70, 90, 110, 130, 150),
                         hrf = HRF_SPMG1,
                         span = 25
                       ),
                       others = list(
                         list(onsets = c(20, 40, 60, 80, 100, 120, 140))
                       )
                     )
                   ))

  beta_multibasis <- lss(Y, X = NULL, method = "oasis",
                        oasis = list(
                          design_spec = list(
                            sframe = sframe,
                            cond = list(
                              onsets = c(10, 30, 50, 70, 90),
                              hrf = HRF_SPMG3,
                              span = 30
                            )
                          ),
                          K = 3
                        ))
} # }
```
