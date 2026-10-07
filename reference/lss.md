# Least Squares Separate (LSS) Analysis

Computes trial-wise beta estimates using the Least Squares Separate
approach of Mumford et al. (2012). This method fits a separate GLM for
each trial, with the trial of interest plus a single regressor formed by
summing all other trials in a one-basis design. A K-basis OASIS model
uses one summed other-trial regressor per basis function.

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
  prewhiten = NULL,
  trial_groups = NULL,
  ridge = NULL
)
```

## Arguments

- Y:

  A numeric matrix of size n × V where n is the number of timepoints and
  V is the number of voxels/variables

- X:

  A numeric matrix of size n × T for a one-basis design, with one column
  per trial. A raw K-basis OASIS design has n × (T K) columns and
  additionally requires `oasis$K`, `oasis$ntrials`, and an explicit
  `oasis$trial_basis_map`. When `X` is an unmodified multi-basis
  [`fmridesign::design_matrix()`](https://bbuchsbaum.github.io/fmridesign/reference/design_matrix.html)
  with one event term, its column metadata is used to infer that
  identity contract and canonicalize rows to trial-major,
  basis-within-trial order.

- Z:

  A numeric matrix of size n × F representing common fixed regressors
  included in every trial-wise model (e.g., intercept, condition
  effects, block effects). Their coefficients are not returned. If NULL,
  an intercept-only design is used. Defaults to NULL.

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
  for the full list. **Note:** `oasis$whiten` is deprecated and ignored
  with a warning. Use the `prewhiten` parameter instead for all temporal
  whitening.

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

- trial_groups:

  Optional vector (character, factor, or integer) with one condition
  label per trial, i.e. per column of `X`. When supplied, each
  trial-wise model uses one summed "other trials" regressor per
  condition (the LSS-N variant of Turner et al., 2012), with the trial
  of interest removed from its own condition's regressor, instead of a
  single regressor pooling every other trial. This is the model used by
  Nilearn's and NiBetaSeries' LSS beta series and is more accurate when
  conditions evoke different responses. Supported by methods
  `"r_optimized"`, `"cpp_optimized"`, `"cpp"`, and `"naive"`. Defaults
  to `NULL` (classic LSS).

- ridge:

  Optional fractional ridge penalty: one number, or two numbers
  `c(trial, others)` for the trial-of-interest and the other-trials
  coefficients. Each is a fraction of the mean design energy (the
  `ridge_mode = "fractional"` convention of OASIS), i.e. the penalty
  added to the trial diagonal is `ridge[1] * mean(c_i'c_i)`. Ridge
  shrinks trial estimates toward zero and can greatly reduce their
  variance in rapid designs where neighbouring trials overlap; it
  composes with `trial_groups` and `prewhiten`. Supported by methods
  `"r_optimized"`, `"cpp_optimized"`, `"cpp"`, and `"naive"`. Defaults
  to `NULL` (no ridge).

## Value

Normally, a numeric matrix of trial-wise beta estimates: T × V for a
one-basis design or (T K) × V for OASIS with K basis functions. With
OASIS `return_diag = TRUE` or `return_se = TRUE`, returns
`list(beta, diag?, se?)`; `beta` and `se` have the same row-by-voxel
shape. One-basis diagnostics contain length-T vectors `d`, `alpha`, and
`s`; multi-basis diagnostics contain K × K × T arrays `D`, `C`, and `E`.
Coefficients for common regressors `Z` are not returned. Multi-basis
beta and SE matrices carry the canonical `trial_basis_map` attribute.
When estimated prewhitening is active, the actual fitted `fmriAR_plan`
is available as `attr(result, "whiten_plan")` (and on `result$beta` for
structured returns).

## Details

The LSS approach fits a separate GLM for each trial, where each model
includes:

- The trial of interest (from column i of X)

- For a one-basis design, all other trials combined into one summed
  regressor. A K-basis OASIS model uses a K-column summed block.

- Common fixed regressors (Z matrix), whose coefficients are not
  returned

**Computation.** Each LSS estimate is a linear functional of the data,
\\\hat\beta_i = w_i^\top y\\. The optimized methods build the \\n \times
T\\ weight matrix \\W\\ from the (small) trial design and compute all
trial betas with one matrix product \\W^\top Y\\. Because the weights
lie in the residual space of the confounds, the data matrix is never
residualized or copied, so the cost is a single \\O(nTV)\\ BLAS call
regardless of the number of confounds.

If Nuisance regressors are provided, the rank-revealed combined span
`cbind(Z, Nuisance)` is projected from both Y and X before fitting.
Without a separate Nuisance matrix, Z remains explicitly in every
trial-wise model.

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
  SPMG3). A raw multi-basis `X` also requires `ntrials` and
  `trial_basis_map`. An unmodified multi-basis
  [`fmridesign::design_matrix()`](https://bbuchsbaum.github.io/fmridesign/reference/design_matrix.html)
  with one event term is recognized from its metadata; an ordinary raw
  `X` is otherwise interpreted as K=1.

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

- `return_se`: Logical, whether to return model-based standard errors
  (default FALSE). This is available only for unpenalized OASIS without
  estimated prewhitening.

- `return_diag`: Logical, whether to return design diagnostics (default
  FALSE). When TRUE, includes diagnostic information about the design
  matrix structure.

- `block_cols`: Integer, voxel block size for memory-efficient
  processing (default 4096). Larger values use more memory but may be
  faster for systems with sufficient RAM.

- `ntrials`: Required number of trials when a raw K \> 1 design is
  supplied.

- `trial_basis_map`: Required data frame for a raw K \> 1 design, with
  one row per X column and fields `column`, `trial`, and `basis`.

- `design_spec$hrf_grid`: Candidate HRFs for grid-based selection within
  an event-built design. A top-level `oasis$hrf_grid` field is invalid
  and rejected.

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

    Voxel-adaptive noise model. Per-voxel residual autocorrelations are
    estimated, voxels are grouped into `voxel_bins` bins (default 50) of
    similar autocorrelation, and an AR model is refitted per bin; each
    bin gets its own filtered design, as in Nilearn's AR(1) GLM.
    Supported by methods `"r_optimized"`, `"cpp_optimized"` and `"cpp"`;
    other methods reject it.

  - `"run"`:

    Fit one AR model per run (requires `runs`). Useful when noise
    structure differs between runs.

  - `"parcel"`:

    Fit one AR model per parcel (requires `parcels`); each parcel is
    fitted with its own filtered design. Supported by methods
    `"r_optimized"`, `"cpp_optimized"` and `"cpp"`.

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
  (see `residual_model`) are computed before fitting the noise model.
  Set to FALSE only if Y is already residualized.

- `design`: Optional numeric design matrix whose projection produced
  those residuals. Supplying it opts in to fmriAR's correction for
  downward bias in residual autocovariance. When
  `compute_residuals = TRUE`, it must span the same columns as the full
  `X`/`Z`/`Nuisance` design, including the intercept that fmrilss adds
  when none is already represented.

- `acvf_correction`: Optional correction matrix or list of matrices from
  [`fmriAR::acvf_bias_matrix()`](https://bbuchsbaum.github.io/fmriAR/reference/acvf_bias_matrix.html),
  used instead of `design` when reusing a correction across datasets.
  The two fields are mutually exclusive.

- `voxel_bins`: Positive integer number of autocorrelation bins for
  `pooling = "voxel"` (default 50).

- `residual_model`: which model's residuals the noise model is estimated
  from. `"aggregate"` (default) uses the confounds plus one summed trial
  regressor per `trial_groups` level (or per basis function); `"full"`
  uses the full trial-wise design, whose many columns bias the
  autocorrelation downward in rapid designs; `"corrected"` uses the full
  design with fmriAR's bias correction (least biased, slower, global/run
  pooling only). `"full"` is implied by `design`/`acvf_correction`. See
  [`vignette("prewhitening")`](https://bbuchsbaum.github.io/fmrilss/articles/prewhitening.md)
  for when to use each.

- `correction_max_lag`: Positive integer lag budget used when `design`
  is supplied (default 25). The correction is intended for
  high-pass-filtered designs; without high-pass filtering, the required
  lag budget can become impractically large. See
  [`fmriAR::fit_noise()`](https://bbuchsbaum.github.io/fmriAR/reference/fit_noise.html)
  for details.

**Typical prewhiten recipes:**


      # Simple AR(1) — good default for most fMRI data
      prewhiten = list(method = "ar", p = 1)

      # Auto-select AR order (AIC), global pooling
      prewhiten = list(method = "ar", p = "auto")

      # Per-run AR(1) for multi-run data
      prewhiten = list(method = "ar", p = 1, pooling = "run",
                       runs = blockids)

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

Turner, B. O., Mumford, J. A., Poldrack, R. A., & Ashby, F. G. (2012).
Spatiotemporal activity estimation for multivoxel pattern analysis with
rapid event-related designs. NeuroImage, 62(3), 1429-1438.

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

# \donttest{
beta_oasis <- lss(Y, X, method = "oasis",
                  oasis = list(ridge_x = 0.1, ridge_b = 0.1,
                              ridge_mode = "fractional"))

result_with_se <- lss(Y, X, method = "oasis",
                     oasis = list(return_se = TRUE,
                                  ridge_mode = "absolute",
                                  ridge_x = 0, ridge_b = 0))
beta_estimates <- result_with_se$beta
standard_errors <- result_with_se$se

  sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1.0)

  beta_auto <- lss(Y, X = NULL, method = "oasis",
                   oasis = list(
                     design_spec = list(
                       sframe = sframe,
                       cond = list(
                         onsets = c(10, 30, 50, 70),
                         hrf = fmrihrf::HRF_SPMG1,
                         span = 25
                       ),
                       others = list(
                         list(onsets = c(20, 40, 60, 80))
                       )
                     )
                   ))

  beta_multibasis <- lss(Y, X = NULL, method = "oasis",
                        oasis = list(
                          design_spec = list(
                            sframe = sframe,
                            cond = list(
                              onsets = c(10, 30, 50, 70),
                              hrf = fmrihrf::HRF_SPMG3,
                              span = 30
                            )
                          ),
                          K = 3
                        ))
# }
```
