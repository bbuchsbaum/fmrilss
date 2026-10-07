# fmrilss News

## fmrilss (development version)

### Performance
- The optimized LSS backends (`r_optimized`, `cpp_optimized`, `cpp`) now
  build the n x T LSS weight matrix from the residualized trial design and
  compute every trial beta with one matrix product. The data matrix is no
  longer residualized or copied, the per-voxel R loop is gone, and
  `cpp_optimized` only splits voxels across OpenMP threads when the linked
  BLAS is single-threaded. Default `lss()` on 600 scans x 20,000 voxels with
  150 trials: 1.9 s -> 0.07 s.
- Prewhitening fits AR models with global or run pooling from an n-column
  factor of the residual Gram matrix instead of the n x V residuals (exact),
  and computes noise residuals with BLAS-3 projections. AR(1) `lss()` on the
  same problem: 4.9 s -> 0.7 s.
- Non-allocating finiteness check for `Y`.

### New features
- `lss(trial_groups = )` fits LSS-N (Turner et al., 2012): one summed
  "other trials" regressor per condition, the model used by Nilearn's and
  NiBetaSeries' beta series.
- `lss(ridge = )` adds a fractional ridge penalty to each trial model (the
  OASIS `ridge_mode = "fractional"` convention), composable with
  `trial_groups` and prewhitening. In rapid designs with overlapping trials
  it lowers beta RMSE substantially (about a third in the benchmark) without
  changing pattern correlations.
- `prewhiten = list(pooling = "voxel")` and `pooling = "parcel"` now work
  with a shared design for `r_optimized`, `cpp_optimized` and `cpp`: each
  whitening operator gets its own filtered design. Voxel pooling bins voxels
  by residual autocorrelation (`voxel_bins`, default 50) and refits an AR
  model per bin, as in Nilearn's AR(1) GLM.

### Bug fixes
- The noise model was estimated from residuals of the full trial-wise (LSA)
  design. In rapid designs with many trials this biased the AR estimate
  strongly downward (e.g. -0.32 for data with AR(1) ~ 0.35), making
  prewhitened betas less accurate than OLS. The new
  `prewhiten$residual_model` defaults to `"aggregate"` (confounds plus one
  summed regressor per trial group or basis function); `"full"` restores the
  previous behaviour and is implied by a user-supplied residual-bias
  correction; `"corrected"` fits the full model with fmriAR's bias
  correction, built automatically (least biased; global/run pooling only).
  The new `vignette("prewhitening")` explains the trade-offs.

### Benchmarks
- `bench/python_comparison/` compares fmrilss with Nilearn (per-trial
  `run_glm`, OLS and AR(1)) and NumPy LSS implementations on a shared
  simulation: fmrilss reproduces the Python OLS estimates to < 5e-13, is
  ~200x faster than the Nilearn per-trial loop at 20k voxels, and its
  voxel-adaptive AR(1) matches Nilearn's AR(1) accuracy.

## fmrilss 0.2.0

### Major Enhancements

#### fmriAR Integration for Advanced Prewhitening
- **New `prewhiten` parameter** in `lss()` function provides comprehensive AR/ARMA noise modeling
  - Automatic AR order selection: `p = "auto"`
  - Voxel-specific parameters: `pooling = "voxel"`
  - Run-aware estimation: `pooling = "run"` with `runs` parameter
  - Parcel-based pooling: `pooling = "parcel"` with `parcels` parameter
  - ARMA models: `method = "arma"` for complex noise structures
- Works with all LSS methods (r_optimized, cpp_optimized, oasis, etc.)
- Leverages fmriAR's optimized C++ implementations with OpenMP
- `prewhiten_options()` now forwards fmriAR 0.3.3's opt-in residual-
  autocovariance bias-correction controls: `design`, `acvf_correction`, and
  `correction_max_lag`.
- Applied prewhitening now records the fitted `fmriAR_plan` in the result's
  `whiten_plan` attribute.
- `lss()` recognizes an unmodified multi-basis fmridesign design matrix and
  returns the same canonical trial-major rows as `lss_design()`.
- `lss_design()` supports every `lss()` estimator for one-basis designs and
  gives an actionable error when a non-OASIS estimator is used with a
  multi-basis model.
- `create_lwu_grid()` output now composes directly with
  `sbhm_build(library_spec = list(..., pgrid = grid))`.

#### API changes
- The legacy `oasis$whiten` option is deprecated and ignored. Use the top-level
  `prewhiten` argument and `prewhiten_options()` instead.
- Matrix responses are now supported by `mixed_solve()` as documented; each
  response column is fitted independently.
- Invalid block sizes and basis dimensions now fail before entering native
  blocked loops.

### Documentation Updates
- Enhanced vignettes with prewhitening examples:
  - `getting_started.Rmd`: New section on temporal autocorrelation
  - `oasis_method.Rmd`: Advanced prewhitening demonstrations
- Comprehensive examples in `examples/prewhitening_examples.R`
- Updated function documentation with detailed parameter descriptions

### Testing
- New test suite for fmriAR integration (`test-fmriAR-integration.R`)
- Updated existing tests to use new API
- Added regression tests for rank-deficient confounds, matrix responses, and
  blocked-loop contracts.

### Dependencies
- Added `fmriAR (>= 0.3.3)` to Imports.
- Raised the minimum R version to 4.0, matching the imported `fmriAR` package.

### Maintenance
- Consolidated the OASIS backend into one implementation owner and removed
  duplicate helper definitions.

## Version 0.1.0

### Added
- Initial support for voxel-wise HRF estimation and LSS using voxel-specific HRFs via `estimate_voxel_hrf()` and `lss_with_hrf()`.
