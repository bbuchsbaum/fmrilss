# fmrilss News

## fmrilss (development version)

### GLMsingle fixes

- Match event run IDs to matrix data in run order, including nonascending IDs.
- Guard fractional-ridge interpolation against zero/nonfinite norms and
  retain regularization when the design has zero eigenvalues.
- Handle constant or single-voxel R^2 distributions in automatic thresholding,
  and skip threshold estimation when denoising is disabled.
- Limit noise PCs to the available rank across runs; empty and rank-zero
  pools use zero PCs instead of arbitrary eigenvectors.

### New: `glmsingle()`

- `glmsingle()` implements GLMsingle (Prince et al., 2022): an ON-OFF model,
  a per-voxel HRF chosen from GLMsingle's 20-HRF library, GLMdenoise noise
  regressors chosen by cross-validation, and voxel-wise fractional ridge
  regression chosen by cross-validation over repeated conditions.
- It computes the same estimator as the reference implementation but solves
  each run separately, applies nuisance projections as low-rank products,
  scores models from sufficient statistics, and compiles the repeated-trial
  cross-validation into fixed per-trial weights.
  On simulated data with 8–12 runs and 20,000 voxels it runs 18–24 times
  faster than pinned Python GLMsingle on a single thread.
- Agreement with pinned Python GLMsingle (commit `1ab54a6`) is tested on 11
  scenarios: all HRF, noise-component and ridge-fraction choices match, and
  betas agree to single-precision accuracy.
- Defaults differ from GLMsingle only where GLMsingle is internally
  inconsistent (`extras_in_denoise`, `zero_sd_cv`); the alternative argument
  values reproduce GLMsingle.
- `glmsingle_design()` fits from an fmridesign event model;
  `glmsingle_hrf()` and `glmsingle_hrf_library()` return GLMsingle's HRFs.
- New vignette: `vignette("glmsingle")`.

### Vignettes

- Use CRAN albersdown (>= 2.1.0) and its `albers_vignette()` output
  format; theme assets are embedded from the installed package.

- Recalibrated `fmrilss`, `oasis_method` and `voxel-wise-hrf` for fmrihrf's
  corrected SPMG HRFs (smaller raw scale and a realistic undershoot). Designs
  now use unit-peak HRFs, and checks that were tied to the old kernel are
  relative or computed.

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
