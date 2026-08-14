# fmrilss News

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
