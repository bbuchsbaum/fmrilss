# fmrilss News

## Version 0.2.0 (Development)

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

#### Backward Compatibility
- Old `oasis$whiten = "ar1"` syntax remains functional with deprecation notice
- Automatic conversion to new API format
- Smooth migration path for existing code

### Documentation Updates
- Enhanced vignettes with prewhitening examples:
  - `getting_started.Rmd`: New section on temporal autocorrelation
  - `oasis_method.Rmd`: Advanced prewhitening demonstrations
- Comprehensive examples in `examples/prewhitening_examples.R`
- Updated function documentation with detailed parameter descriptions

### Testing
- New test suite for fmriAR integration (`test-fmriAR-integration.R`)
- Updated existing tests to use new API
- Backward compatibility tests included

### Dependencies
- Added `fmriAR (>= 0.0.0.9000)` to Imports
- Added `bbuchsbaum/fmriAR` to Remotes

## Version 0.1.0

### Added
- Initial support for voxel-wise HRF estimation and LSS using voxel-specific HRFs via `estimate_voxel_hrf()` and `lss_with_hrf()`.
