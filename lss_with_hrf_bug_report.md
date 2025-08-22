# Bug Report: lss_with_hrf Function Issues

## Date: 2025-08-22
## Package: fmrilss
## Function: lss_with_hrf()
## Location: R/voxel_hrf.R:132-216

## Executive Summary
The `lss_with_hrf()` function has multiple critical issues that prevent it from working correctly with voxel-wise HRF estimates. The function appears to be incomplete or outdated, with fundamental design problems in how it handles HRF basis coefficients and constructs design matrices.

## Test Failures
4 out of 8 tests in `test-voxel-hrf.R` fail:
- Lines 53-90: "lss_with_hrf recovers trial betas"
- Lines 92-126: "lss_with_hrf equivalent to lss when HRF identical"  
- Lines 128-162: "lss_with_hrf basis aliasing"
- One additional test (not shown in current view)

## Critical Issues

### 1. Incorrect HRF Coefficient Dimensions
**Location:** test-voxel-hrf.R:76-82
**Problem:** Tests create incorrect coefficient matrices
```r
# Current (incorrect):
hcoef <- runif(n_vox, 0.5, 1.5)  # Single value per voxel
hrf_est <- list(coefficients = matrix(hcoef, 1, n_vox))  # 1 x n_vox matrix

# Expected:
# Should be nbasis x n_vox matrix for multi-basis HRFs
```

**Impact:** The test setup doesn't match what `estimate_voxel_hrf()` actually produces, which is an `nbasis x n_vox` matrix where `nbasis` is the number of basis functions.

### 2. Missing C++ Implementation
**Location:** R/voxel_hrf.R:208-211
**Problem:** Calls undefined C++ function
```r
lss_engine_vox_hrf(Y, hrf_estimates$coefficients, hrf_basis_kernels,
                   onset_idx, durations, nuisance_regs,
                   betas@address, update_progress,
                   as.integer(chunk_size), verbose)
```

**Impact:** The function `lss_engine_vox_hrf` is not implemented in the C++ codebase, causing immediate failure.

### 3. Incorrect HRF Basis Evaluation
**Location:** R/voxel_hrf.R:176-186
**Problem:** Incorrect method for evaluating HRF basis functions
```r
# Current approach tries to evaluate basis directly:
if (is.function(hrf_func)) {
    hrf_basis_kernels <- hrf_func(fine_grid)
} else {
    r <- fmrihrf::regressor(onsets = 0, hrf = hrf_func, duration = 0, span = hrf_span)
    hrf_basis_kernels <- fmrihrf::evaluate(r, fine_grid, precision = fine_dt, method = "conv")
}
```

**Issues:**
- Doesn't handle multi-basis HRFs correctly
- Should construct kernels from coefficients, not evaluate basis directly
- Mixing up basis functions with actual HRF kernels

### 4. Conceptual Design Flaw
**Core Problem:** The function attempts to use voxel-specific HRFs for LSS, but the implementation strategy is fundamentally flawed:

1. **HRF Reconstruction:** Should reconstruct voxel-specific HRF kernels from basis coefficients:
   ```r
   # For each voxel v:
   hrf_kernel_v <- hrf_basis %*% coefficients[, v]
   ```

2. **Design Matrix Construction:** Should build voxel-specific design matrices:
   ```r
   # For each voxel v:
   X_v <- convolve_design(events, hrf_kernel_v)
   ```

3. **LSS Application:** Should apply standard LSS to each voxel's custom design:
   ```r
   # For each voxel v:
   betas[, v] <- lss(Y[, v], X_v)
   ```

### 5. Test Logic Errors
**Location:** test-voxel-hrf.R:155-161
**Problem:** Incorrect assumptions about basis aliasing
```r
# Test expects:
est_kernel <- fmrihrf::hrf_from_coefficients(basis_fit, hrf_est$coefficients[,1])
expect_equal(est_kernel, hrf_kernel, tolerance = 1e-6)
```

**Issue:** This test assumes the fitted basis can exactly reproduce an HRF created with a different basis, which is not generally true.

## Recommended Solution

### Option 1: Complete Rewrite
Implement proper voxel-wise HRF LSS:

```r
lss_with_hrf <- function(Y, events, hrf_estimates, ...) {
  n_vox <- ncol(Y)
  n_trials <- nrow(events)
  betas <- matrix(NA, n_trials, n_vox)
  
  for (v in 1:n_vox) {
    # Reconstruct voxel-specific HRF
    hrf_kernel <- reconstruct_hrf(hrf_estimates$basis, 
                                  hrf_estimates$coefficients[, v])
    
    # Build design matrix with voxel HRF
    X_v <- build_design_matrix(events, hrf_kernel)
    
    # Apply LSS
    betas[, v] <- lss(Y[, v], X_v)
  }
  
  return(betas)
}
```

### Option 2: Deprecate Function
If voxel-wise HRF for LSS is not a core requirement:
1. Mark function as deprecated
2. Remove or skip failing tests
3. Document that users should use `estimate_voxel_hrf()` followed by standard `lss()`

### Option 3: Simplified Implementation
Implement a basic version that assumes single-coefficient HRFs:
```r
# Only support scaling of canonical HRF
# coefficients would be 1 x n_vox (scaling factors)
```

## Required Actions

1. **Immediate:** Skip failing tests with informative messages
2. **Short-term:** Decide on solution approach (rewrite vs deprecate)
3. **Implementation:** 
   - If rewriting: Implement C++ kernel `lss_engine_vox_hrf`
   - If deprecating: Update documentation and tests
4. **Testing:** Create proper test cases with correct dimensions
5. **Documentation:** Update function documentation to reflect actual behavior

## Test Fix Examples

### Corrected Test Setup
```r
test_that("lss_with_hrf with proper multi-basis HRF", {
  # Create multi-basis HRF estimates
  basis <- fmrihrf::hrf_fir_generator(nbasis = 5)
  n_basis <- 5
  n_vox <- 3
  
  # Coefficients should be nbasis x n_vox
  hrf_coef <- matrix(rnorm(n_basis * n_vox), n_basis, n_vox)
  
  hrf_est <- list(
    coefficients = hrf_coef,
    basis = basis,
    conditions = "A"
  )
  class(hrf_est) <- "VoxelHRF"
  
  # ... rest of test
})
```

## Dependencies
- Requires fmrihrf package updates for proper basis handling
- Needs C++ implementation of `lss_engine_vox_hrf`
- May need additional helper functions for HRF reconstruction

## Impact Assessment
- **Severity:** High - Function is completely non-functional
- **Scope:** Limited - Only affects voxel-wise HRF LSS functionality
- **Users:** Unknown - Function may not be in production use
- **Priority:** Medium - Core OASIS functionality works without this

## Recommendation
Given the extensive issues and the fact that the OASIS implementation provides the core LSS functionality, I recommend:
1. **Temporarily disable** the failing tests with skip() and explanatory messages
2. **Document** the function as experimental/under development
3. **Focus** on ensuring OASIS method is fully functional
4. **Revisit** lss_with_hrf implementation if/when there's a clear use case

## Files Affected
- R/voxel_hrf.R (lines 132-216)
- tests/testthat/test-voxel-hrf.R (lines 53-162)
- src/ directory (missing lss_engine_vox_hrf implementation)