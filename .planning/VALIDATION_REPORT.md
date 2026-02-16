# Validation Report: fmridesign Function Usage

## Overview

This report documents the validation of correct fmridesign function usage in the `lss_design()` implementation, following concerns about:
1. Whether `terms()` is the correct method
2. Whether term names ("drift", "block", "nuisance") are guaranteed

---

## Investigation Method

Three sub-agents performed comprehensive exploration of the fmridesign codebase at `~/code/fmridesign`:
1. Agent 1: baseline_model structure and methods
2. Agent 2: event_model structure and methods
3. Agent 3: Real-world usage examples from tests and vignettes

---

## Key Findings

### ✅ Finding 1: `terms()` is Correct for Both Models

**baseline_model:**
- Has explicit `terms.baseline_model()` method (R/baseline_model.R:250-256)
- Returns `x$terms` (a named list of term objects)
- Also has alias: `baseline_terms()`

**event_model:**
- Has explicit `terms.event_model()` method (R/event_model.R:172-175)
- Returns `x$terms` (a named list of term objects)
- Also has alias: `event_terms()`

**Conclusion:** Using `fmridesign::terms()` is **correct and idiomatic**.

---

### ✅ Finding 2: baseline_model Term Names are GUARANTEED

Term names in baseline_model are **hardcoded** in the constructor (R/baseline_model.R:159-170):

```r
terms_list <- list(
  drift = construct(drift_spec, sframe),                    # ALWAYS present
  block = if (intercept != "none" && basis != "constant"),  # CONDITIONAL
  nuisance = if (!is.null(nuisance_list))                   # CONDITIONAL
)
```

**Guaranteed term names:**
- `"drift"` - Always present (scanner drift basis functions)
- `"block"` - Present when `intercept != "none"` AND `basis != "constant"`
- `"nuisance"` - Present when `nuisance_list` is provided

**Evidence from tests** (test_baseline.R):
- Line 45: `expect_equal(names(terms(bmodel)), c("drift", "block"))`
- Line 63: `expect_true(setequal(names(terms(bmodel)), c("block", "drift", "nuisance")))`
- Lines 107, 112, 125: Confirm "drift" always present

**Conclusion:** Checking for `"drift"`, `"block"`, `"nuisance"` by name is **correct and safe**.

---

### ✅ Finding 3: Better Method Available - `term_matrices()`

While `design_matrix(terms(bmodel)[[i]])` is valid, fmridesign provides a cleaner method:

```r
# Original approach (valid but verbose)
bm_terms <- terms(bmodel)
drift_dm <- design_matrix(bm_terms$drift)
block_dm <- design_matrix(bm_terms$block)

# Improved approach (more idiomatic)
bm_term_mats <- term_matrices(bmodel)
drift_dm <- bm_term_mats$drift
block_dm <- bm_term_mats$block
```

**`term_matrices()` method** (R/baseline_model.R:265-301):
- Extracts design matrices for all terms in one call
- Returns a named list where keys match term names
- Avoids intermediate steps
- Used in multiple vignettes and tests

**Conclusion:** Using `term_matrices()` is **more idiomatic and cleaner**.

---

### ✅ Finding 4: Real-World Patterns Confirm Approach

**Pattern from vignette a_01_introduction.Rmd (lines 84-98):**
```r
# Extract design matrices
X_task <- design_matrix(emodel)
X_baseline <- design_matrix(bmodel)

# Combine
X_full <- cbind(X_task, X_baseline)
```

**Pattern from vignette a_04_event_models.Rmd (lines 458-467):**
```r
bmodel <- baseline_model(basis = "poly", degree = 5, sframe = sframe)
DM_full <- dplyr::bind_cols(design_matrix(emodel), design_matrix(bmodel))
```

**Conclusion:** Extracting design matrices separately then combining is **standard practice**.

---

## Changes Made to Implementation

### Original Code (R/lss_design.R:185-205)
```r
if (!is.null(baseline_model)) {
  bm_terms <- fmridesign::terms(baseline_model)

  # Combine drift and block terms into Z
  z_terms <- list()
  if ("drift" %in% names(bm_terms)) {
    z_terms$drift <- fmridesign::design_matrix(bm_terms$drift)
  }
  if ("block" %in% names(bm_terms)) {
    z_terms$block <- fmridesign::design_matrix(bm_terms$block)
  }

  if (length(z_terms) > 0) {
    Z <- as.matrix(do.call(cbind, z_terms))
  }

  # Extract nuisance term
  if ("nuisance" %in% names(bm_terms)) {
    Nuisance <- as.matrix(fmridesign::design_matrix(bm_terms$nuisance))
  }
}
```

**Issues:**
- ❌ None! Code was functionally correct
- ⚠️ Minor: Used intermediate `design_matrix()` calls instead of `term_matrices()`

### Updated Code
```r
if (!is.null(baseline_model)) {
  # Get term matrices using term_matrices() for cleaner extraction
  # term_matrices() returns a named list where keys are guaranteed to be:
  # "drift" (always), "block" (conditional), "nuisance" (conditional)
  bm_term_mats <- fmridesign::term_matrices(baseline_model)

  # Combine drift and block terms into Z
  # These are fixed effect regressors (scanner drift and run intercepts)
  z_terms <- list()
  if ("drift" %in% names(bm_term_mats)) {
    z_terms$drift <- bm_term_mats$drift
  }
  if ("block" %in% names(bm_term_mats)) {
    z_terms$block <- bm_term_mats$block
  }

  if (length(z_terms) > 0) {
    Z <- as.matrix(do.call(cbind, z_terms))
  }

  # Extract nuisance term (motion, physiology, etc.)
  # These are confounds to be projected out
  if ("nuisance" %in% names(bm_term_mats)) {
    Nuisance <- as.matrix(bm_term_mats$nuisance)
  }
}
```

**Improvements:**
- ✅ Uses `term_matrices()` instead of `terms()` + `design_matrix()`
- ✅ More idiomatic based on fmridesign patterns
- ✅ Better comments explaining guaranteed term structure
- ✅ Clearer separation between fixed effects (Z) and confounds (Nuisance)

---

## Test Updates

**Updated test** (tests/testthat/test-lss-design.R:304-338):
```r
# Method 2: Manual extraction + lss (using term_matrices)
X <- as.matrix(fmridesign::design_matrix(emod))
bm_term_mats <- fmridesign::term_matrices(bmodel)
Z <- as.matrix(cbind(bm_term_mats$drift, bm_term_mats$block))
beta2 <- lss(Y, X, Z, method = "oasis")
```

Now matches the updated implementation approach.

---

## Validation Summary

| Aspect | Status | Details |
|--------|--------|---------|
| `terms()` usage | ✅ Correct | Documented S3 method, widely used |
| Term name assumptions | ✅ Correct | Names are hardcoded, guaranteed |
| `design_matrix()` calls | ✅ Correct | Valid but improved to `term_matrices()` |
| Overall approach | ✅ Correct | Matches fmridesign patterns |

---

## Supporting Evidence

### File Locations Verified

**Method definitions:**
- `/Users/bbuchsbaum/code/fmridesign/R/baseline_model.R`
  - Lines 250-256: `terms.baseline_model()`
  - Lines 258-263: `baseline_terms.baseline_model()`
  - Lines 265-301: `term_matrices.baseline_model()`
  - Lines 159-170: Hardcoded term names in constructor

**Tests confirming term structure:**
- `/Users/bbuchsbaum/code/fmridesign/tests/testthat/test_baseline.R`
  - Lines 2-128: Comprehensive term structure tests
  - Lines 45, 63, 77: Explicit term name checks

**Vignette examples:**
- `/Users/bbuchsbaum/code/fmridesign/vignettes/a_01_introduction.Rmd`
  - Lines 84-98: Combining event and baseline models
- `/Users/bbuchsbaum/code/fmridesign/vignettes/a_03_baseline_model.Rmd`
  - Lines 145-147: Term extraction examples
- `/Users/bbuchsbaum/code/fmridesign/vignettes/a_04_event_models.Rmd`
  - Lines 458-467: Using `bind_cols()` for combination

---

## Conclusion

The original implementation was **functionally correct**, with all assumptions about fmridesign function usage validated:

1. ✅ `terms()` is the correct S3 method
2. ✅ Term names ("drift", "block", "nuisance") are guaranteed
3. ✅ `design_matrix()` calls on individual terms are valid

The updated implementation is **more idiomatic**, using:
- `term_matrices()` instead of `terms()` + `design_matrix()`
- Better comments documenting guaranteed structure
- Clearer semantic separation of term types

All changes are **non-breaking** and improve code quality without affecting functionality.
