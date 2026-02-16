# fmridesign Integration Plan: event_model + baseline_model Support

## Executive Summary

This document outlines a plan to integrate fmridesign's `event_model` and `baseline_model` objects into fmrilss, particularly for the OASIS method. The plan introduces a **new top-level interface** (`lss_design()`) while maintaining full backward compatibility with the existing `lss()` function.

---

## 1. Current Architecture vs. Proposed Architecture

### Current LSS Architecture

```r
lss(Y, X, Z, Nuisance, method = "oasis", oasis = list(design_spec = ...))
```

**Parameters:**
- `Y`: Data matrix (timepoints × voxels)
- `X`: Trial design matrix (manually constructed or via design_spec)
- `Z`: Fixed effects (intercepts, trends, blocks) - *typically NULL*
- `Nuisance`: Confounds (motion, physiology) - *typically NULL*
- `oasis$design_spec`: Simple list for automatic X construction

**Limitations:**
- Manual matrix construction required for complex designs
- No unified handling of baseline/nuisance structure
- `design_spec` lacks multi-condition, parametric modulator support
- Multi-run onset convention issues (global vs. run-relative)

### Proposed Architecture

```r
# New interface - accepts fmridesign objects
lss_design(Y, event_model, baseline_model = NULL, method = "oasis", ...)

# Existing interface - unchanged
lss(Y, X, Z, Nuisance, method = "oasis", oasis = list(), ...)
```

**Mapping:**

| fmridesign Object | LSS Parameter | Contains |
|-------------------|---------------|----------|
| `event_model` | → `X` | Task-related regressors (HRF-convolved trials/conditions) |
| `baseline_model$drift` | → `Z` | Scanner drift basis functions + block intercepts |
| `baseline_model$nuisance` | → `Nuisance` | Motion, physiology, etc. (non-convolved) |

---

## 2. Why a New Interface?

### Advantages of Separate Function

1. **Clear Intent**: `lss_design()` signals fmridesign-based workflow
2. **Simplified Signature**: Fewer parameters, more semantic
3. **Automatic Mapping**: Handles conversion from fmridesign objects internally
4. **Validation**: Can validate compatibility before processing
5. **Backward Compatibility**: Existing `lss()` code unchanged
6. **Documentation**: Clearer documentation per use case

### Alternative Considered: Extend lss()

```r
# Option rejected: adding parameters to lss()
lss(Y, X = NULL, event_model = NULL, baseline_model = NULL, ...)
```

**Why rejected:**
- Too many parameter combinations to validate
- Confusing for users (X vs. event_model?)
- Pollutes main lss() signature
- Harder to maintain separate documentation

---

## 3. Detailed API Design

### 3.1 Core Function: lss_design()

**File:** `R/lss_design.R`

```r
#' LSS Analysis with fmridesign Objects
#'
#' Perform Least Squares Separate (LSS) analysis using event_model and
#' baseline_model objects from the fmridesign package. This provides a
#' streamlined interface for complex designs with multi-condition, parametric
#' modulators, and structured nuisance handling.
#'
#' @param Y Numeric matrix of fMRI data (timepoints × voxels).
#' @param event_model An event_model object from \code{fmridesign::event_model()}.
#'   This defines the trial-wise or condition-wise task design. For LSS, typically
#'   created with \code{trialwise()} to generate one regressor per trial.
#' @param baseline_model Optional baseline_model object from
#'   \code{fmridesign::baseline_model()}. Defines drift correction, block intercepts,
#'   and nuisance regressors. If NULL, no baseline correction is applied.
#' @param method LSS method to use. Currently only "oasis" is supported for
#'   event_model integration.
#' @param oasis List of OASIS-specific options (see \code{?lss} for details).
#'   Note: \code{design_spec} is not used when providing event_model.
#' @param prewhiten Optional prewhitening specification (see \code{?lss}).
#' @param blockids Optional block/run identifiers for event_model. If NULL,
#'   extracted from event_model$blockids.
#' @param validate Logical. If TRUE (default), performs validation checks on
#'   design compatibility, collinearity, and temporal alignment.
#' @param ... Additional arguments passed to the underlying LSS method.
#'
#' @return Matrix of trial-wise beta estimates (trials × voxels), or
#'   (trials × basis_functions) × voxels for multi-basis HRFs.
#'
#' @details
#' \strong{Design Specification:}
#'
#' The \code{event_model} should typically use \code{trialwise()} for LSS:
#' \preformatted{
#'   emod <- event_model(onset ~ trialwise(basis = "spmg1"),
#'                       data = events,
#'                       block = ~run,
#'                       sampling_frame = sframe)
#' }
#'
#' For factorial designs (e.g., estimating condition-level betas separately):
#' \preformatted{
#'   emod <- event_model(onset ~ hrf(condition),
#'                       data = events,
#'                       block = ~run,
#'                       sampling_frame = sframe)
#' }
#'
#' \strong{Baseline Model:}
#'
#' If provided, baseline_model components are mapped as follows:
#' \itemize{
#'   \item \code{drift} and \code{block} terms → Z parameter (fixed effects)
#'   \item \code{nuisance} term → Nuisance parameter (confounds)
#' }
#'
#' \strong{Multi-Run Handling:}
#'
#' Both event_model and baseline_model must use the same \code{sampling_frame}.
#' Run structure is automatically respected. Event onsets should be run-relative
#' (resetting to 0 each run) as per fmridesign convention.
#'
#' \strong{Validation:}
#'
#' When \code{validate = TRUE}, the function checks:
#' \itemize{
#'   \item Temporal alignment: nrow(Y) matches total scans in sampling_frame
#'   \item Collinearity: Design matrix condition number < 30
#'   \item Compatibility: event_model and baseline_model use same sampling_frame
#' }
#'
#' @examples
#' \dontrun{
#' library(fmridesign)
#' library(fmrihrf)
#'
#' # 1. Define temporal structure
#' sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)
#'
#' # 2. Create trial data (run-relative onsets)
#' trials <- data.frame(
#'   onset = c(10, 30, 50, 70, 90, 110,   # Run 1
#'             10, 30, 50, 70, 90, 110),   # Run 2
#'   run = rep(1:2, each = 6)
#' )
#'
#' # 3. Build event model with trialwise
#' emod <- event_model(
#'   onset ~ trialwise(basis = "spmg1"),
#'   data = trials,
#'   block = ~run,
#'   sampling_frame = sframe
#' )
#'
#' # 4. Build baseline model
#' motion <- list(
#'   matrix(rnorm(150 * 6), 150, 6),  # Run 1: 6 motion params
#'   matrix(rnorm(150 * 6), 150, 6)   # Run 2: 6 motion params
#' )
#' bmodel <- baseline_model(
#'   basis = "bs",
#'   degree = 5,
#'   sframe = sframe,
#'   nuisance_list = motion
#' )
#'
#' # 5. Run LSS
#' Y <- matrix(rnorm(300 * 1000), 300, 1000)  # 300 scans, 1000 voxels
#' beta <- lss_design(Y, emod, bmodel, method = "oasis")
#'
#' # Output: 12 × 1000 (12 trials, 1000 voxels)
#' dim(beta)
#' }
#'
#' @seealso
#' \code{\link{lss}} for the traditional matrix-based interface,
#' \code{fmridesign::event_model} for event model creation,
#' \code{fmridesign::baseline_model} for baseline model creation
#'
#' @export
lss_design <- function(Y,
                       event_model,
                       baseline_model = NULL,
                       method = "oasis",
                       oasis = list(),
                       prewhiten = NULL,
                       blockids = NULL,
                       validate = TRUE,
                       ...) {

  # ---- Input Validation ----

  # Check fmridesign availability
  if (!requireNamespace("fmridesign", quietly = TRUE)) {
    stop("Package 'fmridesign' is required for lss_design(). ",
         "Install with: remotes::install_github('bbuchsbaum/fmridesign')",
         call. = FALSE)
  }

  # Check event_model
  if (missing(event_model) || is.null(event_model)) {
    stop("event_model is required. Create with fmridesign::event_model()",
         call. = FALSE)
  }
  if (!inherits(event_model, "event_model")) {
    stop("event_model must be an 'event_model' object from fmridesign",
         call. = FALSE)
  }

  # Check baseline_model if provided
  if (!is.null(baseline_model) && !inherits(baseline_model, "baseline_model")) {
    stop("baseline_model must be a 'baseline_model' object from fmridesign",
         call. = FALSE)
  }

  # Currently only OASIS method supported
  method <- match.arg(method, "oasis")

  # ---- Extract Components ----

  # Get sampling frame
  sframe <- event_model$sampling_frame

  # Validate baseline_model uses same sampling frame
  if (!is.null(baseline_model)) {
    if (!identical(baseline_model$sampling_frame, sframe)) {
      stop("event_model and baseline_model must use the same sampling_frame",
           call. = FALSE)
    }
  }

  # Extract blockids
  if (is.null(blockids)) {
    blockids <- event_model$blockids
  }

  # ---- Build Design Matrices ----

  # Extract X from event_model
  X <- as.matrix(fmridesign::design_matrix(event_model))

  # Extract Z and Nuisance from baseline_model
  Z <- NULL
  Nuisance <- NULL

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

  # ---- Validation ----

  if (validate) {
    # Check temporal alignment
    expected_scans <- sum(fmrihrf::blocklens(sframe))
    if (nrow(Y) != expected_scans) {
      stop(sprintf("Y has %d rows but sampling_frame expects %d scans",
                   nrow(Y), expected_scans),
           call. = FALSE)
    }

    if (nrow(X) != expected_scans) {
      stop(sprintf("Design matrix has %d rows but sampling_frame expects %d scans",
                   nrow(X), expected_scans),
           call. = FALSE)
    }

    # Check collinearity (combined design)
    X_full <- X
    if (!is.null(Z)) X_full <- cbind(X_full, Z)
    if (!is.null(Nuisance)) X_full <- cbind(X_full, Nuisance)

    cond_num <- kappa(crossprod(X_full), exact = TRUE)
    if (cond_num > 30) {
      warning(sprintf("High collinearity detected (condition number = %.1f). ",
                      "Consider using ridge regularization via oasis$ridge_*",
                      cond_num),
              call. = FALSE)
    }

    # Detect multi-basis and extract K
    K <- .detect_K_from_event_model(event_model)
    if (!is.null(K) && K > 1 && is.null(oasis$K)) {
      message(sprintf("Detected multi-basis HRF with K = %d basis functions", K))
      oasis$K <- K
    }
  }

  # ---- Call lss() ----

  result <- lss(
    Y = Y,
    X = X,
    Z = Z,
    Nuisance = Nuisance,
    method = method,
    oasis = oasis,
    prewhiten = prewhiten,
    ...
  )

  # ---- Attach Metadata ----

  attr(result, "event_model") <- event_model
  attr(result, "baseline_model") <- baseline_model
  attr(result, "sampling_frame") <- sframe
  attr(result, "method") <- "lss_design"

  return(result)
}
```

### 3.2 Helper Function: Detect K from event_model

```r
#' Detect number of basis functions from event_model
#' @keywords internal
.detect_K_from_event_model <- function(emod) {
  terms <- fmridesign::terms(emod)

  # Look for nbasis attribute in first term
  if (length(terms) > 0) {
    first_term <- terms[[1]]

    # Try to extract from HRF spec
    if (!is.null(first_term$hrf)) {
      K <- tryCatch(
        fmrihrf::nbasis(first_term$hrf),
        error = function(e) 1L
      )
      return(K)
    }

    # Try to infer from design matrix columns
    dm <- fmridesign::design_matrix(first_term)
    n_events <- length(unique(first_term$blockids))  # Rough estimate
    if (ncol(dm) %% n_events == 0) {
      return(as.integer(ncol(dm) / n_events))
    }
  }

  return(1L)  # Default: single basis
}
```

### 3.3 Validation Helper

```r
#' Validate design compatibility
#' @keywords internal
.validate_design_models <- function(event_model, baseline_model = NULL) {
  issues <- character(0)

  # Check sampling frame consistency
  if (!is.null(baseline_model)) {
    if (!identical(event_model$sampling_frame, baseline_model$sampling_frame)) {
      issues <- c(issues, "event_model and baseline_model have different sampling_frames")
    }
  }

  # Check for empty design
  dm <- fmridesign::design_matrix(event_model)
  if (ncol(dm) == 0) {
    issues <- c(issues, "event_model produces empty design matrix")
  }

  # Check for collinearity using fmridesign utilities
  if (requireNamespace("fmridesign", quietly = TRUE)) {
    tryCatch({
      fmridesign::check_collinearity(event_model, threshold = 0.95)
    }, error = function(e) {
      issues <- c(issues, paste("Collinearity detected:", e$message))
    })
  }

  if (length(issues) > 0) {
    stop("Design validation failed:\n  ", paste(issues, collapse = "\n  "),
         call. = FALSE)
  }

  invisible(TRUE)
}
```

---

## 4. Usage Examples

### Example 1: Basic Trial-Wise LSS

```r
library(fmrilss)
library(fmridesign)
library(fmrihrf)

# Setup: 2 runs, 150 scans each
sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)

# Trial data with run-relative onsets
trials <- data.frame(
  onset = c(10, 30, 50, 70, 90, 110,   # Run 1
            10, 30, 50, 70, 90, 110),   # Run 2
  run = rep(1:2, each = 6)
)

# Create trial-wise event model
emod <- event_model(
  onset ~ trialwise(basis = "spmg1"),
  data = trials,
  block = ~run,
  sampling_frame = sframe
)

# Run LSS (no baseline model)
Y <- matrix(rnorm(300 * 1000), 300, 1000)
beta <- lss_design(Y, emod, method = "oasis")

dim(beta)  # [1] 12 1000
```

### Example 2: LSS with Baseline Correction

```r
# Same setup as Example 1
sframe <- sampling_frame(blocklens = c(150, 150), TR = 2)
trials <- data.frame(
  onset = c(10, 30, 50, 70, 90, 110,
            10, 30, 50, 70, 90, 110),
  run = rep(1:2, each = 6)
)

emod <- event_model(
  onset ~ trialwise(basis = "spmg1"),
  data = trials,
  block = ~run,
  sampling_frame = sframe
)

# Add baseline model with motion correction
motion_run1 <- matrix(rnorm(150 * 6), 150, 6)
motion_run2 <- matrix(rnorm(150 * 6), 150, 6)

bmodel <- baseline_model(
  basis = "bs",
  degree = 5,
  sframe = sframe,
  intercept = "runwise",
  nuisance_list = list(motion_run1, motion_run2)
)

# Run LSS with baseline correction
beta <- lss_design(Y, emod, bmodel, method = "oasis")
```

### Example 3: Multi-Basis HRF with Ridge Regularization

```r
# Use SPMG3 (canonical + temporal + dispersion derivatives)
emod <- event_model(
  onset ~ trialwise(basis = "spmg3", nbasis = 3),
  data = trials,
  block = ~run,
  sampling_frame = sframe
)

bmodel <- baseline_model(
  basis = "bs",
  degree = 5,
  sframe = sframe
)

# LSS with ridge regularization
beta <- lss_design(
  Y, emod, bmodel,
  method = "oasis",
  oasis = list(
    K = 3,  # Auto-detected but can be explicit
    ridge_mode = "fractional",
    ridge_x = 0.02,
    ridge_b = 0.02
  )
)

# Output: (12 trials × 3 basis) × 1000 voxels = 36 × 1000
dim(beta)  # [1] 36 1000

# Extract canonical basis estimates
beta_canonical <- beta[seq(1, 36, by = 3), ]  # Every 3rd row starting at 1
```

### Example 4: Factorial LSS (Multi-Condition)

```r
# Events with condition labels
trials <- data.frame(
  onset = c(10, 30, 50, 70, 90, 110,
            10, 30, 50, 70, 90, 110),
  condition = factor(rep(c("face", "scene", "object"), 4)),
  run = rep(1:2, each = 6)
)

# Create condition-wise event model (not trialwise!)
emod <- event_model(
  onset ~ hrf(condition),
  data = trials,
  block = ~run,
  sampling_frame = sframe
)

# LSS will estimate separate betas per condition
# But OASIS method may need adjustment for this use case
# (Currently optimized for trial-wise)
```

### Example 5: Parametric Modulation

```r
# Events with reaction time modulator
trials <- data.frame(
  onset = c(10, 30, 50, 70, 90, 110,
            10, 30, 50, 70, 90, 110),
  RT = c(0.5, 0.7, 0.6, 0.8, 0.5, 0.9,
         0.6, 0.7, 0.8, 0.5, 0.6, 0.7),
  run = rep(1:2, each = 6)
)

# Center RT
trials$RT_c <- scale(trials$RT, center = TRUE, scale = FALSE)[, 1]

# Model: trial effects + RT modulation
emod <- event_model(
  onset ~ trialwise() + hrf(RT_c),
  data = trials,
  block = ~run,
  sampling_frame = sframe
)

# This creates trial-wise regressors PLUS an RT amplitude modulator
# Note: May need special handling in OASIS
```

---

## 5. Implementation Phases

### Phase 1: Core Infrastructure (Week 1-2)

**Tasks:**
1. Create `R/lss_design.R` with main function
2. Implement helper functions:
   - `.detect_K_from_event_model()`
   - `.validate_design_models()`
3. Add unit tests (`tests/testthat/test-lss-design.R`):
   - Basic trial-wise LSS
   - Multi-run handling
   - Baseline model integration
   - Matrix dimension validation
4. Update NAMESPACE and dependencies

**Deliverables:**
- Working `lss_design()` function
- Comprehensive test suite (>90% coverage)
- Documentation (roxygen2)

### Phase 2: Advanced Features (Week 3)

**Tasks:**
1. Multi-basis HRF support with automatic K detection
2. Ridge regularization integration
3. Prewhitening with run-aware pooling
4. Enhanced validation:
   - Collinearity checks
   - Design diagnostics
   - Temporal alignment
5. Add convenience wrappers:
   ```r
   lss_trialwise()  # Shortcut for common trialwise pattern
   ```

**Deliverables:**
- Multi-basis examples
- Ridge regularization tests
- Validation utilities

### Phase 3: Documentation & Vignettes (Week 4)

**Tasks:**
1. Create vignette: `vignettes/lss_with_fmridesign.Rmd`
   - Introduction to event_model and baseline_model
   - Step-by-step LSS workflow
   - Multi-run examples
   - Comparison with design_spec approach
2. Update existing vignettes:
   - `getting_started.Rmd`: Add lss_design() examples
   - `oasis_method.Rmd`: Show both interfaces
3. Create migration guide:
   - When to use lss_design() vs. lss()
   - Converting existing design_spec to event_model
   - Troubleshooting

**Deliverables:**
- Comprehensive vignette (15-20 pages)
- Updated documentation
- Migration guide

### Phase 4: Edge Cases & Polish (Week 5)

**Tasks:**
1. Handle edge cases:
   - Single-run designs
   - Variable-length runs
   - Missing trials/conditions
   - Factorial designs with LSS
2. Performance optimization:
   - Matrix pre-allocation
   - Efficient extraction
3. User experience:
   - Informative error messages
   - Progress indicators for large datasets
   - Warnings for common pitfalls

**Deliverables:**
- Edge case tests
- Optimized code
- User-friendly error handling

---

## 6. Backward Compatibility Strategy

### Existing Code: UNCHANGED

```r
# All existing lss() calls continue to work exactly as before

# Example 1: Manual X matrix
beta1 <- lss(Y, X, method = "oasis")

# Example 2: design_spec
beta2 <- lss(Y, method = "oasis",
            oasis = list(
              design_spec = list(
                sframe = sframe,
                cond = list(onsets = onsets, hrf = HRF_SPMG1)
              )
            ))

# Example 3: With baseline
beta3 <- lss(Y, X, Z, Nuisance, method = "oasis")
```

### New Code: ADDITIVE

```r
# New interface - opt-in only
beta4 <- lss_design(Y, event_model, baseline_model, method = "oasis")
```

### Internal Changes: MINIMAL

- No changes to `lss()` signature or behavior
- No changes to `.lss_oasis()` internals
- All new code in separate files
- Shared utilities in new `R/fmridesign_utils.R`

### Migration Path

Users can migrate gradually:

**Step 1:** Start using `event_model` for X construction
```r
emod <- event_model(onset ~ trialwise(), ...)
X <- design_matrix(emod)
beta <- lss(Y, X, method = "oasis")  # Old interface
```

**Step 2:** Adopt `lss_design()` when comfortable
```r
emod <- event_model(onset ~ trialwise(), ...)
beta <- lss_design(Y, emod, method = "oasis")  # New interface
```

**Step 3:** Add baseline_model for full integration
```r
emod <- event_model(onset ~ trialwise(), ...)
bmodel <- baseline_model(basis = "bs", ...)
beta <- lss_design(Y, emod, bmodel, method = "oasis")
```

---

## 7. Testing Strategy

### Unit Tests (`tests/testthat/test-lss-design.R`)

```r
test_that("lss_design basic functionality", {
  skip_if_not_installed("fmridesign")

  # Setup
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)
  emod <- event_model(onset ~ trialwise(), data = trials,
                     block = ~run, sampling_frame = sframe)
  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Test basic call
  beta <- lss_design(Y, emod, method = "oasis")

  expect_equal(dim(beta), c(5, 50))  # 5 trials, 50 voxels
  expect_true(!any(is.na(beta)))
})

test_that("lss_design with baseline_model", {
  skip_if_not_installed("fmridesign")

  sframe <- sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- event_model(onset ~ trialwise(), data = trials,
                     block = ~run, sampling_frame = sframe)
  bmodel <- baseline_model(basis = "bs", degree = 5, sframe = sframe)

  Y <- matrix(rnorm(100 * 50), 100, 50)
  beta <- lss_design(Y, emod, bmodel, method = "oasis")

  expect_equal(dim(beta), c(5, 50))
})

test_that("lss_design multi-run handling", {
  skip_if_not_installed("fmridesign")

  sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)
  trials <- data.frame(
    onset = c(10, 30, 50, 10, 30, 50),
    run = rep(1:2, each = 3)
  )

  emod <- event_model(onset ~ trialwise(), data = trials,
                     block = ~run, sampling_frame = sframe)
  Y <- matrix(rnorm(200 * 50), 200, 50)
  beta <- lss_design(Y, emod, method = "oasis")

  expect_equal(nrow(beta), 6)  # 6 trials across 2 runs
})

test_that("lss_design validates temporal alignment", {
  skip_if_not_installed("fmridesign")

  sframe <- sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)
  emod <- event_model(onset ~ trialwise(), data = trials,
                     block = ~run, sampling_frame = sframe)

  # Wrong number of timepoints
  Y_wrong <- matrix(rnorm(150 * 50), 150, 50)

  expect_error(
    lss_design(Y_wrong, emod, method = "oasis"),
    "150 rows but sampling_frame expects 100"
  )
})

test_that("lss_design detects multi-basis K", {
  skip_if_not_installed("fmridesign")

  sframe <- sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- event_model(
    onset ~ trialwise(basis = "spmg3", nbasis = 3),
    data = trials,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Should auto-detect K = 3
  expect_message(
    beta <- lss_design(Y, emod, method = "oasis"),
    "Detected multi-basis HRF with K = 3"
  )

  expect_equal(nrow(beta), 15)  # 5 trials × 3 basis = 15
})
```

### Integration Tests

```r
test_that("lss_design equivalent to lss with manual extraction", {
  skip_if_not_installed("fmridesign")

  # Setup
  sframe <- sampling_frame(blocklens = 100, TR = 2)
  trials <- data.frame(onset = seq(10, 90, by = 20), run = 1)

  emod <- event_model(onset ~ trialwise(), data = trials,
                     block = ~run, sampling_frame = sframe)
  bmodel <- baseline_model(basis = "bs", degree = 5, sframe = sframe)

  Y <- matrix(rnorm(100 * 50), 100, 50)

  # Method 1: lss_design
  beta1 <- lss_design(Y, emod, bmodel, method = "oasis")

  # Method 2: Manual extraction + lss
  X <- as.matrix(design_matrix(emod))
  bm_terms <- terms(bmodel)
  Z <- as.matrix(cbind(
    design_matrix(bm_terms$drift),
    design_matrix(bm_terms$block)
  ))
  beta2 <- lss(Y, X, Z, method = "oasis")

  expect_equal(beta1, beta2, tolerance = 1e-10)
})
```

---

## 8. Documentation Structure

### Function Documentation (`man/lss_design.Rd`)

Already outlined in Section 3.1 above.

### Vignette: `vignettes/lss_with_fmridesign.Rmd`

**Outline:**

1. **Introduction**
   - Why use fmridesign with fmrilss?
   - When to use `lss_design()` vs. `lss()`

2. **Quick Start**
   - Minimal working example
   - Comparison with design_spec

3. **Event Models for LSS**
   - Creating trial-wise designs with `trialwise()`
   - Multi-run experiments
   - Multi-basis HRFs
   - Parametric modulators (advanced)

4. **Baseline Models**
   - Drift correction options
   - Adding motion parameters
   - Block intercepts

5. **Complete Workflow**
   - BIDS dataset example
   - Multi-run, multi-subject
   - Result visualization

6. **Advanced Topics**
   - Ridge regularization
   - Prewhitening
   - Design diagnostics
   - Troubleshooting

### Migration Guide

Separate document or appendix explaining:
- Converting `design_spec` → `event_model`
- Onset convention handling
- Common pitfalls
- Performance considerations

---

## 9. Dependencies & NAMESPACE Updates

### DESCRIPTION Changes

```r
Suggests:
    fmridesign (>= 0.1.0),  # Already present, keep version requirement
    testthat (>= 3.0.0),
    ...
```

**No change to dependency level** - fmridesign stays in Suggests.

### NAMESPACE Updates

```r
# R/lss_design.R
export(lss_design)

# Imports from fmridesign (conditional, via requireNamespace)
# No explicit imports needed since we use :: notation
```

### Conditional Loading Pattern

```r
# In lss_design.R
if (!requireNamespace("fmridesign", quietly = TRUE)) {
  stop("Package 'fmridesign' required. Install with: ...", call. = FALSE)
}
```

This keeps fmridesign optional for users who don't need it.

---

## 10. Future Enhancements (Post-MVP)

### 10.1 Support for Other LSS Methods

Currently targeting OASIS only. Future:

```r
lss_design(Y, emod, bmodel, method = "r_optimized")
lss_design(Y, emod, bmodel, method = "cpp_optimized")
```

### 10.2 SBHM Integration

```r
lss_sbhm_design(Y, event_model, baseline_model, sbhm_spec, ...)
```

### 10.3 Result Object Enhancement

Return rich S3 object instead of plain matrix:

```r
result <- lss_design(...)
class(result)  # "lss_result"

# Methods
plot(result)
coef(result, trial = 5)
residuals(result)
```

### 10.4 Design Optimization Utilities

```r
optimize_design(event_model, baseline_model, criterion = "efficiency")
suggest_regularization(event_model, baseline_model, Y)
```

---

## 11. Success Criteria

### Functionality
- ✅ `lss_design()` produces identical results to manual `lss()` call
- ✅ Handles single-run and multi-run designs correctly
- ✅ Baseline model integration works for all term types
- ✅ Multi-basis HRF auto-detection works
- ✅ Validation catches common errors

### Performance
- ✅ No performance regression vs. direct `lss()` calls (<5% overhead)
- ✅ Memory efficient for large datasets (no extra copies)

### Usability
- ✅ Clear error messages with actionable guidance
- ✅ Comprehensive documentation with examples
- ✅ Vignette covers 90% of use cases

### Compatibility
- ✅ All existing tests pass
- ✅ Backward compatibility maintained
- ✅ Works with fmridesign >= 0.1.0

### Testing
- ✅ >95% code coverage for new functions
- ✅ Integration tests with real-world scenarios
- ✅ Edge cases handled gracefully

---

## 12. Timeline & Milestones

| Week | Milestone | Deliverables |
|------|-----------|--------------|
| 1-2 | Core Implementation | `lss_design()`, helpers, basic tests |
| 3 | Advanced Features | Multi-basis, ridge, validation |
| 4 | Documentation | Vignette, migration guide |
| 5 | Polish & Edge Cases | Error handling, optimization |
| 6 | Review & Release | Code review, CRAN prep |

**Total Estimate: 6 weeks (1 developer, part-time)**

---

## 13. Open Questions & Decisions Needed

1. **Naming:** Is `lss_design()` the best name? Alternatives:
   - `lss_fmridesign()`
   - `lss_model()`
   - `lss_from_design()`

2. **Validation strictness:** Should validation be opt-out (default TRUE) or opt-in?

3. **Method support:** Should we implement non-OASIS methods in Phase 1?

4. **Return object:** Plain matrix (simple) or rich S3 object (future-proof)?

5. **Factorial LSS:** How to handle multi-condition designs in LSS context?
   - Currently LSS is trial-wise
   - Multi-condition might need different output structure

---

## 14. Example Use Case: BIDS Dataset

```r
library(fmrilss)
library(fmridesign)
library(fmrihrf)

# 1. Load BIDS event files
events_run1 <- read.delim("sub-01_task-faces_run-1_events.tsv")
events_run2 <- read.delim("sub-01_task-faces_run-2_events.tsv")

# Combine runs (add run column)
events_run1$run <- 1
events_run2$run <- 2
events <- rbind(events_run1, events_run2)

# 2. Define sampling frame
TR <- 2.0
nscans_run1 <- 200
nscans_run2 <- 200
sframe <- sampling_frame(blocklens = c(nscans_run1, nscans_run2), TR = TR)

# 3. Create event model (trial-wise)
emod <- event_model(
  onset ~ trialwise(basis = "spmg1"),
  data = events,
  block = ~run,
  sampling_frame = sframe,
  durations = events$duration
)

# 4. Load motion parameters
motion_run1 <- read.table("sub-01_task-faces_run-1_motion.txt")
motion_run2 <- read.table("sub-01_task-faces_run-2_motion.txt")

# 5. Create baseline model
bmodel <- baseline_model(
  basis = "bs",
  degree = 5,
  sframe = sframe,
  intercept = "runwise",
  nuisance_list = list(as.matrix(motion_run1), as.matrix(motion_run2))
)

# 6. Load fMRI data
Y <- read_nifti_matrix("sub-01_task-faces_runs-1-2_bold.nii.gz")

# 7. Run LSS
beta <- lss_design(
  Y, emod, bmodel,
  method = "oasis",
  oasis = list(
    ridge_mode = "fractional",
    ridge_x = 0.02
  )
)

# 8. Save results
write_nifti_matrix(beta, "sub-01_task-faces_lss-betas.nii.gz")
```

---

## Conclusion

This plan provides a **comprehensive, backward-compatible integration** of fmridesign's powerful event_model and baseline_model into fmrilss. The new `lss_design()` interface:

- ✅ **Simplifies** complex design specifications
- ✅ **Leverages** fmridesign's validation and diagnostics
- ✅ **Maintains** full backward compatibility
- ✅ **Enables** advanced features (parametric modulators, multi-condition)
- ✅ **Fixes** onset convention issues
- ✅ **Documents** best practices clearly

**Recommended Next Steps:**
1. Review and approve this plan
2. Implement Phase 1 (core infrastructure)
3. Gather user feedback from beta testers
4. Iterate on Phases 2-4 based on feedback
