# fmridesign Integration Analysis

## Executive Summary

This document analyzes the integration between `fmrilss` and `fmridesign`, particularly focusing on design specification systems (`design_spec` vs `event_model`) and multi-run experiment handling.

## Critical Issue: Onset Timing Convention Mismatch

### Problem
- **fmridesign convention**: Run-relative onsets (reset to 0 each run) → automatically converted to global
- **fmrilss design_spec**: Expects absolute/global onsets → NO automatic conversion

### Impact
Silent failures in multi-run LSS analyses if run-relative onsets are passed to `design_spec`.

### Example of Problematic Code
```r
# Multi-run design with run-relative onsets (STANDARD CONVENTION)
events <- data.frame(
  onset = c(10, 30, 50,   # Run 1: relative to run start
            10, 30, 50),  # Run 2: relative to run start (PROBLEM!)
  run = c(1, 1, 1, 2, 2, 2)
)

sframe <- sampling_frame(blocklens = c(100, 100), TR = 2)

# WRONG - will model run 2 trials at same times as run 1
beta <- lss(Y, method = "oasis",
           oasis = list(
             design_spec = list(
               sframe = sframe,
               cond = list(onsets = events$onset)  # INCORRECT!
             )
           ))

# CORRECT - convert to global first
events$global_onset <- fmrihrf::global_onsets(sframe, events$onset, events$run)
beta <- lss(Y, method = "oasis",
           oasis = list(
             design_spec = list(
               sframe = sframe,
               cond = list(onsets = events$global_onset)  # CORRECT
             )
           ))
```

## Capability Comparison: design_spec vs event_model

| Feature | design_spec | event_model |
|---------|-------------|-------------|
| Trial-wise designs | ✅ | ✅ |
| Multi-condition designs | ❌ | ✅ |
| Parametric modulators | ⚠️ Limited | ✅ Full |
| Formula-based DSL | ❌ | ✅ |
| Explicit run handling | ⚠️ Implicit | ✅ Explicit |
| Onset convention | Global/absolute | Run-relative |
| Onset conversion | ❌ Manual | ✅ Automatic |
| Non-convolved covariates | ❌ | ✅ |
| Contrast specifications | ❌ | ✅ |
| Design validation | ❌ | ✅ |
| Multi-basis HRF | ✅ | ✅ |

## Current Integration Architecture

### Dependency Status
- `fmridesign` is in **Suggests** (not Imports)
- Conditional usage in 3 files (scripts, vignettes)
- No core integration in LSS/OASIS/SBHM internals

### Current Usage Pattern
```r
if (requireNamespace("fmridesign", quietly = TRUE)) {
  emod <- fmridesign::event_model(onset ~ hrf(trial, basis = basis), ...)
  X <- fmridesign::design_matrix(emod)
} else {
  X <- fmrihrf::evaluate(fmrihrf::regressor_set(...), ...)
}
beta <- lss(Y, X, ...)
```

## Recommended Integration Strategy: Enhanced Conditional Integration

### Goals
1. Keep backward compatibility with `design_spec`
2. Enable optional `event_model` usage for advanced features
3. Avoid forcing `fmridesign` dependency
4. Fix onset convention mismatch issues

### Implementation Plan

#### Phase 1: Core Integration (1-2 weeks)

**1. Add converter function** (`R/event_model_helpers.R`):
```r
#' Convert event_model to design_spec
#' @keywords internal
.event_model_to_design_spec <- function(emod, trial_term = NULL) {
  # Extract sampling frame
  sframe <- emod$sampling_frame

  # Find trial-wise term (or use specified term)
  if (is.null(trial_term)) {
    # Auto-detect trialwise term
    trial_term <- which(sapply(emod$terms, function(t) {
      inherits(t, "event_term") && isTRUE(attr(t, "trialwise"))
    }))[1]
    if (is.na(trial_term)) {
      stop("No trialwise term found. Specify trial_term explicitly.")
    }
  }

  term <- emod$terms[[trial_term]]

  # Extract onsets (already in global time from event_model)
  # BUT: need to convert BACK to run-relative for design_spec if we fix it
  # OR: document that we expect global onsets

  list(
    sframe = sframe,
    cond = list(
      onsets = term$onsets,  # These are already global from event_model
      hrf = term$hrf,
      duration = term$durations[1],
      span = attr(term, "span") %||% 30
    )
  )
}
```

**2. Update `lss()` signature** (`R/lss.R`):
```r
lss <- function(Y, X = NULL, Z = NULL, Nuisance = NULL,
                event_model = NULL,  # NEW
                method = c("r_optimized", "cpp_optimized", "oasis", ...),
                oasis = list(),
                prewhiten = NULL,
                ...) {

  # Handle event_model input
  if (!is.null(event_model)) {
    if (!requireNamespace("fmridesign", quietly = TRUE)) {
      stop("fmridesign package required for event_model parameter")
    }

    if (!inherits(event_model, "event_model")) {
      stop("event_model must be created with fmridesign::event_model()")
    }

    # Extract design matrix
    X <- as.matrix(fmridesign::design_matrix(event_model))

    # If using OASIS and no design_spec provided, create one
    if (method == "oasis" && is.null(oasis$design_spec)) {
      oasis$design_spec <- .event_model_to_design_spec(event_model)
    }
  }

  # Existing validation and dispatch...
}
```

**3. Fix onset convention in design_spec** (`R/oasis_glue.R`):
```r
.oasis_build_X_from_events <- function(spec) {
  sframe <- spec$sframe

  # NEW: Check if onsets need global conversion
  # If blockids provided, assume run-relative and convert
  if (!is.null(spec$cond$blockids)) {
    onsets <- fmrihrf::global_onsets(
      sframe,
      spec$cond$onsets,
      spec$cond$blockids
    )
  } else {
    # Assume already global
    onsets <- spec$cond$onsets
  }

  # Rest of existing code...
}
```

**4. Add tests** (`tests/testthat/test-event-model-integration.R`):
```r
test_that("event_model integration works", {
  skip_if_not_installed("fmridesign")

  # Create multi-run event_model
  events <- data.frame(
    onset = c(10, 30, 50, 10, 30, 50),
    run = c(1, 1, 1, 2, 2, 2)
  )
  sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2)

  emod <- fmridesign::event_model(
    onset ~ fmridesign::trialwise(basis = "spmg1"),
    data = events,
    block = ~run,
    sampling_frame = sframe
  )

  Y <- matrix(rnorm(200 * 50), 200, 50)

  # Test 1: event_model parameter
  beta1 <- lss(Y, event_model = emod, method = "oasis")

  # Test 2: Manual X extraction (should be equivalent)
  X <- as.matrix(fmridesign::design_matrix(emod))
  beta2 <- lss(Y, X, method = "oasis")

  expect_equal(beta1, beta2)
})

test_that("design_spec with blockids converts onsets", {
  # Test that run-relative onsets get converted properly
  sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2)

  design_spec <- list(
    sframe = sframe,
    cond = list(
      onsets = c(10, 30, 50, 10, 30, 50),
      blockids = c(1, 1, 1, 2, 2, 2),  # NEW: explicit block IDs
      hrf = fmrihrf::HRF_SPMG1
    )
  )

  Y <- matrix(rnorm(200 * 50), 200, 50)
  beta <- lss(Y, method = "oasis", oasis = list(design_spec = design_spec))

  expect_equal(nrow(beta), 6)
})
```

#### Phase 2: Documentation (1 week)

**1. Create integration vignette** (`vignettes/event_model_integration.Rmd`):
```rmd
---
title: "Using fmridesign with fmrilss"
---

## Overview

The `fmridesign` package provides a powerful formula-based interface for
creating fMRI design matrices. While `fmrilss` has its own `design_spec`
format, you can also use `event_model` objects directly.

## When to Use event_model

Use `fmridesign::event_model()` when:

- You have **multi-condition factorial designs**
- You need **parametric modulators** (e.g., RT, difficulty)
- You want **design validation** and diagnostic plots
- You prefer **formula-based specification**

Use `design_spec` when:

- You have simple **trial-wise designs**
- You want **minimal dependencies**
- You're using **SBHM internal pipelines**

## Multi-Run Onset Convention: IMPORTANT

**fmridesign uses run-relative onsets** (resetting to 0 each run):
# ... rest of vignette
```

**2. Update `?lss` documentation**:
```r
#' @param event_model Optional. An event_model object from
#'   \code{fmridesign::event_model()}. If provided, X will be extracted
#'   automatically. For OASIS method, design_spec will also be generated
#'   if not provided. Requires the fmridesign package.
```

#### Phase 3: Migration Guide (1 week)

**Create migration document** (`vignettes/design_spec_migration.Rmd`):
- How to convert existing design_spec to event_model
- Multi-run onset timing guide
- Troubleshooting common issues

## Benefits of This Approach

✅ **Backward compatible** - existing code keeps working
✅ **Lightweight** - fmridesign stays in Suggests
✅ **Fixes onset issue** - adds blockids support to design_spec
✅ **Enables advanced features** - users can opt into event_model power
✅ **Clear migration path** - users can migrate gradually

## Risks & Mitigation

| Risk | Mitigation |
|------|------------|
| Breaking changes to design_spec | Add blockids as optional, maintain backward compatibility |
| Converter limitations | Document which event_model features are unsupported |
| Test coverage | Comprehensive test suite for both pathways |
| Documentation burden | Clear vignettes with use case guidance |

## Timeline Estimate

- Phase 1 (Core): 1-2 weeks
- Phase 2 (Docs): 1 week
- Phase 3 (Migration): 1 week
- **Total**: 3-4 weeks

## Alternative: Status Quo + Documentation

If full integration is too ambitious:

**Minimal effort option** (1-2 days):
1. Add warning to `design_spec` documentation about onset conventions
2. Add example showing how to use `event_model` → extract X → pass to lss()
3. Document `fmrihrf::global_onsets()` for manual conversion

This avoids code changes but leaves the integration burden on users.

## Conclusion

**Recommended**: Proceed with Enhanced Conditional Integration (Phases 1-2)

This provides maximum value with minimal risk, fixing the critical onset
convention issue while enabling optional advanced features for users who
need them.
