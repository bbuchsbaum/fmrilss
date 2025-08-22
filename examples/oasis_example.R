# OASIS Method Examples for fmrilss
# Demonstrates the OASIS implementation with various HRF types and options

library(fmrilss)
library(fmrihrf)

# Set up example data
set.seed(123)
sframe <- sampling_frame(blocklens = c(150, 150), TR = 1.0)
T <- sum(blocklens(sframe))
V <- 100  # number of voxels

# Create synthetic fMRI data
Y <- matrix(rnorm(T * V), T, V)

# Example 1: Single-basis HRF (SPMG1) with event-based design
# ------------------------------------------------------------
cat("\n=== Example 1: Single-basis HRF (SPMG1) ===\n")

beta_spmg1 <- lss(
  Y, X = NULL, 
  method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(
        onsets = c(10, 30, 50, 70, 90, 110, 130, 150, 170, 190),
        hrf = HRF_SPMG1,
        span = 30
      )
    )
  )
)

cat("Single-basis result dimensions:", dim(beta_spmg1), "\n")
cat("  - Rows:", nrow(beta_spmg1), "trials\n")
cat("  - Cols:", ncol(beta_spmg1), "voxels\n")

# Example 2: Multi-basis HRF (SPMG3) with canonical + derivatives
# ----------------------------------------------------------------
cat("\n=== Example 2: Multi-basis HRF (SPMG3) ===\n")

beta_spmg3 <- lss(
  Y, X = NULL,
  method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(
        onsets = c(10, 30, 50, 70, 90, 110, 130, 150, 170, 190),
        hrf = HRF_SPMG3,  # 3 basis functions
        span = 30
      )
    )
  )
)

cat("Multi-basis result dimensions:", dim(beta_spmg3), "\n")
cat("  - Rows:", nrow(beta_spmg3), "(10 trials × 3 basis functions)\n")
cat("  - Cols:", ncol(beta_spmg3), "voxels\n")

# Example 3: Fractional ridge regularization
# -------------------------------------------
cat("\n=== Example 3: Fractional Ridge Regularization ===\n")

beta_ridge <- lss(
  Y, X = NULL,
  method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(
        onsets = c(10, 30, 50, 70, 90, 110, 130, 150, 170, 190),
        hrf = HRF_SPMG1
      )
    ),
    ridge_mode = "fractional",
    ridge_x = 0.05,  # 5% of mean design energy
    ridge_b = 0.05   # 5% of mean aggregator energy
  )
)

cat("Ridge-regularized result dimensions:", dim(beta_ridge), "\n")
cat("Ridge parameters: 5% fractional (scales with design energy)\n")

# Example 4: With other conditions as nuisances
# ----------------------------------------------
cat("\n=== Example 4: Multiple Conditions ===\n")

beta_multi_cond <- lss(
  Y, X = NULL,
  method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(
        onsets = c(10, 50, 90, 130, 170),  # Condition A
        hrf = HRF_SPMG1
      ),
      others = list(
        list(onsets = c(30, 70, 110, 150, 190)),  # Condition B
        list(onsets = c(20, 60, 100, 140, 180))   # Condition C
      )
    )
  )
)

cat("Multi-condition result dimensions:", dim(beta_multi_cond), "\n")
cat("  - Analyzing condition A (5 trials)\n")
cat("  - Conditions B and C included as nuisance regressors\n")

# Example 5: FIR basis (Finite Impulse Response)
# -----------------------------------------------
cat("\n=== Example 5: FIR Basis ===\n")

# Create FIR HRF with 10 time bins
fir_hrf <- hrf_fir_generator(nbasis = 10, span = 20)

beta_fir <- lss(
  Y, X = NULL,
  method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(
        onsets = c(10, 50, 90, 130, 170),
        hrf = fir_hrf
      )
    )
  )
)

cat("FIR result dimensions:", dim(beta_fir), "\n")
cat("  - Rows:", nrow(beta_fir), "(5 trials × 10 time bins)\n")
cat("  - Cols:", ncol(beta_fir), "voxels\n")

# Example 6: Standard errors
# ---------------------------
cat("\n=== Example 6: Standard Errors ===\n")

result_with_se <- lss(
  Y, X = NULL,
  method = "oasis",
  oasis = list(
    design_spec = list(
      sframe = sframe,
      cond = list(
        onsets = c(10, 50, 90, 130, 170),
        hrf = HRF_SPMG1
      )
    ),
    return_se = TRUE
  )
)

cat("Results with standard errors:\n")
cat("  - Beta dimensions:", dim(result_with_se$beta), "\n")
cat("  - SE dimensions:", dim(result_with_se$se), "\n")
cat("  - Mean SE across voxels:", mean(result_with_se$se), "\n")

# Example 7: Using prebuilt design matrix
# ----------------------------------------
cat("\n=== Example 7: Prebuilt Design Matrix ===\n")

# Build your own trial-wise design matrix
times <- samples(sframe, global = TRUE)
n_trials <- 8
X_manual <- matrix(0, T, n_trials)

# Simple boxcar design for demonstration
for (i in 1:n_trials) {
  onset_time <- i * 30
  X_manual[onset_time:(onset_time + 10), i] <- 1
}

# Use OASIS with prebuilt X
beta_manual <- lss(
  Y, X = X_manual,
  method = "oasis",
  oasis = list(
    ridge_mode = "fractional",
    ridge_x = 0.02,
    ridge_b = 0.02
  )
)

cat("Prebuilt design result dimensions:", dim(beta_manual), "\n")

cat("\n=== OASIS Examples Complete ===\n")
cat("The OASIS method provides fast, exact LSS estimates with:\n")
cat("  - Support for any fmrihrf HRF type\n")
cat("  - Automatic handling of single and multi-basis HRFs\n")
cat("  - Fractional ridge regularization\n")
cat("  - Standard error calculation\n")
cat("  - Efficient blocked computation\n")