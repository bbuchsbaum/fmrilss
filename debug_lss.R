library(fmrilss)

# Simple test case
set.seed(123)
n_timepoints <- 20
n_trials <- 3
n_voxels <- 2

# Create design matrices
dmat_base <- cbind(intercept = 1, trend = 1:n_timepoints)
dmat_ran <- matrix(0, n_timepoints, n_trials)

# Simple non-overlapping trials
dmat_ran[1:5, 1] <- 1
dmat_ran[8:12, 2] <- 1  
dmat_ran[15:19, 3] <- 1

# Simple data
Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)

bdes <- list(
  dmat_base = dmat_base,
  dmat_ran = dmat_ran,
  dmat_fixed = NULL,
  fixed_ind = NULL
)

cat("=== DEBUGGING LSS IMPLEMENTATIONS ===\n")

# Let's manually implement what the naive SHOULD be doing
cat("\n--- Manual Implementation for Validation ---\n")

# Step 1: Project out confounds 
X_base_fixed <- dmat_base
qrX <- qr(X_base_fixed)
Q <- diag(n_timepoints) - tcrossprod(qr.Q(qrX))

# Apply projection
data_projected <- Q %*% Y
dmat_ran_projected <- Q %*% dmat_ran

cat("Projected data shape:", dim(data_projected), "\n")
cat("Projected design shape:", dim(dmat_ran_projected), "\n")

# Manual trial-by-trial fitting
manual_result <- matrix(NA, n_trials, n_voxels)

for (i in 1:n_trials) {
  cat(sprintf("\nTrial %d:\n", i))
  
  # Current trial regressor
  ci <- dmat_ran_projected[, i, drop = FALSE]
  
  # Other trials regressor  
  other_indices <- setdiff(1:n_trials, i)
  bi <- rowSums(dmat_ran_projected[, other_indices, drop = FALSE])
  
  cat("ci sum:", sum(ci^2), "\n")
  cat("bi sum:", sum(bi^2), "\n")
  cat("ci'bi:", sum(ci * bi), "\n")
  
  # Design matrix for this trial
  X_trial <- cbind(ci, bi)
  
  # Fit model
  beta_all <- MASS::ginv(X_trial) %*% data_projected
  manual_result[i, ] <- beta_all[1, ]  # First column is the trial of interest
  
  cat("Beta for trial", i, ":", beta_all[1, ], "\n")
}

cat("\nManual result:\n")
print(manual_result)

# Method 1: Naive 
cat("\n--- Naive Method ---\n")

# Let me step through the naive implementation manually to debug
n_timepoints <- nrow(Y)
n_voxels <- ncol(Y)
dmat_base_local <- as.matrix(bdes$dmat_base)
dmat_ran_local <- as.matrix(bdes$dmat_ran)
n_events <- ncol(dmat_ran_local)
X_base_fixed_local <- dmat_base_local

cat("n_timepoints:", n_timepoints, "\n")
cat("n_events:", n_events, "\n")
cat("X_base_fixed shape:", dim(X_base_fixed_local), "\n")

# Project out confounds
qrX_local <- qr(X_base_fixed_local)
Q_naive <- diag(n_timepoints) - tcrossprod(qr.Q(qrX_local))
data_projected_naive <- Q_naive %*% Y
dmat_ran_projected_naive <- Q_naive %*% dmat_ran_local

cat("Q matrix difference:", max(abs(Q_naive - Q)), "\n")
cat("Data projection difference:", max(abs(data_projected_naive - data_projected)), "\n")
cat("Design projection difference:", max(abs(dmat_ran_projected_naive - dmat_ran_projected)), "\n")

# Initialize results
beta_matrix_debug <- matrix(NA, nrow = n_events, ncol = n_voxels)

for (i in seq_len(n_events)) {
  trial_regressor <- dmat_ran_projected_naive[, i, drop = FALSE]
  other_trials_indices <- setdiff(seq_len(n_events), i)
  other_trials_regressor <- rowSums(dmat_ran_projected_naive[, other_trials_indices, drop = FALSE])
  X_trial <- cbind(trial_regressor, other_trials_regressor)
  
  cat("Trial", i, "X_trial shape:", dim(X_trial), "\n")
  cat("Trial", i, "ginv(X_trial) shape:", dim(MASS::ginv(X_trial)), "\n")
  
  beta_all <- MASS::ginv(X_trial) %*% data_projected_naive
  cat("Trial", i, "beta_all shape:", dim(beta_all), "\n")
  cat("Trial", i, "beta_all[1,]:", beta_all[1, ], "\n")
  
  beta_matrix_debug[i, ] <- beta_all[1, ]
}

cat("Debug naive result:\n")
print(beta_matrix_debug)

result_naive <- lss_naive(Y = Y, bdes = bdes)
cat("Actual naive result:\n")
print(result_naive)

# Method 2: R vectorized
cat("\n--- R Vectorized Method ---\n")
result_r <- lss(Y = Y, bdes = bdes, use_cpp = FALSE)
cat("R vectorized result:\n")
print(result_r)

# Method 3: C++
cat("\n--- C++ Method ---\n") 
result_cpp <- lss(Y = Y, bdes = bdes, use_cpp = TRUE)
cat("C++ result:\n")
print(result_cpp)

cat("\n--- Differences ---\n")
cat("Max diff (naive vs R):", max(abs(result_naive - result_r)), "\n")
cat("Max diff (naive vs C++):", max(abs(result_naive - result_cpp)), "\n")
cat("Max diff (R vs C++):", max(abs(result_r - result_cpp)), "\n")

cat("\nDifference manual vs naive:", max(abs(manual_result - result_naive)), "\n")
cat("Difference manual vs R:", max(abs(manual_result - result_r)), "\n") 