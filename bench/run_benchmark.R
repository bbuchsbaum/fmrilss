# Load necessary packages
library(fmrilss)
library(bench)

# --- 1. Data Simulation ---
# Generate data with realistic dimensions for fMRI analysis
set.seed(123)
n_scans    <- 200 # Number of time points
n_voxels   <- 50000 # Number of brain voxels
n_trials   <- 50  # Number of experimental trials/events
n_confounds<- 5   # Number of nuisance regressors

# Y: The fMRI data matrix (time x voxels)
Y <- matrix(rnorm(n_scans * n_voxels), n_scans, n_voxels)

# X: The trial-specific design matrix (time x trials)
# Each column represents a single trial
X <- matrix(0, n_scans, n_trials)
trial_onsets <- sample(1:(n_scans - 10), n_trials, replace = FALSE)
for (i in 1:n_trials) {
  X[trial_onsets[i]:(trial_onsets[i] + 5), i] <- 1
}

# Z: The "base" or shared design matrix (time x confounds)
# These are regressors common to all trials, like motion parameters
Z <- matrix(rnorm(n_scans * n_confounds), n_scans, n_confounds)

# Add dimnames for clarity
colnames(Y) <- paste0("V", 1:n_voxels)
colnames(X) <- paste0("Trial", 1:n_trials)
colnames(Z) <- paste0("Confound", 1:n_confounds)


# --- 2. Run Benchmark ---
# We will compare LSS methods against each other and against LSA
cat("--- Benchmarking LSS vs LSA Implementations ---\n")
print(paste("Dimensions: Voxel x Scans x Trials =", n_voxels, "x", n_scans, "x", n_trials))

# Use bench::mark to get detailed and reliable performance metrics
benchmark_results <- bench::mark(
  `lss_cpp` = fmrilss::lss(Y, X, Z, method = "cpp"),
  `lss_cpp_optimized` = fmrilss::lss(Y, X, Z, method = "cpp_optimized"),
  `lsa_r` = fmrilss::lsa(Y, X, Z, method = "r"),
  iterations = 5,
  check = FALSE # LSS and LSA produce different results, so don't check for equality
)

# --- 3. Display Results ---
cat("\n--- Benchmark Results ---\n")
print(benchmark_results)

# --- 4. Quick Comparison of Results ---
cat("\n--- Quick Results Comparison ---\n")
lss_result <- fmrilss::lss(Y, X, Z, method = "cpp_optimized")
lsa_result <- fmrilss::lsa(Y, X, Z, method = "r")

cat("LSS result dimensions:", paste(dim(lss_result), collapse = " x "), "\n")
cat("LSA result dimensions:", paste(dim(lsa_result), collapse = " x "), "\n")
cat("Results are identical:", identical(lss_result, lsa_result), "\n")
cat("Mean absolute difference:", mean(abs(lss_result - lsa_result)), "\n")

# Optional: Plot the results for a visual comparison
# library(ggplot2)
# plot(benchmark_results) + theme_bw() 