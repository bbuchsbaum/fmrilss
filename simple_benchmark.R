library(fmrilss)
library(microbenchmark)
library(dplyr)

# Function to create test data
create_data <- function(n_timepoints, n_trials, n_voxels) {
  set.seed(123)
  dmat_base <- cbind(1, 1:n_timepoints)
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  for(i in 1:n_trials) {
    start <- (i-1)*8 + 1
    if(start + 5 <= n_timepoints) dmat_ran[start:(start+5), i] <- 1
  }
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  bdes <- list(dmat_base = dmat_base, dmat_ran = dmat_ran, dmat_fixed = NULL, fixed_ind = NULL)
  list(Y = Y, bdes = bdes)
}

# Test configurations
configs <- list(
  list(100, 5, 500),
  list(200, 10, 1000),
  list(300, 15, 2000),
  list(400, 20, 5000),
  list(500, 25, 10000)
)

cat('LSS Implementation Benchmark Results\n')
cat('====================================\n\n')

all_results <- list()

for(i in seq_along(configs)) {
  config <- configs[[i]]
  n_timepoints <- config[[1]]
  n_trials <- config[[2]] 
  n_voxels <- config[[3]]
  
  cat(sprintf('Configuration %d: %d timepoints, %d trials, %d voxels\n', 
              i, n_timepoints, n_trials, n_voxels))
  
  data <- create_data(n_timepoints, n_trials, n_voxels)
  
  # Reduce repetitions for larger datasets
  times <- if(n_voxels > 5000) 3 else 5
  
  mb_result <- microbenchmark(
    naive = lss_naive(Y = data$Y, bdes = data$bdes),
    r_vectorized_old = lss(Y = data$Y, bdes = data$bdes, use_cpp = FALSE),
    r_vectorized_new = lss_optimized(Y = data$Y, bdes = data$bdes, use_cpp = FALSE),
    cpp = lss(Y = data$Y, bdes = data$bdes, use_cpp = TRUE),
    times = times
  )
  
  # Add configuration info
  mb_result$n_timepoints <- n_timepoints
  mb_result$n_trials <- n_trials
  mb_result$n_voxels <- n_voxels
  
  all_results[[i]] <- mb_result
  
  summary_stats <- mb_result %>%
    group_by(expr) %>%
    summarise(median_ms = median(time) / 1e6, .groups = 'drop')
  
  print(summary_stats)
  
  # Calculate speedups
  naive_time <- summary_stats$median_ms[summary_stats$expr == 'naive']
  r_old_time <- summary_stats$median_ms[summary_stats$expr == 'r_vectorized_old']
  r_new_time <- summary_stats$median_ms[summary_stats$expr == 'r_vectorized_new']
  cpp_time <- summary_stats$median_ms[summary_stats$expr == 'cpp']
  
  cat(sprintf('New R speedup vs Old R: %.1fx\n', r_old_time / r_new_time))
  cat(sprintf('New R speedup vs Naive: %.1fx\n', naive_time / r_new_time))
  cat(sprintf('C++ speedup vs New R: %.1fx\n', r_new_time / cpp_time))
  cat('\n')
}

# Combine all results
combined_results <- do.call(rbind, all_results)

# Create overall summary
cat('=== OVERALL SUMMARY ===\n')
overall_summary <- combined_results %>%
  group_by(n_timepoints, n_trials, n_voxels, expr) %>%
  summarise(
    median_ms = median(time) / 1e6,
    .groups = "drop"
  ) %>%
  arrange(n_timepoints, n_trials, n_voxels, expr)

print(overall_summary)

# Calculate average speedups
cat('\n=== AVERAGE SPEEDUPS ===\n')
speedup_summary <- overall_summary %>%
  group_by(n_timepoints, n_trials, n_voxels) %>%
  summarise(
    config = paste0(n_timepoints, "×", n_trials, "×", n_voxels),
    naive_time = median_ms[expr == "naive"],
    r_old_time = median_ms[expr == "r_vectorized_old"],
    r_new_time = median_ms[expr == "r_vectorized_new"],
    cpp_time = median_ms[expr == "cpp"],
    .groups = "drop"
  ) %>%
  mutate(
    r_new_vs_old = r_old_time / r_new_time,
    r_new_vs_naive = naive_time / r_new_time,
    cpp_vs_r_new = r_new_time / cpp_time
  )

print(speedup_summary %>% select(config, r_new_vs_old, r_new_vs_naive, cpp_vs_r_new))

cat('\nAverage speedup factors:\n')
cat(sprintf('New R vs Old R: %.1fx\n', mean(speedup_summary$r_new_vs_old, na.rm = TRUE)))
cat(sprintf('New R vs Naive: %.1fx\n', mean(speedup_summary$r_new_vs_naive, na.rm = TRUE)))
cat(sprintf('C++ vs New R: %.1fx\n', mean(speedup_summary$cpp_vs_r_new, na.rm = TRUE)))

# Save results
saveRDS(combined_results, "simple_benchmark_results.rds")
write.csv(overall_summary, "simple_benchmark_summary.csv", row.names = FALSE)
cat('\nResults saved to simple_benchmark_results.rds and simple_benchmark_summary.csv\n') 