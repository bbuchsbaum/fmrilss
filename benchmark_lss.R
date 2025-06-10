#!/usr/bin/env Rscript

# LSS Method Benchmarking Script
# Tests speed of naive, R vectorized, and C++ implementations

#' USAGE GUIDE:
#' 
#' 1. Quick benchmark (interactive use):
#'    source("benchmark_lss.R")
#'    quick_benchmark()
#' 
#' 2. Full benchmark suite:
#'    source("benchmark_lss.R")
#'    results <- main_benchmark()
#' 
#' 3. Run from command line:
#'    Rscript benchmark_lss.R
#' 
#' 4. Plot existing results:
#'    plot_benchmark_results()
#' 
#' 5. Custom benchmark:
#'    run_benchmark(n_timepoints = 300, n_trials = 25, n_voxels = 10000)

library(microbenchmark)
library(ggplot2)
library(dplyr)
library(fmrilss)

# Helper function to create synthetic fMRI data and design matrices
create_test_data <- function(n_timepoints, n_trials, n_voxels, n_confounds = 2) {
  set.seed(123)  # For reproducible results
  
  # Create base design matrix (intercept + confounds)
  dmat_base <- cbind(
    intercept = 1,
    matrix(rnorm(n_timepoints * (n_confounds - 1)), n_timepoints, n_confounds - 1)
  )
  
  # Create trial design matrix (sparse, realistic timing)
  dmat_ran <- matrix(0, n_timepoints, n_trials)
  trial_length <- 6  # trials last 6 timepoints
  min_gap <- 8       # minimum gap between trials
  
  for (i in seq_len(n_trials)) {
    # Find valid start position
    attempts <- 0
    repeat {
      start_pos <- sample(1:(n_timepoints - trial_length), 1)
      end_pos <- start_pos + trial_length - 1
      
      # Check if this overlaps with existing trials
      if (sum(dmat_ran[start_pos:end_pos, ]) == 0) {
        dmat_ran[start_pos:end_pos, i] <- 1
        break
      }
      
      attempts <- attempts + 1
      if (attempts > 100) {
        # Fallback: just place it somewhere
        start_pos <- ((i - 1) * min_gap) %% (n_timepoints - trial_length) + 1
        end_pos <- start_pos + trial_length - 1
        dmat_ran[start_pos:end_pos, i] <- 1
        break
      }
    }
  }
  
  # Create synthetic fMRI data with some signal
  Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
  
  # Add some realistic signal
  true_betas <- matrix(rnorm(n_trials * n_voxels, 0, 0.5), n_trials, n_voxels)
  for (i in seq_len(n_trials)) {
    Y <- Y + dmat_ran[, i] %*% matrix(true_betas[i, ], 1, n_voxels)
  }
  
  # Create design list
  bdes <- list(
    dmat_base = dmat_base,
    dmat_ran = dmat_ran,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  
  list(Y = Y, bdes = bdes)
}

# Function to run benchmark for a specific configuration
run_benchmark <- function(n_timepoints, n_trials, n_voxels, n_confounds = 2, 
                         times = 10, verbose = TRUE) {
  
  if (verbose) {
    cat(sprintf("Benchmarking: %d timepoints, %d trials, %d voxels, %d confounds\n",
                n_timepoints, n_trials, n_voxels, n_confounds))
  }
  
  # Create test data
  test_data <- create_test_data(n_timepoints, n_trials, n_voxels, n_confounds)
  Y <- test_data$Y
  bdes <- test_data$bdes
  
  # Run benchmark
  if (n_trials == 1) {
    # For single trial, naive and vectorized should be very similar
    mb_result <- microbenchmark(
      naive = lss_naive(Y = Y, bdes = bdes),
      cpp = lss(Y = Y, bdes = bdes, use_cpp = TRUE),
      times = times,
      unit = "ms"
    )
  } else {
    mb_result <- microbenchmark(
      naive = lss_naive(Y = Y, bdes = bdes),
      r_vectorized = lss(Y = Y, bdes = bdes, use_cpp = FALSE),
      cpp = lss(Y = Y, bdes = bdes, use_cpp = TRUE),
      times = times,
      unit = "ms"
    )
  }
  
  # Add configuration info
  mb_result$n_timepoints <- n_timepoints
  mb_result$n_trials <- n_trials
  mb_result$n_voxels <- n_voxels
  mb_result$n_confounds <- n_confounds
  
  return(mb_result)
}

# Function to validate that all methods give the same results
validate_equivalence <- function(n_timepoints = 50, n_trials = 10, n_voxels = 100) {
  cat("Validating equivalence of all methods...\n")
  
  test_data <- create_test_data(n_timepoints, n_trials, n_voxels)
  Y <- test_data$Y
  bdes <- test_data$bdes
  
  # Run all methods
  result_naive <- lss_naive(Y = Y, bdes = bdes)
  result_r <- lss(Y = Y, bdes = bdes, use_cpp = FALSE)
  result_cpp <- lss(Y = Y, bdes = bdes, use_cpp = TRUE)
  
  # Check equivalence
  tol <- 1e-10
  
  diff_naive_r <- max(abs(result_naive - result_r))
  diff_naive_cpp <- max(abs(result_naive - result_cpp))
  diff_r_cpp <- max(abs(result_r - result_cpp))
  
  cat(sprintf("Max difference (naive vs R vectorized): %.2e\n", diff_naive_r))
  cat(sprintf("Max difference (naive vs C++): %.2e\n", diff_naive_cpp))
  cat(sprintf("Max difference (R vectorized vs C++): %.2e\n", diff_r_cpp))
  
  if (diff_naive_r < tol && diff_naive_cpp < tol && diff_r_cpp < tol) {
    cat("✓ All methods are numerically equivalent\n\n")
    return(TRUE)
  } else {
    cat("✗ Methods are NOT equivalent - there may be implementation errors\n\n")
    return(FALSE)
  }
}

# Main benchmarking function
main_benchmark <- function() {
  cat("=== LSS Implementation Benchmark ===\n\n")
  
  # First validate equivalence
  if (!validate_equivalence()) {
    stop("Equivalence validation failed. Fix implementations before benchmarking.")
  }
  
  # Define test configurations
  configs <- list(
    # Small datasets
    list(n_timepoints = 100, n_trials = 5, n_voxels = 100),
    list(n_timepoints = 100, n_trials = 10, n_voxels = 500),
    list(n_timepoints = 200, n_trials = 15, n_voxels = 1000),
    
    # Medium datasets  
    list(n_timepoints = 300, n_trials = 20, n_voxels = 5000),
    list(n_timepoints = 400, n_trials = 30, n_voxels = 10000),
    
    # Large datasets (more realistic for fMRI)
    list(n_timepoints = 500, n_trials = 40, n_voxels = 20000),
    list(n_timepoints = 600, n_trials = 50, n_voxels = 50000),
    
    # Edge cases
    list(n_timepoints = 100, n_trials = 1, n_voxels = 1000),   # Single trial
    list(n_timepoints = 1000, n_trials = 100, n_voxels = 1000), # Many trials
    list(n_timepoints = 200, n_trials = 20, n_voxels = 100000)  # Many voxels
  )
  
  # Run benchmarks
  all_results <- list()
  
  for (i in seq_along(configs)) {
    config <- configs[[i]]
    
    # Adjust number of benchmark repetitions based on expected runtime
    times <- if (config$n_voxels > 20000 || config$n_trials > 50) 3 else 10
    
    result <- do.call(run_benchmark, c(config, list(times = times)))
    all_results[[i]] <- result
  }
  
  # Combine all results
  combined_results <- do.call(rbind, all_results)
  
  # Create summary table
  summary_stats <- combined_results %>%
    group_by(n_timepoints, n_trials, n_voxels, expr) %>%
    summarise(
      median_ms = median(time) / 1e6,
      mean_ms = mean(time) / 1e6,
      min_ms = min(time) / 1e6,
      max_ms = max(time) / 1e6,
      .groups = "drop"
    ) %>%
    arrange(n_timepoints, n_trials, n_voxels, expr)
  
  # Print summary
  cat("\n=== BENCHMARK SUMMARY ===\n")
  print(summary_stats, n = Inf)
  
  # Calculate speedup ratios
  cat("\n=== SPEEDUP ANALYSIS ===\n")
  speedup_analysis <- summary_stats %>%
    group_by(n_timepoints, n_trials, n_voxels) %>%
    mutate(
      config = paste0("T:", n_timepoints, " Tr:", n_trials, " V:", n_voxels)
    ) %>%
    select(config, expr, median_ms) %>%
    tidyr::pivot_wider(names_from = expr, values_from = median_ms) %>%
    mutate(
      r_vs_naive = ifelse(!is.na(r_vectorized) & !is.na(naive), 
                         naive / r_vectorized, NA),
      cpp_vs_naive = ifelse(!is.na(cpp) & !is.na(naive), 
                           naive / cpp, NA),
      cpp_vs_r = ifelse(!is.na(cpp) & !is.na(r_vectorized), 
                       r_vectorized / cpp, NA)
    )
  
  print(speedup_analysis %>% 
        select(config, r_vs_naive, cpp_vs_naive, cpp_vs_r), n = Inf)
  
  # Save results
  saveRDS(combined_results, "lss_benchmark_results.rds")
  write.csv(summary_stats, "lss_benchmark_summary.csv", row.names = FALSE)
  
  cat("\n=== BENCHMARK COMPLETE ===\n")
  cat("Results saved to: lss_benchmark_results.rds and lss_benchmark_summary.csv\n")
  
  return(list(
    raw_results = combined_results,
    summary = summary_stats,
    speedup = speedup_analysis
  ))
}

# Run the benchmark if this script is executed directly
if (!interactive()) {
  main_benchmark()
}

# Function to create visualization of benchmark results
plot_benchmark_results <- function(benchmark_results = NULL, save_plot = TRUE) {
  if (is.null(benchmark_results)) {
    if (file.exists("lss_benchmark_results.rds")) {
      cat("Loading saved benchmark results...\n")
      benchmark_results <- readRDS("lss_benchmark_results.rds")
    } else {
      stop("No benchmark results provided and no saved results found. Run main_benchmark() first.")
    }
  }
  
  # Convert to data frame if it's a list
  if (is.list(benchmark_results) && "raw_results" %in% names(benchmark_results)) {
    benchmark_results <- benchmark_results$raw_results
  }
  
  # Create summary data for plotting
  plot_data <- benchmark_results %>%
    mutate(
      time_ms = time / 1e6,
      config = paste0(n_timepoints, "×", n_trials, "×", n_voxels),
      data_size = n_timepoints * n_trials * n_voxels
    )
  
  # Create speed comparison plot
  p1 <- ggplot(plot_data, aes(x = reorder(config, data_size), y = time_ms, fill = expr)) +
    geom_boxplot() +
    scale_y_log10() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "LSS Implementation Speed Comparison",
      subtitle = "Configuration format: timepoints×trials×voxels",
      x = "Data Configuration",
      y = "Time (ms, log scale)",
      fill = "Method"
    ) +
    scale_fill_manual(
      values = c("naive" = "#FF6B6B", "r_vectorized" = "#4ECDC4", "cpp" = "#45B7D1"),
      labels = c("naive" = "Naive", "r_vectorized" = "R Vectorized", "cpp" = "C++")
    )
  
  # Create speedup plot
  speedup_data <- plot_data %>%
    group_by(config, n_timepoints, n_trials, n_voxels, data_size) %>%
    summarise(
      median_time = median(time_ms),
      method = expr,
      .groups = "keep"
    ) %>%
    ungroup() %>%
    tidyr::pivot_wider(names_from = method, values_from = median_time) %>%
    mutate(
      r_speedup = naive / r_vectorized,
      cpp_speedup = naive / cpp,
      cpp_vs_r = r_vectorized / cpp
    ) %>%
    select(config, data_size, r_speedup, cpp_speedup, cpp_vs_r) %>%
    tidyr::pivot_longer(cols = c(r_speedup, cpp_speedup, cpp_vs_r),
                       names_to = "comparison", values_to = "speedup") %>%
    filter(!is.na(speedup))
  
  p2 <- ggplot(speedup_data, aes(x = reorder(config, data_size), y = speedup, 
                                fill = comparison)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(
      title = "Speedup Ratios",
      subtitle = "Higher values indicate better performance of the optimized method",
      x = "Data Configuration",
      y = "Speedup Factor",
      fill = "Comparison"
    ) +
    scale_fill_manual(
      values = c("r_speedup" = "#4ECDC4", "cpp_speedup" = "#45B7D1", "cpp_vs_r" = "#96CEB4"),
      labels = c("r_speedup" = "R vs Naive", "cpp_speedup" = "C++ vs Naive", 
                "cpp_vs_r" = "C++ vs R")
    ) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.7)
  
  # Print plots
  print(p1)
  cat("\nPress Enter to see speedup plot...")
  readline()
  print(p2)
  
  if (save_plot) {
    ggsave("lss_benchmark_comparison.png", p1, width = 12, height = 6, dpi = 300)
    ggsave("lss_benchmark_speedup.png", p2, width = 12, height = 6, dpi = 300)
    cat("\nPlots saved to lss_benchmark_comparison.png and lss_benchmark_speedup.png\n")
  }
  
  return(list(comparison_plot = p1, speedup_plot = p2))
}

# Simple benchmark function for quick testing
quick_benchmark <- function(n_timepoints = 200, n_trials = 20, n_voxels = 5000, times = 5) {
  cat("Running quick benchmark...\n")
  
  # Validate equivalence first
  if (!validate_equivalence(n_timepoints = 50, n_trials = 10, n_voxels = 100)) {
    stop("Methods are not equivalent!")
  }
  
  result <- run_benchmark(n_timepoints, n_trials, n_voxels, times = times)
  
  # Print summary
  summary <- result %>%
    group_by(expr) %>%
    summarise(
      median_ms = median(time) / 1e6,
      mean_ms = mean(time) / 1e6,
      .groups = "drop"
    )
  
  cat("\nQuick Benchmark Results:\n")
  print(summary)
  
  # Calculate speedups
  if (nrow(summary) > 1) {
    naive_time <- summary$median_ms[summary$expr == "naive"]
    r_time <- summary$median_ms[summary$expr == "r_vectorized"]
    cpp_time <- summary$median_ms[summary$expr == "cpp"]
    
    if (length(naive_time) > 0 && length(r_time) > 0) {
      cat(sprintf("\nR vectorized speedup vs naive: %.1fx\n", naive_time / r_time))
    }
    if (length(naive_time) > 0 && length(cpp_time) > 0) {
      cat(sprintf("C++ speedup vs naive: %.1fx\n", naive_time / cpp_time))
    }
    if (length(r_time) > 0 && length(cpp_time) > 0) {
      cat(sprintf("C++ speedup vs R vectorized: %.1fx\n", r_time / cpp_time))
    }
  }
  
  return(result)
} 