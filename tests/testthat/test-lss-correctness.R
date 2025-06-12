# Test the corrected LSS implementation
library(testthat)
library(fmrilss)

test_that("fmrilss implementation matches ground truth LSS", {
  
  set.seed(789)
  
  # Test parameters
  n <- 100
  T_trials <- 7
  p <- 15
  q <- 4
  lambda <- 1e-6
  
  # Create trial matrices
  X_trials <- list()
  for (t in 1:T_trials) {
    X_t <- matrix(0, n, p)
    onset <- 5 + (t-1) * 12
    if (onset + p <= n) {
      X_t[onset:(onset + p - 1), ] <- diag(p)
    }
    X_trials[[t]] <- X_t
  }
  
  # HRF
  h <- dgamma(0:(p-1), shape = 5, rate = 1.5)
  h <- h / sum(h)
  
  # Confounds
  Z <- cbind(
    1,
    (1:n) / n,
    sin(2 * pi * (1:n) / n),
    cos(2 * pi * (1:n) / n)
  )
  
  # True parameters
  true_betas <- c(1, -0.5, 2, 0, 1.5, -1, 0.8)
  true_conf <- c(5, -2, 1, 0.5)
  
  # Generate data
  y <- Z %*% true_conf
  for (t in 1:T_trials) {
    y <- y + X_trials[[t]] %*% h * true_betas[t]
  }
  y <- y + rnorm(n, sd = 0.2)
  
  # Prepare projection manually
  ZtZ <- crossprod(Z)
  ZtZ_reg <- ZtZ + lambda * diag(ncol(Z))
  P_proj <- diag(n) - Z %*% solve(ZtZ_reg) %*% t(Z)
  y_proj <- P_proj %*% y
  
  # METHOD 1: Direct LSS (ground truth)
  betas_direct <- numeric(T_trials)
  for (t in 1:T_trials) {
    # The regressor for all other trials is their sum.
    X_other <- rowSums(do.call(cbind, lapply(X_trials[-t], function(X) X %*% h)))
    X_full <- cbind(X_trials[[t]] %*% h, X_other, Z)
    XtX <- crossprod(X_full) + lambda * diag(ncol(X_full))
    Xty <- crossprod(X_full, y)
    betas_direct[t] <- solve(XtX, Xty)[1]
  }
  
  # METHOD 2: fmrilss implementations
  # Create convolved regressors
  C <- matrix(0, n, T_trials)
  for (t in 1:T_trials) {
    C[, t] <- X_trials[[t]] %*% h
  }
  
  Y_mat <- matrix(y, ncol = 1)
  
  # Run all lss methods
  betas_r_opt <- fmrilss::lss(Y = Y_mat, X = C, Z = Z, method = "r_optimized")
  betas_cpp_opt <- fmrilss::lss(Y = Y_mat, X = C, Z = Z, method = "cpp_optimized")
  betas_r_vec <- fmrilss::lss(Y = Y_mat, X = C, Z = Z, method = "r_vectorized")
  betas_cpp <- fmrilss::lss(Y = Y_mat, X = C, Z = Z, method = "cpp")
  betas_naive <- fmrilss::lss(Y = Y_mat, X = C, Z = Z, method = "naive")
  
  # Compare
  cat("\n=== LSS Implementation Comparison ===\n")
  cat("True betas:      ", round(true_betas, 3), "\n")
  cat("Direct LSS:      ", round(betas_direct, 3), "\n")
  cat("R Optimized:     ", round(as.vector(betas_r_opt), 3), "\n")
  cat("C++ Optimized:   ", round(as.vector(betas_cpp_opt), 3), "\n")
  cat("R Vectorized:    ", round(as.vector(betas_r_vec), 3), "\n")
  cat("C++ Standard:    ", round(as.vector(betas_cpp), 3), "\n")
  cat("R Naive:         ", round(as.vector(betas_naive), 3), "\n")
  
  cat("\n=== Differences from Direct LSS ===\n")
  cat("R Optimized diff:    ", round(betas_r_opt - betas_direct, 6), "\n")
  cat("C++ Optimized diff:  ", round(betas_cpp_opt - betas_direct, 6), "\n")
  cat("R Vectorized diff:   ", round(betas_r_vec - betas_direct, 6), "\n")
  cat("C++ Standard diff:   ", round(betas_cpp - betas_direct, 6), "\n")
  cat("R Naive diff:        ", round(betas_naive - betas_direct, 6), "\n")
  
  # The fmrilss implementation should do proper trial-wise LSS
  # It should match the direct LSS, not simultaneous estimation
  
  # Note: fmrilss uses a different numerical approach
  # so we allow for some tolerance in the comparison
  # Also handle NaN values that may occur for trials that don't fit in the time series
  
  expect_lss_equivalent <- function(betas, method_name) {
    valid_idx <- !is.na(betas) & !is.na(betas_direct)
    if (sum(valid_idx) > 0) {
      expect_lt(
        max(abs(betas[valid_idx] - betas_direct[valid_idx])), 1e-4,
        paste("fmrilss", method_name, "implementation should match direct LSS")
      )
    }
  }
  
  expect_lss_equivalent(betas_r_opt, "R Optimized")
  expect_lss_equivalent(betas_cpp_opt, "C++ Optimized")
  expect_lss_equivalent(betas_r_vec, "R Vectorized")
  expect_lss_equivalent(betas_cpp, "C++ Standard")
  expect_lss_equivalent(betas_naive, "R Naive")
  
}) 