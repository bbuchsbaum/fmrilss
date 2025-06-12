library(testthat)
library(fmrilss)

test_that("zero regressor detection works correctly", {
  
  # Test data setup
  set.seed(456)
  n <- 50
  Y <- matrix(rnorm(n * 3), n, 3)
  
  # Test 1: Normal regressors (should not warn)
  cat("\n=== Test 1: Normal regressors ===\n")
  X_normal <- matrix(0, n, 3)
  X_normal[1:15, 1] <- 1
  X_normal[16:30, 2] <- 1  
  X_normal[31:45, 3] <- 1
  colnames(X_normal) <- c("Trial1", "Trial2", "Trial3")
  
  # Should not produce warnings
  expect_silent(lss(Y, X_normal, method = "naive"))
  
  # Test 2: Zero regressor (should warn)
  cat("\n=== Test 2: Zero regressor ===\n")
  X_zero <- X_normal
  X_zero[, 3] <- 0  # Make third regressor zero
  
  # Should produce warning about zero regressor
  expect_warning(
    lss(Y, X_zero, method = "naive"),
    "Trial regressor 'Trial3' appears to be zero.*This may cause numerical issues"
  )
  
  # Test 3: Near-zero variance regressor (should warn)
  cat("\n=== Test 3: Near-zero variance regressor ===\n")
  X_near_zero <- X_normal
  X_near_zero[, 2] <- rep(1e-6, n)  # Small constant value (non-zero but low variance)
  
  # Should produce warning about low variance
  expect_warning(
    lss(Y, X_near_zero, method = "naive"),
    "Trial regressor 'Trial2' has very low variance.*may cause numerical instability"
  )
  
  # Test 4: Multiple problematic regressors
  cat("\n=== Test 4: Multiple problematic regressors ===\n")
  X_multiple <- matrix(0, n, 4)
  X_multiple[1:10, 1] <- 1   # Normal regressor
  X_multiple[, 2] <- 0       # Zero regressor
  X_multiple[, 3] <- 1e-6    # Near-zero regressor  
  X_multiple[11:20, 4] <- 1  # Another normal regressor
  colnames(X_multiple) <- c("Good1", "Zero", "NearZero", "Good2")
  
  # Should produce multiple warnings
  expect_warning(
    lss(Y, X_multiple, method = "naive"),
    "Trial regressor 'Zero' appears to be zero"
  )
  
  # Test 5: Verify the function still works despite warnings
  cat("\n=== Test 5: Function still works with problematic regressors ===\n")
  suppressWarnings({
    result_with_issues <- lss(Y, X_zero, method = "naive")
  })
  
  # Should still return a result matrix of correct dimensions
  expect_true(is.matrix(result_with_issues))
  expect_equal(nrow(result_with_issues), ncol(X_zero))
  expect_equal(ncol(result_with_issues), ncol(Y))
  
  cat("\n=== All zero regressor tests completed ===\n")
})

test_that("zero regressor detection works across methods", {
  
  # Test that warnings work for all LSS methods
  set.seed(789)
  n <- 30
  Y <- matrix(rnorm(n * 2), n, 2)
  X <- matrix(0, n, 3)
  X[1:10, 1] <- 1
  X[11:20, 2] <- 1
  # X[, 3] remains zero
  colnames(X) <- c("Good1", "Good2", "Zero")
  
  methods_to_test <- c("naive", "r_vectorized", "r_optimized", "cpp", "cpp_optimized")
  
  for (method in methods_to_test) {
    cat(sprintf("\nTesting method: %s\n", method))
    expect_warning(
      lss(Y, X, method = method),
      "Trial regressor 'Zero' appears to be zero",
      info = sprintf("Method %s should warn about zero regressors", method)
    )
  }
}) 