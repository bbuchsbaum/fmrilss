test_that("benchmark_mixed_solve runs without error", {
  set.seed(123)

  # Small test case for quick benchmarking
  n <- 40
  p <- 2
  q <- 3
  n_voxels <- 5

  X <- cbind(1, rnorm(n))
  Z <- matrix(rnorm(n * q), n, q)
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)

  # Use minimal reps for testing
  expect_output(
    result <- benchmark_mixed_solve(X, Z, Y = Y, n_reps = 2),
    "Benchmarking mixed model implementations"
  )

  # Check result structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_true("method" %in% names(result))
  expect_true("mean_time" %in% names(result))
  expect_true(all(c("standard", "optimized") %in% result$method))
})

test_that("benchmark_mixed_solve accepts custom K matrix", {
  set.seed(234)

  n <- 30
  p <- 2
  q <- 3
  n_voxels <- 4

  X <- cbind(1, rnorm(n))
  Z <- matrix(rnorm(n * q), n, q)
  K <- diag(q)  # Identity kinship matrix
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)

  expect_output(
    result <- benchmark_mixed_solve(X, Z, K = K, Y = Y, n_reps = 2),
    "Benchmarking"
  )

  expect_s3_class(result, "data.frame")
})

test_that("benchmark_mixed_solve defaults K to identity when NULL", {
  set.seed(345)

  n <- 30
  q <- 3
  n_voxels <- 4

  X <- cbind(1, rnorm(n))
  Z <- matrix(rnorm(n * q), n, q)
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)

  expect_output(
    result <- benchmark_mixed_solve(X, Z, K = NULL, Y = Y, n_reps = 2),
    "Benchmarking"
  )

  expect_s3_class(result, "data.frame")
})

test_that("benchmark_mixed_solve reports timing information", {
  set.seed(456)

  n <- 35
  q <- 2
  n_voxels <- 3

  X <- cbind(1, rnorm(n))
  Z <- matrix(rnorm(n * q), n, q)
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)

  # Capture output
  output <- capture.output(
    result <- benchmark_mixed_solve(X, Z, Y = Y, n_reps = 2)
  )

  # Check output contains expected information
  output_text <- paste(output, collapse = "\n")
  expect_true(grepl("Data:", output_text))
  expect_true(grepl("voxels", output_text))
  expect_true(grepl("Results:", output_text))
  expect_true(grepl("Speedup:", output_text))
})
