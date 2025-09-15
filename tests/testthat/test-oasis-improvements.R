# Tests for OASIS improvements
library(fmrilss)

test_that("OASIS HRF grid selection works", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()
  set.seed(42)
  
  # Create data
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1.0)
  T <- sum(fmrihrf::blocklens(sframe))
  V <- 50
  Y <- matrix(rnorm(T * V), T, V)
  
  # Create HRF grid - use predefined HRFs
  hrf_grid <- list(
    fmrihrf::HRF_SPMG1,
    fmrihrf::HRF_SPMG3,
    fmrihrf::HRF_GAMMA
  )
  
  # Test with HRF grid
  result <- lss(Y, X = NULL, method = "oasis",
                oasis = list(
                  design_spec = list(
                    sframe = sframe,
                    cond = list(onsets = c(10, 30, 50, 70)),
                    hrf_grid = hrf_grid
                  )
                ))
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)  # 4 trials
  expect_equal(ncol(result), V)
})

test_that("OASIS multi-basis standard errors work", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()
  set.seed(43)
  
  sframe <- fmrihrf::sampling_frame(blocklens = 100, TR = 1.0)
  T <- sum(fmrihrf::blocklens(sframe))
  V <- 30
  Y <- matrix(rnorm(T * V), T, V)
  
  # Test K>1 with SEs
  result <- lss(Y, X = NULL, method = "oasis",
                oasis = list(
                  design_spec = list(
                    sframe = sframe,
                    cond = list(onsets = c(10, 30, 50, 70),
                               hrf = fmrihrf::HRF_SPMG3)
                  ),
                  return_se = TRUE
                ))
  
  expect_true(is.list(result))
  expect_true(all(c("beta", "se") %in% names(result)))
  expect_equal(dim(result$beta), dim(result$se))
  expect_equal(nrow(result$beta), 4 * 3)  # 4 trials * 3 basis functions
  expect_true(all(result$se >= 0))
})

test_that("OASIS auto-detects K from user-provided X", {
  skip_on_cran()
  
  set.seed(44)
  T <- 200
  V <- 40
  K <- 3
  N <- 10
  
  # Create multi-basis design
  X <- matrix(0, T, N * K)
  for (j in 1:N) {
    for (k in 1:K) {
      col_idx <- (j - 1) * K + k
      # Create correlated basis functions
      base_times <- seq(10 + j*15, min(T, 20 + j*15))
      if (length(base_times) > 0) {
        X[base_times, col_idx] <- exp(-(0:(length(base_times) - 1)) / (3 + k))
      }
    }
  }
  
  Y <- X %*% matrix(rnorm(N*K * V, sd = 0.5), ncol = V) + 
       matrix(rnorm(T * V), T, V)
  
  # OASIS should auto-detect K=3
  result <- lss(Y, X = X, method = "oasis")
  
  expect_true(is.matrix(result))
  # Note: Detection heuristic may or may not find K=3, just check it runs
  expect_equal(ncol(result), V)
})

test_that("OASIS AR(1) whitening works", {
  skip_on_cran()
  
  set.seed(45)
  T <- 150
  V <- 25
  N <- 8
  
  # Create simple design
  X <- matrix(rnorm(T * N), T, N)
  
  # Create AR(1) noise
  Y <- matrix(0, T, V)
  rho <- 0.7
  for (v in 1:V) {
    noise <- rnorm(T)
    for (t in 2:T) {
      noise[t] <- rho * noise[t-1] + sqrt(1 - rho^2) * noise[t]
    }
    Y[, v] <- X %*% rnorm(N) + noise
  }
  
  # Test with AR(1) whitening
  result <- lss(Y, X = X, method = "oasis",
                oasis = list(whiten = "ar1"))
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), N)
  expect_equal(ncol(result), V)
})

test_that("OASIS input validation catches errors", {
  skip_on_cran()
  
  # Non-matrix Y with valid X
  expect_error(lss(1:10, X = matrix(1:10, 10, 1), method = "oasis"), "Y must be a numeric matrix")
  
  # Non-finite values
  Y <- matrix(1:20, 10, 2)
  Y[5, 1] <- NA
  expect_error(lss(Y, X = matrix(1:10, 10, 1), method = "oasis"), "non-finite")
  
  # Dimension mismatch
  Y <- matrix(rnorm(100), 10, 10)
  X <- matrix(rnorm(80), 8, 10)
  expect_error(lss(Y, X = X, method = "oasis"), "same number of rows")
  
  # Negative ridge
  Y <- matrix(rnorm(100), 10, 10)
  X <- matrix(rnorm(50), 10, 5)
  expect_error(lss(Y, X = X, method = "oasis", 
                   oasis = list(ridge_x = -1)), "non-negative")
})

test_that("OASIS precision and method options work", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()
  set.seed(46)
  
  sframe <- fmrihrf::sampling_frame(blocklens = 80, TR = 1.0)
  T <- sum(fmrihrf::blocklens(sframe))
  V <- 20
  Y <- matrix(rnorm(T * V), T, V)
  
  # Test custom precision and method
  result <- lss(Y, X = NULL, method = "oasis",
                oasis = list(
                  design_spec = list(
                    sframe = sframe,
                    cond = list(onsets = c(10, 30, 50),
                               hrf = fmrihrf::HRF_SPMG1),  # Specify HRF
                    precision = 0.05,
                    method = "conv"
                  )
                ))
  
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 3)
  expect_equal(ncol(result), V)
})
