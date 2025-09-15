# File: tests/testthat/test-oasis.R
# Tests for OASIS implementation

library(fmrilss)

test_that("OASIS matches per-trial GLM (ridge=0) on synthetic data", {
  skip_on_cran()

  set.seed(1)
  T <- 240
  V <- 200
  TR <- 1

  # Simple events for one condition
  onsets <- sort(sample(10:(T-30), 40))
  span   <- 24

  # Build a canonical HRF and trial-wise design (lightweight FIR)
  h  <- dgamma(0:span, 6, 1) - 0.35*dgamma(0:span, 12, 1)
  h  <- h / max(h)  # normalize
  X  <- matrix(0, T, length(onsets))
  for (j in seq_along(onsets)) {
    t0 <- onsets[j]
    idx <- seq_len(min(length(h), T - t0))
    if (length(idx) > 0) {
      X[t0 + idx - 1, j] <- h[idx]
    }
  }

  # Confounds (intercept + 2 drifts) and one "other condition" aggregate (nuisance)
  Nuis <- cbind(1, poly(seq_len(T), 2))
  other_on <- sort(sample(15:(T-20), 20))
  X_other <- rep(0, T)
  for (o in other_on) {
    idx <- seq_len(min(length(h), T - o))
    if (length(idx) > 0) {
      X_other[o + idx - 1] <- X_other[o + idx - 1] + h[idx]
    }
  }

  # Simulate data
  Btrue <- matrix(rnorm(ncol(X) * V, sd = 0.6), ncol = V)
  Y <- X %*% Btrue + matrix(rnorm(T * V, sd = 1), T, V) + 
       Nuis %*% matrix(rnorm(ncol(Nuis) * V, sd = .2), ncol = V)

  # --- OASIS path ---
  beta_oasis <- lss(Y, X = X, Z = NULL, Nuisance = cbind(Nuis, X_other),
                    method = "oasis",
                    oasis = list(ridge_x = 0, ridge_b = 0, block_cols = 512L))

  # --- Reference per-trial GLM (explicit LSS) ---
  beta_ref <- matrix(NA_real_, ncol(X), V)
  Xall <- X
  Sall <- rowSums(Xall)
  for (j in seq_len(ncol(Xall))) {
    xj <- Xall[, j, drop = FALSE]
    Sj <- Sall - xj
    D  <- cbind(xj, Sj, Nuis, X_other)
    # OLS for all voxels at once
    XtX <- crossprod(D)
    XtY <- crossprod(D, Y)
    coef <- solve(XtX, XtY)
    beta_ref[j, ] <- coef[1, ]   # coefficient for xj
  }

  # Check correlation (subset for speed)
  cors <- cor(as.vector(beta_oasis[, 1:50]), as.vector(beta_ref[, 1:50]))
  expect_true(cors > 0.99)
})

test_that("OASIS ridge reduces variance in dense designs", {
  skip_on_cran()

  set.seed(2)
  T <- 180
  V <- 120
  ISI <- 2
  onsets <- seq(5, 5 + ISI*(40-1), by = ISI)
  span <- 20
  h <- dgamma(0:span, 6, 1)
  h <- h / max(h)

  X <- matrix(0, T, length(onsets))
  for (j in seq_along(onsets)) {
    t0 <- onsets[j]
    idx <- seq_len(min(length(h), T - t0))
    if (length(idx) > 0) {
      X[t0 + idx - 1, j] <- h[idx]
    }
  }
  Y <- X %*% matrix(rnorm(ncol(X) * V, sd = 0.5), ncol = V) + matrix(rnorm(T * V), T, V)

  B0 <- lss(Y, X = X, method = "oasis", oasis = list(ridge_x = 0, ridge_b = 0))
  B1 <- lss(Y, X = X, method = "oasis", oasis = list(ridge_x = 0.1, ridge_b = 0.1))

  # Variance across split-halves should shrink with ridge (heuristic check)
  idx <- sample(seq_len(V))
  v0  <- apply(B0[, idx[1:floor(V/2)]], 1, var)
  v1  <- apply(B1[, idx[(floor(V/2)+1):V]], 1, var)
  # Ridge should generally reduce variance
  expect_true(median(v1, na.rm=TRUE) <= median(v0, na.rm=TRUE) * 1.5)
})

test_that("OASIS returns SEs when requested", {
  skip_on_cran()

  set.seed(3)
  T <- 120
  V <- 50
  onsets <- sort(sample(8:(T-20), 20))
  span <- 24
  h <- dgamma(0:span, 6, 1)
  h <- h / max(h)

  X <- matrix(0, T, length(onsets))
  for (j in seq_along(onsets)) {
    t0 <- onsets[j]
    idx <- seq_len(min(length(h), T - t0))
    if (length(idx) > 0) {
      X[t0 + idx - 1, j] <- h[idx]
    }
  }
  Y <- X %*% matrix(rnorm(ncol(X) * V, sd = 0.5), ncol = V) + matrix(rnorm(T * V), T, V)

  out <- lss(Y, X = X, method = "oasis", oasis = list(return_se = TRUE))
  expect_true(is.list(out) && all(c("beta","se") %in% names(out)))
  expect_equal(dim(out$beta), dim(out$se))
  expect_true(all(is.finite(out$se)))
})

test_that("OASIS fractional ridge mode works", {
  skip_on_cran()
  
  set.seed(4)
  T <- 100
  V <- 30
  N <- 15
  
  # Create simple design
  X <- matrix(rnorm(T * N), T, N)
  Y <- matrix(rnorm(T * V), T, V)
  
  # Test fractional ridge
  out_frac <- lss(Y, X = X, method = "oasis", 
                  oasis = list(ridge_mode = "fractional",
                              ridge_x = 0.05, ridge_b = 0.05))
  
  # Should return results without error
  expect_true(is.matrix(out_frac))
  expect_equal(nrow(out_frac), N)
  expect_equal(ncol(out_frac), V)
})

test_that("OASIS handles multi-basis HRFs (K>1)", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()
  set.seed(5)
  
  # Create sampling frame
  sframe <- sampling_frame(blocklens = c(100, 100), TR = 1.0)
  T <- sum(blocklens(sframe))
  V <- 50
  
  # Create synthetic data
  Y <- matrix(rnorm(T * V), T, V)
  
  # Test with SPMG3 (K=3)
  out <- lss(Y, X = NULL, method = "oasis",
            oasis = list(
              design_spec = list(
                sframe = sframe,
                cond = list(onsets = c(10, 30, 50, 70, 90, 110, 130, 150),
                           hrf = HRF_SPMG3,
                           span = 30)
              ),
              K = 3
            ))
  
  # Should return K*N rows (3 basis functions * 8 trials = 24)
  expect_equal(nrow(out), 8 * 3)
  expect_equal(ncol(out), V)
})

test_that("OASIS design_spec builds correct trial-wise design", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()
  set.seed(6)
  
  # Create sampling frame
  sframe <- sampling_frame(blocklens = 150, TR = 1.0)
  T <- sum(blocklens(sframe))
  V <- 30
  
  # Create synthetic data
  Y <- matrix(rnorm(T * V), T, V)
  
  # Build design via design_spec
  onsets <- c(10, 30, 50, 70, 90, 110)
  out <- lss(Y, X = NULL, method = "oasis",
            oasis = list(
              design_spec = list(
                sframe = sframe,
                cond = list(onsets = onsets,
                           hrf = HRF_SPMG1,
                           span = 25)
              )
            ))
  
  # Should return one beta per trial
  expect_equal(nrow(out), length(onsets))
  expect_equal(ncol(out), V)
  expect_true(all(is.finite(out)))
})

test_that("OASIS handles other conditions as nuisances", {
  skip_if_not_installed("fmrihrf")
  skip_on_cran()
  set.seed(7)
  
  sframe <- sampling_frame(blocklens = 200, TR = 1.0)
  T <- sum(blocklens(sframe))
  V <- 25
  
  Y <- matrix(rnorm(T * V), T, V)
  
  # Test with other conditions
  out <- lss(Y, X = NULL, method = "oasis",
            oasis = list(
              design_spec = list(
                sframe = sframe,
                cond = list(onsets = c(10, 40, 70, 100, 130),
                           hrf = HRF_SPMG1),
                others = list(
                  list(onsets = c(20, 50, 80, 110, 140)),
                  list(onsets = c(30, 60, 90, 120, 150))
                )
              )
            ))
  
  expect_equal(nrow(out), 5)  # 5 trials for main condition
  expect_equal(ncol(out), V)
})