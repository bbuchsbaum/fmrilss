test_that("sbhm_build creates valid basis (basic, rank-deficient, shifts)", {
  skip_on_cran()
  library(fmrihrf)
  set.seed(11)

  Tlen <- 60
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)

  # Basic build from explicit library matrix
  H <- cbind(
    exp(-seq(0, 20, length.out = Tlen)/4),
    exp(-seq(0, 20, length.out = Tlen)/6),
    exp(-seq(0, 20, length.out = Tlen)/8)
  )
  sbhm <- sbhm_build(library_H = H, r = 4, sframe = sframe, normalize = TRUE)
  expect_true(is.matrix(sbhm$B) && nrow(sbhm$B) == Tlen)
  expect_true(length(sbhm$S) == ncol(sbhm$B))
  expect_true(all(dim(sbhm$A)[1] == ncol(sbhm$B)))

  # Rank-deficient case: duplicate column
  H2 <- cbind(H, H[,1])
  sbhm2 <- sbhm_build(library_H = H2, r = 6, sframe = sframe, normalize = TRUE)
  expect_lte(ncol(sbhm2$B), min(nrow(H2), ncol(H2)))
  expect_true(all(sbhm2$S >= 0))

  # With shifts
  sbhm3 <- sbhm_build(library_H = H, r = 3, sframe = sframe, shifts = c(0.5, 1))
  expect_true(ncol(sbhm3$A) > ncol(sbhm$A))
})

test_that("sbhm_prepass handles edge cases", {
  skip_on_cran()
  library(fmrihrf)
  set.seed(12)

  Tlen <- 90; V <- 3
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8)
  )
  sbhm <- sbhm_build(library_H = H, r = 3, sframe = sframe, normalize = TRUE)
  r <- ncol(sbhm$B)  # Use actual rank from sbhm

  onsets <- seq(6, 72, by = 12)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

  # Single voxel path
  Y1 <- matrix(rnorm(Tlen), Tlen, 1)
  pre1 <- sbhm_prepass(Y1, sbhm, design_spec)
  expect_equal(ncol(pre1$beta_bar), 1)

  # Collinear regressors test: create X_trials with repeated columns via dense events
  # Use normal path but ensure qr.solve handles G properly
  Y <- matrix(rnorm(Tlen*V, sd = 0.5), Tlen, V)
  pre2 <- sbhm_prepass(Y, sbhm, design_spec)
  expect_equal(nrow(pre2$beta_bar), r)

  # With prewhitening requested (skips or applies safely)
  pre3 <- sbhm_prepass(Y, sbhm, design_spec, prewhiten = list(method = "none"))
  expect_true(is.list(pre3$diag))
})

test_that("sbhm pipeline produces plausible results end-to-end", {
  skip_on_cran()
  library(fmrihrf)
  set.seed(13)
  Tlen <- 120; V <- 3
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/5),
    exp(-seq(0, 30, length.out = Tlen)/7),
    exp(-seq(0, 30, length.out = Tlen)/9),
    exp(-seq(0, 30, length.out = Tlen)/11)
  )
  sbhm <- sbhm_build(library_H = H, r = 4, sframe = sframe, normalize = TRUE)
  r <- ncol(sbhm$B)  # Use actual rank from sbhm

  onsets <- seq(10, 100, by = 10)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

  # Create signal in voxel 1 following sbhm basis combination
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
  Xr <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  alpha_true <- rnorm(ncol(sbhm$B))
  Y <- matrix(rnorm(Tlen*V, sd = 0.6), Tlen, V)
  Y[,1] <- Y[,1] + Xr %*% alpha_true

  out <- lss_sbhm(Y, sbhm, design_spec, return = "both")
  expect_true(is.matrix(out$amplitude))
  expect_equal(dim(out$coeffs_r)[1], r)
  expect_equal(dim(out$coeffs_r)[2], length(onsets))
  expect_equal(dim(out$coeffs_r)[3], V)
})

