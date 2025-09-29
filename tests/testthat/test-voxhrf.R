test_that("OASIS-VOXHRF returns finite betas with expected dims", {
  skip_on_cran()
  set.seed(123)
  n_time <- 160
  n_vox  <- 20
  TR <- 1
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)

  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = TR)
  onsets <- seq(10, n_time - 20, by = 10)

  beta <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG3, span = 30)
      ),
      hrf_mode = "voxel_ridge",
      lambda_shape = 1,
      mu_rough = 0,
      shrink_global = 0,
      orient_ref = TRUE
    )
  )

  expect_true(is.matrix(beta))
  expect_equal(dim(beta), c(length(onsets), n_vox))
  expect_true(all(is.finite(beta)))
})

test_that("VOXHRF approximates canonical OASIS with strong shape prior", {
  skip_on_cran()
  set.seed(42)
  n_time <- 180
  n_vox  <- 30
  TR <- 1
  Y <- matrix(rnorm(n_time * n_vox), n_time, n_vox)
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = TR)
  onsets <- seq(12, n_time - 24, by = 12)

  # Baseline OASIS with SPMG1 (K=1)
  beta_base <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG1, span = 25)
      ),
      ridge_x = 0, ridge_b = 0
    )
  )

  # VOXHRF with SPMG1 basis and strong pull to reference
  beta_vox <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(onsets = onsets, hrf = fmrihrf::HRF_SPMG1, span = 25)
      ),
      hrf_mode = "voxel_ridge",
      lambda_shape = 1e6,  # effectively force canonical shape
      mu_rough = 0,
      shrink_global = 0,
      orient_ref = TRUE
    )
  )

  # Should match closely up to numeric tolerance
  expect_lt(max(abs(beta_base - beta_vox)), 1e-6)
})

