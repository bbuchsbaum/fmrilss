test_that("Soft blending brings coordinates closer to mixed truth", {
  set.seed(2201)
  library(fmrihrf)

  Tlen <- 160; V <- 1; r <- 4
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 40, length.out = Tlen)/4),
    exp(-seq(0, 40, length.out = Tlen)/6),
    exp(-seq(0, 40, length.out = Tlen)/8),
    exp(-seq(0, 40, length.out = Tlen)/10)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

  # Design: spaced events within window
  onsets <- seq(12, 120, by = 12)
  ntrials <- length(onsets)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

  # Build per-trial design blocks (T×r each), and simulate HRF as blend of two library members
  regressors_by_trial <- vector("list", ntrials)
  for (t in seq_len(ntrials)) {
    rr_t <- fmrihrf::regressor(onsets = onsets[t], hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
    regressors_by_trial[[t]] <- fmrihrf::evaluate(rr_t, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  }

  k1 <- 2L; k2 <- 4L
  alpha_true <- 0.6 * sbhm$A[, k1] + 0.4 * sbhm$A[, k2]
  amps_true <- rnorm(ntrials, mean = 1.8, sd = 0.3)

  # Simulate Y
  Y <- matrix(rnorm(Tlen * V, sd = 0.25), Tlen, V)
  for (t in seq_len(ntrials)) {
    reg_t <- regressors_by_trial[[t]] %*% alpha_true
    Y[, 1] <- Y[, 1] + amps_true[t] * as.numeric(reg_t)
  }

  # Hard assignment
  hard <- lss_sbhm(
    Y = Y,
    sbhm = sbhm,
    design_spec = design_spec,
    prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
    match = list(topK = 1),
    return = "both"
  )

  # Soft blend across top-3
  soft <- lss_sbhm(
    Y = Y,
    sbhm = sbhm,
    design_spec = design_spec,
    prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
    match = list(topK = 3, soft_blend = TRUE),
    return = "both"
  )

  # Cosine similarity to mixed ground truth in coefficient space
  cos_sim <- function(u, v) {
    sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))
  }
  cs_hard <- cos_sim(hard$alpha_coords[, 1], alpha_true)
  cs_soft <- cos_sim(soft$alpha_coords[, 1], alpha_true)

  expect_gte(cs_soft, cs_hard)

  # Amplitude recovery correlation: soft should be at least as good
  cor_hard <- suppressWarnings(cor(hard$amplitude[, 1], amps_true))
  cor_soft <- suppressWarnings(cor(soft$amplitude[, 1], amps_true))
  expect_gte(cor_soft, cor_hard - 1e-6)
})

test_that("Blend margin controls whether blending occurs", {
  set.seed(2202)
  library(fmrihrf)

  Tlen <- 120; V <- 1; r <- 3
  sframe <- sampling_frame(blocklens = Tlen, TR = 1)
  H <- cbind(
    exp(-seq(0, 30, length.out = Tlen)/4),
    exp(-seq(0, 30, length.out = Tlen)/6),
    exp(-seq(0, 30, length.out = Tlen)/8)
  )
  sbhm <- sbhm_build(library_H = H, r = r, sframe = sframe, normalize = TRUE)
  onsets <- seq(10, 90, by = 10)
  design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

  # Simple signal from a single library HRF (should be unambiguous after ridge)
  hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
  rr <- fmrihrf::regressor(onsets = onsets, hrf = hrf_B, duration = 0, span = 30, summate = TRUE)
  Xagg <- fmrihrf::evaluate(rr, grid = sbhm$tgrid, precision = 0.1, method = "conv")
  alpha_true <- sbhm$A[, 2]
  Y <- matrix(rnorm(Tlen*V, sd = 0.2), Tlen, V)
  Y[, 1] <- Y[, 1] + as.numeric(Xagg %*% alpha_true)

  hard <- lss_sbhm(Y, sbhm, design_spec, prepass = list(ridge = list(mode = "fractional", lambda = 0.01)), return = "both")

  # soft_blend TRUE, but blend_margin = -Inf forces no blending (identical to hard)
  soft_no <- lss_sbhm(
    Y, sbhm, design_spec,
    prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
    match = list(topK = 3, soft_blend = TRUE, blend_margin = -Inf),
    return = "both"
  )

  expect_equal(soft_no$alpha_coords[, 1], hard$alpha_coords[, 1])
  expect_equal(as.numeric(soft_no$amplitude[, 1]), as.numeric(hard$amplitude[, 1]))

  # soft_blend TRUE, high blend_margin -> forces blending (may differ from hard)
  soft_yes <- lss_sbhm(
    Y, sbhm, design_spec,
    prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
    match = list(topK = 3, soft_blend = TRUE, blend_margin = 1e6),
    return = "both"
  )
  # Presence of topK outputs and mode flag
  expect_true(!is.null(soft_yes$topK_idx) && !is.null(soft_yes$weights))
  expect_equal(soft_yes$alpha_mode, "soft")
})

