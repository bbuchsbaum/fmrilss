# Unit tests for glmsingle() building blocks against GLMsingle semantics.

test_that("HRF library and canonical HRF match GLMsingle", {
  ref <- utils::read.csv(testthat::test_path("fixtures", "glmsingle_hrf_reference.csv"))
  for (key in unique(paste(ref$stimdur, ref$tr))) {
    sub <- ref[paste(ref$stimdur, ref$tr) == key, ]
    sd <- sub$stimdur[1]; tr <- sub$tr[1]
    lib <- sub[sub$kind == "library", ]
    L <- glmsingle_hrf_library(sd, tr)
    expect_equal(dim(L), c(max(lib$time), max(lib$hrf)))
    expect_equal(L[cbind(lib$time, lib$hrf)], lib$value, tolerance = 1e-6)
    can <- sub[sub$kind == "canonical", ]
    expect_equal(glmsingle_hrf(sd, tr), can$value, tolerance = 1e-6)
  }
})

test_that("nuisance span follows GLMsingle's normalised-Gram rank rule", {
  set.seed(3)
  q <- qr.Q(qr(matrix(rnorm(40 * 3), 40)))
  Z <- cbind(q[, 1], q[, 2], q[, 2] + 1e-9 * q[, 3])
  U <- .glms_span(Z)
  expect_equal(ncol(U), 2L)
  P_ref <- ref_projection(Z)
  expect_equal(diag(40) - tcrossprod(U), P_ref, tolerance = 1e-6)
  # full-rank case: projector equals the dense reference
  M <- cbind(1, matrix(rnorm(40 * 4), 40))
  expect_equal(diag(40) - tcrossprod(.glms_span(M)), ref_projection(M), tolerance = 1e-10)
})

test_that("compiled cross-validation equals literal calcbadness", {
  set.seed(11)
  runs <- 6; per <- 9; ncond <- 7; V <- 15
  validcolumns <- lapply(seq_len(runs), function(r) (r - 1) * per + seq_len(per))
  stimix <- lapply(seq_len(runs), function(r) sample(ncond, per, replace = TRUE))
  session <- c(1, 1, 1, 2, 2, 2)
  res <- lapply(1:5, function(k) matrix(rnorm(V * runs * per, 0.5, 1 + k / 5), V))
  res[[1]][1, validcolumns[[1]]] <- res[[1]][1, validcolumns[[2]]] <- res[[1]][1, validcolumns[[3]]] <- 2
  geom <- list(n_trials = runs * per, validcolumns = validcolumns,
               session = session, stimorder = unlist(stimix),
               trial_run = rep(seq_len(runs), each = per))
  schemes <- list(as.list(seq_len(runs)), list(1:2, 3:4, 5:6))
  for (xv in schemes) {
    geom$xval_scheme <- xv
    for (mode in c("python", "zero")) {
      ref <- ref_calcbadness(xv, validcolumns, stimix, res, session, python = mode == "python")
      cv <- .glms_cv_compile(geom, t(res[[1]]), mode)
      got <- cbind(.glms_cv_loss_ref(cv, t(res[[1]])[cv$used, , drop = FALSE]),
                   vapply(res[-1], function(r) .glms_cv_loss(cv, t(r)[cv$used, , drop = FALSE]),
                          numeric(V)))
      expect_equal(got, ref, tolerance = 1e-12)
    }
  }
})

test_that("PC-count rule matches select_noise_regressors", {
  expect_equal(.glms_select_pcs(c(-10, -8, -5, -4.9, -6), 1.05), 2L)
  expect_equal(.glms_select_pcs(c(-10, -11, -12), 1.05), 0L)
  expect_equal(.glms_select_pcs(c(-10, -9, -8, -7), 1.05), 3L)
})

test_that("autoscale reproduces olsmatrix, including constant candidates", {
  set.seed(5)
  cand <- cbind(matrix(rnorm(30 * 4), 30), 2, 0)
  ref <- matrix(rnorm(30 * 6), 30)
  got <- .glms_autoscale(cand, ref)
  for (v in seq_len(ncol(cand))) {
    X <- cbind(cand[, v], 1)
    h <- drop(ref_olsmatrix(X) %*% ref[, v])
    if (h[1] < 0) h <- c(1, 0)
    expect_equal(unname(got$h[v, ]), h, tolerance = 1e-10)
    expect_equal(got$fitted[, v], drop(X %*% h), tolerance = 1e-10)
  }
})

test_that("fracridge alpha mapping reproduces fracridge on a stacked design", {
  set.seed(9)
  n1 <- 60; n2 <- 50
  X1 <- matrix(rnorm(n1 * 5), n1); X2 <- matrix(rnorm(n2 * 4), n2)
  y1 <- matrix(rnorm(n1 * 3), n1); y2 <- matrix(rnorm(n2 * 3), n2)
  nuis <- function(n) list(Qp = matrix(1 / sqrt(n), n, 1), E = list(k0 = matrix(0, n, 0)))
  ctr <- function(m) sweep(m, 2L, colMeans(m))
  st <- list(.glms_design_stats_x(X1, nuis(n1), "k0", solver = FALSE, spectral = TRUE),
             .glms_design_stats_x(X2, nuis(n2), "k0", solver = FALSE, spectral = TRUE))
  dt <- list(.glms_data_stats(y1, nuis(n1)), .glms_data_stats(y2, nuis(n2)))
  sp <- .glms_spectral(st, dt)
  Xs <- rbind(cbind(ctr(X1), matrix(0, n1, 4)), cbind(matrix(0, n2, 5), ctr(X2)))
  ys <- rbind(ctr(y1), ctr(y2))
  for (f in c(1, 0.7, 0.3, 0.05)) {
    a <- .glms_frac_alphas(sp, f)
    expect_equal(.glms_ridge_coef(sp, .glms_shrink(sp, a[1, ])), ref_fracridge(Xs, ys, f),
                 tolerance = 1e-9)
  }
})
