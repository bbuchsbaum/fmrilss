# Run-blocked least squares and fractional ridge for glmsingle().
#
# After nuisance projection GLMsingle's single-trial design is block diagonal
# across runs, so every fit factorises per run. Fractional ridge keeps one
# penalty per voxel, calibrated on the pooled spectrum of all runs, exactly as
# fracridge() does on the stacked design.
#
# Notation for one run: X raw trial design, Qp orthonormal polynomial basis,
# E orthonormal basis of the remaining nuisance (extras, PCs) orthogonal to
# Qp, Ap = (I - Qp Qp') X, A = (I - E E') Ap, G = A'A, Gp = Ap'Ap
# = G + (Ap'E)(Ap'E)'. Fits use A (all nuisance removed); reported R^2 uses
# the poly-only residual of the raw-design prediction, as in GLMsingle:
#   SSE = s - 2 beta'c + beta' Gp beta,  c = Ap'y = G beta_ols + (Ap'E)(E'y).

# Voxel-independent statistics of one run's trial design for one HRF and one
# nuisance model. `solver` precomputes the whitened OLS operator; `spectral`
# the eigendecomposition used by fractional ridge.
.glms_design_stats <- function(onsets, n_time, hrf, nuis_run, k_name,
                               singular = "error", solver = TRUE,
                               spectral = FALSE) {
  .glms_design_stats_x(.glms_trial_design(onsets, n_time, hrf), nuis_run,
                       k_name, singular, solver, spectral)
}

.glms_design_stats_x <- function(X, nuis_run, k_name, singular = "error",
                                 solver = TRUE, spectral = FALSE) {
  Ap <- .glms_resid(X, nuis_run$Qp)
  E <- nuis_run$E[[k_name]]
  ApE <- if (ncol(E)) crossprod(Ap, E) else NULL
  A <- .glms_resid(Ap, E)
  G <- crossprod(A)
  out <- list(Ap = Ap, ApE = ApE, G = G, zero = colSums(X != 0) == 0,
              k = k_name, n = ncol(X))
  if (!out$n) return(out)
  if (solver) {
    # W'W = G^-1 (or G^+): z = W A'y gives beta = W'z and beta'G beta = |z|^2
    out$W <- .glms_whiten(G, out$zero, singular)
    out$P <- out$W %*% t(Ap)
    out$p1 <- rowSums(out$P)
    out$PE <- if (is.null(ApE)) NULL else out$W %*% ApE
  }
  if (spectral) {
    e <- eigen(G, symmetric = TRUE)
    out$V <- e$vectors
    out$s2 <- pmax(e$values, 0)
    out$VAp <- crossprod(e$vectors, t(Ap))
    out$v1 <- rowSums(out$VAp)
    out$VApE <- if (is.null(ApE)) NULL else crossprod(e$vectors, ApE)
  }
  out
}

# Whitening operator W (rank x n) with W'W = G^-1, following olsmatrix2():
# exactly-zero design columns get zero coefficients and a singular remaining
# system is an error unless singular = "pinv" (then W'W = G^+).
.glms_whiten <- function(G, zero, singular = "error") {
  n <- nrow(G)
  good <- !zero
  W <- matrix(0, 0L, n)
  if (!any(good)) return(W)
  Gg <- G[good, good, drop = FALSE]
  R <- tryCatch(chol(Gg), error = function(e) NULL)
  if (!is.null(R) && min(diag(R)) > sqrt(.Machine$double.eps) * max(diag(R))) {
    W <- matrix(0, sum(good), n)
    W[, good] <- t(backsolve(R, diag(sum(good))))
    return(W)
  }
  if (!identical(singular, "pinv")) {
    stop("Singular trial design: two or more trial regressors are linearly dependent. ",
         "Use singular = \"pinv\" for a minimum-norm solution.", call. = FALSE)
  }
  # counted by glmsingle(), which warns once per call
  signalCondition(structure(class = c("glms_singular", "condition"),
                            list(message = "trial design", call = NULL)))
  e <- eigen(Gg, symmetric = TRUE)
  keep <- e$values > max(e$values) * n * .Machine$double.eps
  W <- matrix(0, sum(keep), n)
  W[, good] <- t(e$vectors[, keep, drop = FALSE]) / sqrt(e$values[keep])
  W
}

# Data statistics of one run for a voxel tile, without copying or centring
# the data: column means, the poly-residual sum of squares s, and E'y for
# each nuisance model (E is orthogonal to the constant; the mean correction
# keeps products exact).
.glms_data_stats <- function(Yr, nuis_run, k_names = names(nuis_run$E)) {
  ms <- glms_col_mean_ss(Yr, .glms_nt())
  mean <- drop(ms$mean)
  QpY <- crossprod(nuis_run$Qp, Yr) - outer(colSums(nuis_run$Qp), mean)
  EY <- lapply(nuis_run$E[k_names], function(E) {
    if (ncol(E)) crossprod(E, Yr) - outer(colSums(E), mean) else NULL
  })
  list(Y = Yr, mean = mean, s = drop(ms$ss) - colSums(QpY^2), EY = EY)
}

# Data statistics restricted to voxel columns `sel`.
.glms_dats_subset <- function(dats, sel) {
  lapply(dats, function(dt) list(
    Y = dt$Y[, sel, drop = FALSE], mean = dt$mean[sel], s = dt$s[sel],
    EY = lapply(dt$EY, function(x) if (is.null(x)) NULL else x[, sel, drop = FALSE])
  ))
}

# B y for an operator B (rows orthogonal to the constant) with precomputed
# row sums b1, applied to uncentred data with the mean correction.
.glms_apply <- function(B, b1, dt) B %*% dt$Y - outer(b1, dt$mean)

# Poly-side SSE for explicit c = Ap'y (single-regressor type A fits).
.glms_sse <- function(beta, c, Gp, s) {
  s - 2 * colSums(beta * c) + colSums(beta * (Gp %*% beta))
}

.glms_r2 <- function(sse, s) {
  out <- 100 * (1 - sse / s)
  out[s == 0] <- NaN
  out
}

# Right-hand side b = A'y for a single-regressor design (type A).
.glms_rhs <- function(ds, dat) {
  c <- crossprod(ds$Ap, dat$Y) - outer(colSums(ds$Ap), dat$mean)
  EY <- dat$EY[[ds$k]]
  b <- if (is.null(ds$ApE)) c else c - ds$ApE %*% EY
  list(b = b, c = c)
}

# OLS scores for all runs of a voxel tile, from whitened coordinates
# z = W (A'y): beta'G beta = |z|^2 and u = (Ap'E)'beta = (W Ap'E)' z, so the
# poly-side SSE = s - |z|^2 - 2 u'(E'y) + u'u needs no betas. With
# keep_z = TRUE the z blocks are returned for beta reconstruction.
.glms_ols_score <- function(stats, dats, keep_z = FALSE) {
  R <- length(stats)
  V <- length(dats[[1]]$s)
  zs <- if (keep_z) vector("list", R) else NULL
  sse <- s <- matrix(0, R, V)
  for (r in seq_len(R)) {
    s[r, ] <- dats[[r]]$s
    sc <- .glms_ols_score_run(stats[[r]], dats[[r]], keep_z)
    sse[r, ] <- sc$sse
    if (keep_z) zs[r] <- list(sc$z)
  }
  list(sse = sse, s = s, z = zs)
}

# One run of .glms_ols_score().
.glms_ols_score_run <- function(st, dt, keep_z = FALSE) {
  if (!st$n) return(list(sse = dt$s, z = NULL))
  EY <- dt$EY[[st$k]]
  z <- .glms_apply(st$P, st$p1, dt)
  if (!is.null(st$PE)) z <- z - st$PE %*% EY
  q <- colSums(z^2)
  if (!is.null(st$PE)) {
    u <- crossprod(st$PE, z)
    q <- q + 2 * colSums(u * EY) - colSums(u^2)
  }
  list(sse = dt$s - q, z = if (keep_z) z)
}

# OLS betas (trials x voxels) and scores for all runs of a voxel tile.
.glms_fit_ols <- function(stats, dats) {
  sc <- .glms_ols_score(stats, dats, keep_z = TRUE)
  n <- vapply(stats, `[[`, integer(1), "n")
  beta <- matrix(0, sum(n), length(dats[[1]]$s))
  off <- 0L
  for (r in seq_along(stats)) {
    if (n[r]) beta[off + seq_len(n[r]), ] <- crossprod(stats[[r]]$W, sc$z[[r]])
    off <- off + n[r]
  }
  list(beta = beta, sse = sc$sse, s = sc$s)
}

# Spectral form of the stacked problem used by fracridge: per-run eigen
# coordinates vb = V'b, OLS coefficients a = vb / s2 (components with
# singular value < 1e-10 zeroed, as fracridge does) and the pooled spectrum.
.glms_spectral <- function(stats, dats) {
  R <- length(stats)
  vb <- a <- EYs <- vector("list", R)
  for (r in seq_len(R)) {
    st <- stats[[r]]
    if (!st$n) next
    EYs[r] <- list(dats[[r]]$EY[[st$k]])
    x <- .glms_apply(st$VAp, st$v1, dats[[r]])
    if (!is.null(st$VApE)) x <- x - st$VApE %*% EYs[[r]]
    vb[[r]] <- x
    coef <- x / st$s2
    coef[sqrt(st$s2) < 1e-10, ] <- 0
    a[[r]] <- coef
  }
  list(stats = stats, vb = vb, a = a, EY = EYs,
       s = lapply(dats, `[[`, "s"),
       s2_all = unlist(lapply(stats, `[[`, "s2")),
       a_all = do.call(rbind, a))
}

# fracridge's alpha grid for the pooled spectrum.
.glms_alpha_grid <- function(s2_all) {
  sv <- sqrt(s2_all)
  val1 <- 10e3 * max(sv)^2
  val2 <- 10e-3 * min(sv)^2
  if (val2 == 0) val2 <- 10e-3
  lo <- floor(log10(val2))
  hi <- ceiling(log10(val1))
  n <- max(0L, as.integer(ceiling((hi - lo) / 0.2)))
  c(0, 10^(lo + 0.2 * seq_len(n) - 0.2))
}

# Per-voxel alphas (fractions x voxels) for the requested fractions.
.glms_frac_alphas <- function(sp, fracs, method = "fracridge") {
  a2 <- sp$a_all^2
  if (identical(method, "exact")) {
    return(glms_frac_alpha_exact(a2, sp$s2_all, fracs, .glms_nt()))
  }
  if (!any(sp$s2_all > 0)) return(matrix(NaN, length(fracs), ncol(a2)))
  grid <- .glms_alpha_grid(sp$s2_all)
  s2 <- sqrt(sp$s2_all)^2
  sclg_sq <- (outer(grid, s2, function(g, s) s / (s + g)))^2
  # Null-space coefficients are zero, including at alpha = 0. Do not let
  # their 0/0 shrinkage contaminate every voxel's coefficient norm.
  sclg_sq[, s2 == 0] <- 0
  newlen <- sqrt(sclg_sq %*% a2)
  glms_frac_alpha_grid(newlen, grid, fracs, .glms_nt())
}

# Shrunk eigen coefficients for one alpha per voxel, per run.
.glms_shrink <- function(sp, alpha) {
  alpha[!is.finite(alpha)] <- 0
  lapply(seq_along(sp$a), function(r) {
    if (is.null(sp$a[[r]])) return(NULL)
    s2 <- sp$stats[[r]]$s2
    sh <- s2 / outer(s2, alpha, `+`)
    if (any(s2 == 0)) sh[s2 == 0, ] <- 0
    sh * sp$a[[r]]
  })
}

# Ridge coefficients (trials x voxels) from shrunk eigen coefficients;
# optionally only the trial rows in `rows`.
.glms_ridge_coef <- function(sp, coef, rows = NULL) {
  out <- vector("list", length(coef))
  off <- 0L
  for (r in seq_along(coef)) {
    if (is.null(coef[[r]])) next
    V <- sp$stats[[r]]$V
    n <- nrow(V)
    if (is.null(rows)) {
      out[[r]] <- V %*% coef[[r]]
    } else {
      keep <- rows[rows > off & rows <= off + n] - off
      if (length(keep)) out[[r]] <- V[keep, , drop = FALSE] %*% coef[[r]]
    }
    off <- off + n
  }
  do.call(rbind, out)
}

# Per-run poly-side SSE of ridge coefficients, computed in eigen
# coordinates: with beta = V w and u = (Ap'E)' beta,
#   SSE = s - 2 (w'vb + u'(E'y)) + w' diag(s2) w + u'u.
.glms_ridge_sse <- function(sp, coef) {
  R <- length(coef)
  V <- length(sp$s[[1]])
  sse <- s <- matrix(0, R, V)
  for (r in seq_len(R)) {
    s[r, ] <- sp$s[[r]]
    w <- coef[[r]]
    if (is.null(w)) { sse[r, ] <- s[r, ]; next }
    st <- sp$stats[[r]]
    q <- 2 * colSums(w * sp$vb[[r]]) - colSums(st$s2 * w^2)
    if (!is.null(st$VApE)) {
      u <- crossprod(st$VApE, w)
      q <- q + 2 * colSums(u * sp$EY[[r]]) - colSums(u^2)
    }
    sse[r, ] <- s[r, ] - q
  }
  list(sse = sse, s = s)
}
