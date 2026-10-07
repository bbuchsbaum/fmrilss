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
# nuisance model. `solver` precomputes the OLS operator; `spectral` the
# eigendecomposition used by fractional ridge.
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
    Ginv <- .glms_solve(G, diag(out$n), out$zero, singular)
    out$M <- Ginv %*% t(Ap)
    out$K <- if (is.null(ApE)) NULL else Ginv %*% ApE
  }
  if (spectral) {
    e <- eigen(G, symmetric = TRUE)
    out$V <- e$vectors
    out$s2 <- pmax(e$values, 0)
    out$Wt <- crossprod(e$vectors, t(Ap))
    out$VApE <- if (is.null(ApE)) NULL else crossprod(e$vectors, ApE)
  }
  out
}

# Data statistics of one run for a voxel tile: run-centred data (exact under
# the polynomial projection, and it limits cancellation), the poly-residual
# sum of squares s, and E'y for each nuisance model.
.glms_data_stats <- function(Yr, nuis_run, k_names = names(nuis_run$E)) {
  Yc <- Yr - rep(colMeans(Yr), each = nrow(Yr))
  QpY <- crossprod(nuis_run$Qp, Yc)
  EY <- lapply(nuis_run$E[k_names], function(E) if (ncol(E)) crossprod(E, Yc) else NULL)
  list(Yc = Yc, s = colSums(Yc^2) - colSums(QpY^2), EY = EY)
}

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
  c <- crossprod(ds$Ap, dat$Yc)
  EY <- dat$EY[[ds$k]]
  b <- if (is.null(ds$ApE)) c else c - ds$ApE %*% EY
  list(b = b, c = c)
}

# OLS for all runs of a voxel tile: beta = M y - K (E'y), and the poly-side
# SSE = s - beta'G beta - 2 u'(E'y) + u'u with u = (Ap'E)' beta.
.glms_fit_ols <- function(stats, dats) {
  R <- length(stats)
  V <- ncol(dats[[1]]$Yc)
  betas <- vector("list", R)
  sse <- s <- matrix(0, R, V)
  for (r in seq_len(R)) {
    st <- stats[[r]]; dt <- dats[[r]]
    s[r, ] <- dt$s
    if (!st$n) { sse[r, ] <- dt$s; next }
    EY <- dt$EY[[st$k]]
    b <- st$M %*% dt$Yc
    if (!is.null(st$K)) b <- b - st$K %*% EY
    q <- colSums(b * (st$G %*% b))
    if (!is.null(st$ApE)) {
      u <- crossprod(st$ApE, b)
      q <- q + 2 * colSums(u * EY) - colSums(u^2)
    }
    sse[r, ] <- dt$s - q
    betas[[r]] <- b
  }
  list(beta = do.call(rbind, betas), sse = sse, s = s)
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
    x <- st$Wt %*% dats[[r]]$Yc
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
    return(glms_frac_alpha_exact(a2, sp$s2_all, fracs))
  }
  grid <- .glms_alpha_grid(sp$s2_all)
  s2 <- sqrt(sp$s2_all)^2
  sclg_sq <- (outer(grid, s2, function(g, s) s / (s + g)))^2
  newlen <- sqrt(sclg_sq %*% a2)
  glms_frac_alpha_grid(newlen, grid, fracs)
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
