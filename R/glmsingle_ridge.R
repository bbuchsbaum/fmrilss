# Run-blocked least squares and fractional ridge for glmsingle().
#
# After nuisance projection GLMsingle's single-trial design is block diagonal
# across runs, so every fit factorises per run. Fractional ridge keeps one
# penalty per voxel, calibrated on the pooled spectrum of all runs, exactly as
# fracridge() does on the stacked design.

# Voxel-independent statistics of one run's trial design for one HRF and one
# nuisance model: poly-residualised design Ap, fully residualised Gram G,
# poly-only Gram Gp, and Ap' E for the low-rank correction of the data.
.glms_design_stats <- function(onsets, n_time, hrf, nuis_run, k_name) {
  .glms_design_stats_x(.glms_trial_design(onsets, n_time, hrf), nuis_run, k_name)
}

.glms_design_stats_x <- function(X, nuis_run, k_name) {
  Ap <- .glms_resid(X, nuis_run$Qp)
  E <- nuis_run$E[[k_name]]
  A <- .glms_resid(Ap, E)
  list(
    Ap = Ap,
    ApE = if (ncol(E)) crossprod(Ap, E) else NULL,
    G = crossprod(A),
    Gp = crossprod(Ap),
    zero = colSums(X != 0) == 0,
    k = k_name
  )
}

# Data statistics of one run for a voxel tile: run-centred data, poly-residual
# sum of squares s and the nuisance products E'Y for each nuisance model.
.glms_data_stats <- function(Yr, nuis_run, k_names = names(nuis_run$E)) {
  Yc <- sweep(Yr, 2L, colMeans(Yr), "-")
  QpY <- crossprod(nuis_run$Qp, Yc)
  EY <- lapply(nuis_run$E[k_names], function(E) if (ncol(E)) crossprod(E, Yc) else NULL)
  list(Yc = Yc, s = colSums(Yc^2) - colSums(QpY^2), EY = EY)
}

# Fit-side right-hand side b = A'y and report-side c = Ap'y.
.glms_rhs <- function(ds, dat) {
  c <- crossprod(ds$Ap, dat$Yc)
  EY <- dat$EY[[ds$k]]
  b <- if (is.null(ds$ApE)) c else c - ds$ApE %*% EY
  list(b = b, c = c)
}

# Poly-only residual sum of squares of a fit (reported R^2 convention).
.glms_sse <- function(beta, c, Gp, s) {
  s - 2 * colSums(beta * c) + colSums(beta * (Gp %*% beta))
}

.glms_r2 <- function(sse, s) {
  out <- 100 * (1 - sse / s)
  out[s == 0] <- NaN
  out
}

# OLS for all runs of a voxel tile. `stats` is a list (one per run) of design
# statistics, `dats` the matching data statistics. Returns betas (trials x
# voxels), per-run SSE and s, and the poly-side products for R^2.
.glms_fit_ols <- function(stats, dats, singular) {
  R <- length(stats)
  betas <- vector("list", R)
  sse <- s <- matrix(0, R, ncol(dats[[1]]$Yc))
  for (r in seq_len(R)) {
    if (!length(stats[[r]]$zero)) {
      sse[r, ] <- s[r, ] <- dats[[r]]$s
      betas[[r]] <- matrix(0, 0L, ncol(sse))
      next
    }
    rhs <- .glms_rhs(stats[[r]], dats[[r]])
    b <- .glms_solve(stats[[r]]$G, rhs$b, stats[[r]]$zero, singular)
    betas[[r]] <- b
    sse[r, ] <- .glms_sse(b, rhs$c, stats[[r]]$Gp, dats[[r]]$s)
    s[r, ] <- dats[[r]]$s
  }
  list(beta = do.call(rbind, betas), sse = sse, s = s)
}

# Spectral form of the stacked problem used by fracridge: per-run eigen
# decompositions of G, the pooled singular values, and OLS coefficients in
# eigen coordinates (components with singular value < 1e-10 are zeroed).
.glms_spectral <- function(stats, dats) {
  R <- length(stats)
  vecs <- s2 <- a <- c_list <- vector("list", R)
  for (r in seq_len(R)) {
    if (!length(stats[[r]]$zero)) next
    e <- eigen(stats[[r]]$G, symmetric = TRUE)
    vecs[[r]] <- e$vectors
    s2[[r]] <- pmax(e$values, 0)
    rhs <- .glms_rhs(stats[[r]], dats[[r]])
    c_list[[r]] <- rhs$c
    coef <- crossprod(e$vectors, rhs$b) / s2[[r]]
    coef[sqrt(s2[[r]]) < 1e-10, ] <- 0
    a[[r]] <- coef
  }
  list(vecs = vecs, s2 = s2, a = a, c = c_list,
       s2_all = unlist(s2), a_all = do.call(rbind, a))
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

# Coefficients (trials x voxels) for one alpha per voxel; optionally only the
# trial rows in `rows`.
.glms_ridge_coef <- function(sp, alpha, rows = NULL) {
  R <- length(sp$vecs)
  out <- vector("list", R)
  off <- 0L
  for (r in seq_len(R)) {
    if (is.null(sp$vecs[[r]])) next
    s2 <- sp$s2[[r]]
    shrink <- outer(s2, alpha, function(s, al) s / (s + al))
    coef <- shrink * sp$a[[r]]
    coef[!is.finite(coef)] <- 0
    V <- sp$vecs[[r]]
    n <- nrow(V)
    if (!is.null(rows)) {
      keep <- rows[rows > off & rows <= off + n] - off
      out[[r]] <- V[keep, , drop = FALSE] %*% coef
    } else {
      out[[r]] <- V %*% coef
    }
    off <- off + n
  }
  do.call(rbind, out)
}

# Per-run SSE of ridge coefficients under the poly-only reporting convention.
.glms_ridge_sse <- function(sp, beta, stats, dats) {
  R <- length(stats)
  sse <- s <- matrix(0, R, ncol(beta))
  off <- 0L
  for (r in seq_len(R)) {
    s[r, ] <- dats[[r]]$s
    if (is.null(sp$vecs[[r]])) { sse[r, ] <- s[r, ]; next }
    n <- nrow(sp$vecs[[r]])
    b <- beta[off + seq_len(n), , drop = FALSE]
    sse[r, ] <- .glms_sse(b, sp$c[[r]], stats[[r]]$Gp, dats[[r]]$s)
    off <- off + n
  }
  list(sse = sse, s = s)
}
