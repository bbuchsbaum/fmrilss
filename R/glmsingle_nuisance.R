# Nuisance bases for glmsingle(): polynomial drift, user extra regressors and
# GLMdenoise noise PCs, kept per run.

# Orthonormal polynomial basis of degrees 0..deg on linspace(-1, 1, n), as
# make_polynomial_matrix() (only the span matters for projection).
.glms_poly_basis <- function(n, deg) {
  t <- seq(-1, 1, length.out = n)
  P <- outer(t, 0:deg, `^`)
  .glms_span(P)
}

# Per-run nuisance bases.
#
# For a run, `Qp` spans the polynomials. For each requested model `k`
# (number of noise PCs), `E[[k + 1]]` is an orthonormal basis of the part of
# span([poly, extras?, pcs[, 1:k]]) orthogonal to the polynomials, where the
# span follows GLMsingle's projection rule. Extras are included when
# `extras_for(k)` is TRUE.
.glms_nuisance_run <- function(n_time, deg, extras, pcs, ks, extras_for) {
  Qp <- .glms_poly_basis(n_time, deg)
  E <- lapply(ks, function(k) {
    cols <- list(Qp)
    if (!is.null(extras) && ncol(extras) && extras_for(k)) cols <- c(cols, list(extras))
    if (k > 0L) cols <- c(cols, list(pcs[, seq_len(k), drop = FALSE]))
    if (length(cols) == 1L) return(matrix(0, n_time, 0L))
    U <- .glms_span(do.call(cbind, cols))
    Ue <- .glms_resid(U, Qp)
    if (!ncol(Ue)) return(Ue)
    s <- svd(Ue, nv = 0L)
    s$u[, s$d > 0.5, drop = FALSE]
  })
  names(E) <- paste0("k", ks)
  list(Qp = Qp, E = E, ks = ks)
}

# GLMdenoise noise regressors for one run: poly-projected noise-pool time
# series, unit-normalised per voxel, covariance accumulated over voxel tiles,
# leading eigenvectors scaled to unit standard deviation (ddof = 1).
.glms_noise_pcs <- function(Yr, pool, Qp, max_pcs, chunk_size) {
  n_time <- nrow(Yr)
  idx <- which(pool)
  if (!length(idx)) return(matrix(0, n_time, 0L))
  C <- matrix(0, n_time, n_time)
  for (blk in split(idx, ceiling(seq_along(idx) / chunk_size))) {
    Z <- Yr[, blk, drop = FALSE]
    Z <- sweep(Z, 2L, colMeans(Z))
    scale <- sqrt(colSums(Z^2))
    Z <- .glms_resid(Z, Qp)
    nrm <- sqrt(colSums(Z^2))
    keep <- nrm > sqrt(.Machine$double.eps) * scale
    if (!any(keep)) next
    Z <- sweep(Z[, keep, drop = FALSE], 2L, nrm[keep], "/")
    C <- C + tcrossprod(Z)
  }
  e <- eigen(C, symmetric = TRUE)
  rank <- sum(e$values > n_time * .Machine$double.eps * max(e$values))
  if (!rank) return(matrix(0, n_time, 0L))
  U <- e$vectors[, seq_len(min(max_pcs + 1L, rank)), drop = FALSE]
  sweep(U, 2L, apply(U, 2L, stats::sd), "/")
}
