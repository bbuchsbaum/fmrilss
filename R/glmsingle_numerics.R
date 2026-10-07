# Fixed numerical rules shared by the glmsingle() stages.
#
# These reproduce the arithmetic conventions of GLMsingle (cvnlab/GLMsingle,
# Python port at commit 1ab54a6) where they affect estimates. They are not
# user options; variants that users may choose are arguments of glmsingle().

# Per-call settings shared by the C++ kernels (set and restored by glmsingle()).
.glms_state <- new.env(parent = emptyenv())
.glms_state$n_threads <- 1L
.glms_nt <- function() .glms_state$n_threads

# Columns `vox` of Y, without a copy when they are all columns in order.
.glms_cols <- function(Y, vox) {
  if (length(vox) == ncol(Y) && identical(vox, seq_len(ncol(Y)))) Y else Y[, vox, drop = FALSE]
}

# GLMsingle's alt_round(): round half away from zero, to integer.
.glms_alt_round <- function(x) {
  as.integer(sign(x) * ceiling(floor(abs(x) * 2) / 2))
}

# Orthonormal basis for the span of `M` under GLMsingle's projection rule.
#
# make_projection_matrix() calls olsmatrix(): exactly-zero columns are dropped,
# the remaining columns are unit-normalised, and np.linalg.pinv() of the
# normalised Gram is applied (default rcond = 1e-15 on Gram eigenvalues). That
# is equivalent to keeping singular values of the normalised design above
# sqrt(1e-15) times the largest one.
.glms_span <- function(M) {
  M <- as.matrix(M)
  if (!ncol(M)) return(matrix(0, nrow(M), 0L))
  keep <- colSums(M != 0) > 0
  M <- M[, keep, drop = FALSE]
  if (!ncol(M)) return(matrix(0, nrow(M), 0L))
  M <- sweep(M, 2L, sqrt(colSums(M^2)), "/")
  s <- svd(M, nv = 0L)
  rank <- sum(s$d > sqrt(1e-15) * s$d[1L])
  s$u[, seq_len(rank), drop = FALSE]
}

# Remove the component of `X` lying in the span of orthonormal `Q`.
.glms_resid <- function(X, Q) {
  if (is.null(Q) || !ncol(Q)) return(X)
  X - Q %*% crossprod(Q, X)
}

# Autoscale operator: for each voxel, regress the reference betas `ref`
# (trials x voxels) on [cand, 1] with GLMsingle's olsmatrix() operator
# (column-normalised Gram pseudoinverse) and return list(h, fitted).
# A negative scale is replaced by (1, 0), as in GLMsingle.
.glms_autoscale <- function(cand, ref) {
  cand[!is.finite(cand)] <- 0
  N <- nrow(cand)
  len1 <- sqrt(colSums(cand^2))
  len2 <- sqrt(N)
  ok1 <- len1 > 0
  x1r <- ifelse(ok1, colSums(cand * ref) / ifelse(ok1, len1, 1), 0)
  x2r <- colSums(ref) / len2
  c12 <- ifelse(ok1, colSums(cand) / (ifelse(ok1, len1, 1) * len2), 0)
  # pinv of [[1, c], [c, 1]] via its eigenpairs (1 + c, (1,1)/sqrt2), (1 - c, (1,-1)/sqrt2)
  lp <- 1 + c12
  lm <- 1 - c12
  cut <- 1e-15 * pmax(lp, lm)
  ip <- ifelse(lp > cut, 1 / lp, 0)
  im <- ifelse(lm > cut, 1 / lm, 0)
  a11 <- (ip + im) / 2
  a12 <- (ip - im) / 2
  h1 <- ifelse(ok1, (a11 * x1r + a12 * x2r) / ifelse(ok1, len1, 1), 0)
  h2 <- ifelse(ok1, (a12 * x1r + a11 * x2r), x2r) / len2
  neg <- h1 < 0
  h1[neg] <- 1
  h2[neg] <- 0
  fitted <- sweep(cand, 2L, h1, "*") + rep(h2, each = N)
  list(h = cbind(scale = h1, offset = h2), fitted = fitted)
}

# numpy.percentile (linear interpolation) == R quantile type 7.
.glms_percentile <- function(x, p) {
  unname(stats::quantile(x, probs = p / 100, type = 7, names = FALSE))
}
