# Fixed numerical rules shared by the glmsingle() stages.
#
# These reproduce the arithmetic conventions of GLMsingle (cvnlab/GLMsingle,
# Python port at commit 1ab54a6) where they affect estimates. They are not
# user options; variants that users may choose are arguments of glmsingle().

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

# Solve G beta = B for a symmetric positive (semi)definite G.
#
# Mirrors olsmatrix2(): columns whose design column is exactly zero get a zero
# coefficient; the remaining system is solved directly. A singular remaining
# system is an error (as in GLMsingle) unless singular = "pinv".
.glms_solve <- function(G, B, zero_cols = NULL, singular = "error",
                        context = "trial design") {
  B <- as.matrix(B)
  n <- nrow(G)
  out <- matrix(0, n, ncol(B))
  good <- if (is.null(zero_cols)) rep(TRUE, n) else !zero_cols
  if (!any(good)) return(out)
  Gg <- G[good, good, drop = FALSE]
  R <- tryCatch(chol(Gg), error = function(e) NULL)
  if (!is.null(R)) {
    d <- diag(R)
    if (min(d) > sqrt(.Machine$double.eps) * max(d)) {
      out[good, ] <- backsolve(R, forwardsolve(t(R), B[good, , drop = FALSE]))
      return(out)
    }
  }
  if (identical(singular, "pinv")) {
    # counted by glmsingle(), which warns once per call
    signalCondition(structure(class = c("glms_singular", "condition"),
                              list(message = context, call = NULL)))
    e <- eigen(Gg, symmetric = TRUE)
    keep <- e$values > max(e$values) * n * .Machine$double.eps
    V <- e$vectors[, keep, drop = FALSE]
    out[good, ] <- V %*% (crossprod(V, B[good, , drop = FALSE]) / e$values[keep])
    return(out)
  }
  stop(sprintf(
    "Singular %s: two or more trial regressors are linearly dependent. ",
    context),
    "Use singular = \"pinv\" for a minimum-norm solution.", call. = FALSE)
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
