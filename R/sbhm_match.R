#' Match Voxels to Library HRFs in Shared Basis (SBHM)
#'
#' Given per-voxel aggregate coefficients `beta_bar` in the shared basis `B`,
#' and library coordinates `A`, perform shape-only matching in a whitened,
#' L2-normalized coefficient space (cosine similarity). Optionally apply a
#' simple shrinkage of `beta_bar` towards a reference coordinate before
#' matching, and an orientation fix.
#'
#' @param beta_bar Numeric matrix (r×V): per-voxel coefficients from a prepass
#'   GLM in the SBHM basis. Columns correspond to voxels.
#' @param S Numeric vector (length r): singular values of the library SVD.
#' @param A Numeric matrix (r×K): coordinates of library HRFs in SBHM basis.
#' @param shrink List with optional shrinkage options:
#'   - `tau` numeric >=0: global strength (default 0, i.e., no shrinkage).
#'   - `ref` numeric length-r vector (alpha_ref) or NULL. If NULL, uses
#'     the mean of A columns. Shrinkage is: beta_bar <- (1-lambda) beta_bar + lambda ref.
#'   - `snr` optional numeric length-V estimates; if provided, per-voxel
#'     lambda_v = tau/(snr_v + tau). Otherwise lambda = tau.
#' @param topK Integer, return top-K scores/weights if >1 (default 1).
#' @param whiten Logical, divide coefficients by S before normalization (default TRUE).
#' @param orient_ref Logical, flip beta_bar columns when their dot with `ref` is
#'   negative before matching (default TRUE).
#'
#' @return A list with:
#'   - `idx` length-V integer indices of best-matching library HRF (1..K)
#'   - `margin` length-V numeric: score(top1) - score(top2)
#'   - `alpha_hat` r×V matrix: the selected library coordinates (unwhitened, unnormalized)
#'   - `scores` optional K×V cosine score matrix (returned when topK > 1)
#'   - `weights` optional top-K weights per voxel (when topK > 1)
#'
#' @examples
#' \dontrun{
#'   set.seed(42)
#'   r <- 4; K <- 12; V <- 3
#'   A <- matrix(rnorm(r*K), r, K)
#'   S <- seq(1, 0.2, length.out = r)
#'   alpha2 <- A[,2]
#'   beta_bar <- cbind(alpha2 + rnorm(r, sd = 0.1),
#'                     A[,7] + rnorm(r, sd = 0.1),
#'                     A[,10] + rnorm(r, sd = 0.1))
#'   m <- sbhm_match(beta_bar, S, A)
#'   m$idx
#' }
#'
#' @export
sbhm_match <- function(beta_bar, S, A,
                       shrink = list(tau = 0, ref = NULL, snr = NULL),
                       topK = 1,
                       whiten = TRUE,
                       orient_ref = TRUE) {
  stopifnot(is.matrix(beta_bar), is.numeric(S), is.matrix(A))
  r <- nrow(A)
  if (nrow(beta_bar) != r) stop("nrow(beta_bar) must equal nrow(A)")
  if (length(S) != r) stop("length(S) must equal nrow(A)")

  V <- ncol(beta_bar)
  K <- ncol(A)

  # Shrinkage toward reference (optional)
  tau <- shrink$tau %||% 0
  if (!is.null(shrink$ref)) {
    alpha_ref <- as.numeric(shrink$ref)
  } else {
    alpha_ref <- rowMeans(A)
  }
  if (isTRUE(orient_ref)) {
    # Flip beta_bar columns to align with reference
    dots <- as.numeric(crossprod(alpha_ref, beta_bar))
    flip <- which(dots < 0)
    if (length(flip)) beta_bar[, flip] <- -beta_bar[, flip, drop = FALSE]
  }
  if (tau > 0) {
    if (!is.null(shrink$snr)) {
      snr <- as.numeric(shrink$snr)
      if (length(snr) != V) stop("shrink$snr must have length ncol(beta_bar)")
      lam <- tau / (snr + tau)
    } else {
      lam <- rep.int(tau, V)
    }
    beta_bar <- (1 - rep(lam, each = r)) * beta_bar + alpha_ref %*% t(lam)
  }

  # Shape-only whitening and L2 normalization
  if (isTRUE(whiten)) {
    S_safe <- pmax(S, .Machine$double.eps)
    beta_w <- beta_bar / S_safe
    A_w    <- A / S_safe
  } else {
    beta_w <- beta_bar
    A_w    <- A
  }
  beta_w <- sweep(beta_w, 2L, sqrt(colSums(beta_w^2)) + 1e-8, "/")
  A_w    <- sweep(A_w,    2L, sqrt(colSums(A_w^2))    + 1e-8, "/")

  # Cosine scores and assignment
  Scores <- crossprod(A_w, beta_w)  # K×V
  k_star <- max.col(t(Scores), ties.method = "first")
  # margin = top1 - top2 per voxel
  margin <- apply(Scores, 2L, function(s) {
    o <- sort(s, decreasing = TRUE)
    if (length(o) >= 2) o[1] - o[2] else o[1]
  })

  alpha_hat <- A[, k_star, drop = FALSE]

  out <- list(idx = k_star, margin = as.numeric(margin), alpha_hat = alpha_hat)

  if (topK > 1) {
    Kret <- min(topK, K)
    top_idx <- matrix(NA_integer_, nrow = Kret, ncol = V)
    weights <- matrix(NA_real_,  nrow = Kret, ncol = V)
    for (v in seq_len(V)) {
      s <- Scores[, v]
      ord <- order(s, decreasing = TRUE)[seq_len(Kret)]
      top_idx[, v] <- ord
      z <- s[ord]; z <- z - max(z)
      w <- exp(z); w <- w / sum(w)
      weights[, v] <- w
    }
    out$scores   <- Scores
    out$topK_idx <- top_idx
    out$weights  <- weights
  }
  out
}

# local utility for defaulting
`%||%` <- function(a, b) if (is.null(a)) b else a
