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
#'   - `tau` numeric in `[0,1]`: global strength (default 0, i.e., no shrinkage).
#'   - `ref` numeric length-r vector (alpha_ref) or NULL. If NULL, uses
#'     the mean of A columns. Shrinkage is: beta_bar <- (1-lambda) beta_bar + lambda ref.
#'   - `snr` optional numeric length-V estimates; if provided, per-voxel
#'     lambda_v = tau/(snr_v + tau). Otherwise lambda = tau.
#' @param topK Positive integer no larger than the library size; return top-K
#'   scores/weights if greater than one (default 1).
#' @param whiten Logical, divide coefficients by S before normalization (default TRUE).
#' @param sv_floor_rel Relative singular-value floor used when `whiten=TRUE`
#'   (default `1e-6`).
#' @param whiten_power Numeric in `[0, 1]` controlling whitening strength when
#'   `whiten=TRUE`. Uses division by `S^whiten_power` (`1` = full whitening,
#'   `0.5` = partial whitening, `0` = no whitening). Default `1`.
#' @param orient_ref Logical, flip beta_bar columns when their dot with `ref` is
#'   negative before matching (default TRUE).
#'
#' @return A list with:
#'   - `idx` length-V integer indices of best-matching library HRF (1..K)
#'   - `margin` named length-V numeric score difference: top1 - top2. It is not
#'     a calibrated confidence measure.
#'   - `alpha_hat` r×V matrix: the selected library coordinates (unwhitened, unnormalized)
#'   - `scores` optional K×V cosine score matrix (returned when topK > 1)
#'   - `weights` optional cosine-softmax top-K weights per voxel (when topK > 1)
#'
#' Non-finite inputs, zero-norm voxel summaries or library candidates, and
#' malformed matching controls fail explicitly. Named basis axes must identify
#' the same complete set and are aligned before matching.
#'
#' @examples
#' \donttest{
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
                       sv_floor_rel = 1e-6,
                       whiten_power = 1,
                       orient_ref = TRUE) {
  if (!is.matrix(beta_bar) || !is.numeric(beta_bar) || !is.numeric(S) ||
      !is.matrix(A) || !is.numeric(A)) {
    stop("beta_bar and A must be numeric matrices and S a numeric vector", call. = FALSE)
  }
  if (!nrow(beta_bar) || !ncol(beta_bar) || !nrow(A) || !ncol(A) ||
      any(!is.finite(beta_bar)) || any(!is.finite(S)) || any(!is.finite(A))) {
    stop("SBHM matching inputs must be non-empty and finite", call. = FALSE)
  }
  r <- nrow(A)
  if (nrow(beta_bar) != r) stop("nrow(beta_bar) must equal nrow(A)")
  if (length(S) != r) stop("length(S) must equal nrow(A)")
  if (any(S <= 0)) stop("S must contain positive singular values", call. = FALSE)

  beta_basis <- rownames(beta_bar)
  library_basis <- rownames(A)
  if (xor(is.null(beta_basis), is.null(library_basis))) {
    stop("basis names must be supplied on both beta_bar and A or neither",
         call. = FALSE)
  }
  if (!is.null(beta_basis)) {
    if (anyNA(beta_basis) || anyNA(library_basis) || any(!nzchar(beta_basis)) ||
        any(!nzchar(library_basis)) || anyDuplicated(beta_basis) ||
        anyDuplicated(library_basis) || !setequal(beta_basis, library_basis)) {
      stop("beta_bar and A basis names must be complete, unique, and identical",
           call. = FALSE)
    }
    beta_bar <- beta_bar[match(library_basis, beta_basis), , drop = FALSE]
  }
  if (!is.null(names(S))) {
    if (anyNA(names(S)) || any(!nzchar(names(S))) || anyDuplicated(names(S)) ||
        is.null(library_basis) || !setequal(names(S), library_basis)) {
      stop("named S must identify the same complete basis as A", call. = FALSE)
    }
    S <- S[match(library_basis, names(S))]
  }

  V <- ncol(beta_bar)
  K <- ncol(A)
  if (!is.null(colnames(beta_bar)) &&
      (anyNA(colnames(beta_bar)) || any(!nzchar(colnames(beta_bar))) ||
       anyDuplicated(colnames(beta_bar)))) {
    stop("beta_bar voxel names must be complete and unique", call. = FALSE)
  }
  if (!is.null(colnames(A)) &&
      (anyNA(colnames(A)) || any(!nzchar(colnames(A))) || anyDuplicated(colnames(A)))) {
    stop("A candidate names must be complete and unique", call. = FALSE)
  }
  topK_num <- suppressWarnings(as.numeric(topK))
  if (length(topK_num) != 1L || !is.finite(topK_num) || topK_num < 1 ||
      topK_num != round(topK_num) || topK_num > K) {
    stop("topK must be one integer between 1 and ncol(A)", call. = FALSE)
  }
  topK <- as.integer(topK_num)
  if (!is.logical(whiten) || length(whiten) != 1L || is.na(whiten) ||
      !is.logical(orient_ref) || length(orient_ref) != 1L || is.na(orient_ref)) {
    stop("whiten and orient_ref must be TRUE or FALSE", call. = FALSE)
  }
  p <- suppressWarnings(as.numeric(whiten_power))
  if (length(p) != 1L || !is.finite(p) || p < 0 || p > 1) {
    stop("whiten_power must be one finite value in [0, 1]", call. = FALSE)
  }
  rel_floor_input <- suppressWarnings(as.numeric(sv_floor_rel))
  if (length(rel_floor_input) != 1L || !is.finite(rel_floor_input) || rel_floor_input < 0) {
    stop("sv_floor_rel must be one nonnegative finite value", call. = FALSE)
  }
  if (!is.list(shrink)) stop("shrink must be a list", call. = FALSE)

  # Shrinkage toward reference (optional)
  tau <- suppressWarnings(as.numeric(shrink$tau %||% 0))
  if (length(tau) != 1L || !is.finite(tau) || tau < 0 || tau > 1) {
    stop("shrink$tau must be one finite value in [0, 1]", call. = FALSE)
  }
  if (!is.null(shrink$ref)) {
    alpha_ref <- shrink$ref
    if (!is.null(names(alpha_ref)) && !is.null(library_basis)) {
      if (anyNA(names(alpha_ref)) || any(!nzchar(names(alpha_ref))) ||
          anyDuplicated(names(alpha_ref)) || !setequal(names(alpha_ref), library_basis)) {
        stop("named shrink$ref must identify the same complete basis as A",
             call. = FALSE)
      }
      alpha_ref <- alpha_ref[match(library_basis, names(alpha_ref))]
    }
    alpha_ref <- as.numeric(alpha_ref)
  } else {
    alpha_ref <- rowMeans(A)
  }
  if (length(alpha_ref) != r || any(!is.finite(alpha_ref))) {
    stop("shrink$ref must contain one finite value per basis", call. = FALSE)
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
      if (length(snr) != V || any(!is.finite(snr)) || any(snr < 0)) {
        stop("shrink$snr must contain one nonnegative finite value per voxel", call. = FALSE)
      }
      lam <- tau / (snr + tau)
    } else {
      lam <- rep.int(tau, V)
    }
    beta_bar <- (1 - rep(lam, each = r)) * beta_bar + alpha_ref %*% t(lam)
  }

  # Shape-only whitening and L2 normalization
  if (isTRUE(whiten)) {
    smax <- suppressWarnings(max(S, na.rm = TRUE))
    rel_floor <- if (is.finite(smax) && smax > 0) {
      rel_floor_input * smax
    } else {
      .Machine$double.eps
    }
    S_safe <- pmax(S, rel_floor, .Machine$double.eps)
    S_scale <- S_safe^p
    beta_w <- beta_bar / S_scale
    A_w    <- A / S_scale
  } else {
    beta_w <- beta_bar
    A_w    <- A
  }
  beta_norm <- sqrt(colSums(beta_w^2))
  library_norm <- sqrt(colSums(A_w^2))
  if (any(beta_norm <= .Machine$double.eps * max(beta_norm, 1))) {
    stop("beta_bar contains a zero-norm voxel shape and cannot be matched", call. = FALSE)
  }
  if (any(library_norm <= .Machine$double.eps * max(library_norm, 1))) {
    stop("A contains a zero-norm library candidate", call. = FALSE)
  }
  beta_w <- sweep(beta_w, 2L, beta_norm, "/")
  A_w    <- sweep(A_w,    2L, library_norm, "/")

  # Cosine scores and assignment
  Scores <- crossprod(A_w, beta_w)  # K×V
  k_star <- max.col(t(Scores), ties.method = "first")
  # margin = top1 - top2 per voxel
  margin <- apply(Scores, 2L, function(s) {
    o <- sort(s, decreasing = TRUE)
    if (length(o) >= 2) o[1] - o[2] else o[1]
  })

  alpha_hat <- A[, k_star, drop = FALSE]
  voxel_names <- colnames(beta_bar) %||% paste0("voxel_", seq_len(V))
  candidate_names <- colnames(A) %||% paste0("candidate_", seq_len(K))
  names(k_star) <- names(margin) <- voxel_names
  dimnames(alpha_hat) <- list(rownames(A), voxel_names)

  out <- list(idx = k_star, margin = as.numeric(margin), alpha_hat = alpha_hat)
  names(out$margin) <- voxel_names

  if (topK > 1) {
    Kret <- topK
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
    dimnames(out$scores) <- list(candidate_names, voxel_names)
    dimnames(out$topK_idx) <- list(paste0("rank_", seq_len(Kret)), voxel_names)
    dimnames(out$weights) <- dimnames(out$topK_idx)
  }
  out
}
