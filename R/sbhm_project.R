#' Project Trial-wise SBHM Coefficients to Scalar Amplitudes
#'
#' Given trial-wise coefficients in the shared basis (rxntrialsxV) and the
#' voxel-specific matched library coordinates `alpha_hat` (rxV), compute scalar
#' amplitudes per trial and voxel via least-squares projection:
#' `a = (alpha' beta) / (alpha' alpha)`.
#'
#' @param beta_rt 3D array of shape r x ntrials x V containing per-trial
#'   coefficients in the SBHM basis (as returned by OASIS with K=r, reshaped).
#' @param alpha_hat Numeric matrix rxV of matched library coordinates per voxel
#'   (e.g., `sbhm_match()$alpha_hat`). These should be in the same coordinate
#'   system as `beta_rt` (unwhitened, not L2-normalized) for interpretable amplitudes.
#'
#' @return Numeric matrix ntrials x V of scalar amplitudes.
#'   Zero-norm or non-finite coordinates are unidentified and fail explicitly.
#'   When basis or voxel names are present, both inputs must provide the same
#'   complete unique identities; `alpha_hat` is reordered to `beta_rt`.
#'
#' @examples
#' \donttest{
#'   set.seed(1)
#'   r <- 2; ntrials <- 3; nvox <- 2
#'   alpha_hat <- matrix(rnorm(r * nvox), r, nvox)
#'   beta_rt <- array(rnorm(r * ntrials * nvox), c(r, ntrials, nvox))
#'   amps <- sbhm_project(beta_rt, alpha_hat)
#'   dim(amps)
#' }
#'
#' @export
sbhm_project <- function(beta_rt, alpha_hat) {
  if (length(dim(beta_rt)) != 3L || !is.numeric(beta_rt) ||
      !is.matrix(alpha_hat) || !is.numeric(alpha_hat) ||
      any(!is.finite(beta_rt)) || any(!is.finite(alpha_hat))) {
    stop("beta_rt and alpha_hat must be finite numeric arrays with the documented shapes",
         call. = FALSE)
  }
  r <- dim(beta_rt)[1]
  ntrials <- dim(beta_rt)[2]
  V <- dim(beta_rt)[3]
  if (nrow(alpha_hat) != r || ncol(alpha_hat) != V) {
    stop("alpha_hat must be rxV to match beta_rt dims")
  }
  align_axis <- function(reference, supplied, label) {
    if (is.null(reference) && is.null(supplied)) return(NULL)
    if (is.null(reference) || is.null(supplied)) {
      stop(label, " names must be supplied on both inputs or neither", call. = FALSE)
    }
    if (anyNA(reference) || anyNA(supplied) || any(!nzchar(reference)) ||
        any(!nzchar(supplied)) || anyDuplicated(reference) || anyDuplicated(supplied) ||
        !setequal(reference, supplied)) {
      stop(label, " names must be complete, unique, and identical", call. = FALSE)
    }
    match(reference, supplied)
  }
  basis_order <- align_axis(dimnames(beta_rt)[[1L]], rownames(alpha_hat), "basis")
  if (!is.null(basis_order)) alpha_hat <- alpha_hat[basis_order, , drop = FALSE]
  voxel_order <- align_axis(dimnames(beta_rt)[[3L]], colnames(alpha_hat), "voxel")
  if (!is.null(voxel_order)) alpha_hat <- alpha_hat[, voxel_order, drop = FALSE]
  denom <- colSums(alpha_hat^2)
  if (any(denom <= .Machine$double.eps * max(denom, 1))) {
    stop("alpha_hat contains a zero-norm voxel coordinate", call. = FALSE)
  }

  amps <- matrix(
    NA_real_, nrow = ntrials, ncol = V,
    dimnames = list(dimnames(beta_rt)[[2L]], dimnames(beta_rt)[[3L]])
  )
  for (v in seq_len(V)) {
    Bv <- matrix(beta_rt[, , v, drop = FALSE], nrow = r, ncol = ntrials)
    num <- as.numeric(crossprod(alpha_hat[, v], Bv))
    amps[, v] <- num / denom[v]
  }
  amps
}
