#' Project Trial-wise SBHM Coefficients to Scalar Amplitudes
#'
#' Given trial-wise coefficients in the shared basis (rxntrialsxV) and the
#' voxel-specific matched library coordinates `alpha_hat` (rxV), compute scalar
#' amplitudes per trial and voxel via inner products.
#'
#' @param beta_rt 3D array of shape r x ntrials x V containing per-trial
#'   coefficients in the SBHM basis (as returned by OASIS with K=r, reshaped).
#' @param alpha_hat Numeric matrix rxV of matched library coordinates per voxel
#'   (e.g., `sbhm_match()$alpha_hat`). These should be in the same coordinate
#'   system as `beta_rt` (unwhitened, not L2-normalized) for interpretable amplitudes.
#'
#' @return Numeric matrix ntrials x V of scalar amplitudes.
#'
#' @examples
#' \dontrun{
#'   r <- nrow(alpha_hat)
#'   ntrials <- nrow(beta_mat) / r
#'   beta_rt <- array(beta_mat, dim = c(r, ntrials, ncol(beta_mat)))
#'   amps <- sbhm_project(beta_rt, alpha_hat)
#' }
#'
#' @export
sbhm_project <- function(beta_rt, alpha_hat) {
  stopifnot(length(dim(beta_rt)) == 3L, is.matrix(alpha_hat))
  r <- dim(beta_rt)[1]
  ntrials <- dim(beta_rt)[2]
  V <- dim(beta_rt)[3]
  if (nrow(alpha_hat) != r || ncol(alpha_hat) != V) {
    stop("alpha_hat must be rxV to match beta_rt dims")
  }
  amps <- vapply(seq_len(V), function(v) {
    colSums(beta_rt[, , v] * alpha_hat[, v])
  }, numeric(ntrials))
  amps
}
