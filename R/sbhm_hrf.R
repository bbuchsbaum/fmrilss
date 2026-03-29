#' Wrap a Learned Basis as an HRF (SBHM HRF)
#'
#' Convert a shared time basis matrix `B` (T×r) into an `fmrihrf::HRF` object
#' so it can be used directly in OASIS design construction. The HRF returns the
#' r basis columns evaluated at arbitrary times by piecewise-linear
#' interpolation on the provided time grid.
#'
#' @param B A numeric matrix (T×r) with orthonormal columns (shared basis).
#' @param tgrid Numeric vector of length T giving the global times (seconds)
#'   corresponding to the rows of `B`.
#' @param span Numeric HRF span passed to `fmrihrf::HRF()` metadata.
#'
#' @return An `fmrihrf::HRF` object with `nbasis = ncol(B)`.
#'
#' @examples
#' \dontrun{
#'   hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
#' }
#'
#' @export
sbhm_hrf <- function(B, tgrid, span) {
  stopifnot(is.matrix(B), is.numeric(B), is.numeric(tgrid))
  if (nrow(B) != length(tgrid)) stop("nrow(B) must equal length(tgrid)")
  r <- ncol(B)
  tgrid <- as.numeric(tgrid)

  # Ensure monotone grid for stable interpolation.
  ord <- order(tgrid)
  if (!all(ord == seq_along(tgrid))) {
    tgrid <- tgrid[ord]
    B <- B[ord, , drop = FALSE]
  }
  nT <- length(tgrid)

  # Closure capturing B and tgrid; returns len(t) × r matrix
  fun <- function(t) {
    t <- as.numeric(t)
    n <- length(t)
    out <- matrix(0, n, r)

    # HRF support: explicitly zero outside [0, span] and outside interpolation grid.
    inside <- (t >= 0) & (t <= span) & (t >= tgrid[1L]) & (t <= tgrid[nT])
    if (!any(inside)) return(out)

    ti <- t[inside]
    idx <- findInterval(ti, tgrid)
    idx <- pmax(1L, pmin(idx, nT - 1L))
    idx2 <- idx + 1L
    dt <- tgrid[idx2] - tgrid[idx]
    w <- (ti - tgrid[idx]) / ifelse(dt == 0, 1, dt)

    out[inside, ] <- B[idx, , drop = FALSE] * (1 - w) + B[idx2, , drop = FALSE] * w
    out
  }
  fmrihrf::HRF(fun = fun, name = "SBHM", nbasis = r, span = span)
}
