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
#'   # Given sbhm <- sbhm_build(...)
#'   hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
#'   # Use in lss(..., method="oasis", oasis=list(design_spec=list(hrf=hrf_B, ...)))
#' }
#'
#' @export
sbhm_hrf <- function(B, tgrid, span) {
  stopifnot(is.matrix(B), is.numeric(B), is.numeric(tgrid))
  if (nrow(B) != length(tgrid)) stop("nrow(B) must equal length(tgrid)")
  r <- ncol(B)

  # Closure capturing B and tgrid; returns len(t) × r matrix
  fun <- function(t) {
    t <- as.numeric(t)
    n <- length(t)
    out <- matrix(0, n, r)
    # Precompute interval indices on tgrid
    idx <- findInterval(t, tgrid, all.inside = TRUE)
    idx2 <- pmin(idx + 1L, length(tgrid))
    dt  <- tgrid[idx2] - tgrid[idx]
    w   <- (t - tgrid[idx]) / ifelse(dt == 0, 1, dt)
    for (j in seq_len(r)) {
      b1 <- B[cbind(idx,  j)]
      b2 <- B[cbind(idx2, j)]
      out[, j] <- (1 - w) * b1 + w * b2
    }
    out
  }
  fmrihrf::HRF(fun = fun, name = "SBHM", nbasis = r, span = span)
}

