#' Build a Shared-Basis HRF Library (SBHM)
#'
#' Learn a low-rank shared time basis from a parameterized HRF library using
#' `fmrihrf::hrf_library()`. The library is evaluated on the TR grid, optionally
#' baseline-removed and L2-normalized, and decomposed via SVD into
#' `B = U_r` (shared basis), singular values `S`, and library coordinates
#' `A = diag(S) %*% t(V_r)`.
#'
#' @param library_spec Either NULL (when `library_H` is provided) or a list with:
#'   - `fun`: a function compatible with `fmrihrf::hrf_library(fun, pgrid, ...)`
#'            that returns an `fmrihrf` HRF object when called with parameters.
#'   - `pgrid`: a data.frame of parameter combinations (see examples).
#'   - `span`: numeric, HRF span in seconds (default `span`).
#'   - `precision`: numeric, evaluation precision (default 0.1 sec).
#'   - `method`: evaluation method for `fmrihrf::evaluate()` (default "conv").
#'   - `extras`: optional list of additional arguments passed to `hrf_library`.
#' @param library_H Optional precomputed T×K matrix of candidate HRFs, already
#'   aligned to the TR grid `tgrid` (or `sframe`). Mutually exclusive with `library_spec`.
#' @param r Target rank for the shared basis (default 6). Clipped to `min(T, K)`.
#' @param sframe Optional `fmrihrf::sampling_frame`, used to derive the global time
#'   grid when `tgrid` is not provided.
#' @param tgrid Optional numeric vector of global times (in seconds). If provided,
#'   takes precedence over `sframe`.
#' @param span HRF span in seconds (default 32). Used for reference HRF when needed.
#' @param normalize Logical, L2-normalize library columns (default TRUE).
#' @param baseline Numeric length-2 vector specifying a time window (in seconds)
#'   used for baseline removal (column-wise mean subtraction within this window).
#'   Set to NULL to skip.
#' @param shifts Optional numeric vector of time shifts (in seconds). When provided,
#'   shifted copies of the library are added by linear interpolation on `tgrid`.
#' @param ref Reference for coefficient-space shrinkage and orientation.
#'   One of `"mean"` (default) or `"spmg1"`. If `"spmg1"`, the SPMG1 HRF is
#'   projected onto the learned basis to form `alpha_ref`.
#'
#' @return A list with components:
#'   - `B` (T×r): shared orthonormal time basis
#'   - `S` (length r): singular values
#'   - `A` (r×K): coordinates of library HRFs in the shared basis
#'   - `tgrid`: the global time grid used (seconds)
#'   - `span`: span used for reference HRF
#'   - `ref`: list with `alpha_ref` (length r) and `name`
#'   - `meta`: list with `r`, `K`, `normalize`, `baseline`
#'
#' @examples
#' \dontrun{
#'   library(fmrihrf)
#'   # Build a gamma HRF library across a small parameter grid
#'   param_grid <- expand.grid(shape = c(6, 8, 10), rate = c(0.9, 1.0, 1.1))
#'   gamma_fun  <- function(shape, rate) fmrihrf::as_hrf(
#'     fmrihrf::hrf_gamma, params = list(shape = shape, rate = rate)
#'   )
#'
#'   sframe <- fmrihrf::sampling_frame(blocklens = 200, TR = 1)
#'   sbhm <- sbhm_build(
#'     library_spec = list(fun = gamma_fun, pgrid = param_grid, span = 32),
#'     r = 6, sframe = sframe, baseline = c(0, 0.5)
#'   )
#'
#'   # Use the learned basis as an HRF in OASIS designs
#'   hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)
#' }
#'
#' @export
sbhm_build <- function(library_spec = NULL,
                       library_H    = NULL,
                       r = 6,
                       sframe = NULL,
                       tgrid  = NULL,
                       span   = 32,
                       normalize = TRUE,
                       baseline  = c(0, 0.5),
                       shifts = NULL,
                       ref = c("mean","spmg1")) {
  ref <- match.arg(ref)
  if (is.null(library_spec) == is.null(library_H)) {
    stop("Provide exactly one of library_spec or library_H")
  }

  # Establish evaluation grid (seconds)
  if (!is.null(tgrid)) {
    times <- as.numeric(tgrid)
  } else {
    if (is.null(sframe)) stop("Provide sframe or tgrid")
    times <- fmrihrf::samples(sframe, global = TRUE)
  }

  # 1) Build/evaluate library matrix H (T×K)
  if (!is.null(library_H)) {
    H <- as.matrix(library_H)
    if (nrow(H) != length(times)) stop("library_H rows must match length(tgrid)")
  } else {
    if (!is.list(library_spec) || !is.function(library_spec$fun) || !is.data.frame(library_spec$pgrid)) {
      stop("library_spec must be a list with elements fun (function) and pgrid (data.frame)")
    }
    extras <- library_spec$extras
    if (is.null(extras)) extras <- list()
    lib_hrf <- do.call(fmrihrf::hrf_library,
                       c(list(fun = library_spec$fun, pgrid = library_spec$pgrid), extras))

    rr <- fmrihrf::regressor(
      onsets   = 0,
      hrf      = lib_hrf,
      duration = 0,
      span     = if (!is.null(library_spec$span)) library_spec$span else span,
      summate  = FALSE
    )
    H <- fmrihrf::evaluate(
      rr,
      grid      = times,
      precision = if (!is.null(library_spec$precision)) library_spec$precision else 0.1,
      method    = if (!is.null(library_spec$method)) library_spec$method else "conv"
    )
    if (inherits(H, "Matrix")) H <- as.matrix(H)
    if (!is.matrix(H)) H <- cbind(H)
  }

  # 2) Optional time-shifts: augment columns by linear interpolation on tgrid
  if (!is.null(shifts) && length(shifts) > 0) {
    add_shift_block <- function(mat, s) {
      # Fast column-wise linear interpolation using approx for each column
      out <- matrix(NA_real_, nrow(mat), ncol(mat))
      for (j in seq_len(ncol(mat))) {
        out[, j] <- stats::approx(times, mat[, j], xout = times - s, rule = 2, ties = "ordered")$y
      }
      out
    }
    shifted_list <- lapply(shifts, function(s) add_shift_block(H, s))
    H <- cbind(H, do.call(cbind, shifted_list))
  }

  # 3) Baseline removal within window and L2 normalization per column
  if (!is.null(baseline) && length(baseline) == 2) {
    idx <- which(times >= baseline[1] & times <= baseline[2])
    if (length(idx) > 0) {
      H <- sweep(H, 2L, colMeans(H[idx, , drop = FALSE]), "-")
    }
  }
  if (isTRUE(normalize)) {
    H <- sweep(H, 2L, sqrt(colSums(H^2)) + 1e-8, "/")
  }

  # 4) SVD → shared basis (rank r)
  r_eff <- min(r, nrow(H), ncol(H))
  sv <- La.svd(H, nu = r_eff, nv = r_eff)
  # Sanity check: detect degenerate libraries (e.g., identical columns)
  if (r_eff >= 2) {
    if (!is.null(sv$d) && length(sv$d) >= 2) {
      if (sv$d[2] <= (1e-8 * sv$d[1])) {
        warning("HRF library appears nearly rank-1 (second singular value ≈ 0). " ,
                "Check that `library_spec$fun` closes over its parameters. For example:\n",
                "  fun <- function(shape, rate) as_hrf(function(t) hrf_gamma(t, shape=shape, rate=rate), span=32)")
      }
    }
  }
  B  <- sv$u[, seq_len(r_eff), drop = FALSE]
  S  <- sv$d[seq_len(r_eff)]
  Vt <- sv$vt[seq_len(r_eff), , drop = FALSE]
  A  <- diag(S, r_eff, r_eff) %*% Vt  # r x K

  # 5) Reference alpha for shrinkage/matching
  alpha_ref <- switch(ref,
    mean  = rowMeans(A),
    spmg1 = {
      rr1 <- fmrihrf::regressor(onsets = 0,
                                hrf = fmrihrf::make_hrf("spmg1"),
                                duration = 0,
                                span = span,
                                summate = TRUE)
      x1 <- fmrihrf::evaluate(rr1, grid = times, precision = 0.1, method = "conv")
      x1 <- as.numeric(if (is.matrix(x1)) x1[, 1] else x1)
      drop(crossprod(B, x1))
    }
  )

  list(
    B = B,
    S = S,
    A = A,
    tgrid = times,
    span = span,
    ref = list(alpha_ref = as.numeric(alpha_ref), name = ref),
    meta = list(r = r_eff, K = ncol(A), normalize = normalize, baseline = baseline)
  )
}
