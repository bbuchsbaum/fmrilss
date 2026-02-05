# Fast per-voxel HRF estimation helpers for OASIS VOXHRF mode

# --- Internal: Fast per-voxel HRF estimation via aggregated per-basis regressors ---
#' @keywords internal
#' @noRd
.estimate_voxel_hrf_fast <- function(Y, X_trials, design_spec, N_nuis = NULL, K = NULL,
                                     lambda_shape = 0, mu_rough = 0,
                                     ref_hrf = NULL, shrink_global = 0, orient_ref = TRUE) {
  # Coerce
  Y <- .as_base_matrix(Y)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  X_trials <- .as_base_matrix(X_trials)
  if (!is.matrix(X_trials)) X_trials <- as.matrix(X_trials)
  N_nuis <- .as_base_matrix(N_nuis)
  if (!is.null(N_nuis) && !is.matrix(N_nuis)) N_nuis <- as.matrix(N_nuis)

  # Prefer HRF nbasis for K; fall back to detection when needed
  K_hrf <- tryCatch(fmrihrf::nbasis(design_spec$cond$hrf), error = function(e) NULL)
  if (is.null(K)) {
    K <- if (!is.null(K_hrf)) as.integer(K_hrf) else .detect_basis_dimension(X_trials)
  } else {
    K <- as.integer(K)
  }
  if ((ncol(X_trials) %% K) != 0L) {
    K <- .detect_basis_dimension(X_trials)
  }
  stopifnot(ncol(X_trials) %% K == 0L)

  # 1) Build aggregated per-basis regressors A [T x K] by summing K-blocked columns
  idx_by_basis <- lapply(seq_len(K), function(k) seq.int(k, ncol(X_trials), by = K))
  A <- do.call(cbind, lapply(idx_by_basis, function(idx) rowSums(X_trials[, idx, drop = FALSE])))

  # 2) Residualize Y and A by confounds (Z+Nuisance+others already combined into N_nuis)
  if (!is.null(N_nuis) && ncol(N_nuis) > 0) {
    qrN <- qr(N_nuis)
    Y <- qr.resid(qrN, Y)
    A <- qr.resid(qrN, A)
  }

  # 3) Shape-prior ridge in basis space
  sframe <- design_spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)

  # Evaluate basis kernels on TR grid using a single onset at 0
  span <- design_spec$cond$span %||% (max(times) - min(times)) %||% 40
  r <- fmrihrf::regressor(onsets = 0, hrf = design_spec$cond$hrf, duration = 0, span = span)
  B_time <- fmrihrf::evaluate(r, grid = times, precision = design_spec$precision %||% 0.1,
                              method = design_spec$method %||% "conv")
  if (inherits(B_time, "Matrix")) B_time <- as.matrix(B_time)
  if (is.vector(B_time)) B_time <- cbind(B_time)
  # If evaluated basis dimension disagrees with K, rebuild A accordingly
  K_eval <- ncol(B_time)
  if (K_eval != K) {
    K <- as.integer(K_eval)
    stopifnot(ncol(X_trials) %% K == 0L)
    idx_by_basis <- lapply(seq_len(K), function(k) seq.int(k, ncol(X_trials), by = K))
    A <- do.call(cbind, lapply(idx_by_basis, function(idx) rowSums(X_trials[, idx, drop = FALSE])))
  }

  # Reference HRF projection -> W0 (K-vector)
  if (is.null(ref_hrf)) {
    # Default: use first basis kernel as reference shape
    ref_hrf <- B_time[, 1]
  } else {
    # Coerce to a vector of same length as time grid if needed
    ref_hrf <- as.numeric(ref_hrf)
    if (length(ref_hrf) != nrow(B_time)) {
      # If provided on a different grid, evaluate HRF object if possible; otherwise fallback to first kernel
      ref_hrf <- tryCatch({
        rr <- fmrihrf::regressor(onsets = 0, hrf = ref_hrf, duration = 0, span = span)
        as.numeric(fmrihrf::evaluate(rr, grid = times, precision = design_spec$precision %||% 0.1,
                                     method = design_spec$method %||% "conv"))
      }, error = function(e) B_time[, 1])
    }
  }
  BtB <- crossprod(B_time)                          # K x K
  # Small ridge on BtB to avoid singularity in projection
  W0 <- drop(solve(BtB + 1e-8 * diag(ncol(B_time)), crossprod(B_time, ref_hrf)))

  # Roughness penalty R in basis space: R = (D2 B)^T (D2 B)
  if (mu_rough > 0) {
    D2 <- diff(diag(nrow(B_time)), differences = 2)
    R  <- crossprod(D2 %*% B_time)
  } else {
    R <- matrix(0, ncol(B_time), ncol(B_time))
  }

  # 4) Solve tiny KxK ridge system batched over voxels
  AtA <- crossprod(A)                              # K x K
  AtY <- crossprod(A, Y)                           # K x V
  P   <- AtA + lambda_shape * diag(K) + mu_rough * R
  rhs <- AtY + lambda_shape * W0 %*% matrix(1, 1, ncol(Y))
  W   <- solve(P, rhs)
  W   <- matrix(W, nrow = K, ncol = ncol(Y))        # ensure matrix even when K=1

  # 5) Normalize each voxel's HRF shape to unit energy on the time grid
  # energy_v = sqrt( diag(W^T (B^T B) W) ) computed without forming H = B W
  M <- BtB
  MW <- M %*% W                                   # K x V
  s2 <- sqrt(pmax(colSums(W * MW), .Machine$double.eps))
  W  <- sweep(W, 2L, s2, "/")
  W  <- matrix(W, nrow = K, ncol = ncol(Y))

  ref_norm <- sqrt(sum((B_time %*% W0)^2))

  # 6) Optional global shrinkage for stability
  if (shrink_global > 0) {
    Wbar <- rowMeans(W)
    W <- (1 - shrink_global) * W + shrink_global * matrix(Wbar, nrow = K, ncol = ncol(W))
    W <- matrix(W, nrow = K, ncol = ncol(W))
  }

  # 7) Optional orientation: align with reference HRF so betas are positive-interpretable
  if (isTRUE(orient_ref)) {
    cvec <- as.numeric(crossprod(B_time, ref_hrf))   # length K
    if (length(cvec) != nrow(W)) {
      stop("Orientation reference length mismatch with basis dimension")
    }
    dots <- as.numeric(crossprod(cvec, W))          # length V
    flip <- which(dots < 0)
    if (length(flip)) W[, flip] <- -W[, flip, drop = FALSE]
  }

  structure(list(
    coefficients = W,                          # K x V
    basis       = design_spec$cond$hrf,        # fmrihrf basis object
    conditions  = "cond",
    ref_norm    = as.numeric(ref_norm)
  ), class = "VoxelHRF")
}
