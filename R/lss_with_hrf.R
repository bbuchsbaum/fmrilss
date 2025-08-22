#' Least Squares Separate with voxel-wise HRF (basis-weighted)
#'
#' Compute LSS trial-wise betas when each voxel has its own HRF formed
#' as a linear combination of K basis kernels sampled on the TR grid.
#'
#' **Design & nuisance handling match `lss()`**:
#'   - The trial-of-interest (Xi) and the sum of all other trials (Xother) are
#'     included in each per-trial GLM.
#'   - If `Nuisance` is supplied, it is projected out of **Y** and the trial
#'     regressors before LSS (standard residualization). Experimental regressors
#'     `Z` are *not* residualized, matching `lss()` documentation.
#'   - If `Z` is `NULL`, an intercept-only design is used.
#'
#' @param Y numeric matrix (n_time x n_vox)
#' @param onset_idx integer vector (length n_trials), 1-based TR indices
#' @param durations numeric vector (length n_trials), in TRs; 0 means an impulse.
#'                  Uses inclusive indexing: if duration = d, samples \code{o:(o+d)} are 1.
#' @param hrf_basis_kernels numeric matrix (L x K), K basis kernels on TR grid
#' @param coefficients numeric matrix (K x n_vox), voxel-wise HRF weights
#' @param Z optional numeric matrix (n_time x F) of experimental regressors;
#'          if NULL, an intercept (column of 1s) is used.
#' @param Nuisance optional numeric matrix (n_time x q) of confounds to project out
#' @param verbose logical; print progress every 1000 voxels
#' @param method character: "r" (default, pure R), "cpp" (C++ backend), 
#'   "cpp_arma" (Armadillo backend), or "cpp_omp" (OpenMP parallel backend).
#'   Falls back automatically: cpp_omp -> cpp_arma -> cpp -> r.
#'
#' @return numeric matrix (n_trials x n_vox) of trial-wise beta estimates
#' @examples
#' \dontrun{
#' # Minimal use (R backend):
#' betas <- lss_with_hrf_pure_r(Y, onset_idx, durations, basis, coeffs, Z = NULL, Nuisance = NULL)
#' # Or with C++ backend:
#' betas <- lss_with_hrf_pure_r(Y, onset_idx, durations, basis, coeffs, method = "cpp")
#' }
#' @keywords internal
lss_with_hrf_pure_r <- function(
  Y,
  onset_idx,
  durations = NULL,
  hrf_basis_kernels,
  coefficients,
  Z = NULL,
  Nuisance = NULL,
  verbose = FALSE,
  method = c("r", "cpp", "cpp_arma", "cpp_omp")
) {
  method <- match.arg(method)
  # ---- basic checks ----
  if (!is.matrix(Y)) stop("Y must be a matrix [n_time x n_vox].")
  if (!is.matrix(hrf_basis_kernels)) stop("hrf_basis_kernels must be a matrix [L x K].")
  if (!is.matrix(coefficients)) stop("coefficients must be a matrix [K x n_vox].")

  n_time <- nrow(Y)
  n_vox  <- ncol(Y)

  if (!is.integer(onset_idx)) onset_idx <- as.integer(onset_idx)
  n_trials <- length(onset_idx)
  if (n_trials < 1L) stop("Need at least one trial (length(onset_idx) >= 1).")

  if (is.null(durations)) durations <- rep(0L, n_trials)
  if (length(durations) != n_trials) stop("durations must match length(onset_idx).")
  durations <- as.integer(round(pmax(0, durations)))

  K <- ncol(hrf_basis_kernels)
  if (nrow(coefficients) != K) stop("nrow(coefficients) must equal ncol(hrf_basis_kernels).")
  if (ncol(coefficients) != n_vox) stop("ncol(coefficients) must equal ncol(Y).")

  if (!is.null(Z) && nrow(Z) != n_time) stop("Z must have nrow == n_time.")
  if (!is.null(Nuisance) && nrow(Nuisance) != n_time) stop("Nuisance must have nrow == n_time.")

  # ---- helper: safe open convolution, then truncate to n_time ----
  conv_open_trim <- function(x, k) {
    # convolve() expects rev(k) to do standard causal filtering
    as.numeric(stats::convolve(x, rev(as.numeric(k)), type = "open"))[seq_len(length(x))]
  }

  # ---- 1) Event "stick" design on TR grid: n_time x n_trials ----
  Xev <- matrix(0, n_time, n_trials)
  for (i in seq_len(n_trials)) {
    o <- onset_idx[i]
    if (is.na(o) || o < 1L || o > n_time) next
    d <- durations[i]
    i1 <- o
    i2 <- min(n_time, o + d)
    Xev[i1:i2, i] <- 1
  }

  # ---- 2) Convolve events with each basis kernel -> list of n_time x n_trials ----
  basis_convolved <- vector("list", K)
  for (k in seq_len(K)) {
    bk <- as.numeric(hrf_basis_kernels[, k])
    out <- vapply(seq_len(n_trials), function(j) conv_open_trim(Xev[, j], bk),
                  numeric(n_time))
    storage.mode(out) <- "double"
    basis_convolved[[k]] <- out
  }

  # ---- 3) Residualize Y and trial designs by Nuisance (match lss()) ----
  # Use QR-based residualization without forming a dense projection matrix.
  if (!is.null(Nuisance)) {
    qrN <- qr(Nuisance)
    # residualize Y
    Y <- Y - Nuisance %*% qr.coef(qrN, Y)
    # residualize each basis-convolved design
    for (k in seq_len(K)) {
      Dk <- basis_convolved[[k]]
      basis_convolved[[k]] <- Dk - Nuisance %*% qr.coef(qrN, Dk)
    }
    # NOTE: Do NOT residualize Z (matches lss() semantics).
  }

  # ---- 4) Prepare Z (experimental regressors); default intercept-only ----
  if (is.null(Z)) {
    Z_use <- matrix(1, n_time, 1L)
  } else {
    Z_use <- as.matrix(Z)
  }
  pz <- ncol(Z_use)

  # ---- output ----
  betas <- matrix(NA_real_, n_trials, n_vox)
  colnames(betas) <- colnames(Y)
  rownames(betas) <- paste0("trial_", seq_len(n_trials))

  # ---- optional: C++/Armadillo/OpenMP backends ----
  if (method != "r") {
    # Fallback chain: cpp_omp -> cpp_arma -> cpp -> r
    methods_to_try <- switch(method,
      cpp_omp = c("cpp_omp", "cpp_arma", "cpp", "r"),
      cpp_arma = c("cpp_arma", "cpp", "r"),
      cpp = c("cpp", "r"),
      "r"
    )
    
    for (try_method in methods_to_try) {
      if (try_method == "r") break  # Will use R implementation below
      
      backend_fn <- switch(try_method,
        cpp_omp = "lss_engine_vox_hrf_omp",
        cpp_arma = "lss_engine_vox_hrf_arma",
        cpp = "lss_engine_vox_hrf_cpp"
      )
      
      have_backend <- FALSE
      try({
        get(backend_fn, envir = asNamespace("fmrilss"))
        have_backend <- TRUE
      }, silent = TRUE)
      
      if (have_backend) {
        if (verbose && try_method != method) {
          message("Using ", try_method, " backend (", method, " not available)")
        }
        
        # Call the appropriate backend
        betas_cpp <- switch(try_method,
          cpp_omp = fmrilss:::lss_engine_vox_hrf_omp(
            Y, coefficients, basis_convolved, Z_use
          ),
          cpp_arma = fmrilss:::lss_engine_vox_hrf_arma(
            Y, coefficients, basis_convolved, Z_use
          ),
          cpp = fmrilss:::lss_engine_vox_hrf_cpp(
            Y, coefficients, basis_convolved, Z_use
          )
        )
        
        dimnames(betas_cpp) <- dimnames(betas)
        return(betas_cpp)
      }
    }
    
    if (verbose) message("All C++ backends unavailable, falling back to R")
  }

  # ---- 5) For each voxel, combine basis designs with that voxel's HRF weights ----
  for (v in seq_len(n_vox)) {
    wv <- as.numeric(coefficients[, v])  # length K

    # X_v = sum_k wv[k] * Dk  (n_time x n_trials)
    X_v <- matrix(0.0, n_time, n_trials)
    for (k in seq_len(K)) {
      wk <- wv[k]
      if (wk == 0.0) next
      X_v <- X_v + basis_convolved[[k]] * wk
    }

    # precompute the "all others" column once
    x_all <- rowSums(X_v)

    yv <- Y[, v]

    # ---- 6) Fit per-trial GLMs ----
    for (i in seq_len(n_trials)) {
      Xi <- X_v[, i]
      if (n_trials == 1L) {
        Xother <- NULL
        Xdesign <- cbind(Z_use, Xi)
      } else {
        Xother <- x_all - Xi
        Xdesign <- cbind(Z_use, Xi, Xother)
      }

      fit <- try(stats::lm.fit(Xdesign, yv), silent = TRUE)
      if (inherits(fit, "try-error") || is.null(fit$coefficients)) {
        betas[i, v] <- NA_real_
      } else {
        coef <- fit$coefficients
        # Xi is the (pz + 1)-th coefficient
        if (length(coef) >= pz + 1L && !is.na(coef[pz + 1L])) {
          betas[i, v] <- coef[pz + 1L]
        } else {
          betas[i, v] <- NA_real_
        }
      }
    }

    if (verbose && (v %% 1000L == 0L)) message("Processed voxel ", v, " / ", n_vox)
  }

  betas
}