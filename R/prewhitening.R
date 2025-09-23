#' Prewhiten fMRI data using AR/ARMA models
#'
#' Internal function providing a unified interface for prewhitening fMRI data
#' using the fmriAR package. Supports various AR/ARMA models with flexible
#' parameter estimation strategies.
#' @importFrom fmriAR fit_noise whiten_apply
#'
#' @param Y Numeric matrix (timepoints x voxels) of fMRI data
#' @param X Optional design matrix for trials (timepoints x trials)
#' @param Z Optional experimental design matrix (timepoints x regressors)
#' @param Nuisance Optional nuisance regressors (timepoints x nuisance)
#' @param prewhiten List of prewhitening options:
#'   \describe{
#'     \item{method}{Character: "ar" (default), "arma", or "none"}
#'     \item{p}{Integer or "auto": AR order (default "auto")}
#'     \item{q}{Integer: MA order for ARMA (default 0)}
#'     \item{p_max}{Integer: Maximum AR order when p="auto" (default 6)}
#'     \item{pooling}{Character: "global" (default), "voxel", "run", or "parcel"}
#'     \item{runs}{Integer vector: Run identifiers for run-aware estimation}
#'     \item{parcels}{Integer vector: Parcel memberships for parcel-based pooling}
#'     \item{exact_first}{Character: "ar1" or "none" for exact AR(1) scaling}
#'     \item{compute_residuals}{Logical: Whether to compute residuals first (default TRUE)}
#'   }
#' @return List containing:
#'   \describe{
#'     \item{Y_whitened}{Whitened data matrix}
#'     \item{X_whitened}{Whitened trial design matrix (if provided)}
#'     \item{Z_whitened}{Whitened experimental design (if provided)}
#'     \item{Nuisance_whitened}{Whitened nuisance regressors (if provided)}
#'     \item{whiten_plan}{fmriAR plan object for diagnostics}
#'     \item{applied}{Logical: Whether whitening was applied}
#'   }
#' @importFrom utils modifyList
#' @keywords internal
.prewhiten_data <- function(Y, X = NULL, Z = NULL, Nuisance = NULL,
                           prewhiten = list()) {

  # Set defaults
  defaults <- list(
    method = "ar",
    p = "auto",
    q = 0L,
    p_max = 6L,
    pooling = "global",
    runs = NULL,
    parcels = NULL,
    exact_first = "ar1",
    compute_residuals = TRUE
  )

  # Merge with user options
  opts <- modifyList(defaults, prewhiten)

  # Early exit if no whitening requested
  if (opts$method == "none" || is.null(opts$method)) {
    return(list(
      Y_whitened = Y,
      X_whitened = X,
      Z_whitened = Z,
      Nuisance_whitened = Nuisance,
      whiten_plan = NULL,
      applied = FALSE
    ))
  }

  # Validate inputs
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  n_time <- nrow(Y)
  n_vox <- ncol(Y)

  pool_mode <- opts$pooling
  if (length(pool_mode) == 0L || is.na(pool_mode[1])) {
    pool_mode <- "global"
  } else {
    pool_mode <- tolower(pool_mode[1])
  }
  pool_mode <- match.arg(pool_mode, c("global", "voxel", "run", "parcel"))
  if (pool_mode == "voxel") {
    if (!is.null(opts$parcels)) {
      warning("prewhiten$parcels is ignored when pooling='voxel'; using one parcel per voxel.")
    }
    opts$pooling <- "parcel"
    opts$parcels <- seq_len(n_vox)
  }
  plan_parcels <- opts$parcels
  if (is.null(plan_parcels)) {
    plan_parcels <- rep(1L, n_vox)
  } else {
    plan_parcels <- as.integer(plan_parcels)
    if (length(plan_parcels) != n_vox) {
      plan_parcels <- rep(plan_parcels, length.out = n_vox)
    }
  }
  parcel_ids <- function(n_cols) {
    if (n_cols == 0L) {
      integer(0)
    } else {
      rep(plan_parcels, length.out = n_cols)
    }
  }

  # Prepare residuals for AR estimation
  if (opts$compute_residuals) {
    # Combine all design matrices
    design_full <- NULL

    if (!is.null(Z)) {
      if (!is.matrix(Z)) Z <- as.matrix(Z)
      design_full <- cbind(design_full, Z)
    }

    if (!is.null(X)) {
      if (!is.matrix(X)) X <- as.matrix(X)
      design_full <- cbind(design_full, X)
    }

    if (!is.null(Nuisance)) {
      if (!is.matrix(Nuisance)) Nuisance <- as.matrix(Nuisance)
      design_full <- cbind(design_full, Nuisance)
    }

    # If we have design matrices, compute residuals
    if (!is.null(design_full)) {
      # Add intercept if not present
      if (!any(apply(design_full, 2, function(x) all(x == x[1])))) {
        design_full <- cbind(1, design_full)
      }

      # Compute OLS residuals
      resid <- Y - design_full %*% qr.solve(design_full, Y)
    } else {
      # No design matrices, use demeaned data
      resid <- sweep(Y, 2, colMeans(Y))
    }
  } else {
    # Use Y directly as "residuals" (e.g., if already residualized)
    resid <- Y
  }

  # Fit noise model using fmriAR
  whiten_plan <- fmriAR::fit_noise(
    resid = resid,
    runs = opts$runs,
    method = opts$method,
    p = opts$p,
    q = opts$q,
    p_max = opts$p_max,
    exact_first = opts$exact_first,
    pooling = opts$pooling,
    parcels = opts$parcels
  )

  # Apply whitening to all matrices
  # fmriAR requires both X and Y, so we use a dummy X if needed
  dummy_X <- matrix(1, nrow(Y), 1)

  whitened_Y <- fmriAR::whiten_apply(whiten_plan, X = dummy_X, Y = Y,
                                     parcels = parcel_ids(ncol(Y)))
  Y_whitened <- whitened_Y$Y

  X_whitened <- if (!is.null(X)) {
    whitened_X <- fmriAR::whiten_apply(whiten_plan, X = dummy_X, Y = X,
                                     parcels = parcel_ids(ncol(X)))
    whitened_X$Y
  } else NULL

  Z_whitened <- if (!is.null(Z)) {
    whitened_Z <- fmriAR::whiten_apply(whiten_plan, X = dummy_X, Y = Z,
                                     parcels = parcel_ids(ncol(Z)))
    whitened_Z$Y
  } else NULL

  Nuisance_whitened <- if (!is.null(Nuisance)) {
    whitened_N <- fmriAR::whiten_apply(whiten_plan, X = dummy_X, Y = Nuisance,
                                     parcels = parcel_ids(ncol(Nuisance)))
    whitened_N$Y
  } else NULL

  # Return whitened data and plan
  list(
    Y_whitened = Y_whitened,
    X_whitened = X_whitened,
    Z_whitened = Z_whitened,
    Nuisance_whitened = Nuisance_whitened,
    whiten_plan = whiten_plan,
    applied = TRUE
  )
}

#' Convert old-style whitening options to new format
#'
#' Internal helper to maintain backward compatibility with the old
#' oasis$whiten = "ar1" syntax by converting to new prewhiten format.
#'
#' @param oasis_opts List of OASIS options potentially containing whiten field
#' @return List of prewhiten options or NULL if no whitening
#' @keywords internal
.convert_legacy_whiten <- function(oasis_opts) {
  if (is.null(oasis_opts$whiten) || oasis_opts$whiten == "none") {
    return(NULL)
  }

  if (tolower(oasis_opts$whiten) == "ar1") {
    # Convert old AR(1) syntax to new format
    list(
      method = "ar",
      p = 1L,
      pooling = "global",  # Match old behavior of single rho
      exact_first = "ar1"
    )
  } else {
    warning("Unknown whiten option '", oasis_opts$whiten,
            "'. Use prewhiten parameter instead.")
    NULL
  }
}

#' Check if advanced prewhitening is needed
#'
#' Determines whether to use fmriAR (advanced) or keep simple AR(1)
#' based on the requested features.
#'
#' @param prewhiten List of prewhitening options
#' @return Logical: TRUE if fmriAR is needed, FALSE for simple AR(1)
#' @keywords internal
.needs_advanced_prewhitening <- function(prewhiten) {
  if (is.null(prewhiten)) return(FALSE)

  # Check for features that require fmriAR
  needs_fmriAR <- any(c(
    !is.null(prewhiten$method) && prewhiten$method == "arma",
    !is.null(prewhiten$p) && (prewhiten$p != 1L && prewhiten$p != "auto"),
    !is.null(prewhiten$q) && prewhiten$q > 0L,
    !is.null(prewhiten$pooling) && prewhiten$pooling != "global",
    !is.null(prewhiten$parcels),
    !is.null(prewhiten$runs)
  ))

  needs_fmriAR
}
