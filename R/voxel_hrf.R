#' VoxelHRF object
#'
#' Simple list-based S3 class returned by \code{estimate_voxel_hrf} containing
#' voxel-wise HRF basis coefficients and related metadata.
#' @name VoxelHRF
NULL

#' LSSBeta object
#'
#' Simple list-based S3 class returned by \code{lss_with_hrf} containing
#' trial-wise beta estimates.
#' @name LSSBeta
NULL

#' Estimate Voxel-wise HRF Basis Coefficients
#'
#' Fits a GLM to estimate HRF basis coefficients for every voxel.
#'
#' @param Y Numeric matrix of BOLD data (time \times voxels).
#' @param events Data frame with \code{onset}, \code{duration} and \code{condition} columns.
#' @param basis HRF object from the \code{fmrihrf} package.
#' @param nuisance_regs Optional numeric matrix of nuisance regressors.
#'
#' @return A \link{VoxelHRF} object containing at least:
#'   \item{coefficients}{Matrix of HRF basis coefficients.}
#'   \item{basis}{The HRF basis object used.}
#'   \item{conditions}{Character vector of modeled conditions.}
#'
#' @examples
#' if (requireNamespace("fmrihrf", quietly = TRUE)) {
#'   set.seed(1)
#'   Y <- matrix(rnorm(100), 50, 2)
#'   events <- data.frame(onset = c(5, 25), duration = 1,
#'                        condition = "A")
#'   basis <- fmrihrf::HRF_GAMMA()
#'   X <- fmrihrf::regressor_set(events, basis, n = nrow(Y))$X
#'   coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
#'   Y <- X %*% coef + Y * 0.1
#'   est <- estimate_voxel_hrf(Y, events, basis)
#'   str(est)
#' }
#' @export
estimate_voxel_hrf <- function(Y, events, basis, nuisance_regs = NULL) {
  # Input validation
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix")
  }

  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition")
  }

  if (!inherits(basis, "HRF")) {
    stop("basis must be an 'HRF' object from fmrihrf")
  }

  if (!is.null(nuisance_regs)) {
    if (!is.matrix(nuisance_regs) || !is.numeric(nuisance_regs)) {
      stop("nuisance_regs must be a numeric matrix")
    }
    if (nrow(nuisance_regs) != nrow(Y)) {
      stop("nuisance_regs must have same number of rows as Y")
    }
  }

  # Construct HRF regressor matrix using fmrihrf
  # Create sampling frame and time grid
  sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1)
  times <- fmrihrf::samples(sframe, global = TRUE)
  
  # Build regressor set with correct API
  rset <- fmrihrf::regressor_set(
    onsets = events$onset,
    fac = factor(seq_len(nrow(events))),
    hrf = basis,
    duration = events$duration,
    span = if (!is.null(attr(basis, "span"))) attr(basis, "span") else 30
  )
  X_basis <- fmrihrf::evaluate(rset, grid = times, precision = 0.1, method = "conv")
  if (inherits(X_basis, "Matrix")) X_basis <- as.matrix(X_basis)

  X_full <- X_basis
  if (!is.null(nuisance_regs)) {
    X_full <- cbind(X_full, nuisance_regs)
  }

  coef_mat <- estimate_hrf_cpp(X_full, Y)

  coef_basis <- coef_mat[seq_len(ncol(X_basis)), , drop = FALSE]

  result <- list(
    coefficients = coef_basis,
    basis = basis,
    conditions = unique(as.character(events$condition))
  )
  class(result) <- "VoxelHRF"
  result
}

#' Perform LSS using Voxel-wise HRFs
#'
#' Computes trial-wise beta estimates using voxel-specific HRFs.
#'
#' @param Y Numeric matrix of BOLD data (time \times voxels).
#' @param events Data frame with \code{onset}, \code{duration} and \code{condition} columns.
#' @param hrf_estimates A \link{VoxelHRF} object returned by \code{estimate_voxel_hrf}.
#' @param nuisance_regs Optional numeric matrix of nuisance regressors.
#' @param engine Computational engine: "R" for pure R implementation (default), 
#'   "C++" for optimized C++ (experimental).
#' @param chunk_size Number of voxels to process per batch (C++ engine only).
#' @param verbose Logical; display progress bar.
#' @param backing_dir Directory for bigmemory backing files. If NULL, a
#'   temporary directory is used (C++ engine only).
#'
#' @return An object of class \link{LSSBeta} for C++ engine, or a numeric matrix
#'   (n_trials x n_vox) for R engine.
#'
#' @examples
#' if (requireNamespace("fmrihrf", quietly = TRUE)) {
#'   set.seed(1)
#'   Y <- matrix(rnorm(100), 50, 2)
#'   events <- data.frame(onset = c(5, 25), duration = 1,
#'                        condition = "A")
#'   basis <- fmrihrf::HRF_GAMMA()
#'   X <- fmrihrf::regressor_set(events, basis, n = nrow(Y))$X
#'   coef <- matrix(rnorm(ncol(X) * ncol(Y)), ncol(X), ncol(Y))
#'   Y <- X %*% coef + Y * 0.1
#'   est <- estimate_voxel_hrf(Y, events, basis)
#'   betas <- lss_with_hrf(Y, events, est, verbose = FALSE, engine = "R")
#'   dim(betas)
#' }
#' @export
lss_with_hrf <- function(Y, events, hrf_estimates, nuisance_regs = NULL,
                         engine = "R", chunk_size = 5000, verbose = TRUE,
                         backing_dir = NULL) {
  engine <- match.arg(engine, c("R", "C++"))

  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix")
  }

  required_cols <- c("onset", "duration", "condition")
  if (!is.data.frame(events) || !all(required_cols %in% names(events))) {
    stop("events must be a data.frame with columns onset, duration, condition")
  }

  if (!inherits(hrf_estimates, "VoxelHRF") ||
      !is.matrix(hrf_estimates$coefficients)) {
    stop("hrf_estimates must be a VoxelHRF object with coefficient matrix")
  }

  if (!is.null(nuisance_regs)) {
    if (!is.matrix(nuisance_regs) || !is.numeric(nuisance_regs) ||
        nrow(nuisance_regs) != nrow(Y)) {
      stop("nuisance_regs must be a numeric matrix with same number of rows as Y")
    }
  }

  if (!is.numeric(chunk_size) || length(chunk_size) != 1 || chunk_size <= 0) {
    stop("chunk_size must be a positive integer")
  }
  
  # Use R implementation if requested
  if (engine == "R") {
    # Convert events to onset indices (1-based TR indices)
    onset_idx <- as.integer(round(events$onset))
    durations <- as.integer(round(events$duration))
    
    # Get HRF basis kernels on TR grid
    # Evaluate basis at TR timepoints
    n_time <- nrow(Y)
    hrf_span <- if (!is.null(attr(hrf_estimates$basis, "span"))) {
      attr(hrf_estimates$basis, "span") 
    } else {
      30
    }
    
    # Create TR grid for HRF evaluation
    tr_grid <- seq(0, hrf_span, by = 1)  # Assuming TR = 1
    
    # Evaluate basis to get kernels
    if (inherits(hrf_estimates$basis, "HRF")) {
      # Use fmrihrf to evaluate basis
      r <- fmrihrf::regressor(onsets = 0, hrf = hrf_estimates$basis, duration = 0, span = hrf_span)
      hrf_basis_kernels <- fmrihrf::evaluate(r, tr_grid, precision = 1, method = "conv")
      if (inherits(hrf_basis_kernels, "Matrix")) {
        hrf_basis_kernels <- as.matrix(hrf_basis_kernels)
      }
      # Ensure it's a matrix even if single column
      if (!is.matrix(hrf_basis_kernels)) {
        hrf_basis_kernels <- as.matrix(hrf_basis_kernels)
      }
    } else {
      stop("Cannot evaluate HRF basis - unsupported type")
    }
    
    # Call the pure R implementation
    betas <- lss_with_hrf_pure_r(
      Y = Y,
      onset_idx = onset_idx,
      durations = durations,
      hrf_basis_kernels = hrf_basis_kernels,
      coefficients = hrf_estimates$coefficients,
      Z = NULL,
      Nuisance = nuisance_regs,
      verbose = verbose
    )
    
    return(betas)
  }

  onsets <- as.numeric(events$onset)
  durations <- as.numeric(events$duration)
  conditions <- as.character(events$condition)

  fine_dt <- 0.1
  # Get HRF span from the basis object attribute
  hrf_span <- if (!is.null(attr(hrf_estimates$basis, "span"))) {
    attr(hrf_estimates$basis, "span") 
  } else {
    30
  }
  max_time <- max(onsets + durations) + hrf_span
  fine_grid <- seq(0, max_time, by = fine_dt)
  # Evaluate HRF at fine grid points
  hrf_func <- hrf_estimates$basis
  if (is.function(hrf_func)) {
    hrf_basis_kernels <- hrf_func(fine_grid)
  } else {
    # For HRF objects, create a simple regressor and evaluate
    r <- fmrihrf::regressor(onsets = 0, hrf = hrf_func, duration = 0, span = hrf_span)
    hrf_basis_kernels <- fmrihrf::evaluate(r, fine_grid, precision = fine_dt, method = "conv")
  }
  if (inherits(hrf_basis_kernels, "Matrix")) {
    hrf_basis_kernels <- as.matrix(hrf_basis_kernels)
  }

  onset_idx <- as.integer(round(onsets / fine_dt)) + 1L

  if (is.null(backing_dir)) {
    backing_dir <- tempdir()
  }
  bfile <- tempfile("betas", tmpdir = backing_dir, fileext = ".bin")
  descfile <- tempfile("betas", tmpdir = backing_dir, fileext = ".desc")
  betas <- bigmemory::filebacked.big.matrix(
    nrow = length(onsets), ncol = ncol(Y), type = "double",
    backingfile = basename(bfile), descriptorfile = basename(descfile),
    backingpath = dirname(bfile)
  )

  pb <- NULL
  update_progress <- function(step) {}
  if (verbose) {
    pb <- progress::progress_bar$new(total = ncol(Y))
    update_progress <- function(step) pb$tick(step)
  }

  lss_engine_vox_hrf(Y, hrf_estimates$coefficients, hrf_basis_kernels,
                     onset_idx, durations, nuisance_regs,
                     betas@address, update_progress,
                     as.integer(chunk_size), verbose)

  result <- list(betas = betas)
  class(result) <- "LSSBeta"
  result
}
