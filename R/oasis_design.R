# OASIS design construction + HRF/grid helpers

# --- Helper: build X (trial-wise) and X_other (aggregates for other conditions) via fmrihrf ---
#' @keywords internal
#' @noRd
.oasis_build_X_from_events <- function(spec) {
  if (is.null(spec)) stop("design_spec must be provided when X is NULL.")
  
  sframe <- spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)   # global seconds grid

  # Minimal safety for multi-run: warn if onsets look run‑relative
  # Heuristic: in multi-run, if max(onsets) is no larger than a single run's
  # duration in seconds, users likely supplied run-relative onsets.
  bl <- tryCatch(fmrihrf::blocklens(sframe), error = function(e) NULL)
  if (!is.null(bl) && length(bl) > 1L && !is.null(spec$cond$onsets)) {
    TR <- as.numeric(median(diff(times)))
    run_sec <- max(bl) * TR
    onv <- tryCatch(as.numeric(unlist(spec$cond$onsets)), error = function(e) NA_real_)
    max_onset <- suppressWarnings(max(onv, na.rm = TRUE))
    if (is.finite(max_onset) && max_onset <= run_sec + 1e-6) {
      warning(
        paste0(
          "Onsets appear run-relative but design_spec expects global time.\n",
          "For multi-run designs, prefer lss_design() with fmridesign::event_model(),\n",
          "or convert onsets to global seconds (offset each run)."
        ), call. = FALSE
      )
    }
  }
  
  # Get HRF object
  hrf_obj <- spec$cond$hrf %||% fmrihrf::make_hrf("spmg1")
  
  # Detect K from HRF
  K <- tryCatch(fmrihrf::nbasis(hrf_obj), error = function(e) 1L)
  
  # Build trial-wise design
  fac    <- factor(seq_along(spec$cond$onsets))
  rset   <- fmrihrf::regressor_set(onsets   = spec$cond$onsets,
                                   fac      = fac,
                                   hrf      = hrf_obj,
                                   duration = spec$cond$duration %||% 0,
                                   amplitude= spec$cond$amplitude %||% 1,
                                   span     = spec$cond$span %||% 40,
                                   summate  = TRUE)
  X_trials <- fmrihrf::evaluate(rset, grid = times, 
                                precision = spec$precision %||% 0.1,
                                method = spec$method %||% "conv")
  X_trials <- if (inherits(X_trials, "Matrix")) as.matrix(X_trials) else X_trials

  ntrials <- length(spec$cond$onsets)
  output_names <- if (K == 1L) {
    .default_trial_names(ntrials)
  } else {
    as.vector(t(outer(
      .default_trial_names(ntrials),
      sprintf("Basis_%d", seq_len(K)), paste, sep = ":"
    )))
  }
  if (length(output_names) != ncol(X_trials)) {
    stop("event-built trial design dimensions disagree with its HRF basis", call. = FALSE)
  }
  colnames(X_trials) <- output_names
  trial_basis_map <- data.frame(
    column = output_names,
    trial = rep(seq_len(ntrials), each = K),
    basis = rep(seq_len(K), times = ntrials),
    output_name = output_names,
    stringsAsFactors = FALSE
  )

  X_other <- .oasis_build_X_other(spec, hrf_obj = hrf_obj, times = times)

  return(list(
    X_trials = X_trials, X_other = X_other, K = K,
    trial_basis_map = trial_basis_map
  ))
}

# Build the aggregate basis columns for modeled non-target conditions. The
# target HRF is the fallback for an `others` entry that does not specify one.
# Keeping this separate lets HRF-grid selection construct the exact common
# span associated with each candidate.
.oasis_build_X_other <- function(spec, hrf_obj, times = NULL) {
  if (!length(spec$others)) return(NULL)
  if (is.null(times)) times <- fmrihrf::samples(spec$sframe, global = TRUE)

  X_other <- do.call(cbind, lapply(spec$others, function(oc) {
    rr <- fmrihrf::regressor(
      onsets = oc$onsets,
      hrf = oc$hrf %||% hrf_obj,
      duration = oc$duration %||% 0,
      amplitude = oc$amplitude %||% 1,
      span = oc$span %||% 40,
      summate = TRUE
    )
    x_eval <- fmrihrf::evaluate(
      rr, times,
      precision = spec$precision %||% 0.1,
      method = spec$method %||% "conv"
    )
    if (inherits(x_eval, "Matrix")) x_eval <- as.matrix(x_eval)
    if (is.matrix(x_eval)) x_eval else matrix(as.numeric(x_eval), ncol = 1L)
  }))
  if (is.null(dim(X_other))) X_other <- matrix(X_other, ncol = 1L)
  X_other
}

# --- Helper: LWU HRF grid selection via matched-filter score ---
#' @keywords internal
#' @noRd
.oasis_pick_hrf_lwu_fast <- function(Y, design_spec, hrf_grid, confounds = NULL, 
                                     block_cols = 4096L) {
  # Build aggregate regressor for all trials with each HRF candidate
  sframe <- design_spec$sframe
  times  <- fmrihrf::samples(sframe, global = TRUE)
  onsets <- design_spec$cond$onsets
  
  scores <- numeric(length(hrf_grid))
  
  for (i in seq_along(hrf_grid)) {
    # Build aggregate regressor with this HRF
    r <- fmrihrf::regressor(
      onsets = onsets,
      hrf = hrf_grid[[i]],
      duration = design_spec$cond$duration %||% 0,
      amplitude = design_spec$cond$amplitude %||% 1,
      span = design_spec$cond$span %||% 40,
      summate = TRUE
    )
    x_eval <- fmrihrf::evaluate(r, times, 
                               precision = design_spec$precision %||% 0.1,
                               method = design_spec$method %||% "conv")
    
    # Handle multi-basis HRFs - sum across basis functions for aggregate
    if (is.matrix(x_eval)) {
      x_raw <- rowSums(x_eval)
    } else {
      x_raw <- as.numeric(x_eval)
    }

    # Candidate selection must use the same common span as the subsequent
    # model. Other conditions without their own HRF inherit this candidate.
    X_other <- .oasis_build_X_other(
      design_spec, hrf_obj = hrf_grid[[i]], times = times
    )
    common <- cbind(confounds, X_other)
    Y_use <- Y
    x <- x_raw
    if (!is.null(common) && ncol(common) > 0L) {
      qr_common <- qr(common)
      if (qr_common$rank > 0L) {
        Q <- qr.Q(qr_common)[, seq_len(qr_common$rank), drop = FALSE]
        Y_use <- Y_use - Q %*% crossprod(Q, Y_use)
        x <- x - Q %*% crossprod(Q, x)
      }
    }

    # Matched-filter score: sum of correlations
    x_energy <- sum(x^2)
    raw_energy <- sum(x_raw^2)
    if (!is.finite(x_energy) ||
        x_energy <= .Machine$double.eps * max(raw_energy, 1)) {
      stop(
        sprintf(
          "HRF-grid candidate %d has zero residual norm after projection onto the common design",
          i
        ),
        call. = FALSE
      )
    }
    x_norm <- matrix(x / sqrt(x_energy), ncol = 1L)
    
    # Process in blocks for memory efficiency
    score <- 0
    for (v_start in seq(1, ncol(Y), by = block_cols)) {
      v_end <- min(v_start + block_cols - 1, ncol(Y))
      Y_block <- Y_use[, v_start:v_end, drop = FALSE]
      
      # Correlation with each voxel
      y_norm <- sqrt(colSums(Y_block^2))
      valid <- is.finite(y_norm) & y_norm > 0
      cors <- numeric(ncol(Y_block))
      cors[valid] <- as.numeric(t(x_norm) %*% Y_block[, valid, drop = FALSE]) /
        y_norm[valid]
      score <- score + sum(abs(cors), na.rm = TRUE)
    }
    scores[i] <- score
  }
  
  # Return HRF with highest score
  best_idx <- which.max(scores)
  hrf_grid[[best_idx]]
}

# --- Helper: detect basis dimension from design matrix ---
#' @keywords internal
#' @noRd
.detect_basis_dimension <- function(X) {
  N <- ncol(X)
  if (N <= 1) return(1L)
  
  # Check for common basis dimensions
  for (K in c(2, 3, 4, 5, 6, 8, 10, 12)) {
    if (N %% K == 0) {
      # Check if columns group into K-sized blocks with high within-block correlation
      n_trials <- N / K
      block_cors <- numeric(n_trials)
      
      for (i in seq_len(n_trials)) {
        block_start <- (i - 1) * K + 1
        block_end <- i * K
        block <- X[, block_start:block_end, drop = FALSE]
        
        if (K > 1) {
          cors <- cor(block)
          diag(cors) <- NA
          block_cors[i] <- mean(abs(cors), na.rm = TRUE)
        }
      }
      
      # If blocks show high internal correlation, likely found K
      if (mean(block_cors) > 0.5) {
        return(as.integer(K))
      }
    }
  }
  
  # Default to 1 if no clear pattern
  return(1L)
}
