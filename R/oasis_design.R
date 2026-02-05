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

  # Aggregates for "other" conditions (one col each) -> nuisances
  X_other <- NULL
  if (length(spec$others)) {
    X_other <- do.call(cbind, lapply(spec$others, function(oc) {
      rr <- fmrihrf::regressor(onsets   = oc$onsets,
                               hrf      = oc$hrf %||% hrf_obj,
                               duration = oc$duration %||% 0,
                               amplitude= oc$amplitude %||% 1,
                               span     = oc$span %||% 40,
                               summate  = TRUE)
      x_eval <- fmrihrf::evaluate(rr, times,
                                  precision = spec$precision %||% 0.1,
                                  method = spec$method %||% "conv")
      if (inherits(x_eval, "Matrix")) x_eval <- as.matrix(x_eval)
      # For multi-basis HRFs, aggregate columns to a single condition column
      if (is.matrix(x_eval)) rowSums(x_eval) else as.numeric(x_eval)
    }))
    if (is.null(dim(X_other))) X_other <- matrix(X_other, ncol = 1)
  }

  return(list(X_trials = X_trials, X_other = X_other, K = K))
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
  
  # Residualize Y against confounds once
  if (!is.null(confounds) && ncol(confounds) > 0) {
    Q <- qr.Q(qr(confounds))
    Y <- Y - Q %*% crossprod(Q, Y)
  }
  
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
      x <- rowSums(x_eval)
    } else {
      x <- as.numeric(x_eval)
    }
    
    # Matched-filter score: sum of correlations
    x_norm <- matrix(x / sqrt(sum(x^2)), ncol = 1)
    
    # Process in blocks for memory efficiency
    score <- 0
    for (v_start in seq(1, ncol(Y), by = block_cols)) {
      v_end <- min(v_start + block_cols - 1, ncol(Y))
      Y_block <- Y[, v_start:v_end, drop = FALSE]
      
      # Correlation with each voxel
      cors <- t(x_norm) %*% Y_block / 
        sqrt(colSums(Y_block^2))
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
