# Run geometry for glmsingle(): trials, runs, conditions, sessions and folds.

# Split Y into a list of time x voxel matrices, one per run.
.glms_split_runs <- function(Y, runs) {
  if (is.list(Y) && !is.data.frame(Y)) {
    Y <- lapply(Y, .as_base_matrix)
    nv <- vapply(Y, ncol, integer(1))
    if (length(unique(nv)) != 1L) stop("All runs of Y must have the same number of voxels", call. = FALSE)
    return(Y)
  }
  Y <- .as_base_matrix(Y)
  if (is.null(runs)) {
    stop("runs must be supplied when Y is a single matrix (one run id per row)", call. = FALSE)
  }
  runs <- .as_integer_ids(runs, "runs")
  if (length(runs) != nrow(Y)) stop("runs must have one entry per row of Y", call. = FALSE)
  if (is.unsorted(match(runs, unique(runs))) ) {
    stop("rows of Y must be grouped by run (runs must be contiguous)", call. = FALSE)
  }
  lapply(unique(runs), function(r) Y[runs == r, , drop = FALSE])
}

# Normalise the design into per-run onset indices (0-based TRs) and condition
# ids. Accepts GLMsingle-style 0/1 matrices (time x conditions, one per run)
# or an events data frame with columns run, onset (seconds) and condition.
.glms_parse_design <- function(design, n_time, tr, run_ids = NULL) {
  R <- length(n_time)
  if (is.data.frame(design)) {
    need <- c("run", "onset", "condition")
    if (!all(need %in% names(design))) {
      stop("design data frame needs columns: run, onset, condition", call. = FALSE)
    }
    if (anyNA(design$run)) stop("design run IDs must not be missing", call. = FALSE)
    if (is.null(run_ids)) {
      run_ids <- sort(unique(design$run))
      if (length(run_ids) != R) stop("design must describe the same number of runs as Y", call. = FALSE)
    }
    run_code <- match(design$run, run_ids)
    if (anyNA(run_code)) stop("design run IDs must match the run IDs of Y", call. = FALSE)
    on_tr <- design$onset / tr
    if (any(abs(on_tr - round(on_tr)) > 1e-6)) {
      stop("glmsingle() requires onsets on the TR grid (onset / tr must be an integer)", call. = FALSE)
    }
    cond <- factor(design$condition)
    levels_out <- levels(cond)
    onsets <- conds <- vector("list", R)
    for (r in seq_len(R)) {
      sel <- which(run_code == r)
      o <- order(on_tr[sel])
      onsets[[r]] <- as.integer(round(on_tr[sel][o]))
      conds[[r]] <- as.integer(cond[sel][o])
    }
  } else {
    if (is.matrix(design)) design <- list(design)
    if (!is.list(design) || length(design) != R) {
      stop("design must be a list with one time x condition matrix per run", call. = FALSE)
    }
    design <- lapply(design, as.matrix)
    nc <- vapply(design, ncol, integer(1))
    if (length(unique(nc)) != 1L) stop("All design matrices need the same number of conditions", call. = FALSE)
    if (any(!unlist(design) %in% c(0, 1))) stop("design matrices must contain only 0 and 1", call. = FALSE)
    levels_out <- colnames(design[[1]]) %||% paste0("cond", seq_len(nc[1]))
    onsets <- conds <- vector("list", R)
    for (r in seq_len(R)) {
      D <- design[[r]]
      if (nrow(D) != n_time[r]) {
        stop(sprintf("design run %d has %d rows but Y run %d has %d time points",
                     r, nrow(D), r, n_time[r]), call. = FALSE)
      }
      hits <- which(D == 1, arr.ind = TRUE)
      if (anyDuplicated(hits[, 1])) {
        stop("two conditions have exactly the same trial onset; this is not allowed", call. = FALSE)
      }
      o <- order(hits[, 1])
      onsets[[r]] <- as.integer(hits[o, 1] - 1L)
      conds[[r]] <- as.integer(hits[o, 2])
    }
  }
  for (r in seq_len(R)) {
    if (any(onsets[[r]] < 0L | onsets[[r]] >= n_time[r])) {
      stop(sprintf("run %d has onsets outside the run", r), call. = FALSE)
    }
    if (anyDuplicated(onsets[[r]])) {
      stop("two trials have exactly the same onset; this is not allowed", call. = FALSE)
    }
  }
  list(onsets = onsets, conds = conds, levels = levels_out)
}

# Build the trial/run/condition geometry used by every stage.
.glms_geometry <- function(parsed, n_time, tr, session_indicator, xval_scheme) {
  R <- length(n_time)
  n_per_run <- lengths(parsed$onsets)
  ends <- cumsum(n_per_run)
  starts <- ends - n_per_run + 1L
  validcolumns <- lapply(seq_len(R), function(r) {
    if (n_per_run[r]) seq.int(starts[r], ends[r]) else integer(0)
  })
  stimorder <- unlist(parsed$conds, use.names = FALSE)
  session_indicator <- if (is.null(session_indicator)) rep(1L, R) else {
    s <- .as_integer_ids(session_indicator, "session_indicator")
    if (length(s) != R) stop("session_indicator needs one entry per run", call. = FALSE)
    s
  }
  if (is.null(xval_scheme)) {
    xval_scheme <- as.list(seq_len(R))
  } else {
    if (!is.list(xval_scheme)) xval_scheme <- as.list(xval_scheme)
    xval_scheme <- lapply(xval_scheme, function(f) {
      f <- as.integer(f)
      if (!length(f) || any(f < 1L | f > R)) stop("xval_scheme entries must be run indices in 1..number of runs", call. = FALSE)
      f
    })
  }
  cond_in_runs <- vapply(seq_along(parsed$levels), function(c) {
    sum(vapply(parsed$conds, function(x) c %in% x, logical(1)))
  }, integer(1))
  list(
    n_runs = R, n_time = n_time, tr = tr,
    onsets = parsed$onsets, conds = parsed$conds, levels = parsed$levels,
    n_trials = sum(n_per_run), n_per_run = n_per_run,
    validcolumns = validcolumns, stimorder = stimorder,
    trial_run = rep(seq_len(R), n_per_run),
    session = session_indicator, xval_scheme = xval_scheme,
    cond_in_runs = cond_in_runs,
    cond_trials = split(seq_along(stimorder) - 1L, stimorder)
  )
}
