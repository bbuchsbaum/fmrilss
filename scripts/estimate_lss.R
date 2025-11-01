#!/usr/bin/env Rscript

# Estimate trial-wise betas using fmrilss (LSS/OASIS), replacing legacy AFNI usage.
# - Reads a 4D BOLD NIfTI and an events.tsv (BIDS-like) file
# - Builds a per-trial HRF-convolved design
# - Runs fmrilss::lss with selected backend (R, cpp_optimized, or oasis)
# - Writes per-trial beta volumes as a 4D NIfTI (X×Y×Z×ntrials) and optional CSV
#
# Notes:
# - AFNI support has been removed from this script by design. Keep AFNI in a separate script if needed.
# - Optional prewhitening is available via fmriAR.
# - If fmridesign is installed, uses event_model to construct trial-wise columns; otherwise falls back to fmrihrf.

suppressPackageStartupMessages({
  requireNamespace("fmrilss")
  requireNamespace("fmrihrf")
})

as_flag <- function(x) tolower(x) %in% c("1","true","t","yes","y")

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- strsplit(args, "=", fixed = TRUE)
  m <- vapply(kv, function(p) if (length(p)==2) p[1] else NA_character_, character(1))
  v <- vapply(kv, function(p) if (length(p)==2) p[2] else NA_character_, character(1))
  aa <- stats::setNames(as.list(v), m)
  list(
    bold    = aa$"--bold",
    mask    = aa$"--mask",
    events  = aa$"--events",
    tr      = as.numeric(aa$"--tr" %||% NA_real_),
    basis   = aa$"--basis" %||% "spmg1",
    method  = aa$"--method" %||% "oasis",         # "r", "cpp_optimized", or "oasis"
    prew    = aa$"--prewhiten" %||% "none",        # "none", "ar1", "auto"
    onset_col    = aa$"--onset_col" %||% NULL,      # optional custom onset column name
    duration_col = aa$"--duration_col" %||% NULL,   # optional custom duration column name
    outdir  = aa$"--outdir" %||% "lss_out",
    outcsv  = as_flag(aa$"--write-csv" %||% "true"),
    outnii  = as_flag(aa$"--write-nifti" %||% "true")
  )
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

stop_if_missing <- function(x, msg) if (is.null(x) || !nzchar(x)) stop(msg, call. = FALSE)

read_events_tsv <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    dat <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    dat <- data.table::fread(path, data.table = FALSE)
  }
  dat
}

# Standardize events: choose onset and duration columns; duration optional (defaults 0)
standardize_events <- function(dat, onset_name = NULL, duration_name = NULL) {
  nms <- names(dat)
  low <- tolower(nms)
  pick_col <- function(target, fallback_patterns) {
    if (!is.null(target)) {
      # case-insensitive match against provided target name
      i <- which(low == tolower(target))
      if (length(i)) return(nms[i[1]])
      stop(sprintf("Column '%s' not found in events (available: %s)", target, paste(nms, collapse=",")))
    }
    for (pat in fallback_patterns) {
      i <- which(low == pat)
      if (length(i)) return(nms[i[1]])
    }
    NULL
  }
  onset_col <- pick_col(onset_name, c("onset","onsets","start"))
  if (is.null(onset_col)) stop("Could not detect onset column. Use --onset_col=<name>.")
  dur_col <- pick_col(duration_name, c("duration","dur","durations"))

  onsets <- suppressWarnings(as.numeric(dat[[onset_col]]))
  if (any(!is.finite(onsets))) stop("Onsets must be numeric (seconds)")
  if (!is.null(dur_col)) {
    durations <- suppressWarnings(as.numeric(dat[[dur_col]]))
    if (any(!is.finite(durations))) stop("Durations must be numeric (seconds)")
  } else {
    durations <- rep(0, length(onsets))
  }
  out <- data.frame(onset = onsets, duration = durations, stringsAsFactors = FALSE)
  # Preserve run column if present (case-insensitive "run")
  run_idx <- which(low == "run")
  out$run <- if (length(run_idx)) dat[[nms[run_idx[1]]]] else 1L
  out
}

load_nifti <- function(path) {
  if (!requireNamespace("RNifti", quietly = TRUE)) stop("Please install RNifti for NIfTI I/O")
  RNifti::readNifti(path)
}

to_matrix <- function(img, mask = NULL) {
  dim4 <- dim(img)
  stopifnot(length(dim4) == 4)
  Tlen <- dim4[4]
  if (is.null(mask)) {
    # simple variance-based mask
    volm <- array(img, dim = dim4)
    Vmat <- matrix(volm, nrow = prod(dim4[1:3]), ncol = Tlen)
    keep <- apply(Vmat, 1L, function(x) sd(x) > 0)
  } else {
    m <- as.array(mask)
    keep <- as.vector(m != 0)
  }
  vox <- which(keep)
  Y <- matrix(0, nrow = Tlen, ncol = length(vox))
  for (t in seq_len(Tlen)) Y[t, ] <- as.vector(img[,,,t])[keep]
  list(Y = Y, idx = vox, mask_keep = keep, dim3 = dim4[1:3])
}

write_betas_nifti <- function(beta, dim3, keep, template, out_path) {
  if (!requireNamespace("RNifti", quietly = TRUE)) stop("Please install RNifti for NIfTI I/O")
  ntrials <- nrow(beta)
  out <- array(0, dim = c(dim3, ntrials))
  for (j in seq_len(ntrials)) {
    vol <- numeric(prod(dim3))
    vol[keep] <- beta[j, ]
    out[,,,j] <- array(vol, dim = dim3)
  }
  img <- RNifti::updateNifti(out, template)
  RNifti::writeNifti(img, out_path)
}

build_design <- function(events, TR, Tlen, basis = c("spmg1","spmg3")) {
  basis <- match.arg(basis)
  sframe <- fmrihrf::sampling_frame(blocklens = Tlen, TR = TR)
  onsets <- events$onset
  durations <- events$duration %||% 0
  # Prefer fmridesign if present for robust hrf() interface
  if (requireNamespace("fmridesign", quietly = TRUE)) {
    hrf <- fmridesign::hrf
    des <- data.frame(onset = onsets, run = events$run %||% 1L, trial = factor(seq_along(onsets)))
    emod <- fmridesign::event_model(onset ~ hrf(trial, basis = basis),
                                    data = des, block = ~run, sampling_frame = sframe,
                                    durations = durations, precision = 0.1)
    X <- fmridesign::design_matrix(emod)
    if (!is.matrix(X)) X <- as.matrix(X)
    list(X = X, sframe = sframe)
  } else {
    # Fallback: fmrihrf regressor_set builds one column per trial
    Rset <- fmrihrf::regressor_set(onsets = onsets,
                                   fac    = factor(seq_along(onsets)),
                                   hrf    = switch(basis, spmg1 = fmrihrf::HRF_SPMG1, fmrihrf::HRF_SPMG3),
                                   duration = durations, span = 30, summate = FALSE)
    X <- fmrihrf::evaluate(Rset, grid = fmrihrf::samples(sframe, global = TRUE), precision = 0.1, method = "conv")
    if (!is.matrix(X)) X <- as.matrix(X)
    list(X = X, sframe = sframe)
  }
}

main <- function() {
  a <- parse_args()
  stop_if_missing(a$bold,   "--bold=/path/to/bold.nii.gz is required")
  stop_if_missing(a$events, "--events=/path/to/events.tsv is required")
  if (is.na(a$tr)) stop("--tr=<seconds> is required", call. = FALSE)

  message("[fmrilss] Loading BOLD and events...")
  bold_img <- load_nifti(a$bold)
  mask_img <- if (!is.null(a$mask) && nzchar(a$mask)) load_nifti(a$mask) else NULL
  Tlen <- dim(bold_img)[4]
  ev_raw <- read_events_tsv(a$events)
  ev <- standardize_events(ev_raw, onset_name = a$onset_col, duration_name = a$duration_col)

  # Build design
  des <- build_design(ev, TR = a$tr, Tlen = Tlen, basis = a$basis)
  X <- des$X
  if (nrow(X) != Tlen) stop(sprintf("Design rows (%d) != BOLD timepoints (%d)", nrow(X), Tlen))

  # Extract time×voxel matrix Y
  mv <- to_matrix(bold_img, mask_img)
  Y <- mv$Y

  # Prewhitening options
  prew <- NULL
  if (a$prew == "ar1") prew <- list(method = "ar", p = 1L)
  if (a$prew == "auto") prew <- list(method = "ar", p = "auto", p_max = 4)

  # Run fmrilss
  message("[fmrilss] Estimating trial-wise betas with method=", a$method)
  if (a$method %in% c("r","cpp_optimized")) {
    beta <- fmrilss::lss(Y, X, Z = cbind(Intercept = 1), Nuisance = NULL,
                         method = a$method, prewhiten = prew)
  } else if (a$method == "oasis") {
    # Let OASIS build its own design_spec to match X (optional)
    beta <- fmrilss::lss(Y = Y, X = X, Z = cbind(Intercept = 1),
                         method = "oasis",
                         oasis = list(ridge_mode = "fractional", ridge_x = 0.02, ridge_b = 0.02),
                         prewhiten = prew)
  } else {
    stop("Unknown --method; use r, cpp_optimized, or oasis")
  }

  # Prepare outputs
  if (!dir.exists(a$outdir)) dir.create(a$outdir, recursive = TRUE)
  ntrials <- nrow(beta); V <- ncol(beta)
  message(sprintf("[fmrilss] Done. Trials=%d Voxels=%d", ntrials, V))

  # Write CSV (long) for quick diagnostics
  if (isTRUE(a$outcsv)) {
    out_csv <- file.path(a$outdir, sprintf("betas_%s.csv", a$method))
    if (!requireNamespace("data.table", quietly = TRUE)) {
      write.csv(beta, out_csv, row.names = FALSE)
    } else {
      dt <- data.table::as.data.table(beta)
      dt[, trial := .I]
      data.table::fwrite(dt, out_csv)
    }
    message("[fmrilss] Wrote ", out_csv)
  }

  # Write 4D NIfTI (X×Y×Z×ntrials) in mask
  if (isTRUE(a$outnii)) {
    out_nii <- file.path(a$outdir, sprintf("betas_%s.nii.gz", a$method))
    write_betas_nifti(beta, mv$dim3, mv$mask_keep, bold_img, out_nii)
    message("[fmrilss] Wrote ", out_nii)
  }
}

if (identical(environment(), globalenv())) {
  tryCatch(main(), error = function(e) {message("ERROR: ", conditionMessage(e)); quit(status = 1)})
}
