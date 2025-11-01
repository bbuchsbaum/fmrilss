#!/usr/bin/env Rscript

# BIDS-integrated trial-wise beta estimation using fmrilss
#
# Combines BIDS dataset navigation from estimate_betas.R with modern fmrilss
# estimation methods from estimate_lss.R.
#
# Usage:
#   estimate_lss_bids.R --bids_path=/data/myproject --subid=01 --task=stroop \
#                       --method=oasis --prewhiten=ar1 --confounds=confounds.txt

suppressPackageStartupMessages({
  requireNamespace("fmrilss")
  requireNamespace("fmrihrf")
})

# =============================================================================
# Utilities
# =============================================================================

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

as_flag <- function(x) tolower(x) %in% c("1", "true", "t", "yes", "y")

stop_if_missing <- function(x, msg) {
  if (is.null(x) || (is.character(x) && !nzchar(x))) {
    stop(msg, call. = FALSE)
  }
}

# =============================================================================
# Argument Parsing
# =============================================================================

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  kv <- strsplit(args, "=", fixed = TRUE)
  m <- vapply(kv, function(p) if (length(p) == 2) p[1] else NA_character_, character(1))
  v <- vapply(kv, function(p) if (length(p) == 2) p[2] else NA_character_, character(1))
  aa <- stats::setNames(as.list(v), m)

  list(
    # BIDS parameters
    bids_path    = aa$"--bids_path" %||% ".",
    bids_session = aa$"--bids_session" %||% "",
    deriv_folder = aa$"--deriv_folder" %||% "derivatives/fmriprep",
    subid        = aa$"--subid",
    task         = aa$"--task",
    space        = aa$"--space" %||% "MNI152NLin2009cAsym",

    # Confounds
    confounds    = aa$"--confounds",
    percvar      = as.numeric(aa$"--percvar" %||% 95),

    # Design parameters
    tr           = as.numeric(aa$"--tr" %||% NA_real_),
    basis        = aa$"--basis" %||% "spmg1",
    duration     = as.numeric(aa$"--duration" %||% 0),
    onset_col    = aa$"--onset_col" %||% "onset",
    duration_col = aa$"--duration_col" %||% NULL,
    milliseconds = as_flag(aa$"--milliseconds" %||% "false"),

    # Estimation
    method       = aa$"--method" %||% "oasis",
    prew         = aa$"--prewhiten" %||% "none",
    polort       = as.numeric(aa$"--polort" %||% 3),

    # Mask
    mask         = aa$"--mask",

    # Output
    outdir       = aa$"--outdir" %||% "betas",
    out_stem     = aa$"--out" %||% "betas",
    concatenate  = as_flag(aa$"--concatenate" %||% "false"),
    outcsv       = as_flag(aa$"--write-csv" %||% "false"),
    outnii       = as_flag(aa$"--write-nifti" %||% "true")
  )
}

# =============================================================================
# BIDS Discovery
# =============================================================================

check_bidser <- function() {
  if (!requireNamespace("bidser", quietly = TRUE)) {
    stop("bidser package required. Install from GitHub: bbuchsbaum/bidser", call. = FALSE)
  }
}

build_bids_project <- function(bids_path) {
  if (bids_path == ".") bids_path <- getwd()
  message("[fmrilss] Building BIDS project from ", bids_path)
  bidser::bids_project(bids_path, fmriprep = TRUE)
}

discover_scans <- function(bids, bids_path, subid, task) {
  message("[fmrilss] Discovering preprocessed scans...")
  scans <- paste0(bids_path, "/", bidser:::preproc_scans(bids, subid = subid, task = task))

  fexists <- sapply(scans, file.exists)
  if (!all(fexists)) {
    stop("Could not find scan(s): ", paste(scans[!fexists], collapse = " "), call. = FALSE)
  }

  scans
}

discover_events <- function(bids, subid, task) {
  message("[fmrilss] Discovering event files...")
  evfiles <- bidser:::search_files(bids, subid = subid, task = task,
                                   kind = "events", full_path = TRUE)

  if (length(evfiles) == 0) {
    stop("No event files found for subid=", subid, " task=", task, call. = FALSE)
  }

  evfiles
}

load_confound_vars <- function(confounds_file) {
  if (is.null(confounds_file) || !nzchar(confounds_file)) {
    return(NULL)
  }

  if (!file.exists(confounds_file)) {
    stop("Confounds file not found: ", confounds_file, call. = FALSE)
  }

  cvars <- scan(confounds_file, what = "", quiet = TRUE)
  message("[fmrilss] Confound variables: ", paste(cvars, collapse = " "))
  cvars
}

extract_confounds <- function(bids, cvars, percvar, subid, task) {
  if (is.null(cvars)) return(NULL)

  message("[fmrilss] Extracting confounds with ", percvar, "% variance retention...")
  cdat <- bidser:::read_confounds(bids, perc_var = percvar, cvars = cvars,
                                  subid = subid, task = task, nest = TRUE)

  # Convert to list of matrices per run
  lapply(cdat$data, function(df) {
    mat <- as.matrix(df)
    # Remove all-NA columns
    keep <- !apply(mat, 2, function(z) all(is.na(z)))
    mat[, keep, drop = FALSE]
  })
}

construct_mask_path <- function(bids_path, deriv_folder, subid, session, mask_name) {
  if (is.null(mask_name) || !nzchar(mask_name)) {
    return(NULL)
  }

  if (session == "") {
    mask_path <- file.path(bids_path, deriv_folder,
                          paste0("sub-", subid), "func",
                          paste0("sub-", subid, "_", mask_name))
  } else {
    # Try with session in filename first
    mask_path <- file.path(bids_path, deriv_folder,
                          paste0("sub-", subid),
                          paste0("ses-", session), "func",
                          paste0("sub-", subid, "_ses-", session, "_", mask_name))

    # Fallback without session in filename
    if (!file.exists(mask_path)) {
      mask_path <- file.path(bids_path, deriv_folder,
                            paste0("sub-", subid),
                            paste0("ses-", session), "func",
                            paste0("sub-", subid, "_", mask_name))
    }
  }

  if (!file.exists(mask_path)) {
    stop("Mask not found: ", mask_path, call. = FALSE)
  }

  message("[fmrilss] Loading mask: ", mask_path)
  mask_path
}

# =============================================================================
# Event Processing
# =============================================================================

read_events_tsv <- function(path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
               na.strings = c("n/a", "NA", "N/A", ""))
  } else {
    data.table::fread(path, data.table = FALSE,
                     na.strings = c("n/a", "NA", "N/A", ""))
  }
}

standardize_events <- function(dat, onset_name = "onset", duration_name = NULL,
                               milliseconds = FALSE, default_duration = 0) {
  nms <- names(dat)
  low <- tolower(nms)

  # Find onset column
  onset_idx <- which(low == tolower(onset_name))
  if (length(onset_idx) == 0) {
    stop(sprintf("Onset column '%s' not found (available: %s)",
                 onset_name, paste(nms, collapse = ",")), call. = FALSE)
  }

  onsets <- suppressWarnings(as.numeric(dat[[nms[onset_idx[1]]]]))
  if (any(!is.finite(onsets))) {
    stop("Onsets must be numeric", call. = FALSE)
  }

  # Convert milliseconds if needed
  if (milliseconds) onsets <- onsets / 1000

  # Get durations
  durations <- default_duration
  if (!is.null(duration_name)) {
    dur_idx <- which(low == tolower(duration_name))
    if (length(dur_idx) > 0) {
      durations <- suppressWarnings(as.numeric(dat[[nms[dur_idx[1]]]]))
      if (any(!is.finite(durations))) {
        stop("Durations must be numeric", call. = FALSE)
      }
    }
  }

  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }

  out <- data.frame(onset = onsets, duration = durations, stringsAsFactors = FALSE)

  # Preserve run column if present
  run_idx <- which(low == "run")
  out$run <- if (length(run_idx)) dat[[nms[run_idx[1]]]] else 1L

  out
}

# =============================================================================
# Design Matrix Construction
# =============================================================================

build_design <- function(events, TR, Tlen, basis = c("spmg1", "spmg3")) {
  basis <- match.arg(basis)
  sframe <- fmrihrf::sampling_frame(blocklens = Tlen, TR = TR)
  onsets <- events$onset
  durations <- events$duration

  # Prefer fmridesign if available
  if (requireNamespace("fmridesign", quietly = TRUE)) {
    hrf <- fmridesign::hrf
    des <- data.frame(onset = onsets, run = events$run %||% 1L,
                     trial = factor(seq_along(onsets)))
    emod <- fmridesign::event_model(onset ~ hrf(trial, basis = basis),
                                    data = des, block = ~run,
                                    sampling_frame = sframe,
                                    durations = durations, precision = 0.1)
    X <- fmridesign::design_matrix(emod)
  } else {
    # Fallback to fmrihrf
    Rset <- fmrihrf::regressor_set(
      onsets = onsets,
      fac = factor(seq_along(onsets)),
      hrf = switch(basis, spmg1 = fmrihrf::HRF_SPMG1, fmrihrf::HRF_SPMG3),
      duration = durations, span = 30, summate = FALSE
    )
    X <- fmrihrf::evaluate(Rset, grid = fmrihrf::samples(sframe, global = TRUE),
                          precision = 0.1, method = "conv")
  }

  if (!is.matrix(X)) X <- as.matrix(X)
  list(X = X, sframe = sframe)
}

build_Z_matrix <- function(Tlen, polort) {
  if (polort > 0) {
    t_norm <- seq(0, 1, length.out = Tlen)
    Z <- cbind(Intercept = 1, sapply(1:polort, function(p) t_norm^p))
    colnames(Z)[-1] <- paste0("poly", 1:polort)
  } else {
    Z <- cbind(Intercept = 1)
  }
  Z
}

# =============================================================================
# NIfTI I/O
# =============================================================================

load_nifti <- function(path) {
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("Please install RNifti for NIfTI I/O", call. = FALSE)
  }
  RNifti::readNifti(path)
}

to_matrix <- function(img, mask = NULL) {
  dim4 <- dim(img)
  stopifnot(length(dim4) == 4)
  Tlen <- dim4[4]

  if (is.null(mask)) {
    # Variance-based mask
    volm <- array(img, dim = dim4)
    Vmat <- matrix(volm, nrow = prod(dim4[1:3]), ncol = Tlen)
    keep <- apply(Vmat, 1L, function(x) sd(x) > 0)
  } else {
    m <- as.array(mask)
    keep <- as.vector(m != 0)
  }

  vox <- which(keep)
  Y <- matrix(0, nrow = Tlen, ncol = length(vox))
  for (t in seq_len(Tlen)) {
    Y[t, ] <- as.vector(img[, , , t])[keep]
  }

  list(Y = Y, idx = vox, mask_keep = keep, dim3 = dim4[1:3])
}

write_betas_nifti <- function(beta, dim3, keep, template, out_path) {
  if (!requireNamespace("RNifti", quietly = TRUE)) {
    stop("Please install RNifti for NIfTI I/O", call. = FALSE)
  }

  ntrials <- nrow(beta)
  out <- array(0, dim = c(dim3, ntrials))

  for (j in seq_len(ntrials)) {
    vol <- numeric(prod(dim3))
    vol[keep] <- beta[j, ]
    out[, , , j] <- array(vol, dim = dim3)
  }

  img <- RNifti::updateNifti(out, template)
  RNifti::writeNifti(img, out_path)
}

# =============================================================================
# LSS Estimation
# =============================================================================

configure_prewhitening <- function(prew_method) {
  if (prew_method == "ar1") {
    list(method = "ar", p = 1L)
  } else if (prew_method == "auto") {
    list(method = "ar", p = "auto", p_max = 4)
  } else {
    NULL
  }
}

estimate_lss <- function(Y, X, Z, Nuisance, method, prewhiten) {
  if (method %in% c("r", "r_optimized", "cpp_optimized", "naive")) {
    fmrilss::lss(Y, X, Z = Z, Nuisance = Nuisance,
                method = method, prewhiten = prewhiten)
  } else if (method == "oasis") {
    fmrilss::lss(Y = Y, X = X, Z = Z, Nuisance = Nuisance,
                method = "oasis",
                oasis = list(ridge_mode = "fractional",
                           ridge_x = 0.02, ridge_b = 0.02),
                prewhiten = prewhiten)
  } else {
    stop("Unknown method: ", method, " (use r, cpp_optimized, or oasis)", call. = FALSE)
  }
}

# =============================================================================
# Run Processing
# =============================================================================

process_run <- function(run_idx, scan_path, event_path, args, mask_img, nuisance_run) {
  message(sprintf("[fmrilss] Processing run %d: %s", run_idx, basename(scan_path)))

  # Load BOLD
  bold_img <- load_nifti(scan_path)
  Tlen <- dim(bold_img)[4]

  # Load and standardize events
  ev_raw <- read_events_tsv(event_path)
  ev <- standardize_events(ev_raw,
                          onset_name = args$onset_col,
                          duration_name = args$duration_col,
                          milliseconds = args$milliseconds,
                          default_duration = args$duration)

  message(sprintf("  - %d events, %d timepoints", nrow(ev), Tlen))

  # Build design matrix
  des <- build_design(ev, TR = args$tr, Tlen = Tlen, basis = args$basis)
  X <- des$X

  if (nrow(X) != Tlen) {
    stop(sprintf("Design rows (%d) != BOLD timepoints (%d)", nrow(X), Tlen),
         call. = FALSE)
  }

  # Extract Y matrix
  mv <- to_matrix(bold_img, mask_img)
  Y <- mv$Y
  message(sprintf("  - %d voxels in mask", ncol(Y)))

  # Nuisance regressors
  if (!is.null(nuisance_run)) {
    message(sprintf("  - %d nuisance regressors", ncol(nuisance_run)))
  }

  # Build Z matrix (trends)
  Z <- build_Z_matrix(Tlen, args$polort)

  # Configure prewhitening
  prew <- configure_prewhitening(args$prew)

  # Estimate
  message("  - Estimating with method=", args$method)
  beta <- estimate_lss(Y, X, Z, nuisance_run, args$method, prew)
  message(sprintf("  - Estimated %d trial betas", nrow(beta)))

  list(beta = beta, mv = mv, template = bold_img)
}

# =============================================================================
# Output Writing
# =============================================================================

write_csv_output <- function(beta, out_path) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    write.csv(beta, out_path, row.names = FALSE)
  } else {
    dt <- data.table::as.data.table(beta)
    dt[, trial := .I]
    data.table::fwrite(dt, out_path)
  }
  message("[fmrilss] Wrote ", out_path)
}

write_outputs <- function(beta_list, args, outpath, outfiles) {
  if (isTRUE(args$concatenate)) {
    # Concatenate all runs
    beta_all <- do.call(rbind, lapply(beta_list, function(x) x$beta))
    message(sprintf("[fmrilss] Concatenating %d total trials", nrow(beta_all)))

    if (isTRUE(args$outnii)) {
      out_nii <- file.path(outpath, paste0(args$out_stem, "_all.nii.gz"))
      write_betas_nifti(beta_all, beta_list[[1]]$mv$dim3,
                       beta_list[[1]]$mv$mask_keep,
                       beta_list[[1]]$template, out_nii)
      message("[fmrilss] Wrote ", out_nii)
    }

    if (isTRUE(args$outcsv)) {
      out_csv <- file.path(outpath, paste0(args$out_stem, "_all.csv"))
      write_csv_output(beta_all, out_csv)
    }

  } else {
    # Write per-run outputs
    for (run in seq_along(beta_list)) {
      bl <- beta_list[[run]]

      if (isTRUE(args$outnii)) {
        out_nii <- file.path(outpath, outfiles[run])
        write_betas_nifti(bl$beta, bl$mv$dim3, bl$mv$mask_keep,
                         bl$template, out_nii)
        message("[fmrilss] Wrote ", out_nii)
      }

      if (isTRUE(args$outcsv)) {
        out_csv <- file.path(outpath, sub("\\.nii(\\.gz)?$", ".csv", outfiles[run]))
        write_csv_output(bl$beta, out_csv)
      }
    }
  }
}

# =============================================================================
# Main Pipeline
# =============================================================================

main <- function() {
  args <- parse_args()

  # Validate required arguments
  stop_if_missing(args$subid, "--subid=<subject_id> is required")
  stop_if_missing(args$task, "--task=<task_name> is required")
  if (is.na(args$tr)) stop("--tr=<seconds> is required", call. = FALSE)

  # Check dependencies
  check_bidser()

  # Build BIDS project
  bids <- build_bids_project(args$bids_path)

  # Handle session
  session <- if (nzchar(args$bids_session)) {
    message("[fmrilss] Using session: ", args$bids_session)
    args$bids_session
  } else {
    ""
  }

  # Discover scans and events
  scans <- discover_scans(bids, args$bids_path, args$subid, args$task)
  evfiles <- discover_events(bids, args$subid, args$task)

  if (length(scans) != length(evfiles)) {
    stop(sprintf("Scan count (%d) != event file count (%d)",
                 length(scans), length(evfiles)), call. = FALSE)
  }

  message("[fmrilss] Found ", length(scans), " run(s)")

  # Load confounds
  cvars <- load_confound_vars(args$confounds)
  nuisanceCovars <- extract_confounds(bids, cvars, args$percvar, args$subid, args$task)

  # Set up output directory
  outpath <- file.path(dirname(dirname(scans[1])), args$outdir)
  if (!dir.exists(outpath)) {
    dir.create(outpath, recursive = TRUE)
  }
  message("[fmrilss] Output directory: ", outpath)

  # Load mask
  mask_path <- construct_mask_path(args$bids_path, args$deriv_folder,
                                   args$subid, session, args$mask)
  mask_img <- if (!is.null(mask_path)) load_nifti(mask_path) else NULL

  # Generate output filenames
  outfiles <- basename(gsub("desc-preproc", "desc-betas", scans))

  # Process each run
  beta_list <- lapply(seq_along(scans), function(run) {
    nuisance_run <- if (!is.null(nuisanceCovars)) nuisanceCovars[[run]] else NULL
    process_run(run, scans[run], evfiles[run], args, mask_img, nuisance_run)
  })

  # Write outputs
  write_outputs(beta_list, args, outpath, outfiles)

  message("[fmrilss] Complete!")
}

# =============================================================================
# Entry Point
# =============================================================================

if (identical(environment(), globalenv())) {
  tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1)
  })
}
