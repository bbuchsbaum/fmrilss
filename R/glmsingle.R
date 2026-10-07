#' GLMsingle single-trial response estimation
#'
#' Estimates single-trial response amplitudes with the GLMsingle procedure
#' (Prince et al., 2022): a library-based HRF per voxel, GLMdenoise noise
#' regressors chosen by cross-validation, and voxel-wise fractional ridge
#' regression chosen by cross-validation over repeated conditions. The
#' implementation reorganises the computation (per-run factorisation, a
#' single pass of data products per stage, compiled cross-validation and
#' sufficient-statistic R^2) so it is much faster than the reference
#' implementation while computing the same estimator.
#'
#' Four models are returned, matching GLMsingle:
#' \describe{
#'   \item{typea}{ON-OFF model: one canonical-HRF regressor for all trials.}
#'   \item{typeb}{Single-trial OLS with the best library HRF per voxel.}
#'   \item{typec}{Type B plus GLMdenoise noise regressors.}
#'   \item{typed}{Type C plus voxel-wise fractional ridge regression.}
#' }
#'
#' @section Argument names: Arguments use snake_case; the GLMsingle names are
#'   `wantlibrary` (`want_library`), `wantglmdenoise` (`want_glmdenoise`),
#'   `wantfracridge` (`want_fracridge`), `xvalscheme` (`xval_scheme`),
#'   `sessionindicator` (`session_indicator`), `maxpolydeg` (`max_poly_deg`),
#'   `brainthresh` (`brain_thresh`), `brainR2` (`brain_r2`), `brainexclude`
#'   (`brain_exclude`), `pcR2cutoff` (`pc_r2_cutoff`), `pcR2cutoffmask`
#'   (`pc_r2_cutoff_mask`), `wantpercentbold` (`want_percent_bold`),
#'   `wantautoscale` (`want_autoscale`) and `chunklen` (`chunk_size`). Output
#'   fields keep GLMsingle's names.
#'
#' @section Differences from GLMsingle: Defaults reproduce GLMsingle (Python,
#'   commit 1ab54a6) except where GLMsingle is internally inconsistent:
#'   \itemize{
#'     \item `extras_in_denoise = "always"` keeps user extra regressors in every
#'       GLMdenoise and ridge fit; GLMsingle (Python) drops them when zero noise
#'       PCs are used, including from the cross-validation reference.
#'       `"with_pcs"` reproduces GLMsingle.
#'     \item `zero_sd_cv = "zero"` gives voxels with zero beta variance no
#'       weight in cross-validation; GLMsingle's in-place division helper
#'       treats them inconsistently across candidates. `"python"` reproduces it.
#'     \item The ON-OFF R^2 threshold uses a deterministic two-Gaussian mixture
#'       fit with a small variance floor (as in GLMsingle's MATLAB code);
#'       GLMsingle's Python fit is unseeded.
#'     \item `want_percent_bold` applies to all model types (GLMsingle always
#'       scales types C and D).
#'     \item Computation is in double precision (GLMsingle uses single).
#'     \item HRF indices are 1-based; betas are trials x voxels.
#'   }
#'   The FIR diagnostic model, `hrfmodel = "optimize"`, `wantlss`, bootstrap
#'   modes, figures and file outputs are not implemented.
#'
#' @param Y Data: a list of time x voxel matrices (one per run), or one time
#'   x voxel matrix with `runs` giving the run of each row.
#' @param design Either a list of time x condition 0/1 matrices, one per run
#'   (GLMsingle's format; a 1 marks a trial onset), or a data frame with
#'   columns `run`, `onset` (seconds; must be on the TR grid) and
#'   `condition`. Repeated conditions drive the cross-validation.
#'   For matrix `Y`, event run IDs are matched to `runs` in data order.
#'   For list `Y`, ascending event run IDs correspond to the list order.
#' @param tr Repetition time in seconds.
#' @param stimdur Trial duration in seconds.
#' @param runs Run identifier per row of `Y` when `Y` is a single matrix.
#' @param hrf_library Optional time x HRF matrix sampled at the TR (first row
#'   at onset). Default: [glmsingle_hrf_library()].
#' @param want_library Fit the HRF library per voxel (otherwise the canonical
#'   HRF is used for every voxel).
#' @param want_glmdenoise Derive and select GLMdenoise noise regressors.
#' @param want_fracridge Fit the fractional ridge model (type D).
#' @param fracs Ridge fractions in (0, 1] to evaluate. A single value skips
#'   cross-validation and uses that fraction.
#' @param n_pcs Maximum number of noise PCs to evaluate. Capped, with a
#'   warning, at the available noise-pool rank across runs. Empty or
#'   rank-zero pools use zero PCs and skip PC-count cross-validation.
#' @param pcstop Stopping factor for choosing the number of PCs. A value
#'   `<= 0` uses `-pcstop` PCs without cross-validation.
#' @param xval_scheme List of integer vectors of runs held out together in
#'   each cross-validation fold. Default: leave one run out.
#' @param session_indicator Integer session per run, used to z-score betas
#'   within session for cross-validation. Default: one session.
#' @param extra_regressors Optional list (one per run) of time x regressor
#'   nuisance matrices (e.g. motion), or `NULL`.
#' @param max_poly_deg Polynomial drift degree, scalar or one per run.
#'   Default: `round(run_seconds / 120)` as in GLMsingle.
#' @param brain_thresh Two numbers: a percentile of the mean volume and a
#'   fraction of it; voxels brighter than their product may enter the noise
#'   pool.
#' @param brain_r2 ON-OFF R^2 (percent) below which bright voxels enter the
#'   noise pool. Default: estimated tail threshold. With fewer than two
#'   distinct finite ON-OFF R^2 values, uses the common value (or zero if
#'   none are finite). Thresholds are only estimated when denoising is used.
#' @param brain_exclude Optional logical vector of voxels to exclude from the
#'   noise pool (`FALSE` excludes).
#' @param pc_r2_cutoff ON-OFF R^2 above which voxels summarise the PC-count
#'   cross-validation. Default: estimated tail threshold.
#' @param pc_r2_cutoff_mask Optional logical vector restricting those voxels.
#' @param want_percent_bold Express betas as percent signal change of the mean
#'   volume.
#' @param want_autoscale Rescale ridge betas to best match the unregularised
#'   betas (type D).
#' @param extras_in_denoise How user extra regressors enter the GLMdenoise and
#'   ridge fits: `"always"` (default, recommended) or `"with_pcs"` (only when
#'   at least one PC is used; GLMsingle's behaviour).
#' @param zero_sd_cv Treatment of voxels whose reference betas have zero
#'   variance within a session during cross-validation: `"zero"` (default) or
#'   `"python"` (reproduce GLMsingle's Python behaviour).
#' @param singular What to do when a run's trial regressors are linearly
#'   dependent: `"error"` (default, as GLMsingle) or `"pinv"`
#'   (minimum-norm solution with a warning).
#' @param frac_alpha How fractions map to ridge penalties: `"fracridge"`
#'   (default; GLMsingle's grid interpolation) or `"exact"` (solve the
#'   fraction equation exactly; experimental).
#' @param full_glmbadness Compute the PC-count cross-validation for every
#'   voxel instead of only the voxels that decide the PC count. Estimates are
#'   identical; only the `glmbadness` diagnostic is filled for all voxels.
#' @param chunk_size Number of voxels processed at a time (GLMsingle's
#'   `chunklen`). Lower it to reduce peak memory.
#' @param n_threads Threads for the per-voxel C++ loops (`0` = OpenMP
#'   default). Results do not depend on the thread count. Matrix products use
#'   the BLAS library's own threads; use `n_threads > 1` only with a
#'   single-threaded BLAS, because the two thread pools compete for cores.
#' @param verbose Print progress messages.
#'
#' @return An object of class `glmsingle_fit`: a list with elements `typea`,
#'   `typeb`, `typec`, `typed` (each a list using GLMsingle's field names;
#'   betas are trials x voxels), `meanvol`, `hrf_library`, `design`
#'   (trial bookkeeping), `settings` and `timing`.
#'
#' @references Prince, J. S., Charest, I., Kurzawski, J. W., Pyles, J. A.,
#'   Tarr, M. J., & Kay, K. N. (2022). Improving the accuracy of single-trial
#'   fMRI response estimates using GLMsingle. eLife, 11, e77599.
#'
#' @examples
#' set.seed(1)
#' tr <- 1; n_time <- 120; n_vox <- 40
#' design <- lapply(1:3, function(r) {
#'   D <- matrix(0, n_time, 6)
#'   onsets <- seq(5, 95, by = 8)
#'   D[cbind(onsets, rep_len(sample(6), length(onsets)))] <- 1
#'   D
#' })
#' Y <- lapply(design, function(D) {
#'   X <- apply(D, 2, function(s) stats::filter(s, glmsingle_hrf(3, tr),
#'                                               sides = 1, circular = FALSE))
#'   X[is.na(X)] <- 0
#'   100 + X %*% matrix(rnorm(6 * n_vox, 2), 6) + matrix(rnorm(n_time * n_vox), n_time)
#' })
#' fit <- glmsingle(Y, design, tr = tr, stimdur = 3, n_pcs = 2, verbose = FALSE)
#' dim(fit$typed$betasmd)
#' @export
glmsingle <- function(Y, design, tr, stimdur,
                      runs = NULL,
                      hrf_library = NULL,
                      want_library = TRUE,
                      want_glmdenoise = TRUE,
                      want_fracridge = TRUE,
                      fracs = seq(1, 0.05, by = -0.05),
                      n_pcs = 10L,
                      pcstop = 1.05,
                      xval_scheme = NULL,
                      session_indicator = NULL,
                      extra_regressors = NULL,
                      max_poly_deg = NULL,
                      brain_thresh = c(99, 0.1),
                      brain_r2 = NULL,
                      brain_exclude = NULL,
                      pc_r2_cutoff = NULL,
                      pc_r2_cutoff_mask = NULL,
                      want_percent_bold = TRUE,
                      want_autoscale = TRUE,
                      extras_in_denoise = c("always", "with_pcs"),
                      zero_sd_cv = c("zero", "python"),
                      singular = c("error", "pinv"),
                      frac_alpha = c("fracridge", "exact"),
                      full_glmbadness = FALSE,
                      chunk_size = 50000L,
                      n_threads = 1L,
                      verbose = TRUE) {
  t_start <- proc.time()[["elapsed"]]
  timing <- list()
  tick <- function(stage) {
    now <- proc.time()[["elapsed"]]
    timing[[stage]] <<- now - t_start - sum(unlist(timing))
  }
  say <- function(...) if (verbose) message(...)

  extras_in_denoise <- match.arg(extras_in_denoise)
  zero_sd_cv <- match.arg(zero_sd_cv)
  singular <- match.arg(singular)
  frac_alpha <- match.arg(frac_alpha)
  want_library <- .as_scalar_logical(want_library, "want_library")
  want_glmdenoise <- .as_scalar_logical(want_glmdenoise, "want_glmdenoise")
  want_fracridge <- .as_scalar_logical(want_fracridge, "want_fracridge")
  want_percent_bold <- .as_scalar_logical(want_percent_bold, "want_percent_bold")
  want_autoscale <- .as_scalar_logical(want_autoscale, "want_autoscale")
  full_glmbadness <- .as_scalar_logical(full_glmbadness, "full_glmbadness")
  verbose <- .as_scalar_logical(verbose, "verbose")
  tr <- .glms_positive_scalar(tr, "tr")
  stimdur <- .as_nonnegative_scalar(stimdur, "stimdur")
  n_pcs <- .as_nonnegative_integer(n_pcs, "n_pcs")
  chunk_size <- .as_positive_integer(chunk_size, "chunk_size")
  n_threads <- .as_nonnegative_integer(n_threads, "n_threads")
  old_threads <- .glms_state$n_threads
  .glms_state$n_threads <- n_threads
  on.exit(.glms_state$n_threads <- old_threads, add = TRUE)
  if (!is.numeric(pcstop) || length(pcstop) != 1L || !is.finite(pcstop)) {
    stop("pcstop must be a single finite number", call. = FALSE)
  }
  if (!is.numeric(fracs) || !length(fracs) || any(!is.finite(fracs)) ||
      any(fracs <= 0 | fracs > 1)) {
    stop("fracs must be numbers in (0, 1]", call. = FALSE)
  }
  fracs <- sort(unique(fracs), decreasing = TRUE)
  if (!is.numeric(brain_thresh) || length(brain_thresh) != 2L) {
    stop("brain_thresh must be two numbers", call. = FALSE)
  }

  # ---- inputs and geometry -------------------------------------------------
  Ylist <- .glms_split_runs(Y, runs)
  run_ids <- if (is.list(Y) && !is.data.frame(Y)) NULL else unique(runs)
  R <- length(Ylist)
  n_time <- vapply(Ylist, nrow, integer(1))
  n_vox <- ncol(Ylist[[1]])
  if (any(vapply(Ylist, function(y) !all(is.finite(range(y))), logical(1)))) {
    stop("Y contains non-finite values", call. = FALSE)
  }
  parsed <- .glms_parse_design(design, n_time, tr, run_ids)
  geom <- .glms_geometry(parsed, n_time, tr, session_indicator, xval_scheme)
  if (!geom$n_trials) stop("design contains no trials", call. = FALSE)

  max_poly_deg <- if (is.null(max_poly_deg)) {
    .glms_alt_round(((n_time * tr) / 60) / 2)
  } else {
    d <- .as_integer_ids(max_poly_deg, "max_poly_deg")
    if (length(d) == 1L) d <- rep(d, R)
    if (length(d) != R || any(d < 0L)) stop("max_poly_deg must be one non-negative integer or one per run", call. = FALSE)
    d
  }
  extras <- .glms_check_extras(extra_regressors, n_time)
  for (nm in c("brain_exclude", "pc_r2_cutoff_mask")) {
    m <- get(nm)
    if (!is.null(m) && (length(m) != n_vox || anyNA(as.logical(m)))) {
      stop(sprintf("%s must be a logical vector with one entry per voxel", nm), call. = FALSE)
    }
  }

  hrf0 <- glmsingle_hrf(stimdur, tr)
  hrf_lib <- if (!want_library) {
    matrix(hrf0, ncol = 1L)
  } else if (is.null(hrf_library)) {
    glmsingle_hrf_library(stimdur, tr)
  } else {
    lib <- .as_base_matrix(hrf_library)
    if (any(!is.finite(lib)) || any(apply(lib, 2L, max) <= 0)) {
      stop("hrf_library columns must be finite with a positive peak", call. = FALSE)
    }
    sweep(lib, 2L, apply(lib, 2L, max), "/")
  }

  if (all(geom$cond_in_runs <= 1L)) {
    warning("No condition occurs in more than one run; cross-validation is not possible.",
            call. = FALSE)
    if (want_glmdenoise && pcstop > 0) {
      warning("Setting want_glmdenoise = FALSE (no repeats).", call. = FALSE)
      want_glmdenoise <- FALSE
    }
    if (want_fracridge && length(fracs) > 1L) {
      warning("Setting want_fracridge = FALSE (no repeats).", call. = FALSE)
      want_fracridge <- FALSE
    }
  }
  tiles <- split(seq_len(n_vox), ceiling(seq_len(n_vox) / chunk_size))

  meanvol <- Reduce(`+`, lapply(Ylist, colSums)) / sum(n_time)
  pb <- if (want_percent_bold) 100 / abs(meanvol) else rep(1, n_vox)
  beta_names <- list(paste0("trial", seq_len(geom$n_trials)), colnames(Ylist[[1]]))
  scale_betas <- function(B) {
    structure(glms_scale_cols(B, pb, .glms_nt()), dimnames = beta_names)
  }
  tick("setup")

  nuis_ab <- lapply(seq_len(R), function(r) {
    .glms_nuisance_run(n_time[r], max_poly_deg[r], extras[[r]], NULL, 0L,
                       function(k) TRUE)
  })
  n_singular <- 0L
  withCallingHandlers({
    say("Fitting type-A (ON-OFF) and type-B (HRF library) models")
    b <- .glms_fit_types_ab(Ylist, geom, hrf0, hrf_lib, nuis_ab, tiles, singular)
    typea <- list(onoffR2 = b$onoffR2, meanvol = meanvol, betasmd = b$beta_a * pb)
    typeb <- c(b[c("FitHRFR2", "FitHRFR2run", "HRFindex", "HRFindexrun", "R2", "R2run")],
               list(betasmd = scale_betas(b$beta), meanvol = meanvol))
    b$beta <- NULL
    tick("typeab")

    thresh <- NULL
    if (want_glmdenoise && (is.null(brain_r2) || is.null(pc_r2_cutoff))) {
      thresh <- .glms_tail_threshold(b$onoffR2)
    }
    brain_r2 <- brain_r2 %||% thresh
    pc_r2_cutoff <- pc_r2_cutoff %||% thresh

    dn <- list(pcregressors = NULL, noisepool = NULL, pcnum = 0L,
               xvaltrend = NULL, glmbadness = NULL, pcvoxels = NULL)
    if (want_glmdenoise) {
      say("Deriving GLMdenoise regressors")
      dn <- .glms_denoise(Ylist, geom, hrf_lib, b$HRFindex, b$onoffR2, meanvol,
                          nuis_ab, max_poly_deg, extras, n_pcs, pcstop,
                          brain_thresh, brain_r2, brain_exclude, pc_r2_cutoff,
                          pc_r2_cutoff_mask, full_glmbadness, extras_in_denoise,
                          zero_sd_cv, singular, chunk_size, say)
    }
    tick("glmdenoise")

    typec <- typed <- NULL
    if (want_glmdenoise || want_fracridge) {
      say(sprintf("Fitting type-C/D models (%d noise PCs)", dn$pcnum))
      fracstouse <- if (want_fracridge) unique(c(1, fracs)) else 1
      cd <- .glms_fit_cd(Ylist, geom, hrf_lib, b$HRFindex, max_poly_deg, extras,
                         dn$pcregressors, dn$pcnum, extras_in_denoise, fracstouse,
                         fracs, want_fracridge, want_autoscale, zero_sd_cv,
                         frac_alpha, chunk_size)
      common <- c(list(HRFindex = b$HRFindex, HRFindexrun = b$HRFindexrun),
                  dn[c("glmbadness", "pcvoxels", "pcnum", "xvaltrend",
                       "noisepool", "pcregressors")],
                  list(meanvol = meanvol))
      if (want_glmdenoise) {
        typec <- c(common, list(betasmd = scale_betas(cd$beta_c),
                                R2 = cd$r2_c, R2run = cd$r2run_c))
      }
      if (want_fracridge) {
        typed <- c(common, list(betasmd = scale_betas(cd$beta_d),
                                R2 = cd$r2_d, R2run = cd$r2run_d,
                                FRACvalue = cd$frac_value,
                                scaleoffset = cd$scaleoffset,
                                rrbadness = cd$rrbadness))
      }
    }
    tick("typecd")
  }, glms_singular = function(cnd) n_singular <<- n_singular + 1L)
  if (n_singular) {
    warning(sprintf(paste(
      "Singular trial design in %d fit(s); used the minimum-norm",
      "(pseudoinverse) solution."), n_singular), call. = FALSE)
  }

  structure(list(
    typea = typea, typeb = typeb, typec = typec, typed = typed,
    meanvol = meanvol, hrf_library = hrf_lib, hrf_assumed = hrf0,
    design = list(stimorder = geom$stimorder, condition_levels = geom$levels,
                  trial_run = geom$trial_run, validcolumns = geom$validcolumns,
                  onsets = geom$onsets, n_time = n_time),
    settings = list(tr = tr, stimdur = stimdur, fracs = fracs, n_pcs = n_pcs,
                    pcstop = pcstop, max_poly_deg = max_poly_deg,
                    brain_r2 = brain_r2, pc_r2_cutoff = pc_r2_cutoff,
                    tail_threshold = thresh,
                    want_percent_bold = want_percent_bold,
                    extras_in_denoise = extras_in_denoise,
                    zero_sd_cv = zero_sd_cv, singular = singular,
                    frac_alpha = frac_alpha),
    timing = unlist(timing)
  ), class = "glmsingle_fit")
}

.glms_check_extras <- function(extra_regressors, n_time) {
  R <- length(n_time)
  if (is.null(extra_regressors)) return(vector("list", R))
  if (is.matrix(extra_regressors)) extra_regressors <- list(extra_regressors)
  if (!is.list(extra_regressors) || length(extra_regressors) != R) {
    stop("extra_regressors must be a list with one matrix (or NULL) per run", call. = FALSE)
  }
  lapply(seq_len(R), function(r) {
    x <- extra_regressors[[r]]
    if (is.null(x)) return(NULL)
    x <- as.matrix(x)
    if (!ncol(x)) return(NULL)
    if (nrow(x) != n_time[r] || any(!is.finite(x))) {
      stop(sprintf("extra_regressors[[%d]] must be finite with %d rows", r, n_time[r]),
           call. = FALSE)
    }
    x
  })
}

# numpy.argmax semantics for one vector: first maximum; an all-NaN (or
# leading NaN) vector gives the first NaN position.
.glms_argmax <- function(x) {
  nan <- which(is.na(x))
  if (length(nan)) return(nan[1L])
  which.max(x)
}

.glms_argmax_rows <- function(M) {
  out <- max.col(replace(M, is.na(M), -Inf), ties.method = "first")
  anyna <- which(rowSums(is.na(M)) > 0)
  for (i in anyna) out[i] <- .glms_argmax(M[i, ])
  out
}
