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
#' Arguments are grouped by stage: HRF (`fit_hrf`, `hrf_library`), noise
#' regressors (`denoise`, `max_pcs`, `n_pcs`, `pc_*`, `noise_pool_*`), ridge
#' (`ridge`, `ridge_*`) and cross-validation (`cv_folds`, `sessions`,
#' `cv_zero_variance`). The defaults are GLMsingle's defaults and are a good
#' choice for typical event-related designs; the arguments most often worth
#' setting are `nuisance` (e.g. motion), `cv_folds`/`sessions` for
#' multi-session data, and `chunk_size` when memory is tight.
#'
#' @section Correspondence with GLMsingle options: `fit_hrf` = `wantlibrary`;
#'   `denoise` = `wantglmdenoise`; `ridge` = `wantfracridge`; `ridge_fracs` =
#'   `fracs`; `ridge_rescale` = `wantautoscale`; `max_pcs` = `n_pcs`;
#'   `pc_stop` = `pcstop` (a fixed count, GLMsingle's `pcstop = -k`, is
#'   `n_pcs = k`); `pc_voxel_r2` = `pcR2cutoff`; `pc_voxel_mask` =
#'   `pcR2cutoffmask`; `noise_pool_r2` = `brainR2`; `noise_pool_brightness` =
#'   `brainthresh`; `noise_pool_mask` = `brainexclude`; `cv_folds` =
#'   `xvalscheme` (1-based run indices); `sessions` = `sessionindicator`;
#'   `nuisance` = `extra_regressors`; `poly_degree` = `maxpolydeg`;
#'   `percent_signal` = `wantpercentbold`; `chunk_size` = `chunklen`;
#'   `stim_dur` = `stimdur`. Output fields keep GLMsingle's names.
#'
#' @section Differences from GLMsingle: Defaults reproduce GLMsingle (Python,
#'   commit 1ab54a6) except where GLMsingle is internally inconsistent:
#'   \itemize{
#'     \item `nuisance_in_denoise = "always"` keeps `nuisance` regressors in
#'       every GLMdenoise and ridge fit; GLMsingle (Python) drops them when
#'       zero noise PCs are used, including from the cross-validation
#'       reference. `"with_pcs"` reproduces GLMsingle.
#'     \item `cv_zero_variance = "ignore"` gives voxels with zero beta variance
#'       no weight in cross-validation; GLMsingle's in-place division helper
#'       treats them inconsistently across candidates. `"glmsingle"`
#'       reproduces it.
#'     \item The ON-OFF R^2 threshold uses a deterministic two-Gaussian mixture
#'       fit with a small variance floor (as in GLMsingle's MATLAB code);
#'       GLMsingle's Python fit is unseeded.
#'     \item `percent_signal` applies to all model types (GLMsingle always
#'       scales types C and D).
#'     \item Computation is in double precision (GLMsingle uses single).
#'     \item HRF indices are 1-based; betas are trials x voxels.
#'   }
#'   The FIR diagnostic model, `hrfmodel = "optimize"`, `wantlss`, bootstrap
#'   modes, figures and file outputs are not implemented.
#'
#' @param Y Data: a list of time x voxel matrices (one per run), or a single
#'   time x voxel matrix with `runs` giving the run of each row. Runs may
#'   differ in length but must share the voxels.
#' @param design Either a list of time x condition 0/1 matrices, one per run
#'   (GLMsingle's format; a 1 marks a trial onset), or a data frame with
#'   columns `run`, `onset` (seconds from the start of the run; must lie on
#'   the TR grid) and `condition`. Conditions that occur in more than one run
#'   drive the cross-validation; with no repeats, noise-PC and ridge
#'   selection are switched off with a warning. For matrix `Y`, event run
#'   IDs are matched to `runs` in data order; for list `Y`, ascending event
#'   run IDs correspond to the list order.
#' @param tr Repetition time in seconds.
#' @param stim_dur Trial duration in seconds (sets the duration of the
#'   canonical HRF and of the HRF library).
#' @param runs Run identifier for each row of `Y` when `Y` is a single
#'   matrix; ignored when `Y` is a list.
#' @param nuisance Optional nuisance regressors (e.g. motion, physiology):
#'   a list with one time x regressor matrix (or `NULL`) per run. They are
#'   included in every model; see `nuisance_in_denoise`.
#' @param poly_degree Degree of the per-run polynomial drift model: a single
#'   integer or one per run. Default (`NULL`): `round(run_minutes / 2)`, as in
#'   GLMsingle.
#' @param fit_hrf If `TRUE` (default), choose the best-fitting HRF from
#'   `hrf_library` for each voxel. If `FALSE`, use the canonical HRF
#'   ([glmsingle_hrf()]) for every voxel.
#' @param hrf_library Optional time x HRF matrix sampled at the TR (first row
#'   at trial onset); columns are scaled to unit peak. Default (`NULL`): the
#'   20-HRF GLMsingle library, [glmsingle_hrf_library()]. Ignored when
#'   `fit_hrf = FALSE`.
#' @param denoise If `TRUE` (default), derive GLMdenoise noise regressors
#'   (principal components of a noise pool of voxels) and return the type C
#'   model.
#' @param max_pcs Largest number of noise PCs considered (default 10).
#'   Capped, with a warning, at the noise pool's rank across runs; an empty
#'   or rank-zero pool uses zero PCs and skips the PC-count cross-validation.
#' @param n_pcs Use exactly this many noise PCs instead of choosing the count
#'   by cross-validation. Default (`NULL`): choose by cross-validation.
#' @param pc_stop Stopping factor for the cross-validated PC count (default
#'   1.05): the smallest count whose improvement over zero PCs is at least
#'   `1 / pc_stop` of the best improvement is chosen. Larger values choose
#'   fewer PCs.
#' @param pc_voxel_r2 ON-OFF R^2 (percent) above which a voxel contributes to
#'   choosing the PC count. Default (`NULL`): a data-driven threshold from a
#'   two-Gaussian fit to the ON-OFF R^2 distribution (also stored in
#'   `settings$tail_threshold`).
#' @param pc_voxel_mask Optional logical vector (one per voxel); only `TRUE`
#'   voxels may contribute to choosing the PC count.
#' @param noise_pool_r2 ON-OFF R^2 (percent) below which a voxel may enter the
#'   noise pool. Default (`NULL`): the same data-driven threshold as
#'   `pc_voxel_r2`. The threshold is only estimated when `denoise = TRUE`;
#'   with fewer than two distinct finite ON-OFF R^2 values it is their common
#'   value (or zero if none is finite).
#' @param noise_pool_brightness Two numbers `c(p, f)`: only voxels whose mean
#'   signal exceeds `f` times the `p`-th percentile of the mean volume may
#'   enter the noise pool (excludes out-of-brain voxels). Default `c(99, 0.1)`.
#' @param noise_pool_mask Optional logical vector (one per voxel); only `TRUE`
#'   voxels may enter the noise pool.
#' @param ridge If `TRUE` (default), fit voxel-wise fractional ridge
#'   regression and return the type D model.
#' @param ridge_fracs Ridge fractions in (0, 1] to choose from by
#'   cross-validation; a fraction is the length of the ridge solution relative
#'   to the OLS solution (1 = no shrinkage). A single value skips
#'   cross-validation and uses that fraction everywhere. Default: 1 to 0.05 in
#'   steps of 0.05.
#' @param ridge_rescale If `TRUE` (default), rescale each voxel's ridge betas
#'   with the linear fit (slope, offset) that best matches its unregularised
#'   betas, undoing the overall shrinkage while keeping the ridge pattern.
#' @param ridge_alpha How fractions are converted to ridge penalties:
#'   `"grid"` (default; fracridge's interpolation on a log-spaced penalty
#'   grid, as GLMsingle) or `"exact"` (solve for the penalty that gives the
#'   fraction exactly; experimental).
#' @param cv_folds Cross-validation folds: a list of integer vectors of runs
#'   held out together. Default (`NULL`): leave one run out.
#' @param sessions Session of each run (one entry per run). Betas are
#'   z-scored within session before cross-validation, removing
#'   session-level gain differences. Default (`NULL`): one session.
#' @param percent_signal If `TRUE` (default), express betas as percent signal
#'   change of each voxel's mean; if `FALSE`, in raw data units.
#' @param nuisance_in_denoise How `nuisance` regressors enter the GLMdenoise
#'   and ridge fits: `"always"` (default, recommended) or `"with_pcs"`, which
#'   drops them when zero PCs are used (GLMsingle's behaviour).
#' @param cv_zero_variance Treatment, in cross-validation, of voxels whose
#'   reference betas have zero variance within a session (e.g. all-zero
#'   voxels): `"ignore"` (default; they carry no weight) or `"glmsingle"`
#'   (reproduce GLMsingle's Python behaviour).
#' @param singular What to do when a run's trial regressors are linearly
#'   dependent (e.g. two trials with the same onset): `"error"` (default, as
#'   GLMsingle) or `"pinv"` (minimum-norm solution with a warning).
#' @param pc_cv_all_voxels If `TRUE`, run the PC-count cross-validation for
#'   every voxel rather than only those that decide the count. Estimates are
#'   identical; only the `glmbadness` diagnostic is filled for all voxels.
#'   Default `FALSE` (faster).
#' @param chunk_size Number of voxels processed at a time (default 50000).
#'   Lower it to reduce peak memory; results do not depend on it.
#' @param n_threads Threads for the per-voxel C++ loops (default 1; `0` =
#'   OpenMP default). Results do not depend on the thread count. Matrix
#'   products use the BLAS library's own threads; use `n_threads > 1` only
#'   with a single-threaded BLAS, because the two thread pools compete for
#'   cores.
#' @param verbose Print progress messages (default `TRUE`).
#'
#' @return An object of class `glmsingle_fit`: a list with elements `typea`,
#'   `typeb`, `typec`, `typed` (each a list using GLMsingle's field names;
#'   betas are trials x voxels; `typec`/`typed` are `NULL` when `denoise`/
#'   `ridge` is `FALSE`), `meanvol`, `hrf_library`, `hrf_assumed`, `design`
#'   (trial bookkeeping), `settings` (the resolved settings, including
#'   data-driven thresholds) and `timing`.
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
#' fit <- glmsingle(Y, design, tr = tr, stim_dur = 3, max_pcs = 2, verbose = FALSE)
#' dim(fit$typed$betasmd)
#' @export
glmsingle <- function(Y, design, tr, stim_dur,
                      runs = NULL,
                      nuisance = NULL,
                      poly_degree = NULL,
                      fit_hrf = TRUE,
                      hrf_library = NULL,
                      denoise = TRUE,
                      max_pcs = 10L,
                      n_pcs = NULL,
                      pc_stop = 1.05,
                      pc_voxel_r2 = NULL,
                      pc_voxel_mask = NULL,
                      noise_pool_r2 = NULL,
                      noise_pool_brightness = c(99, 0.1),
                      noise_pool_mask = NULL,
                      ridge = TRUE,
                      ridge_fracs = seq(1, 0.05, by = -0.05),
                      ridge_rescale = TRUE,
                      ridge_alpha = c("grid", "exact"),
                      cv_folds = NULL,
                      sessions = NULL,
                      percent_signal = TRUE,
                      nuisance_in_denoise = c("always", "with_pcs"),
                      cv_zero_variance = c("ignore", "glmsingle"),
                      singular = c("error", "pinv"),
                      pc_cv_all_voxels = FALSE,
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

  nuisance_in_denoise <- match.arg(nuisance_in_denoise)
  cv_zero_variance <- match.arg(cv_zero_variance)
  singular <- match.arg(singular)
  ridge_alpha <- match.arg(ridge_alpha)
  fit_hrf <- .as_scalar_logical(fit_hrf, "fit_hrf")
  denoise <- .as_scalar_logical(denoise, "denoise")
  ridge <- .as_scalar_logical(ridge, "ridge")
  percent_signal <- .as_scalar_logical(percent_signal, "percent_signal")
  ridge_rescale <- .as_scalar_logical(ridge_rescale, "ridge_rescale")
  pc_cv_all_voxels <- .as_scalar_logical(pc_cv_all_voxels, "pc_cv_all_voxels")
  verbose <- .as_scalar_logical(verbose, "verbose")
  tr <- .glms_positive_scalar(tr, "tr")
  stim_dur <- .as_nonnegative_scalar(stim_dur, "stim_dur")
  max_pcs <- .as_nonnegative_integer(max_pcs, "max_pcs")
  chunk_size <- .as_positive_integer(chunk_size, "chunk_size")
  n_threads <- .as_nonnegative_integer(n_threads, "n_threads")
  old_threads <- .glms_state$n_threads
  .glms_state$n_threads <- n_threads
  on.exit(.glms_state$n_threads <- old_threads, add = TRUE)
  if (!is.numeric(pc_stop) || length(pc_stop) != 1L || !is.finite(pc_stop) || pc_stop <= 0) {
    stop("pc_stop must be a single positive number", call. = FALSE)
  }
  if (!is.null(n_pcs)) {
    # A fixed count is GLMsingle's pcstop = -n_pcs.
    n_pcs <- .as_nonnegative_integer(n_pcs, "n_pcs")
    max_pcs <- max(max_pcs, n_pcs)
  }
  pc_rule <- if (is.null(n_pcs)) pc_stop else -n_pcs
  fracs <- ridge_fracs
  if (!is.numeric(fracs) || !length(fracs) || any(!is.finite(fracs)) ||
      any(fracs <= 0 | fracs > 1)) {
    stop("ridge_fracs must be numbers in (0, 1]", call. = FALSE)
  }
  fracs <- sort(unique(fracs), decreasing = TRUE)
  if (!is.numeric(noise_pool_brightness) || length(noise_pool_brightness) != 2L ||
      any(!is.finite(noise_pool_brightness))) {
    stop("noise_pool_brightness must be two finite numbers", call. = FALSE)
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
  geom <- .glms_geometry(parsed, n_time, tr, sessions, cv_folds)
  if (!geom$n_trials) stop("design contains no trials", call. = FALSE)

  poly_degree <- if (is.null(poly_degree)) {
    .glms_alt_round(((n_time * tr) / 60) / 2)
  } else {
    d <- .as_integer_ids(poly_degree, "poly_degree")
    if (length(d) == 1L) d <- rep(d, R)
    if (length(d) != R || any(d < 0L)) stop("poly_degree must be one non-negative integer or one per run", call. = FALSE)
    d
  }
  extras <- .glms_check_extras(nuisance, n_time)
  for (nm in c("noise_pool_mask", "pc_voxel_mask")) {
    m <- get(nm)
    if (!is.null(m) && (length(m) != n_vox || anyNA(as.logical(m)))) {
      stop(sprintf("%s must be a logical vector with one entry per voxel", nm), call. = FALSE)
    }
  }

  hrf0 <- glmsingle_hrf(stim_dur, tr)
  hrf_lib <- if (!fit_hrf) {
    matrix(hrf0, ncol = 1L)
  } else if (is.null(hrf_library)) {
    glmsingle_hrf_library(stim_dur, tr)
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
    if (denoise && is.null(n_pcs)) {
      warning("Setting denoise = FALSE (no repeats).", call. = FALSE)
      denoise <- FALSE
    }
    if (ridge && length(fracs) > 1L) {
      warning("Setting ridge = FALSE (no repeats).", call. = FALSE)
      ridge <- FALSE
    }
  }
  tiles <- split(seq_len(n_vox), ceiling(seq_len(n_vox) / chunk_size))

  meanvol <- Reduce(`+`, lapply(Ylist, colSums)) / sum(n_time)
  pb <- if (percent_signal) 100 / abs(meanvol) else rep(1, n_vox)
  beta_names <- list(paste0("trial", seq_len(geom$n_trials)), colnames(Ylist[[1]]))
  scale_betas <- function(B) {
    structure(glms_scale_cols(B, pb, .glms_nt()), dimnames = beta_names)
  }
  tick("setup")

  nuis_ab <- lapply(seq_len(R), function(r) {
    .glms_nuisance_run(n_time[r], poly_degree[r], extras[[r]], NULL, 0L,
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
    if (denoise && (is.null(noise_pool_r2) || is.null(pc_voxel_r2))) {
      thresh <- .glms_tail_threshold(b$onoffR2)
    }
    noise_pool_r2 <- noise_pool_r2 %||% thresh
    pc_voxel_r2 <- pc_voxel_r2 %||% thresh

    dn <- list(pcregressors = NULL, noisepool = NULL, pcnum = 0L,
               xvaltrend = NULL, glmbadness = NULL, pcvoxels = NULL)
    if (denoise) {
      say("Deriving GLMdenoise regressors")
      dn <- .glms_denoise(Ylist, geom, hrf_lib, b$HRFindex, b$onoffR2, meanvol,
                          nuis_ab, poly_degree, extras, max_pcs, pc_rule,
                          noise_pool_brightness, noise_pool_r2, noise_pool_mask, pc_voxel_r2,
                          pc_voxel_mask, pc_cv_all_voxels, nuisance_in_denoise,
                          cv_zero_variance, singular, chunk_size, say)
    }
    tick("glmdenoise")

    typec <- typed <- NULL
    if (denoise || ridge) {
      say(sprintf("Fitting type-C/D models (%d noise PCs)", dn$pcnum))
      fracstouse <- if (ridge) unique(c(1, fracs)) else 1
      cd <- .glms_fit_cd(Ylist, geom, hrf_lib, b$HRFindex, poly_degree, extras,
                         dn$pcregressors, dn$pcnum, nuisance_in_denoise, fracstouse,
                         fracs, ridge, ridge_rescale, cv_zero_variance,
                         ridge_alpha, chunk_size)
      common <- c(list(HRFindex = b$HRFindex, HRFindexrun = b$HRFindexrun),
                  dn[c("glmbadness", "pcvoxels", "pcnum", "xvaltrend",
                       "noisepool", "pcregressors")],
                  list(meanvol = meanvol))
      if (denoise) {
        typec <- c(common, list(betasmd = scale_betas(cd$beta_c),
                                R2 = cd$r2_c, R2run = cd$r2run_c))
      }
      if (ridge) {
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
    settings = list(tr = tr, stim_dur = stim_dur, poly_degree = poly_degree,
                    fit_hrf = fit_hrf, denoise = denoise, max_pcs = max_pcs,
                    n_pcs = n_pcs, pc_stop = pc_stop, ridge = ridge,
                    ridge_fracs = fracs, ridge_rescale = ridge_rescale,
                    noise_pool_r2 = noise_pool_r2, pc_voxel_r2 = pc_voxel_r2,
                    tail_threshold = thresh,
                    percent_signal = percent_signal,
                    nuisance_in_denoise = nuisance_in_denoise,
                    cv_zero_variance = cv_zero_variance, singular = singular,
                    ridge_alpha = ridge_alpha),
    timing = unlist(timing)
  ), class = "glmsingle_fit")
}

.glms_check_extras <- function(nuisance, n_time) {
  R <- length(n_time)
  if (is.null(nuisance)) return(vector("list", R))
  if (is.matrix(nuisance)) nuisance <- list(nuisance)
  if (!is.list(nuisance) || length(nuisance) != R) {
    stop("nuisance must be a list with one matrix (or NULL) per run", call. = FALSE)
  }
  lapply(seq_len(R), function(r) {
    x <- nuisance[[r]]
    if (is.null(x)) return(NULL)
    x <- as.matrix(x)
    if (!ncol(x)) return(NULL)
    if (nrow(x) != n_time[r] || any(!is.finite(x))) {
      stop(sprintf("nuisance[[%d]] must be finite with %d rows", r, n_time[r]),
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
