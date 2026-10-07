#' fmridesign front end for glmsingle()
#'
#' Fits [glmsingle()] using an `fmridesign::event_model()` to describe the
#' trials. The event model must contain one event term with a single factor
#' whose levels are the conditions (repeated conditions drive
#' cross-validation). Onsets must lie on the TR grid.
#'
#' A `baseline_model` contributes its nuisance term (e.g. motion) as
#' `extra_regressors`. Its drift and block terms are not used: GLMsingle
#' models drift with its own per-run polynomials (`max_poly_deg`).
#'
#' @param Y Time x voxel data matrix covering all runs of the sampling frame.
#' @param event_model An `event_model` from fmridesign.
#' @param baseline_model Optional `baseline_model` from fmridesign.
#' @param stimdur Trial duration in seconds. Default: the event durations,
#'   which must then be a single positive value.
#' @param ... Further arguments passed to [glmsingle()].
#' @return A `glmsingle_fit`; see [glmsingle()].
#' @examplesIf requireNamespace("fmridesign", quietly = TRUE)
#' set.seed(1)
#' sf <- fmrihrf::sampling_frame(c(80L, 80L), TR = 1)
#' events <- data.frame(onset = rep(c(5, 20, 35, 50), 2),
#'                      run = rep(1:2, each = 4),
#'                      stimulus = factor(rep(c("A", "B"), 4)))
#' em <- fmridesign::event_model(
#'   onset ~ fmridesign::hrf(stimulus), data = events, block = ~run,
#'   sampling_frame = sf, durations = rep(2, nrow(events)))
#' Y <- matrix(100 + rnorm(160 * 4), 160, 4)
#' fit <- glmsingle_design(Y, em, want_glmdenoise = FALSE,
#'                         want_fracridge = FALSE, verbose = FALSE)
#' dim(coef(fit, type = "b"))
#' @export
glmsingle_design <- function(Y, event_model, baseline_model = NULL, stimdur = NULL, ...) {
  if (!requireNamespace("fmridesign", quietly = TRUE)) {
    stop("Package 'fmridesign' is required for glmsingle_design()", call. = FALSE)
  }
  if (!inherits(event_model, "event_model")) {
    stop("event_model must be an fmridesign event_model", call. = FALSE)
  }
  if (length(event_model$terms) != 1L) {
    stop("glmsingle_design() needs an event model with exactly one event term", call. = FALSE)
  }
  term <- event_model$terms[[1L]]
  tab <- as.data.frame(term$event_table)
  if (ncol(tab) != 1L || !is.factor(tab[[1L]])) {
    stop("the event term must have a single factor variable (the condition)", call. = FALSE)
  }
  sframe <- event_model$sampling_frame
  blocklens <- fmrihrf::blocklens(sframe)
  tr <- unique(sframe$TR)
  if (length(tr) != 1L) stop("glmsingle_design() requires a single TR", call. = FALSE)
  Y <- .as_base_matrix(Y)
  if (nrow(Y) != sum(blocklens)) {
    stop(sprintf("Y has %d rows but the sampling frame has %d scans", nrow(Y), sum(blocklens)),
         call. = FALSE)
  }
  if (is.null(stimdur)) {
    stimdur <- unique(term$durations)
    if (length(stimdur) != 1L || stimdur <= 0) {
      stop("supply stimdur: event durations are not a single positive value", call. = FALSE)
    }
  }
  events <- data.frame(run = match(term$blockids, unique(event_model$blockids)),
                       onset = term$onsets, condition = tab[[1L]])
  runs <- rep(seq_along(blocklens), blocklens)
  extras <- NULL
  if (!is.null(baseline_model)) {
    mats <- fmridesign::term_matrices(baseline_model)
    if (!is.null(mats$nuisance)) {
      N <- as.matrix(mats$nuisance)
      extras <- lapply(seq_along(blocklens), function(r) {
        x <- N[runs == r, , drop = FALSE]
        x[, colSums(x != 0) > 0, drop = FALSE]
      })
    }
  }
  glmsingle(Y, events, tr = tr, stimdur = stimdur, runs = runs,
            extra_regressors = extras, ...)
}

#' @export
print.glmsingle_fit <- function(x, ...) {
  d <- x$design
  cat("<glmsingle_fit>\n")
  cat(sprintf("  %d runs, %d trials, %d conditions, %d voxels\n",
              length(d$n_time), length(d$stimorder), length(d$condition_levels),
              length(x$meanvol)))
  types <- c(typea = "A (ON-OFF)", typeb = "B (HRF library)",
             typec = "C (GLMdenoise)", typed = "D (fractional ridge)")
  cat("  models:", paste(types[names(types)[!vapply(x[names(types)], is.null, logical(1))]],
                         collapse = ", "), "\n")
  if (!is.null(x$typec %||% x$typed)) {
    cat(sprintf("  noise PCs selected: %d\n", (x$typec %||% x$typed)$pcnum))
  }
  if (!is.null(x$typed)) {
    cat(sprintf("  median ridge fraction: %.2f\n", stats::median(x$typed$FRACvalue, na.rm = TRUE)))
  }
  cat(sprintf("  elapsed: %.1f s\n", sum(x$timing)))
  invisible(x)
}

#' @export
summary.glmsingle_fit <- function(object, ...) {
  out <- list(
    n_voxels = length(object$meanvol),
    n_trials = length(object$design$stimorder),
    hrf_index = table(object$typeb$HRFindex),
    pcnum = (object$typec %||% object$typed)$pcnum,
    frac = if (!is.null(object$typed)) table(object$typed$FRACvalue) else NULL,
    r2 = vapply(c("typeb", "typec", "typed"), function(t) {
      if (is.null(object[[t]])) NA_real_ else stats::median(object[[t]]$R2, na.rm = TRUE)
    }, numeric(1)),
    timing = object$timing
  )
  class(out) <- "summary.glmsingle_fit"
  out
}

#' @export
print.summary.glmsingle_fit <- function(x, ...) {
  cat(sprintf("glmsingle fit: %d voxels, %d trials\n", x$n_voxels, x$n_trials))
  cat("HRF index counts:\n"); print(x$hrf_index)
  if (!is.null(x$pcnum)) cat(sprintf("Noise PCs: %d\n", x$pcnum))
  if (!is.null(x$frac)) { cat("Ridge fraction counts:\n"); print(x$frac) }
  cat("Median R^2 (%):\n"); print(round(x$r2, 2))
  cat("Stage timing (s):\n"); print(round(x$timing, 2))
  invisible(x)
}

#' Single-trial betas from a glmsingle fit
#'
#' @param object A `glmsingle_fit`.
#' @param type Model type: `"d"` (fractional ridge, default), `"c"`
#'   (GLMdenoise), `"b"` (HRF library) or `"a"` (ON-OFF, one value per voxel).
#' @param ... Unused.
#' @return Trials x voxels matrix (a vector for type `"a"`).
#' @examples
#' example("glmsingle", echo = FALSE)
#' beta <- coef(fit)
#' dim(beta)
#' @export
coef.glmsingle_fit <- function(object, type = c("d", "c", "b", "a"), ...) {
  type <- match.arg(type)
  m <- object[[paste0("type", type)]]
  if (is.null(m)) stop(sprintf("model type %s was not fitted", toupper(type)), call. = FALSE)
  m$betasmd
}
