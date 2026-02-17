#' Build ITEM design metadata
#'
#' Build and validate trial-wise design objects used by the ITEM helper layer.
#' The returned object has class `item_bundle` and carries trial-level metadata
#' needed by `item_cv()` and `item_slice_fold()`.
#'
#' @param X_t Numeric trial-wise design matrix with shape `n_time x n_trials`.
#' @param T_target Optional supervised targets with `n_trials` rows. Accepts:
#'   numeric matrix/data.frame, numeric vector (regression), or
#'   factor/character/logical vector (classification labels).
#' @param run_id Optional run/session identifier of length `n_trials`.
#'   If `NULL`, all trials are assigned to a single run.
#' @param C_transform Optional transformation matrix used for
#'   `X = X_t %*% C_transform`. Must have `n_trials` rows when provided.
#' @param trial_id Optional trial identifier vector of length `n_trials`.
#'   Defaults to `colnames(X_t)` when available, else `Trial_1..Trial_n`.
#' @param trial_hash Optional hash used by alignment guards.
#'   If supplied, `item_cv(..., check_hash = TRUE)` validates it.
#' @param meta Optional metadata list.
#' @param diagnostics Optional diagnostics list to attach to the bundle.
#' @param validate Logical; when `TRUE` run strict structure checks.
#'
#' @return An object of class `item_bundle` with fields:
#'   `Gamma`, `X_t`, `C_transform`, `T_target`, `U`, `U_by_run`, `run_id`,
#'   `trial_id`, `trial_hash`, `trial_info`, `meta`, and `diagnostics`.
#'
#' @export
item_build_design <- function(X_t,
                              T_target = NULL,
                              run_id = NULL,
                              C_transform = NULL,
                              trial_id = NULL,
                              trial_hash = NULL,
                              meta = list(),
                              diagnostics = list(),
                              validate = TRUE) {
  X_t <- .item_as_numeric_matrix(X_t, "X_t")
  n_trials <- ncol(X_t)

  run_id <- .item_as_run_id(run_id, n_trials)
  trial_id <- .item_as_trial_id(trial_id, X_t, n_trials)

  target_info <- .item_prepare_targets(T_target, n_trials)

  C_transform <- .item_as_matrix_optional(C_transform, "C_transform")
  if (!is.null(C_transform) && nrow(C_transform) != n_trials) {
    stop(
      sprintf(
        "C_transform must have %d rows (got %d).",
        n_trials,
        nrow(C_transform)
      ),
      call. = FALSE
    )
  }

  trial_info <- data.frame(
    trial_index = seq_len(n_trials),
    trial_id = as.character(trial_id),
    run_id = run_id,
    stringsAsFactors = FALSE
  )

  meta <- utils::modifyList(
    list(
      target_type = target_info$target_type,
      class_levels = target_info$class_levels
    ),
    meta
  )

  bundle <- .item_new_bundle(
    Gamma = NULL,
    X_t = X_t,
    C_transform = C_transform,
    T_target = target_info$T_target,
    U = NULL,
    U_by_run = NULL,
    run_id = run_id,
    trial_id = trial_id,
    trial_hash = trial_hash,
    trial_info = trial_info,
    meta = meta,
    diagnostics = diagnostics
  )

  if (isTRUE(validate)) {
    .item_validate_bundle(
      bundle,
      require_gamma = FALSE,
      require_u = FALSE,
      check_hash = FALSE
    )
  }

  bundle
}

#' Validate an ITEM bundle
#'
#' @param bundle `item_bundle` object.
#' @param require_gamma Require non-NULL `Gamma`.
#' @param require_u Require either `U` or `U_by_run`.
#' @param check_hash Validate `trial_hash` when present.
#'
#' @return Invisibly returns `TRUE` when validation passes.
#' @keywords internal
#' @noRd
.item_validate_bundle <- function(bundle,
                                  require_gamma = TRUE,
                                  require_u = TRUE,
                                  check_hash = FALSE) {
  if (!inherits(bundle, "item_bundle")) {
    stop("bundle must inherit class 'item_bundle'.", call. = FALSE)
  }

  run_id <- bundle[["run_id"]]
  Gamma <- bundle[["Gamma"]]
  T_target <- bundle[["T_target"]]
  U <- bundle[["U"]]
  U_by_run <- bundle[["U_by_run"]]
  trial_id <- bundle[["trial_id"]]
  trial_hash <- bundle[["trial_hash"]]

  if (is.null(run_id)) {
    stop("bundle$run_id is required.", call. = FALSE)
  }

  n_trials <- length(run_id)
  if (n_trials < 2L) {
    stop("bundle must contain at least 2 trials.", call. = FALSE)
  }

  if (isTRUE(require_gamma) && is.null(Gamma)) {
    stop("bundle$Gamma is required for this operation.", call. = FALSE)
  }

  if (!is.null(Gamma)) {
    if (!is.matrix(Gamma)) {
      stop("bundle$Gamma must be a matrix.", call. = FALSE)
    }
    if (nrow(Gamma) != n_trials) {
      stop(
        sprintf(
          "bundle$Gamma has %d rows but run_id has length %d.",
          nrow(Gamma),
          n_trials
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(T_target)) {
    if (!is.matrix(T_target)) {
      stop("bundle$T_target must be a matrix.", call. = FALSE)
    }
    if (nrow(T_target) != n_trials) {
      stop(
        sprintf(
          "bundle$T_target has %d rows but run_id has length %d.",
          nrow(T_target),
          n_trials
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(U)) {
    if (!is.matrix(U)) {
      stop("bundle$U must be a matrix.", call. = FALSE)
    }
    if (nrow(U) != n_trials || ncol(U) != n_trials) {
      stop(
        sprintf(
          "bundle$U must be %d x %d but is %d x %d.",
          n_trials,
          n_trials,
          nrow(U),
          ncol(U)
        ),
        call. = FALSE
      )
    }
  }

  if (!is.null(U_by_run)) {
    .item_validate_u_by_run(U_by_run, run_id)
  }

  if (isTRUE(require_u) && is.null(U) && is.null(U_by_run)) {
    stop("bundle must include either U (full matrix) or U_by_run (run-block list).", call. = FALSE)
  }

  if (!is.null(trial_id) && length(trial_id) != n_trials) {
    stop("bundle$trial_id must have length equal to number of trials.", call. = FALSE)
  }

  if (isTRUE(check_hash) && !is.null(trial_hash)) {
    .item_check_trial_hash(trial_id, trial_hash)
  }

  invisible(TRUE)
}

#' Update or create an ITEM bundle
#'
#' @param Gamma Trial-wise beta matrix (`n_trials x n_features`).
#' @param X_t Trial-wise design matrix (`n_time x n_trials`).
#' @param C_transform Optional transformation matrix.
#' @param T_target Trial-level targets (`n_trials x p`).
#' @param U Full trial covariance matrix (`n_trials x n_trials`).
#' @param U_by_run Optional run-block trial covariance list.
#' @param run_id Run/session identifiers.
#' @param trial_id Trial identifiers.
#' @param trial_hash Optional trial hash.
#' @param trial_info Optional trial info table.
#' @param meta Metadata list.
#' @param diagnostics Diagnostics list.
#'
#' @return Object of class `item_bundle`.
#' @keywords internal
#' @noRd
.item_new_bundle <- function(Gamma,
                             X_t,
                             C_transform,
                             T_target,
                             U,
                             U_by_run,
                             run_id,
                             trial_id,
                             trial_hash,
                             trial_info,
                             meta,
                             diagnostics) {
  structure(
    list(
      Gamma = Gamma,
      X_t = X_t,
      C_transform = C_transform,
      T_target = T_target,
      U = U,
      U_by_run = U_by_run,
      run_id = run_id,
      trial_id = trial_id,
      trial_hash = trial_hash,
      trial_info = trial_info,
      meta = meta,
      diagnostics = diagnostics
    ),
    class = "item_bundle"
  )
}

#' @keywords internal
#' @noRd
.item_prepare_targets <- function(T_target, n_trials) {
  if (is.null(T_target)) {
    return(list(
      T_target = matrix(numeric(0), nrow = n_trials, ncol = 0),
      target_type = "unspecified",
      class_levels = NULL
    ))
  }

  if (is.matrix(T_target) || is.data.frame(T_target)) {
    tmat <- .as_base_matrix(T_target)
    if (nrow(tmat) != n_trials) {
      stop(
        sprintf("T_target must have %d rows (got %d).", n_trials, nrow(tmat)),
        call. = FALSE
      )
    }
    if (!is.numeric(tmat)) {
      stop(
        paste(
          "T_target matrix/data.frame must be numeric.",
          "For class labels, pass a factor/character vector."
        ),
        call. = FALSE
      )
    }

    storage.mode(tmat) <- "double"
    return(list(
      T_target = tmat,
      target_type = "matrix",
      class_levels = NULL
    ))
  }

  if (!is.atomic(T_target) || is.list(T_target)) {
    stop("Unsupported T_target type.", call. = FALSE)
  }

  if (length(T_target) != n_trials) {
    stop(
      sprintf("T_target vector must have length %d.", n_trials),
      call. = FALSE
    )
  }

  if (is.factor(T_target) || is.character(T_target) || is.logical(T_target)) {
    if (anyNA(T_target)) {
      stop("T_target labels must not contain missing values.", call. = FALSE)
    }

    class_levels <- sort(unique(as.character(T_target)))
    if (length(class_levels) < 2L) {
      stop("Classification targets must include at least 2 classes.", call. = FALSE)
    }

    fac <- factor(as.character(T_target), levels = class_levels)
    tmat <- stats::model.matrix(~ fac - 1)
    colnames(tmat) <- class_levels

    return(list(
      T_target = tmat,
      target_type = "classification",
      class_levels = class_levels
    ))
  }

  if (!is.numeric(T_target)) {
    stop("Numeric target vector expected for regression mode.", call. = FALSE)
  }

  tmat <- matrix(as.numeric(T_target), ncol = 1)
  colnames(tmat) <- "target_1"

  list(
    T_target = tmat,
    target_type = "regression",
    class_levels = NULL
  )
}

#' @keywords internal
#' @noRd
.item_as_numeric_matrix <- function(x, name) {
  x <- .as_base_matrix(x)
  if (!is.matrix(x)) {
    stop(sprintf("%s must be coercible to a matrix.", name), call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

#' @keywords internal
#' @noRd
.item_as_matrix_optional <- function(x, name) {
  if (is.null(x)) return(NULL)
  x <- .as_base_matrix(x)
  if (!is.matrix(x)) {
    stop(sprintf("%s must be coercible to a matrix.", name), call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(sprintf("%s must be numeric.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

#' @keywords internal
#' @noRd
.item_as_run_id <- function(run_id, n_trials) {
  if (is.null(run_id)) {
    return(rep.int(1L, n_trials))
  }
  if (length(run_id) != n_trials) {
    stop(
      sprintf("run_id must have length %d.", n_trials),
      call. = FALSE
    )
  }
  if (anyNA(run_id)) {
    stop("run_id must not contain missing values.", call. = FALSE)
  }
  run_id
}

#' @keywords internal
#' @noRd
.item_as_trial_id <- function(trial_id, X_t, n_trials) {
  if (is.null(trial_id)) {
    trial_id <- .trial_names_from_cols(X_t, n_trials)
  }
  if (length(trial_id) != n_trials) {
    stop(
      sprintf("trial_id must have length %d.", n_trials),
      call. = FALSE
    )
  }
  if (anyNA(trial_id)) {
    stop("trial_id must not contain missing values.", call. = FALSE)
  }
  as.character(trial_id)
}

#' @keywords internal
#' @noRd
.item_validate_u_by_run <- function(U_by_run, run_id) {
  if (!is.list(U_by_run) || length(U_by_run) < 1L) {
    stop("U_by_run must be a non-empty list of square matrices.", call. = FALSE)
  }

  runs <- sort(unique(run_id))

  has_names <- !is.null(names(U_by_run)) && all(names(U_by_run) != "")
  if (!has_names && length(U_by_run) != length(runs)) {
    stop(
      paste(
        "Unnamed U_by_run must have one block per unique run_id",
        sprintf("(expected %d, got %d).", length(runs), length(U_by_run))
      ),
      call. = FALSE
    )
  }

  for (i in seq_along(runs)) {
    run_value <- runs[[i]]
    key <- as.character(run_value)

    block <- if (has_names) {
      U_by_run[[key]]
    } else {
      U_by_run[[i]]
    }

    if (is.null(block)) {
      stop(
        sprintf("Missing U block for run '%s'.", key),
        call. = FALSE
      )
    }

    block <- .as_base_matrix(block)
    if (!is.matrix(block) || !is.numeric(block)) {
      stop(
        sprintf("U_by_run[['%s']] must be a numeric matrix.", key),
        call. = FALSE
      )
    }
    if (nrow(block) != ncol(block)) {
      stop(
        sprintf("U_by_run[['%s']] must be square.", key),
        call. = FALSE
      )
    }

    expected_n <- sum(run_id == run_value)
    if (nrow(block) != expected_n) {
      stop(
        sprintf(
          "U_by_run block for run '%s' must be %d x %d (got %d x %d).",
          key,
          expected_n,
          expected_n,
          nrow(block),
          ncol(block)
        ),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.item_check_trial_hash <- function(trial_id, trial_hash) {
  if (is.null(trial_id)) {
    stop("trial_id required for hash checks.", call. = FALSE)
  }

  expected <- if (is.list(trial_hash)) {
    trial_hash$trial_id %||% trial_hash[[1]]
  } else {
    trial_hash
  }

  if (is.null(expected) || length(expected) != 1L) {
    stop("trial_hash must be a scalar string or list with trial_id hash.", call. = FALSE)
  }

  actual <- .item_simple_hash(trial_id)
  if (!identical(as.character(expected), actual)) {
    stop(
      paste(
        "Trial hash mismatch:",
        sprintf("expected '%s' but computed '%s'.", expected, actual)
      ),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.item_simple_hash <- function(x) {
  raw <- serialize(x, connection = NULL, version = 2)
  ints <- as.integer(raw)
  if (length(ints) == 0L) return("00000000")

  weights <- ((seq_along(ints) - 1L) %% 251L) + 1L
  accum <- sum((ints + 1L) * weights)
  sprintf("%08x", as.integer(accum %% 2147483647))
}
