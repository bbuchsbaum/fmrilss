#' Crossvalidated ITEM decoding
#'
#' Run deterministic leave-one-run-out (LOSO) crossvalidation for classification
#' or regression using ITEM decoder weights.
#'
#' @param Gamma Trial-wise beta matrix (`n_trials x n_features`) or an
#'   `item_bundle` object.
#' @param T_target Supervised targets (`n_trials x p`) or target vector.
#'   Ignored when `Gamma` is an `item_bundle`.
#' @param U Trial covariance as full matrix (`n_trials x n_trials`) or run-block
#'   list. Ignored when `Gamma` is an `item_bundle`.
#' @param run_id Run/session id vector (`n_trials`). Ignored when `Gamma` is an
#'   `item_bundle`.
#' @param mode Decoding mode: `"classification"` or `"regression"`.
#' @param metric Optional metric name.
#'   Classification: `"accuracy"` (default), `"balanced_accuracy"`.
#'   Regression: `"correlation"` (default), `"rmse"`.
#' @param ridge Ridge passed to `item_fit()`.
#' @param method Solver preference for `item_fit()`.
#' @param class_levels Optional fixed class order for classification.
#' @param trial_id Optional trial id vector (used if `Gamma` is not a bundle).
#' @param trial_hash Optional trial hash (used if `Gamma` is not a bundle).
#' @param check_hash Logical; validate stored trial hash before CV.
#'
#' @return Object of class `item_cv_result` with per-fold metrics, aggregate
#'   metric summary, and trial-level predictions.
#'
#' @export
item_cv <- function(Gamma,
                    T_target = NULL,
                    U = NULL,
                    run_id = NULL,
                    mode = c("classification", "regression"),
                    metric = NULL,
                    ridge = 0,
                    method = c("chol", "svd", "pinv"),
                    class_levels = NULL,
                    trial_id = NULL,
                    trial_hash = NULL,
                    check_hash = FALSE) {
  mode <- match.arg(mode)
  method <- match.arg(method)

  bundle <- .item_prepare_cv_bundle(
    Gamma = Gamma,
    T_target = T_target,
    U = U,
    run_id = run_id,
    mode = mode,
    class_levels = class_levels,
    trial_id = trial_id,
    trial_hash = trial_hash
  )

  .item_validate_bundle(
    bundle,
    require_gamma = TRUE,
    require_u = TRUE,
    check_hash = isTRUE(check_hash)
  )

  folds <- sort(unique(bundle$run_id))
  if (length(folds) < 2L) {
    stop("item_cv requires at least 2 unique runs for LOSO CV.", call. = FALSE)
  }

  class_levels <- if (identical(mode, "classification")) {
    .item_class_levels(bundle$T_target, bundle$meta$class_levels)
  } else {
    NULL
  }

  metric_name <- .item_metric_name(mode = mode, metric = metric)

  n_trials <- nrow(bundle$Gamma)
  T_hat_all <- matrix(
    NA_real_,
    nrow = n_trials,
    ncol = ncol(bundle$T_target),
    dimnames = list(bundle$trial_id, colnames(bundle$T_target))
  )

  pred_label_all <- if (identical(mode, "classification")) {
    rep(NA_character_, n_trials)
  } else {
    NULL
  }

  fold_rows <- vector("list", length(folds))

  for (i in seq_along(folds)) {
    test_run <- folds[[i]]
    fold <- item_slice_fold(bundle, test_run = test_run, check_hash = FALSE)

    W_hat <- item_fit(
      Gamma_train = fold$Gamma_train,
      T_train = fold$T_train,
      U_train = fold$U_train,
      ridge = ridge,
      method = method
    )

    T_hat <- item_predict(fold$Gamma_test, W_hat)
    T_hat_all[fold$test_idx, ] <- T_hat

    scored <- .item_score_fold(
      T_true = fold$T_test,
      T_hat = T_hat,
      mode = mode,
      metric = metric_name,
      class_levels = class_levels
    )

    if (!is.null(pred_label_all)) {
      pred_label_all[fold$test_idx] <- scored$pred_labels
    }

    fold_rows[[i]] <- data.frame(
      fold = i,
      test_run = as.character(test_run),
      n_train = nrow(fold$Gamma_train),
      n_test = nrow(fold$Gamma_test),
      metric = scored$metric,
      stringsAsFactors = FALSE
    )
  }

  fold_df <- do.call(rbind, fold_rows)
  aggregate <- list(
    metric = metric_name,
    mean = mean(fold_df$metric, na.rm = TRUE),
    sd = stats::sd(fold_df$metric, na.rm = TRUE),
    n_folds = nrow(fold_df)
  )

  predictions <- list(
    T_hat = T_hat_all,
    T_true = bundle$T_target
  )
  if (!is.null(pred_label_all)) {
    predictions$predicted_class <- pred_label_all
    predictions$true_class <- .item_true_class(bundle$T_target, class_levels)
    predictions$class_levels <- class_levels
  }

  structure(
    list(
      mode = mode,
      metric = metric_name,
      folds = fold_df,
      aggregate = aggregate,
      predictions = predictions,
      diagnostics = list(
        fold_order = as.character(folds),
        solver = method,
        ridge = ridge
      )
    ),
    class = "item_cv_result"
  )
}

#' @keywords internal
#' @noRd
.item_prepare_cv_bundle <- function(Gamma,
                                    T_target,
                                    U,
                                    run_id,
                                    mode,
                                    class_levels,
                                    trial_id,
                                    trial_hash) {
  if (inherits(Gamma, "item_bundle")) {
    bundle <- Gamma
    if (!is.null(T_target) || !is.null(U) || !is.null(run_id)) {
      warning(
        "item_cv: Gamma is an item_bundle; T_target/U/run_id arguments are ignored.",
        call. = FALSE
      )
    }

    if (identical(mode, "classification")) {
      lvl <- .item_class_levels(bundle$T_target, class_levels %||% bundle$meta$class_levels)
      bundle$meta$class_levels <- lvl
    }

    return(bundle)
  }

  Gamma <- .item_as_numeric_matrix(Gamma, "Gamma")
  n_trials <- nrow(Gamma)

  if (is.null(T_target)) {
    stop("T_target is required when Gamma is not an item_bundle.", call. = FALSE)
  }

  target_info <- .item_prepare_targets(T_target, n_trials)
  tmat <- target_info$T_target

  if (identical(mode, "classification")) {
    if (ncol(tmat) < 2L) {
      stop("Classification mode requires at least 2 target columns/classes.", call. = FALSE)
    }

    lvl <- .item_class_levels(tmat, class_levels %||% target_info$class_levels)
    colnames(tmat) <- lvl
    target_type <- "classification"
    class_meta <- lvl
  } else {
    if (ncol(tmat) < 1L) {
      stop("Regression mode requires at least one target column.", call. = FALSE)
    }
    target_type <- "regression"
    class_meta <- NULL
  }

  run_id <- .item_as_run_id(run_id, n_trials)
  trial_id <- .item_as_trial_id(trial_id, X_t = NULL, n_trials = n_trials)

  U_matrix <- NULL
  U_by_run <- NULL
  if (is.null(U)) {
    stop("U is required when Gamma is not an item_bundle.", call. = FALSE)
  }

  if (is.list(U)) {
    U_by_run <- U
  } else {
    U_matrix <- .item_as_numeric_matrix(U, "U")
  }

  trial_info <- data.frame(
    trial_index = seq_len(n_trials),
    trial_id = trial_id,
    run_id = run_id,
    stringsAsFactors = FALSE
  )

  .item_new_bundle(
    Gamma = Gamma,
    X_t = NULL,
    C_transform = NULL,
    T_target = tmat,
    U = U_matrix,
    U_by_run = U_by_run,
    run_id = run_id,
    trial_id = trial_id,
    trial_hash = trial_hash,
    trial_info = trial_info,
    meta = list(target_type = target_type, class_levels = class_meta),
    diagnostics = list()
  )
}

#' @keywords internal
#' @noRd
.item_metric_name <- function(mode, metric) {
  if (is.null(metric)) {
    return(if (identical(mode, "classification")) "accuracy" else "correlation")
  }

  metric <- tolower(as.character(metric)[1])

  if (identical(mode, "classification") && !(metric %in% c("accuracy", "balanced_accuracy"))) {
    stop("Unsupported classification metric.", call. = FALSE)
  }

  if (identical(mode, "regression") && !(metric %in% c("correlation", "rmse"))) {
    stop("Unsupported regression metric.", call. = FALSE)
  }

  metric
}

#' @keywords internal
#' @noRd
.item_class_levels <- function(T_target, class_levels = NULL) {
  if (!is.null(class_levels)) {
    class_levels <- as.character(class_levels)
    if (length(class_levels) != ncol(T_target)) {
      stop(
        sprintf("class_levels must have length %d.", ncol(T_target)),
        call. = FALSE
      )
    }
    return(class_levels)
  }

  cn <- colnames(T_target)
  if (!is.null(cn) && length(cn) == ncol(T_target)) {
    return(as.character(cn))
  }

  as.character(seq_len(ncol(T_target)))
}

#' @keywords internal
#' @noRd
.item_score_fold <- function(T_true, T_hat, mode, metric, class_levels = NULL) {
  if (identical(mode, "classification")) {
    truth_idx <- max.col(T_true, ties.method = "first")
    pred_idx <- max.col(T_hat, ties.method = "first")

    truth_labels <- class_levels[truth_idx]
    pred_labels <- class_levels[pred_idx]

    score <- if (identical(metric, "balanced_accuracy")) {
      tab <- table(
        factor(truth_labels, levels = class_levels),
        factor(pred_labels, levels = class_levels)
      )
      recalls <- diag(tab) / rowSums(tab)
      mean(recalls, na.rm = TRUE)
    } else {
      mean(pred_labels == truth_labels)
    }

    return(list(metric = score, pred_labels = pred_labels))
  }

  # regression
  score <- if (identical(metric, "rmse")) {
    sqrt(mean((T_hat - T_true)^2))
  } else {
    cors <- vapply(seq_len(ncol(T_true)), function(j) {
      x <- T_true[, j]
      y <- T_hat[, j]
      if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
      stats::cor(x, y)
    }, numeric(1))
    mean(cors, na.rm = TRUE)
  }

  list(metric = score, pred_labels = NULL)
}

#' @keywords internal
#' @noRd
.item_true_class <- function(T_target, class_levels) {
  idx <- max.col(T_target, ties.method = "first")
  class_levels[idx]
}
