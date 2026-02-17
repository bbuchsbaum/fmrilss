#' Slice an ITEM bundle into train/test fold objects
#'
#' Create deterministic leave-one-run-out slices for `Gamma`, `T_target`, and
#' trial covariance (`U` or `U_by_run`).
#'
#' @param bundle Object of class `item_bundle`.
#' @param test_run Run/session id to hold out for testing.
#' @param check_hash Logical; if `TRUE`, validate stored `trial_hash`.
#'
#' @return A list with train/test slices:
#' `Gamma_train`, `Gamma_test`, `T_train`, `T_test`,
#' `U_train`, `U_test`, `train_idx`, `test_idx`, `train_runs`, `test_run`.
#'
#' @export
item_slice_fold <- function(bundle, test_run, check_hash = FALSE) {
  .item_validate_bundle(
    bundle,
    require_gamma = TRUE,
    require_u = TRUE,
    check_hash = isTRUE(check_hash)
  )

  run_id <- bundle[["run_id"]]
  Gamma <- bundle[["Gamma"]]
  T_target <- bundle[["T_target"]]
  U <- bundle[["U"]]
  U_by_run <- bundle[["U_by_run"]]
  trial_id <- bundle[["trial_id"]]

  run_values <- sort(unique(run_id))
  if (!(test_run %in% run_values)) {
    stop(
      sprintf("test_run '%s' not found in run_id.", as.character(test_run)),
      call. = FALSE
    )
  }

  test_idx <- which(run_id == test_run)
  train_idx <- which(run_id != test_run)

  if (length(test_idx) == 0L || length(train_idx) == 0L) {
    stop("Fold split must contain both train and test trials.", call. = FALSE)
  }

  Gamma_train <- Gamma[train_idx, , drop = FALSE]
  Gamma_test <- Gamma[test_idx, , drop = FALSE]

  T_train <- T_target[train_idx, , drop = FALSE]
  T_test <- T_target[test_idx, , drop = FALSE]

  if (!is.null(U)) {
    U_train <- U[train_idx, train_idx, drop = FALSE]
    U_test <- U[test_idx, test_idx, drop = FALSE]
  } else {
    train_runs <- run_values[run_values != test_run]
    train_blocks <- lapply(train_runs, function(r) {
      .item_u_block_for_run(U_by_run, r, run_values)
    })

    U_train <- .item_block_diag(train_blocks)
    U_test <- .item_u_block_for_run(U_by_run, test_run, run_values)
  }

  list(
    Gamma_train = Gamma_train,
    Gamma_test = Gamma_test,
    T_train = T_train,
    T_test = T_test,
    U_train = U_train,
    U_test = U_test,
    train_idx = train_idx,
    test_idx = test_idx,
    train_runs = run_values[run_values != test_run],
    test_run = test_run,
    trial_id_train = trial_id[train_idx],
    trial_id_test = trial_id[test_idx]
  )
}

#' @keywords internal
#' @noRd
.item_u_block_for_run <- function(U_by_run, run_value, run_values) {
  has_names <- !is.null(names(U_by_run)) && all(names(U_by_run) != "")

  if (has_names) {
    key <- as.character(run_value)
    block <- U_by_run[[key]]
    if (is.null(block)) {
      stop(sprintf("Missing U_by_run block for run '%s'.", key), call. = FALSE)
    }
    return(.item_as_numeric_matrix(block, sprintf("U_by_run[['%s']]", key)))
  }

  idx <- match(run_value, run_values)
  if (is.na(idx) || idx < 1L || idx > length(U_by_run)) {
    stop(
      sprintf("Could not map run '%s' into unnamed U_by_run list.", as.character(run_value)),
      call. = FALSE
    )
  }

  .item_as_numeric_matrix(U_by_run[[idx]], sprintf("U_by_run[[%d]]", idx))
}
