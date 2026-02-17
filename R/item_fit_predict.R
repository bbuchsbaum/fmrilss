#' Fit ITEM decoder weights
#'
#' Fit ITEM weights with generalized least squares:
#' `W_hat = (Gamma' U^{-1} Gamma + ridge * I)^{-1} Gamma' U^{-1} T_train`.
#'
#' @param Gamma_train Numeric matrix (`n_train x n_features`).
#' @param T_train Numeric target matrix (`n_train x p`).
#' @param U_train Trial covariance (`n_train x n_train`) or run-block list.
#' @param ridge Non-negative ridge term added to the left-hand system.
#' @param method Preferred solver path (`"chol"`, `"svd"`, `"pinv"`).
#' @param tol Numerical tolerance for rank/solver fallbacks.
#'
#' @return Numeric weight matrix `W_hat` (`n_features x p`) with
#'   `item_diagnostics` attribute.
#'
#' @export
item_fit <- function(Gamma_train,
                     T_train,
                     U_train,
                     ridge = 0,
                     method = c("chol", "svd", "pinv"),
                     tol = sqrt(.Machine$double.eps)) {
  method <- match.arg(method)

  if (!is.numeric(ridge) || length(ridge) != 1L || is.na(ridge) || ridge < 0) {
    stop("ridge must be a single non-negative number.", call. = FALSE)
  }

  Gamma_train <- .item_as_numeric_matrix(Gamma_train, "Gamma_train")
  T_train <- .item_as_numeric_matrix(T_train, "T_train")

  U_train <- if (is.list(U_train)) {
    .item_block_diag(U_train)
  } else {
    .item_as_numeric_matrix(U_train, "U_train")
  }

  n_train <- nrow(Gamma_train)
  if (nrow(T_train) != n_train) {
    stop(
      sprintf("T_train must have %d rows to match Gamma_train.", n_train),
      call. = FALSE
    )
  }

  if (nrow(U_train) != n_train || ncol(U_train) != n_train) {
    stop(
      sprintf(
        "U_train must be %d x %d (got %d x %d).",
        n_train, n_train, nrow(U_train), ncol(U_train)
      ),
      call. = FALSE
    )
  }

  solve_u_gamma <- .item_safe_solve(
    A = U_train,
    B = Gamma_train,
    method = method,
    tol = tol,
    context = "item_fit(U^{-1}Gamma)"
  )
  solve_u_t <- .item_safe_solve(
    A = U_train,
    B = T_train,
    method = method,
    tol = tol,
    context = "item_fit(U^{-1}T)"
  )

  lhs <- crossprod(Gamma_train, solve_u_gamma$value)
  if (ridge > 0) {
    lhs <- lhs + diag(ridge, nrow(lhs))
  }
  rhs <- crossprod(Gamma_train, solve_u_t$value)

  fit <- .item_safe_solve(
    A = lhs,
    B = rhs,
    method = method,
    tol = tol,
    context = "item_fit(W solve)"
  )

  diagnostics <- list(
    rank = qr(lhs, tol = tol)$rank,
    condition_number = .item_condition_number(lhs),
    solver_path = fit$method,
    ridge = ridge,
    warnings = unique(c(solve_u_gamma$warnings, solve_u_t$warnings, fit$warnings))
  )

  warn_msgs <- Filter(
    f = Negate(is.null),
    x = list(
      .item_solver_warning("item_fit(U^{-1}Gamma)", method, solve_u_gamma$method, solve_u_gamma$warnings),
      .item_solver_warning("item_fit(U^{-1}T)", method, solve_u_t$method, solve_u_t$warnings),
      .item_solver_warning("item_fit(W solve)", method, fit$method, fit$warnings)
    )
  )
  if (length(warn_msgs) > 0L) {
    warning(paste(unique(warn_msgs), collapse = " "), call. = FALSE)
  }

  W_hat <- fit$value
  attr(W_hat, "item_diagnostics") <- diagnostics
  W_hat
}

#' Predict targets from ITEM weights
#'
#' Compute `T_hat = Gamma_test %*% W_hat`.
#'
#' @param Gamma_test Numeric matrix (`n_test x n_features`).
#' @param W_hat Numeric matrix (`n_features x p`).
#'
#' @return Numeric matrix of predictions (`n_test x p`).
#' @export
item_predict <- function(Gamma_test, W_hat) {
  Gamma_test <- .item_as_numeric_matrix(Gamma_test, "Gamma_test")
  W_hat <- .item_as_numeric_matrix(W_hat, "W_hat")

  if (ncol(Gamma_test) != nrow(W_hat)) {
    stop(
      sprintf(
        "ncol(Gamma_test) must equal nrow(W_hat); got %d and %d.",
        ncol(Gamma_test),
        nrow(W_hat)
      ),
      call. = FALSE
    )
  }

  Gamma_test %*% W_hat
}

#' Build an ITEM bundle from LS-A estimates
#'
#' Convenience wrapper that runs `lsa()` to estimate trial-wise amplitudes,
#' computes `U`, and returns an `item_bundle` ready for crossvalidation.
#'
#' @param Y Numeric data matrix (`n_time x n_features`).
#' @param X_t Numeric trial-wise design matrix (`n_time x n_trials`).
#' @param T_target Supervised targets with `n_trials` rows.
#' @param run_id Trial-level run/session identifier (`n_trials`).
#' @param Z Optional nuisance matrix passed to `lsa()`.
#' @param Nuisance Alias for `Z` in `lsa()`.
#' @param V Optional covariance/precision for `item_compute_u()`.
#' @param v_type Whether `V` is covariance or precision.
#' @param ridge Ridge passed to `item_compute_u()`.
#' @param lsa_method LS-A backend (`"r"` or `"cpp"`).
#' @param solver Solver preference for `item_compute_u()`.
#' @param u_output Return full `U` matrix or `U_by_run` blocks.
#' @param C_transform Optional transform matrix for `X = X_t %*% C_transform`.
#' @param trial_id Optional trial id vector.
#' @param trial_hash Optional trial hash.
#' @param meta Optional metadata list.
#' @param validate Logical; enforce strict checks before returning.
#'
#' @return `item_bundle` with fields `Gamma`, `X_t`, `T_target`, `U`/`U_by_run`,
#'   `run_id`, `meta`, and `diagnostics`.
#'
#' @export
item_from_lsa <- function(Y,
                          X_t,
                          T_target,
                          run_id,
                          Z = NULL,
                          Nuisance = NULL,
                          V = NULL,
                          v_type = c("cov", "precision"),
                          ridge = 0,
                          lsa_method = c("r", "cpp"),
                          solver = c("chol", "svd", "pinv"),
                          u_output = c("matrix", "by_run"),
                          C_transform = NULL,
                          trial_id = NULL,
                          trial_hash = NULL,
                          meta = list(),
                          validate = TRUE) {
  lsa_method <- match.arg(lsa_method)
  solver <- match.arg(solver)
  u_output <- match.arg(u_output)
  v_type <- match.arg(v_type)

  design_bundle <- item_build_design(
    X_t = X_t,
    T_target = T_target,
    run_id = run_id,
    C_transform = C_transform,
    trial_id = trial_id,
    trial_hash = trial_hash,
    meta = meta,
    diagnostics = list(),
    validate = TRUE
  )

  Gamma <- lsa(
    Y = Y,
    X = design_bundle$X_t,
    Z = Z,
    Nuisance = Nuisance,
    method = lsa_method
  )

  if (nrow(Gamma) != length(design_bundle$run_id)) {
    stop(
      paste(
        "Alignment mismatch after lsa():",
        sprintf("Gamma has %d trials but run_id has %d.", nrow(Gamma), length(design_bundle$run_id))
      ),
      call. = FALSE
    )
  }

  u_obj <- item_compute_u(
    X_t = design_bundle$X_t,
    V = V,
    v_type = v_type,
    ridge = ridge,
    method = solver,
    run_id = design_bundle$run_id,
    output = u_output
  )

  if (identical(u_output, "matrix")) {
    design_bundle$U <- u_obj
    design_bundle$U_by_run <- NULL
    u_diag <- attr(u_obj, "item_diagnostics")
  } else {
    design_bundle$U <- NULL
    design_bundle$U_by_run <- u_obj
    u_diag <- attr(u_obj, "item_diagnostics")
  }

  design_bundle$Gamma <- Gamma
  design_bundle$meta <- utils::modifyList(
    design_bundle$meta,
    list(
      lsa_method = lsa_method,
      u_solver = solver,
      v_type = v_type
    )
  )
  design_bundle$diagnostics <- utils::modifyList(
    design_bundle$diagnostics,
    list(
      lsa = list(method = lsa_method),
      u = u_diag
    )
  )

  if (isTRUE(validate)) {
    .item_validate_bundle(
      design_bundle,
      require_gamma = TRUE,
      require_u = TRUE,
      check_hash = FALSE
    )
  }

  design_bundle
}
