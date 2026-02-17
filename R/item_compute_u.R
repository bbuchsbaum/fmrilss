#' Compute ITEM trial covariance matrix
#'
#' Compute the trial covariance term
#' `U = (X_t' V^{-1} X_t + ridge * I)^{-1}`
#' using stable solve paths with fallbacks.
#'
#' @param X_t Numeric trial-wise design matrix (`n_time x n_trials`).
#' @param V Optional temporal covariance/precision object.
#'   Accepted forms:
#'   - `NULL` (identity),
#'   - dense/sparse matrix (`n_time x n_time`),
#'   - run-block list of square matrices whose block sizes sum to `n_time`.
#' @param v_type Whether `V` is covariance (`"cov"`) or precision (`"precision"`).
#' @param ridge Non-negative ridge term added to `X_t' V^{-1} X_t`.
#' @param method Preferred solver path (`"chol"`, `"svd"`, or `"pinv"`).
#' @param run_id Optional trial-level run ids (length `n_trials`). Required only
#'   when `output = "by_run"`.
#' @param output Return full `U` (`"matrix"`) or split run blocks (`"by_run"`).
#' @param tol Numerical tolerance for rank/solver fallbacks.
#'
#' @return Numeric matrix `U` by default, or a named list `U_by_run` when
#'   `output = "by_run"`. Returned object includes `item_diagnostics` attribute
#'   with rank/condition/solver details.
#'
#' @export
item_compute_u <- function(X_t,
                           V = NULL,
                           v_type = c("cov", "precision"),
                           ridge = 0,
                           method = c("chol", "svd", "pinv"),
                           run_id = NULL,
                           output = c("matrix", "by_run"),
                           tol = sqrt(.Machine$double.eps)) {
  v_type <- match.arg(v_type)
  method <- match.arg(method)
  output <- match.arg(output)

  if (!is.numeric(ridge) || length(ridge) != 1L || is.na(ridge) || ridge < 0) {
    stop("ridge must be a single non-negative number.", call. = FALSE)
  }

  X_t <- .item_as_numeric_matrix(X_t, "X_t")

  x_vinv <- .item_apply_vinv(
    V = V,
    X = X_t,
    v_type = v_type,
    method = method,
    tol = tol
  )

  xt_vinv_x <- crossprod(X_t, x_vinv)
  if (ridge > 0) {
    xt_vinv_x <- xt_vinv_x + diag(ridge, nrow(xt_vinv_x))
  }

  inv_fit <- .item_safe_solve(
    A = xt_vinv_x,
    B = NULL,
    method = method,
    tol = tol,
    context = "item_compute_u"
  )

  U <- inv_fit$value
  U <- (U + t(U)) / 2

  diagnostics <- list(
    rank = qr(xt_vinv_x, tol = tol)$rank,
    condition_number = .item_condition_number(xt_vinv_x),
    solver_path = inv_fit$method,
    ridge = ridge,
    warnings = inv_fit$warnings,
    v_type = v_type,
    n_time = nrow(X_t),
    n_trials = ncol(X_t)
  )

  fallback_warning <- .item_solver_warning(
    context = "item_compute_u",
    requested = method,
    used = inv_fit$method,
    warnings = inv_fit$warnings
  )
  if (!is.null(fallback_warning)) {
    warning(fallback_warning, call. = FALSE)
  }

  if (identical(output, "by_run")) {
    run_id <- .item_as_run_id(run_id, ncol(X_t))
    U_by_run <- .item_split_u_by_run(U, run_id)
    attr(U_by_run, "item_diagnostics") <- diagnostics
    return(U_by_run)
  }

  attr(U, "item_diagnostics") <- diagnostics
  U
}

#' @keywords internal
#' @noRd
.item_apply_vinv <- function(V, X, v_type, method, tol) {
  if (is.null(V)) return(X)

  if (is.list(V)) {
    return(.item_apply_vinv_blocks(V, X, v_type, method, tol))
  }

  V <- .item_as_numeric_matrix(V, "V")
  if (nrow(V) != nrow(X) || ncol(V) != nrow(X)) {
    stop(
      sprintf(
        "V must be %d x %d to match n_time, got %d x %d.",
        nrow(X), nrow(X), nrow(V), ncol(V)
      ),
      call. = FALSE
    )
  }

  if (identical(v_type, "precision")) {
    return(V %*% X)
  }

  .item_safe_solve(
    A = V,
    B = X,
    method = method,
    tol = tol,
    context = "item_compute_u(V solve)"
  )$value
}

#' @keywords internal
#' @noRd
.item_apply_vinv_blocks <- function(V_blocks, X, v_type, method, tol) {
  if (length(V_blocks) == 0L) {
    stop("V block list is empty.", call. = FALSE)
  }

  starts <- integer(length(V_blocks))
  ends <- integer(length(V_blocks))
  offset <- 0L

  for (i in seq_along(V_blocks)) {
    block <- .item_as_numeric_matrix(V_blocks[[i]], sprintf("V[[%d]]", i))
    if (nrow(block) != ncol(block)) {
      stop(sprintf("V[[%d]] must be square.", i), call. = FALSE)
    }

    starts[[i]] <- offset + 1L
    ends[[i]] <- offset + nrow(block)
    offset <- ends[[i]]

    V_blocks[[i]] <- block
  }

  if (offset != nrow(X)) {
    stop(
      sprintf(
        "Sum of V block sizes (%d) must equal nrow(X_t) (%d).",
        offset,
        nrow(X)
      ),
      call. = FALSE
    )
  }

  out <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (i in seq_along(V_blocks)) {
    idx <- starts[[i]]:ends[[i]]
    Xi <- X[idx, , drop = FALSE]
    Vi <- V_blocks[[i]]

    out[idx, ] <- if (identical(v_type, "precision")) {
      Vi %*% Xi
    } else {
      .item_safe_solve(
        A = Vi,
        B = Xi,
        method = method,
        tol = tol,
        context = sprintf("item_compute_u(V block %d solve)", i)
      )$value
    }
  }

  out
}

#' @keywords internal
#' @noRd
.item_safe_solve <- function(A,
                             B = NULL,
                             method = c("chol", "svd", "pinv"),
                             tol = sqrt(.Machine$double.eps),
                             context = "solve") {
  methods <- unique(c(method, "chol", "svd", "pinv"))
  warnings <- character(0)

  for (m in methods) {
    attempt <- try(.item_solve_once(A = A, B = B, method = m, tol = tol), silent = TRUE)
    if (!inherits(attempt, "try-error")) {
      return(list(value = attempt, method = m, warnings = warnings))
    }

    warnings <- c(warnings, sprintf("%s failed with method '%s'.", context, m))
  }

  stop(sprintf("All solver paths failed in %s.", context), call. = FALSE)
}

#' @keywords internal
#' @noRd
.item_solve_once <- function(A, B = NULL, method = c("chol", "svd", "pinv"), tol = sqrt(.Machine$double.eps)) {
  method <- match.arg(method)

  if (identical(method, "chol")) {
    R <- chol(A)
    d <- abs(diag(R))
    if (length(d) == 0L || any(!is.finite(d))) {
      stop("Cholesky factor contains invalid diagonal values.")
    }
    if (max(d) == 0 || min(d) <= (tol * max(d))) {
      stop("Cholesky factor indicates rank-deficient or near-singular system.")
    }
    if (is.null(B)) {
      return(chol2inv(R))
    }

    B <- as.matrix(B)
    return(backsolve(R, forwardsolve(t(R), B)))
  }

  if (identical(method, "svd")) {
    s <- svd(A)
    if (length(s$d) == 0L) stop("Empty singular value decomposition.")

    d_max <- max(s$d)
    keep <- s$d > (d_max * tol)
    if (!any(keep)) stop("No stable singular values above tolerance.")

    U_k <- s$u[, keep, drop = FALSE]
    V_k <- s$v[, keep, drop = FALSE]
    inv_d <- 1 / s$d[keep]

    if (is.null(B)) {
      return(V_k %*% (inv_d * t(U_k)))
    }

    B <- as.matrix(B)
    return(V_k %*% (inv_d * (t(U_k) %*% B)))
  }

  # pinv
  A_pinv <- MASS::ginv(A, tol = tol)
  if (is.null(B)) {
    return(A_pinv)
  }

  A_pinv %*% as.matrix(B)
}

#' @keywords internal
#' @noRd
.item_condition_number <- function(A) {
  out <- try(kappa(A, exact = FALSE), silent = TRUE)
  if (inherits(out, "try-error") || is.na(out) || !is.finite(out)) {
    return(Inf)
  }
  as.numeric(out)
}

#' @keywords internal
#' @noRd
.item_solver_warning <- function(context, requested, used, warnings) {
  changed <- !identical(requested, used)
  has_warnings <- length(warnings) > 0L

  if (!changed && !has_warnings) return(NULL)

  details <- character(0)
  if (changed) {
    details <- c(details, sprintf("requested '%s' but used '%s'", requested, used))
  }
  if (has_warnings) {
    details <- c(details, paste(unique(warnings), collapse = " "))
  }

  sprintf("%s: %s", context, paste(details, collapse = "; "))
}

#' @keywords internal
#' @noRd
.item_split_u_by_run <- function(U, run_id) {
  run_values <- sort(unique(run_id))
  out <- vector("list", length(run_values))
  names(out) <- as.character(run_values)

  for (i in seq_along(run_values)) {
    r <- run_values[[i]]
    idx <- which(run_id == r)
    out[[i]] <- U[idx, idx, drop = FALSE]
  }

  out
}

#' @keywords internal
#' @noRd
.item_block_diag <- function(blocks) {
  if (length(blocks) == 0L) {
    return(matrix(0, nrow = 0, ncol = 0))
  }

  dims <- vapply(blocks, nrow, integer(1))
  total <- sum(dims)
  out <- matrix(0, nrow = total, ncol = total)

  offset <- 0L
  for (i in seq_along(blocks)) {
    block <- .item_as_numeric_matrix(blocks[[i]], sprintf("blocks[[%d]]", i))
    n <- nrow(block)
    if (n != ncol(block)) {
      stop(sprintf("blocks[[%d]] must be square.", i), call. = FALSE)
    }

    idx <- (offset + 1L):(offset + n)
    out[idx, idx] <- block
    offset <- offset + n
  }

  out
}
