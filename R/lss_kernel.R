#' LSS weight-matrix kernel
#'
#' Every closed-form LSS estimate is a linear functional of the data: for
#' trial i, beta_i = w_i' y for a weight vector w_i that depends only on the
#' design. Building the n x T weight matrix W once and computing
#' `crossprod(W, Y)` therefore yields all trial betas with a single
#' matrix product, without residualizing Y. Because each w_i lies in the
#' residual space of the confounds, `crossprod(W, Y)` equals
#' `crossprod(W, Q %*% Y)` exactly, so confound projection of the (large)
#' data matrix is never needed.
#'
#' @name lss_kernel
#' @keywords internal
#' @noRd
NULL

#' Residualize trial regressors against confounds
#'
#' Uses a rank-revealing QR basis so duplicated or collinear confound columns
#' do not over-project.
#'
#' @param C Trial design (n x T).
#' @param Xc Confound design (n x p) or NULL.
#' @return Residualized trial design (n x T).
#' @keywords internal
#' @noRd
.lss_residualize_trials <- function(C, Xc) {
  if (is.null(Xc) || ncol(Xc) == 0L) return(C)
  qrX <- qr(Xc)
  if (qrX$rank == 0L) return(C)
  basis <- qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]
  C - basis %*% crossprod(basis, C)
}

#' Normalize trial group labels to integer codes
#'
#' @param trial_groups NULL or a vector of length T.
#' @param n_trials Number of trials.
#' @return NULL (single pooled "other" regressor) or an integer vector of
#'   codes in 1..G.
#' @keywords internal
#' @noRd
.lss_group_codes <- function(trial_groups, n_trials) {
  if (is.null(trial_groups)) return(NULL)
  if (is.factor(trial_groups)) trial_groups <- as.character(trial_groups)
  if (!is.atomic(trial_groups) || length(trial_groups) != n_trials) {
    stop("trial_groups must be an atomic vector with one label per trial (ncol(X))",
         call. = FALSE)
  }
  if (anyNA(trial_groups)) stop("trial_groups must not contain NA", call. = FALSE)
  codes <- match(trial_groups, unique(trial_groups))
  if (max(codes) == 1L) return(NULL)
  as.integer(codes)
}

#' Build the LSS weight matrix
#'
#' For the classic LSS model (one pooled "other trials" regressor) the
#' weights have the closed form
#' `w_i = ((1 + a_i) c_i - a_i t) / d_i`, with `t` the sum of all trial
#' regressors, `a_i = c_i'b_i / b_i'b_i`, `d_i = c_i'c_i - (c_i'b_i)^2 / b_i'b_i`
#' and `b_i = t - c_i`.
#'
#' With `groups` (LSS-N; Turner et al., 2012) the nuisance block for trial i
#' holds one summed regressor per trial group, with trial i removed from its
#' own group. The weights solve a (G + 1)-dimensional normal system per trial
#' and are again assembled into W, so the data are touched only once.
#'
#' @param C Residualized trial design (n x T).
#' @param groups NULL or integer group codes (length T, values 1..G).
#' @param eps Numerical tolerance.
#' @return Weight matrix W (n x T) with `crossprod(W, Y)` giving T x V betas.
#' @keywords internal
#' @noRd
.lss_weight_matrix <- function(C, groups = NULL, eps = 1e-12) {
  n_trials <- ncol(C)
  if (n_trials == 1L) {
    cc <- sum(C^2)
    return(if (cc <= eps) C * 0 else C / cc)
  }
  if (!is.null(groups)) {
    return(.lss_weight_matrix_grouped(C, groups, eps))
  }

  total <- rowSums(C)
  ss_tot <- sum(total^2)
  CtC <- colSums(C^2)
  CtT <- drop(crossprod(C, total))

  bt2 <- ss_tot - 2 * CtT + CtC
  ctb <- CtT - CtC
  bt2[bt2 < eps] <- Inf
  alpha <- ctb / bt2
  den <- pmax(CtC - ctb^2 / bt2, eps)

  # W = C diag((1 + alpha) / den) - t (alpha / den)'
  sweep(C, 2L, (1 + alpha) / den, `*`) - tcrossprod(total, alpha / den)
}

#' @rdname dot-lss_weight_matrix
#' @keywords internal
#' @noRd
.lss_weight_matrix_grouped <- function(C, groups, eps = 1e-12) {
  n_trials <- ncol(C)
  G <- max(groups)
  A <- vapply(seq_len(G), function(g) {
    rowSums(C[, groups == g, drop = FALSE])
  }, numeric(nrow(C)))
  A <- matrix(A, nrow(C), G)

  AtA <- crossprod(A)
  AtC <- crossprod(A, C)              # G x T
  CtC <- colSums(C^2)

  s <- numeric(n_trials)              # coefficient on c_i
  U <- matrix(0, G, n_trials)         # coefficients on group sums
  for (i in seq_len(n_trials)) {
    g <- groups[i]
    # Gram of B_i = A - c_i e_g'
    a_c <- AtC[, i]
    M_bb <- AtA
    M_bb[g, ] <- M_bb[g, ] - a_c
    M_bb[, g] <- M_bb[, g] - a_c
    M_bb[g, g] <- M_bb[g, g] + CtC[i]
    m_cb <- a_c
    m_cb[g] <- m_cb[g] - CtC[i]

    # Drop group columns that are empty once trial i is removed.
    keep <- diag(M_bb) > eps
    if (!any(keep)) {
      s[i] <- if (CtC[i] > eps) 1 / CtC[i] else 0
      next
    }
    Mk <- M_bb[keep, keep, drop = FALSE]
    mk <- m_cb[keep]
    # Schur complement of the nuisance block gives beta_i directly:
    #   beta_i = (c' - m' Mk^-1 B') y / (c'c - m' Mk^-1 m)
    h <- .lss_small_solve(Mk, mk)
    den <- max(CtC[i] - sum(mk * h), eps)
    hk <- numeric(G)
    hk[keep] <- h
    # w_i = (c_i - B_i h) / den,  B_i h = A h - c_i h_g
    s[i] <- (1 + hk[g]) / den
    U[, i] <- -hk / den
  }
  sweep(C, 2L, s, `*`) + A %*% U
}

#' Solve a small symmetric positive semi-definite system robustly
#' @keywords internal
#' @noRd
.lss_small_solve <- function(M, b) {
  R <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(R)) {
    return(backsolve(R, forwardsolve(t(R), b)))
  }
  ev <- eigen(M, symmetric = TRUE)
  tol <- max(dim(M)) * max(abs(ev$values)) * .Machine$double.eps
  inv_vals <- ifelse(ev$values > tol, 1 / ev$values, 0)
  drop(ev$vectors %*% (inv_vals * crossprod(ev$vectors, b)))
}

#' LSS betas via the weight-matrix kernel
#'
#' @param Y Data (n x V); need not be residualized.
#' @param C Trial design (n x T).
#' @param Xc Confounds (n x p) or NULL.
#' @param groups NULL or integer group codes.
#' @return T x V beta matrix.
#' @keywords internal
#' @noRd
.lss_kernel_r <- function(Y, C, Xc, groups = NULL) {
  C_res <- .lss_residualize_trials(C, Xc)
  W <- .lss_weight_matrix(C_res, groups)
  crossprod(W, Y)
}
