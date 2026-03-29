# Internal backend for overlap-aware single-trial estimation with glmnet.

.stg_or <- function(x, y) {
  if (is.null(x)) y else x
}

.stg_validate_numeric_matrix <- function(x, name) {
  if (is.null(x)) {
    return(NULL)
  }
  if (inherits(x, "Matrix") || is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix", call. = FALSE)
  }
  x
}

.stg_design_overlap <- function(X, method = c("corr", "cosine")) {
  method <- match.arg(method)
  X <- as.matrix(X)
  if (!is.numeric(X)) stop("X must be numeric.", call. = FALSE)
  if (ncol(X) < 2L) stop("X must contain at least two columns.", call. = FALSE)

  if (method == "corr") {
    M <- stats::cor(X)
  } else {
    xs <- sqrt(colSums(X * X))
    Xn <- sweep(X, 2, pmax(xs, .Machine$double.eps), "/")
    M <- crossprod(Xn)
  }
  diag(M) <- 0

  absM <- abs(M)
  score <- rowMeans(absM, na.rm = TRUE)
  max_overlap <- apply(absM, 1L, max, na.rm = TRUE)

  list(
    overlap = M,
    score = score,
    max_overlap = max_overlap,
    method = method
  )
}

.stg_overlap_penalty_factor <- function(X = NULL,
                                        overlap_score = NULL,
                                        strategy = c("multiplicative", "additive", "hybrid", "threshold"),
                                        base = 1,
                                        min_pf = 0.25,
                                        max_pf = 4,
                                        exponent = 1,
                                        strength = 1,
                                        mix = 0.5,
                                        threshold = 0.75) {
  strategy <- match.arg(strategy)
  if (is.null(overlap_score)) {
    if (is.null(X)) stop("Provide either X or overlap_score.", call. = FALSE)
    overlap_score <- .stg_design_overlap(X, method = "corr")$score
  }

  s <- as.numeric(overlap_score)
  if (!length(s)) stop("Empty overlap score.", call. = FALSE)

  s <- s - min(s, na.rm = TRUE)
  den <- max(s, na.rm = TRUE)
  if (is.finite(den) && den > 0) {
    s <- s / den
  }

  base <- as.numeric(base)
  exponent <- as.numeric(exponent)
  strength <- as.numeric(strength)
  mix <- min(1, max(0, as.numeric(mix)))
  threshold <- min(1, max(0, as.numeric(threshold)))

  pows <- s^exponent
  pf_mult <- base * (1 + strength * s)^exponent
  pf_add <- base + strength * pows
  pf_thr <- base + strength * (pmax(0, s - threshold) / pmax(1 - threshold, .Machine$double.eps))^exponent

  pf <- switch(
    strategy,
    multiplicative = pf_mult,
    additive = pf_add,
    hybrid = mix * pf_mult + (1 - mix) * pf_add,
    threshold = pf_thr
  )

  pmin(as.numeric(max_pf), pmax(as.numeric(min_pf), pf))
}

.stg_pred_array <- function(pred, nobs, nresp, nlam) {
  if (is.null(dim(pred))) stop("Unexpected prediction shape.", call. = FALSE)
  d <- dim(pred)
  if (length(d) == 2L) {
    out <- array(pred, dim = c(nobs, 1L, d[2]))
  } else if (length(d) == 3L) {
    out <- pred
  } else {
    stop("Unsupported prediction array rank.", call. = FALSE)
  }
  if (dim(out)[1] != nobs) stop("Prediction row mismatch.", call. = FALSE)
  if (dim(out)[3] != nlam) stop("Prediction lambda mismatch.", call. = FALSE)
  if (dim(out)[2] != nresp && !(nresp == 1L && dim(out)[2] == 1L)) {
    stop("Prediction response mismatch.", call. = FALSE)
  }
  out
}

.stg_pool_basis <- function(n_trial) {
  n_trial <- as.integer(n_trial)
  if (n_trial < 1L) stop("n_trial must be >= 1.", call. = FALSE)

  u <- rep(1 / sqrt(n_trial), n_trial)
  if (n_trial == 1L) {
    Q <- matrix(numeric(0), nrow = 1L, ncol = 0L)
  } else {
    Q <- stats::contr.helmert(n_trial)
    Q <- sweep(Q, 2L, sqrt(colSums(Q * Q)), "/")
  }

  list(u = u, Q = Q, T = cbind(u, Q))
}

.stg_backtransform_trial <- function(theta, pool_fit = NULL) {
  theta <- as.matrix(theta)
  if (is.null(pool_fit) || !isTRUE(pool_fit$enabled)) {
    return(theta)
  }

  Tmat <- as.matrix(pool_fit$trial_to_original)
  if (nrow(theta) != ncol(Tmat)) {
    stop("Pooling back-transform dimension mismatch.", call. = FALSE)
  }

  as.matrix(Tmat %*% theta)
}

.stg_trial_beta_matrix <- function(fit, n_trial, pool_fit = NULL) {
  nlam <- length(fit$lambda)
  if (is.list(fit$beta)) {
    nresp <- length(fit$beta)
    out <- matrix(NA_real_, nrow = n_trial * nresp, ncol = nlam)
    for (k in seq_len(nresp)) {
      bk_theta <- as.matrix(fit$beta[[k]][seq_len(n_trial), , drop = FALSE])
      bk <- .stg_backtransform_trial(bk_theta, pool_fit = pool_fit)
      i1 <- (k - 1L) * n_trial + 1L
      i2 <- k * n_trial
      out[i1:i2, ] <- bk
    }
    out
  } else {
    b_theta <- as.matrix(fit$beta[seq_len(n_trial), , drop = FALSE])
    .stg_backtransform_trial(b_theta, pool_fit = pool_fit)
  }
}

.stg_extract_trial_betas <- function(fit, n_trial, s, pool_fit = NULL) {
  if (length(s) != 1L) stop("'s' must be a single lambda value.", call. = FALSE)
  b <- stats::coef(fit, s = s)
  if (is.list(b)) {
    nresp <- length(b)
    theta <- matrix(NA_real_, nrow = n_trial, ncol = nresp)
    for (k in seq_len(nresp)) {
      theta[, k] <- as.numeric(b[[k]][1L + seq_len(n_trial), 1L, drop = TRUE])
    }
    .stg_backtransform_trial(theta, pool_fit = pool_fit)
  } else {
    theta <- matrix(as.numeric(b[1L + seq_len(n_trial), 1L, drop = TRUE]), ncol = 1L)
    .stg_backtransform_trial(theta, pool_fit = pool_fit)
  }
}

.stg_choose_lambda <- function(lambda, cvm, cvsd, maximize = FALSE) {
  if (!any(is.finite(cvm))) {
    iopt <- length(lambda)
    return(list(iopt = iopt, i1se = iopt, lambda.opt = lambda[iopt], lambda.1se = lambda[iopt]))
  }

  score <- cvm
  score[!is.finite(score)] <- if (isTRUE(maximize)) -Inf else Inf
  cvsd2 <- cvsd
  cvsd2[!is.finite(cvsd2)] <- 0

  if (maximize) {
    iopt <- which.max(score)
    thr <- score[iopt] - cvsd2[iopt]
    cand <- which(score >= thr)
    i1se <- if (length(cand)) max(cand) else iopt
    list(iopt = iopt, i1se = i1se, lambda.opt = lambda[iopt], lambda.1se = lambda[i1se])
  } else {
    iopt <- which.min(score)
    thr <- score[iopt] + cvsd2[iopt]
    cand <- which(score <= thr)
    i1se <- if (length(cand)) max(cand) else iopt
    list(iopt = iopt, i1se = i1se, lambda.opt = lambda[iopt], lambda.1se = lambda[i1se])
  }
}

.stg_parse_composite_weights <- function(w) {
  def <- c(mse = 0.4, correlation = 0.3, reliability = 0.3)
  if (is.null(w)) return(def)
  w <- as.numeric(w)
  nms <- names(w)
  if (is.null(nms)) {
    if (length(w) == length(def)) {
      names(w) <- names(def)
    } else {
      stop("'composite_weights' must be named or length 3.", call. = FALSE)
    }
  }
  out <- def
  for (nm in names(def)) {
    if (nm %in% names(w) && length(w[nm]) > 0L && is.finite(w[nm][1L])) {
      out[nm] <- max(0, as.numeric(w[nm][1L]))
    }
  }
  if (!any(out > 0)) stop("'composite_weights' must contain at least one positive weight.", call. = FALSE)
  out / sum(out)
}

.stg_metric_loss <- function(M, maximize = FALSE) {
  M <- as.matrix(M)
  out <- matrix(NA_real_, nrow = nrow(M), ncol = ncol(M))
  fin <- is.finite(M)
  if (!any(fin)) {
    return(out)
  }

  lo <- min(M[fin], na.rm = TRUE)
  hi <- max(M[fin], na.rm = TRUE)
  den <- hi - lo
  if (!is.finite(den) || den <= 0) {
    out[fin] <- 0
    return(out)
  }

  if (isTRUE(maximize)) {
    out[fin] <- (hi - M[fin]) / den
  } else {
    out[fin] <- (M[fin] - lo) / den
  }
  out
}

.stg_residual_lag1 <- function(Y, X) {
  if (nrow(Y) < 3L) {
    return(0)
  }
  qrX <- qr(X)
  resid <- qr.resid(qrX, Y)
  vals <- apply(resid, 2L, function(z) {
    suppressWarnings(stats::cor(z[-1], z[-length(z)], use = "pairwise.complete.obs"))
  })
  mean(vals, na.rm = TRUE)
}

.stg_glmnet_y <- function(Y, family) {
  if (identical(family, "gaussian") && ncol(Y) == 1L) {
    return(as.numeric(Y[, 1L]))
  }
  Y
}

.stg_fit_core <- function(y,
                          x_trial,
                          x_nuisance = NULL,
                          run_id = NULL,
                          alpha = 0.2,
                          lambda = NULL,
                          family = NULL,
                          standardize = FALSE,
                          intercept = FALSE,
                          overlap_strategy = c("none", "multiplicative", "additive", "hybrid", "threshold"),
                          overlap_strength = 1,
                          overlap_mix = 0.5,
                          overlap_threshold = 0.75,
                          overlap_exponent = 1,
                          graph_pool = FALSE,
                          graph_strength = 1,
                          graph_exponent = 1,
                          graph_mean_penalty = 0.2,
                          graph_metric = c("corr", "cosine"),
                          graph_scale_by_overlap = TRUE,
                          pool_to_mean = FALSE,
                          pool_strength = 1,
                          pool_mean_penalty = 0,
                          pool_scale_by_overlap = TRUE,
                          nuisance_penalty = 0) {
  y <- as.matrix(y)
  x_trial <- as.matrix(x_trial)
  x_nuisance <- if (is.null(x_nuisance)) NULL else as.matrix(x_nuisance)

  if (!is.numeric(y)) stop("y must be numeric.", call. = FALSE)
  if (!is.numeric(x_trial)) stop("x_trial must be numeric.", call. = FALSE)
  if (!is.null(x_nuisance) && nrow(x_nuisance) != nrow(y)) {
    stop("x_nuisance and y must have the same number of rows.", call. = FALSE)
  }

  overlap_strategy <- match.arg(as.character(overlap_strategy), c("none", "multiplicative", "additive", "hybrid", "threshold"))
  graph_metric <- match.arg(as.character(graph_metric), c("corr", "cosine"))

  trial_pf <- rep(1, ncol(x_trial))
  overlap <- NULL
  if (!identical(overlap_strategy, "none")) {
    overlap <- .stg_design_overlap(x_trial, method = "corr")
    trial_pf <- .stg_overlap_penalty_factor(
      overlap_score = overlap$score,
      strategy = overlap_strategy,
      exponent = overlap_exponent,
      strength = overlap_strength,
      mix = overlap_mix,
      threshold = overlap_threshold
    )
  }

  x_trial_fit <- x_trial
  trial_pf_fit <- as.numeric(trial_pf)
  pooling <- list(enabled = FALSE)

  if (isTRUE(graph_pool) && isTRUE(pool_to_mean)) {
    stop("Use at most one of 'graph_pool' and 'pool_to_mean'.", call. = FALSE)
  }

  if (isTRUE(graph_pool)) {
    ov_graph <- .stg_design_overlap(x_trial, method = graph_metric)
    A <- abs(as.matrix(ov_graph$overlap))
    diag(A) <- 0
    D <- diag(rowSums(A))
    L <- 0.5 * ((D - A) + t(D - A))

    if (ncol(x_trial) == 1L) {
      eigvals <- 0
      U <- matrix(1, nrow = 1L, ncol = 1L)
    } else {
      ee <- eigen(L, symmetric = TRUE)
      eigvals <- pmax(0, as.numeric(Re(ee$values)))
      U <- as.matrix(Re(ee$vectors))
      ord <- order(eigvals, decreasing = FALSE)
      eigvals <- eigvals[ord]
      U <- U[, ord, drop = FALSE]
    }

    x_trial_fit <- as.matrix(x_trial %*% U)
    colnames(x_trial_fit) <- paste0("trial_graph_", seq_len(ncol(x_trial_fit)))

    s <- eigvals - min(eigvals, na.rm = TRUE)
    den <- max(s, na.rm = TRUE)
    if (is.finite(den) && den > 0) {
      s <- s / den
    }

    overlap_scale <- if (isTRUE(graph_scale_by_overlap)) mean(as.numeric(trial_pf), na.rm = TRUE) else 1
    if (!is.finite(overlap_scale)) overlap_scale <- 1

    trial_pf_fit <- (as.numeric(graph_mean_penalty) + as.numeric(graph_strength) * (s^as.numeric(graph_exponent))) *
      as.numeric(overlap_scale)
    trial_pf_fit <- pmax(0, as.numeric(trial_pf_fit))

    pooling <- list(
      enabled = TRUE,
      method = "graph_spectral",
      trial_to_original = U,
      graph_metric = graph_metric,
      graph_eigenvalues = eigvals,
      graph_strength = as.numeric(graph_strength),
      graph_exponent = as.numeric(graph_exponent),
      graph_mean_penalty = as.numeric(graph_mean_penalty),
      overlap_scale = as.numeric(overlap_scale),
      graph_scale_by_overlap = isTRUE(graph_scale_by_overlap)
    )
  } else if (isTRUE(pool_to_mean)) {
    pbasis <- .stg_pool_basis(ncol(x_trial))
    x_mean <- matrix(as.numeric(x_trial %*% pbasis$u), ncol = 1L)
    x_contrast <- if (ncol(pbasis$Q) > 0L) {
      as.matrix(x_trial %*% pbasis$Q)
    } else {
      matrix(numeric(0), nrow = nrow(x_trial), ncol = 0L)
    }
    x_trial_fit <- cbind(x_mean, x_contrast)
    colnames(x_trial_fit) <- c("trial_mean", paste0("trial_contrast_", seq_len(ncol(x_contrast))))

    overlap_scale <- if (isTRUE(pool_scale_by_overlap)) mean(as.numeric(trial_pf), na.rm = TRUE) else 1
    if (!is.finite(overlap_scale)) overlap_scale <- 1

    trial_pf_fit <- c(as.numeric(pool_mean_penalty), rep(as.numeric(pool_strength), ncol(x_trial_fit) - 1L)) *
      as.numeric(overlap_scale)
    trial_pf_fit <- pmax(0, as.numeric(trial_pf_fit))

    pooling <- list(
      enabled = TRUE,
      method = "mean_contrast",
      trial_to_original = pbasis$T,
      mean_basis = pbasis$u,
      contrast_basis = pbasis$Q,
      overlap_scale = as.numeric(overlap_scale),
      pool_strength = as.numeric(pool_strength),
      pool_mean_penalty = as.numeric(pool_mean_penalty),
      pool_scale_by_overlap = isTRUE(pool_scale_by_overlap)
    )
  }

  x_fit <- if (is.null(x_nuisance)) x_trial_fit else cbind(x_trial_fit, x_nuisance)
  nuisance_pf <- if (is.null(x_nuisance)) numeric(0) else rep(as.numeric(nuisance_penalty), ncol(x_nuisance))
  pf <- c(trial_pf_fit, nuisance_pf)

  fit_family <- family
  if (is.null(fit_family)) {
    fit_family <- if (ncol(y) > 1L) "mgaussian" else "gaussian"
  }

  fit <- glmnet::glmnet(
    x = x_fit,
    y = .stg_glmnet_y(y, fit_family),
    family = fit_family,
    alpha = as.numeric(alpha),
    lambda = lambda,
    penalty.factor = pf,
    standardize = isTRUE(standardize),
    intercept = isTRUE(intercept)
  )

  out <- list(
    fit = fit,
    X_trial = x_trial,
    X_nuisance = x_nuisance,
    x_fit = x_fit,
    y_fit = y,
    run_id = run_id,
    overlap = overlap,
    overlap_strategy = overlap_strategy,
    overlap_strength = as.numeric(overlap_strength),
    overlap_mix = as.numeric(overlap_mix),
    overlap_threshold = as.numeric(overlap_threshold),
    overlap_exponent = as.numeric(overlap_exponent),
    graph_pool = isTRUE(graph_pool),
    graph_strength = as.numeric(graph_strength),
    graph_exponent = as.numeric(graph_exponent),
    graph_mean_penalty = as.numeric(graph_mean_penalty),
    graph_metric = graph_metric,
    graph_scale_by_overlap = isTRUE(graph_scale_by_overlap),
    pool_to_mean = isTRUE(pool_to_mean),
    pool_strength = as.numeric(pool_strength),
    pool_mean_penalty = as.numeric(pool_mean_penalty),
    pool_scale_by_overlap = isTRUE(pool_scale_by_overlap),
    pooling = pooling,
    penalty.factor = pf,
    trial_penalty.factor = trial_pf_fit,
    trial_penalty.factor_raw = as.numeric(trial_pf),
    alpha = as.numeric(alpha),
    family = fit_family,
    standardize = isTRUE(standardize),
    intercept = isTRUE(intercept)
  )
  class(out) <- "fmrilss_stglmnet_fit"
  out
}

.cv_stglmnet <- function(object,
                         foldid = NULL,
                         nfolds = 10,
                         fold_scheme = c("run", "random"),
                         type.measure = c("auto", "mse", "correlation", "reliability", "composite"),
                         overlap_low_threshold = 0.12,
                         composite_weights = c(mse = 0.4, correlation = 0.3, reliability = 0.3),
                         trace.it = 0) {
  if (!inherits(object, "fmrilss_stglmnet_fit")) {
    stop("object must inherit from 'fmrilss_stglmnet_fit'.", call. = FALSE)
  }

  fold_scheme <- match.arg(fold_scheme)
  requested.measure <- match.arg(type.measure)
  overlap_low_threshold <- as.numeric(overlap_low_threshold)
  composite_weights <- .stg_parse_composite_weights(composite_weights)

  overlap_score_mean <- function(obj) {
    if (!is.null(obj$overlap) && !is.null(obj$overlap$score)) {
      return(mean(abs(obj$overlap$score), na.rm = TRUE))
    }
    mean(abs(.stg_design_overlap(obj$X_trial, method = "corr")$score), na.rm = TRUE)
  }

  resolved.measure <- if (requested.measure == "auto") {
    ov <- overlap_score_mean(object)
    if (is.finite(ov) && ov < overlap_low_threshold) "correlation" else "reliability"
  } else {
    requested.measure
  }

  if (!is.null(foldid)) {
    foldid <- as.integer(foldid)
    if (length(foldid) != nrow(object$x_fit)) {
      stop(sprintf("foldid length (%d) must match number of observations (%d).", length(foldid), nrow(object$x_fit)), call. = FALSE)
    }
  }

  if (is.null(foldid)) {
    if (fold_scheme == "run" && !is.null(object$run_id)) {
      foldid <- as.integer(as.factor(object$run_id))
      nfolds <- length(unique(foldid))
    } else {
      n <- nrow(object$x_fit)
      foldid <- sample(rep(seq_len(nfolds), length.out = n))
    }
  }

  folds <- sort(unique(foldid))
  if (length(folds) < 2L) stop("Need at least two folds for cv.stglmnet.", call. = FALSE)

  lambda <- object$fit$lambda
  nlam <- length(lambda)
  y <- as.matrix(object$y_fit)
  nresp <- ncol(y)
  ntrial <- ncol(object$X_trial)
  cvraw <- matrix(NA_real_, nrow = length(folds), ncol = nlam)
  cvraw_mse <- if (resolved.measure == "composite") matrix(NA_real_, nrow = length(folds), ncol = nlam) else NULL
  cvraw_cor <- if (resolved.measure == "composite") matrix(NA_real_, nrow = length(folds), ncol = nlam) else NULL
  cvraw_rel <- if (resolved.measure == "composite") matrix(NA_real_, nrow = length(folds), ncol = nlam) else NULL

  beta_ref <- if (resolved.measure %in% c("reliability", "composite")) {
    .stg_trial_beta_matrix(object$fit, ntrial, pool_fit = object$pooling)
  } else {
    NULL
  }

  for (i in seq_along(folds)) {
    if (isTRUE(trace.it > 0)) cat(sprintf("cv.stglmnet fold %d/%d\n", i, length(folds)))
    test <- foldid == folds[i]
    train <- !test

    fit_i <- glmnet::glmnet(
      x = object$x_fit[train, , drop = FALSE],
      y = .stg_glmnet_y(y[train, , drop = FALSE], object$family),
      family = object$family,
      alpha = as.numeric(object$alpha),
      lambda = lambda,
      penalty.factor = object$penalty.factor,
      standardize = object$standardize,
      intercept = object$intercept
    )

    if (resolved.measure == "reliability") {
      beta_i <- .stg_trial_beta_matrix(fit_i, ntrial, pool_fit = object$pooling)
      for (j in seq_len(nlam)) {
        cvraw[i, j] <- suppressWarnings(stats::cor(beta_i[, j], beta_ref[, j], use = "pairwise.complete.obs"))
      }
      next
    }

    pred <- stats::predict(fit_i, newx = object$x_fit[test, , drop = FALSE], s = lambda, type = "response")
    parr <- .stg_pred_array(pred, sum(test), nresp, nlam)
    yte <- y[test, , drop = FALSE]

    if (resolved.measure == "mse") {
      for (j in seq_len(nlam)) {
        cvraw[i, j] <- mean((parr[, , j] - yte)^2)
      }
    } else if (resolved.measure == "correlation") {
      for (j in seq_len(nlam)) {
        cors <- rep(NA_real_, nresp)
        for (k in seq_len(nresp)) {
          cors[k] <- suppressWarnings(stats::cor(parr[, k, j], yte[, k], use = "pairwise.complete.obs"))
        }
        cvraw[i, j] <- mean(cors, na.rm = TRUE)
      }
    } else if (resolved.measure == "composite") {
      beta_i <- .stg_trial_beta_matrix(fit_i, ntrial, pool_fit = object$pooling)
      for (j in seq_len(nlam)) {
        cvraw_mse[i, j] <- mean((parr[, , j] - yte)^2)
        cors <- rep(NA_real_, nresp)
        for (k in seq_len(nresp)) {
          cors[k] <- suppressWarnings(stats::cor(parr[, k, j], yte[, k], use = "pairwise.complete.obs"))
        }
        cvraw_cor[i, j] <- mean(cors, na.rm = TRUE)
        cvraw_rel[i, j] <- suppressWarnings(stats::cor(beta_i[, j], beta_ref[, j], use = "pairwise.complete.obs"))
      }
    }
  }

  composite_meta <- NULL
  if (resolved.measure == "composite") {
    avail <- c(
      mse = any(is.finite(cvraw_mse)),
      correlation = any(is.finite(cvraw_cor)),
      reliability = any(is.finite(cvraw_rel))
    )
    w <- composite_weights
    w[!avail] <- 0
    if (!any(w > 0)) {
      return(.cv_stglmnet(
        object = object,
        foldid = foldid,
        nfolds = nfolds,
        fold_scheme = fold_scheme,
        type.measure = "correlation",
        overlap_low_threshold = overlap_low_threshold,
        composite_weights = composite_weights,
        trace.it = trace.it
      ))
    }
    w <- w / sum(w)
    L_mse <- .stg_metric_loss(cvraw_mse, maximize = FALSE)
    L_cor <- .stg_metric_loss(cvraw_cor, maximize = TRUE)
    L_rel <- .stg_metric_loss(cvraw_rel, maximize = TRUE)
    cvraw <- w["mse"] * L_mse + w["correlation"] * L_cor + w["reliability"] * L_rel
    composite_meta <- list(
      weights = as.numeric(w),
      weights_named = w,
      available = avail,
      cvm_mse = colMeans(cvraw_mse, na.rm = TRUE),
      cvm_correlation = colMeans(cvraw_cor, na.rm = TRUE),
      cvm_reliability = colMeans(cvraw_rel, na.rm = TRUE)
    )
  }

  n_eff <- colSums(is.finite(cvraw))
  cvm <- colMeans(cvraw, na.rm = TRUE)
  cvsd <- apply(cvraw, 2L, stats::sd, na.rm = TRUE) / sqrt(pmax(1, n_eff))

  maximize <- resolved.measure %in% c("correlation", "reliability")
  lam <- .stg_choose_lambda(lambda, cvm, cvsd, maximize = maximize)

  out <- list(
    lambda = lambda,
    cvm = cvm,
    cvsd = cvsd,
    cvup = cvm + cvsd,
    cvlo = cvm - cvsd,
    cvraw = cvraw,
    foldid = foldid,
    type.measure = resolved.measure,
    requested.measure = requested.measure,
    auto.measure = if (requested.measure == "auto") resolved.measure else NULL,
    fold_scheme = fold_scheme,
    glmnet.fit = object$fit,
    nfolds = length(folds)
  )
  if (!is.null(composite_meta)) out$composite <- composite_meta
  if (maximize) {
    out$lambda.max <- lam$lambda.opt
  } else {
    out$lambda.min <- lam$lambda.opt
  }
  out$lambda.1se <- lam$lambda.1se
  class(out) <- "fmrilss_cv_stglmnet"
  out
}

.stg_resolve_options <- function(stglmnet = list()) {
  defaults <- list(
    mode = "cv",
    run_id = NULL,
    alpha = 0.2,
    lambda = NULL,
    family = NULL,
    standardize = FALSE,
    intercept = FALSE,
    overlap_strategy = "none",
    overlap_strength = 1,
    overlap_mix = 0.5,
    overlap_threshold = 0.75,
    overlap_exponent = 1,
    graph_pool = FALSE,
    graph_strength = 1,
    graph_exponent = 1,
    graph_mean_penalty = 0.2,
    graph_metric = "corr",
    graph_scale_by_overlap = TRUE,
    pool_to_mean = FALSE,
    pool_strength = 1,
    pool_mean_penalty = 0,
    pool_scale_by_overlap = TRUE,
    nuisance_penalty = 0,
    whiten = "inherit",
    whiten_threshold = 0.15,
    prewhiten_args = list(method = "ar", p = 1, pooling = "global"),
    cv_folds = 5L,
    cv_foldid = NULL,
    cv_type.measure = "auto",
    cv_fold_scheme = "run",
    cv_select = "optimal",
    overlap_low_threshold = 0.12,
    composite_weights = c(mse = 0.4, correlation = 0.3, reliability = 0.3),
    return_fit = FALSE
  )

  opts <- utils::modifyList(defaults, stglmnet, keep.null = TRUE)
  opts$mode <- match.arg(as.character(opts$mode), c("cv", "fixed"))
  opts$overlap_strategy <- match.arg(as.character(opts$overlap_strategy), c("none", "multiplicative", "additive", "hybrid", "threshold"))
  opts$graph_metric <- match.arg(as.character(opts$graph_metric), c("corr", "cosine"))
  opts$whiten <- match.arg(as.character(opts$whiten), c("inherit", "auto", "never", "always"))
  opts$cv_type.measure <- match.arg(as.character(opts$cv_type.measure), c("auto", "mse", "correlation", "reliability", "composite"))
  opts$cv_fold_scheme <- match.arg(as.character(opts$cv_fold_scheme), c("run", "random"))
  opts$cv_select <- match.arg(as.character(opts$cv_select), c("optimal", "1se"))
  opts$cv_folds <- as.integer(opts$cv_folds)
  opts$cv_foldid <- if (is.null(opts$cv_foldid)) NULL else as.integer(opts$cv_foldid)
  opts$composite_weights <- .stg_parse_composite_weights(opts$composite_weights)
  opts
}

.lss_stglmnet <- function(Y,
                          X,
                          Z = NULL,
                          Nuisance = NULL,
                          stglmnet = list(),
                          prewhiten = NULL) {
  opts <- .stg_resolve_options(stglmnet)

  Y <- .stg_validate_numeric_matrix(Y, "Y")
  X <- .stg_validate_numeric_matrix(X, "X")
  Z <- .stg_validate_numeric_matrix(Z, "Z")
  Nuisance <- .stg_validate_numeric_matrix(Nuisance, "Nuisance")

  if (is.null(X)) {
    stop("For method='stglmnet', X must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(Y) != nrow(X)) {
    stop("Y and X must have the same number of rows (timepoints)", call. = FALSE)
  }
  if (!is.null(Z) && nrow(Z) != nrow(Y)) {
    stop("Z must be a numeric matrix with the same number of rows as Y", call. = FALSE)
  }
  if (!is.null(Nuisance) && nrow(Nuisance) != nrow(Y)) {
    stop("Nuisance must be a numeric matrix with the same number of rows as Y", call. = FALSE)
  }

  if (is.null(Z)) {
    Z <- matrix(1, nrow(Y), 1L)
    colnames(Z) <- "Intercept"
  }

  run_id <- .stg_or(opts$run_id, if (is.null(prewhiten)) NULL else prewhiten$runs)
  prewhiten_spec <- prewhiten
  if (is.null(prewhiten_spec) || is.null(prewhiten_spec$method) || identical(prewhiten_spec$method, "none")) {
    if (!identical(opts$whiten, "inherit")) {
      do_whiten <- FALSE
      if (identical(opts$whiten, "always")) {
        do_whiten <- TRUE
      } else if (identical(opts$whiten, "auto")) {
        X_full <- cbind(X, Z, Nuisance)
        do_whiten <- is.finite(.stg_residual_lag1(Y, X_full)) &&
          .stg_residual_lag1(Y, X_full) >= as.numeric(opts$whiten_threshold)
      }
      if (isTRUE(do_whiten)) {
        prewhiten_spec <- utils::modifyList(list(method = "ar", p = 1, pooling = "global"), opts$prewhiten_args)
        if (is.null(prewhiten_spec$runs) && !is.null(run_id)) {
          prewhiten_spec$runs <- run_id
        }
      }
    }
  }

  if (!is.null(prewhiten_spec) && !is.null(prewhiten_spec$method) && prewhiten_spec$method != "none") {
    whitened <- .prewhiten_data(Y, X, Z, Nuisance, prewhiten_spec)
    Y <- whitened$Y_whitened
    X <- whitened$X_whitened
    if (!is.null(whitened$Z_whitened)) Z <- whitened$Z_whitened
    if (!is.null(whitened$Nuisance_whitened)) Nuisance <- whitened$Nuisance_whitened
  }

  X_nuisance <- cbind(Z, Nuisance)
  proj <- .project_out_nuisance(Y, X, X_nuisance)
  Y_use <- proj$Y_residual
  X_use <- proj$X_residual

  fit_obj <- .stg_fit_core(
    y = Y_use,
    x_trial = X_use,
    x_nuisance = Z,
    run_id = run_id,
    alpha = opts$alpha,
    lambda = opts$lambda,
    family = opts$family,
    standardize = opts$standardize,
    intercept = opts$intercept,
    overlap_strategy = opts$overlap_strategy,
    overlap_strength = opts$overlap_strength,
    overlap_mix = opts$overlap_mix,
    overlap_threshold = opts$overlap_threshold,
    overlap_exponent = opts$overlap_exponent,
    graph_pool = opts$graph_pool,
    graph_strength = opts$graph_strength,
    graph_exponent = opts$graph_exponent,
    graph_mean_penalty = opts$graph_mean_penalty,
    graph_metric = opts$graph_metric,
    graph_scale_by_overlap = opts$graph_scale_by_overlap,
    pool_to_mean = opts$pool_to_mean,
    pool_strength = opts$pool_strength,
    pool_mean_penalty = opts$pool_mean_penalty,
    pool_scale_by_overlap = opts$pool_scale_by_overlap,
    nuisance_penalty = opts$nuisance_penalty
  )

  cv_obj <- NULL
  s <- NULL
  if (identical(opts$mode, "cv")) {
    cv_obj <- .cv_stglmnet(
      object = fit_obj,
      foldid = opts$cv_foldid,
      nfolds = opts$cv_folds,
      fold_scheme = opts$cv_fold_scheme,
      type.measure = opts$cv_type.measure,
      overlap_low_threshold = opts$overlap_low_threshold,
      composite_weights = opts$composite_weights
    )
    if (identical(opts$cv_select, "1se")) {
      s <- cv_obj$lambda.1se
    } else if (!is.null(cv_obj$lambda.min)) {
      s <- cv_obj$lambda.min
    } else if (!is.null(cv_obj$lambda.max)) {
      s <- cv_obj$lambda.max
    } else {
      s <- cv_obj$lambda.1se
    }
  } else {
    if (!is.null(opts$lambda) && length(opts$lambda) == 1L) {
      s <- as.numeric(opts$lambda)
    } else {
      lam_finite <- fit_obj$fit$lambda[is.finite(fit_obj$fit$lambda)]
      s <- if (length(lam_finite)) utils::tail(lam_finite, 1L) else NA_real_
    }
  }

  if (!(length(s) == 1L && is.finite(s))) {
    lam_finite <- fit_obj$fit$lambda[is.finite(fit_obj$fit$lambda)]
    if (!length(lam_finite)) stop("Could not resolve a finite lambda.", call. = FALSE)
    s <- utils::tail(lam_finite, 1L)
  }

  beta <- .stg_extract_trial_betas(fit_obj$fit, n_trial = ncol(X), s = s, pool_fit = fit_obj$pooling)

  if (!is.null(colnames(X))) {
    rownames(beta) <- colnames(X)
  } else {
    rownames(beta) <- paste0("Trial_", seq_len(nrow(beta)))
  }
  if (!is.null(colnames(Y))) {
    colnames(beta) <- colnames(Y)[seq_len(ncol(beta))]
  } else {
    colnames(beta) <- paste0("Voxel_", seq_len(ncol(beta)))
  }

  if (!isTRUE(opts$return_fit)) {
    return(beta)
  }

  structure(
    list(beta = beta, fit = fit_obj, cv = cv_obj, lambda = s, mode = opts$mode),
    class = "fmrilss_stglmnet_result"
  )
}
