# Helpers for glmsingle() tests: Python reference fixtures, compact dense
# reference implementations of individual GLMsingle stages, and a small
# simulator.

glms_fixture_dir <- function() {
  testthat::test_path("fixtures", "glmsingle")
}

skip_if_no_glms_fixtures <- function() {
  if (!file.exists(file.path(glms_fixture_dir(), "index.json"))) {
    testthat::skip("GLMsingle reference fixtures not available (dev-only)")
  }
}

# Read one scenario written by tools/glmsingle_ref/make_fixtures.py.
# Arrays are C-ordered little-endian float32; returned as R arrays with the
# same dimensions.
read_glms_fixture <- function(name) {
  dir <- file.path(glms_fixture_dir(), name)
  man <- jsonlite::fromJSON(file.path(dir, "manifest.json"), simplifyVector = FALSE)
  arrays <- lapply(man$arrays, function(a) {
    shape <- unlist(a$shape)
    n <- if (length(shape)) prod(shape) else 1L
    con <- gzfile(file.path(dir, a$file), "rb")
    on.exit(close(con))
    x <- readBin(con, "double", n = n, size = 4L, endian = "little")
    if (length(shape) <= 1L) return(x)
    # C order -> R column-major: fill reversed dims, then transpose axes
    arr <- array(x, rev(shape))
    aperm(arr, rev(seq_along(shape)))
  })
  list(tr = man$tr, stimdur = man$stimdur, params = man$params,
       variants = unlist(man$variants), a = arrays)
}

# Inputs of a fixture in glmsingle() form (time x voxel data per run).
glms_fixture_inputs <- function(fx) {
  runs <- sort(as.integer(sub("data_run", "", grep("^data_run", names(fx$a), value = TRUE))))
  Y <- lapply(runs, function(r) t(fx$a[[paste0("data_run", r)]]))
  design <- lapply(runs, function(r) fx$a[[paste0("design_run", r)]])
  extras <- if (!is.null(fx$a$extras_run1)) lapply(runs, function(r) fx$a[[paste0("extras_run", r)]]) else NULL
  list(Y = Y, design = design, extras = extras)
}

# Run glmsingle() with a fixture's settings.
run_glms_fixture <- function(fx, variant = c("upstream", "fmrilss"), ...) {
  variant <- match.arg(variant)
  inp <- glms_fixture_inputs(fx)
  p <- fx$params
  args <- list(
    Y = inp$Y, design = inp$design, tr = fx$tr, stimdur = fx$stimdur,
    extra_regressors = inp$extras,
    brain_r2 = p$brainR2, pc_r2_cutoff = p$pcR2cutoff,
    verbose = FALSE
  )
  if (!is.null(p$pcstop)) args$pcstop <- p$pcstop
  if (!is.null(p$fracs)) args$fracs <- unlist(p$fracs)
  if (!is.null(p$wantlibrary)) args$want_library <- as.logical(p$wantlibrary)
  if (!is.null(p$sessionindicator)) args$session_indicator <- unlist(p$sessionindicator)
  if (!is.null(p$xvalscheme)) args$xval_scheme <- lapply(p$xvalscheme, function(f) unlist(f) + 1L)
  if (variant == "upstream") {
    args$extras_in_denoise <- "with_pcs"
    args$zero_sd_cv <- "python"
  }
  do.call(glmsingle, c(args, list(...)))
}

# Fixture array for a model type in a variant (fmrilss variant only stores
# types C and D; A and B are shared).
glms_ref <- function(fx, variant, type, key) {
  nm <- paste(variant, type, key, sep = "_")
  if (is.null(fx$a[[nm]])) nm <- paste("upstream", type, key, sep = "_")
  fx$a[[nm]]
}

rel_err <- function(a, b) {
  ok <- is.finite(a) & is.finite(b)
  sqrt(sum((a[ok] - b[ok])^2)) / max(sqrt(sum(b[ok]^2)), 1e-300)
}

# ---------------------------------------------------------------------------
# Dense float64 references, transcribed from GLMsingle's Python source. They
# build the full stacked single-trial design (every run carries a column for
# every trial), dense T x T projectors, and literal calcbadness/fracridge, so
# they share no code path with the fast implementation.

ref_olsmatrix <- function(X) {
  good <- colSums(X != 0) > 0
  f <- matrix(0, ncol(X), nrow(X))
  Xg <- X[, good, drop = FALSE]
  len <- sqrt(colSums(Xg^2))
  Xn <- sweep(Xg, 2L, len, "/")
  f[good, ] <- diag(1 / len, length(len)) %*% MASS::ginv(crossprod(Xn), tol = 1e-15) %*% t(Xn)
  f
}

ref_projection <- function(M) diag(nrow(M)) - M %*% ref_olsmatrix(M)

ref_poly <- function(n, deg) {
  t <- seq(-1, 1, length.out = n)
  P <- matrix(0, n, deg + 1L)
  for (i in 0:deg) {
    v <- t^i
    if (i > 0) v <- ref_projection(P[, seq_len(i), drop = FALSE]) %*% v
    P[, i + 1L] <- v / sqrt(sum(v^2))
  }
  P
}

# designSINGLE convolved with one HRF: list of T_r x N matrices.
ref_design_single <- function(onsets, n_time, hrf) {
  N <- sum(lengths(onsets))
  off <- 0L
  lapply(seq_along(onsets), function(r) {
    X <- matrix(0, n_time[r], N)
    for (i in seq_along(onsets[[r]])) {
      s <- numeric(n_time[r]); s[onsets[[r]][i] + 1L] <- 1
      X[, off + i] <- stats::convolve(s, rev(hrf), type = "open")[seq_len(n_time[r])]
    }
    off <<- off + length(onsets[[r]])
    X
  })
}

ref_fracridge <- function(X, y, frac) {
  XtX <- crossprod(X)
  s <- svd(XtX)
  selt <- sqrt(s$d)
  ynew <- diag(1 / selt, length(selt)) %*% t(s$v) %*% crossprod(X, y)
  ols <- ynew / selt
  ols[selt < 1e-10, ] <- 0
  val1 <- 1e4 * selt[1]^2
  val2 <- 1e-2 * selt[length(selt)]^2
  if (val2 == 0) val2 <- 1e-2
  lo <- floor(log10(val2)); hi <- ceiling(log10(val1))
  grid <- c(0, 10^(lo + 0.2 * (seq_len(ceiling((hi - lo) / 0.2)) - 1L)))
  sclg_sq <- (outer(grid, selt^2, function(g, s) s / (s + g)))^2
  coef <- matrix(0, ncol(X), ncol(y))
  for (v in seq_len(ncol(y))) {
    nl <- sqrt(sclg_sq %*% ols[, v]^2); nl <- nl / nl[1]
    a <- exp(stats::approx(rev(nl), rev(log(1 + grid)), xout = frac, rule = 2, ties = "ordered")$y) - 1
    coef[, v] <- s$v %*% ((selt^2 / (selt^2 + a)) * ols[, v])
  }
  coef
}

# glm_estimatemodel('assume') for a single-trial design: OLS (frac = NULL) or
# fracridge, with R^2 reported against poly-only projected data.
ref_fit_assume <- function(Ylist, onsets, hrf, deg, extras = NULL, frac = NULL) {
  R <- length(Ylist)
  n_time <- vapply(Ylist, nrow, integer(1))
  Xraw <- ref_design_single(onsets, n_time, hrf)
  Pp <- lapply(seq_len(R), function(r) ref_projection(ref_poly(n_time[r], deg[r])))
  Pc <- lapply(seq_len(R), function(r) {
    M <- ref_poly(n_time[r], deg[r])
    if (!is.null(extras) && !is.null(extras[[r]])) M <- cbind(M, extras[[r]])
    ref_projection(M)
  })
  Xs <- do.call(rbind, lapply(seq_len(R), function(r) Pc[[r]] %*% Xraw[[r]]))
  Ys <- do.call(rbind, lapply(seq_len(R), function(r) Pc[[r]] %*% Ylist[[r]]))
  beta <- if (is.null(frac)) {
    good <- colSums(Xs != 0) > 0
    b <- matrix(0, ncol(Xs), ncol(Ys))
    b[good, ] <- solve(crossprod(Xs[, good]), crossprod(Xs[, good], Ys))
    b
  } else ref_fracridge(Xs, Ys, frac)
  sse <- s <- matrix(0, R, ncol(Ys))
  for (r in seq_len(R)) {
    fit <- Pp[[r]] %*% (Xraw[[r]] %*% beta)
    dat <- Pp[[r]] %*% Ylist[[r]]
    sse[r, ] <- colSums((fit - dat)^2); s[r, ] <- colSums(dat^2)
  }
  list(beta = beta, R2 = 100 * (1 - colSums(sse) / colSums(s)), R2run = t(100 * (1 - sse / s)))
}

# Literal port of calcbadness(), including zerodiv's in-place divisor update
# when python = TRUE. `results` are lists of voxels x trials matrices.
ref_calcbadness <- function(xvals, validcolumns, stimix, results, session, python = TRUE) {
  V <- nrow(results[[1]])
  bad <- matrix(0, V, length(results))
  dm <- results
  for (s in unique(session)) {
    whcol <- unlist(validcolumns[session == s])
    x0 <- results[[1]][, whcol, drop = FALSE]
    mn <- rowMeans(x0)
    sd <- sqrt(rowSums((x0 - mn)^2) / (length(whcol) - 1))
    for (k in seq_along(results)) {
      d <- results[[k]][, whcol, drop = FALSE] - mn
      z0 <- sd == 0
      div <- ifelse(z0, 1, sd)
      z <- d / div
      z[z0, ] <- 0
      if (python) sd[z0] <- 1
      dm[[k]][, whcol] <- z
    }
  }
  allruns <- seq_along(validcolumns)
  for (test in xvals) {
    train <- setdiff(allruns, test)
    testcols <- unlist(validcolumns[test]); testids <- unlist(stimix[test])
    traincols <- unlist(validcolumns[train]); trainids <- unlist(stimix[train])
    for (k in seq_along(results)) {
      for (j in seq_along(testids)) {
        have <- which(trainids == testids[j])
        if (!length(have)) next
        b1 <- dm[[k]][, traincols[have], drop = FALSE]
        b2 <- dm[[1]][, testcols[j]]
        bad[, k] <- bad[, k] + rowSums((b1 - b2)^2)
      }
    }
  }
  bad
}

# ---------------------------------------------------------------------------
# Small simulator with TR-locked onsets, repeated conditions, voxel-specific
# library HRFs, shared structured noise and drift. Uses its own RNG state.
sim_glms <- function(seed = 1, n_runs = 4, n_time = 90, n_vox = 24, n_cond = 8,
                     tr = 1, stimdur = 3, isi = c(3, 5), n_extras = 0) {
  old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv()))
  set.seed(seed)
  lib <- glmsingle_hrf_library(stimdur, tr)
  hidx <- sample(ncol(lib), n_vox, replace = TRUE)
  mu <- matrix(rnorm(n_cond * n_vox, 1, 0.8), n_cond)
  load <- matrix(rnorm(3 * n_vox), 3)
  base <- runif(n_vox, 500, 1500)
  design <- Y <- extras <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    on <- integer(0); t <- sample(2:4, 1)
    while (t < n_time - 12) { on <- c(on, t); t <- t + sample(isi[1]:isi[2], 1) }
    cond <- sample(n_cond, length(on), replace = TRUE)
    D <- matrix(0, n_time, n_cond); D[cbind(on + 1L, cond)] <- 1
    S <- matrix(0, n_time, n_vox)
    for (v in seq_len(n_vox)) {
      st <- numeric(n_time); st[on + 1L] <- mu[cond, v] + rnorm(length(on), 0, 0.4)
      S[, v] <- stats::convolve(st, rev(lib[, hidx[v]]), type = "open")[seq_len(n_time)]
    }
    latent <- apply(matrix(rnorm(3 * n_time), n_time), 2, cumsum) * 0.3 + matrix(rnorm(3 * n_time), n_time)
    tt <- seq(-1, 1, length.out = n_time)
    noise <- matrix(rnorm(n_time * n_vox), n_time) + latent %*% load + outer(tt, rnorm(n_vox, 0, 2))
    Yr <- sweep(1 + 0.01 * (S + noise), 2L, base, "*")
    if (n_extras) {
      E <- apply(matrix(rnorm(n_time * n_extras), n_time), 2, cumsum) * 0.1
      Yr <- Yr + sweep(E %*% matrix(rnorm(n_extras * n_vox), n_extras), 2L, base * 0.01, "*")
      extras[[r]] <- E
    }
    design[[r]] <- D; Y[[r]] <- Yr
  }
  list(Y = Y, design = design, extras = if (n_extras) extras else NULL, tr = tr, stimdur = stimdur)
}

# Condition number of each library HRF's residualised design (max over runs).
glms_kappa <- function(fit, extras = NULL) {
  g <- fit$design
  vapply(seq_len(ncol(fit$hrf_library)), function(h) {
    max(vapply(seq_along(g$onsets), function(r) {
      nu <- .glms_nuisance_run(g$n_time[r], fit$settings$max_poly_deg[r], extras[[r]], NULL, 0L, function(k) TRUE)
      e <- eigen(.glms_design_stats(g$onsets[[r]], g$n_time[r], fit$hrf_library[, h], nu, "k0")$G,
                 symmetric = TRUE, only.values = TRUE)$values
      sqrt(max(e) / min(e))
    }, numeric(1)))
  }, numeric(1))
}
