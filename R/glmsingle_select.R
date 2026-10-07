# Model selection pieces of glmsingle(): compiled repeated-trial
# cross-validation, the GLMdenoise PC-count rule, and the ON-OFF R^2 tail
# threshold.

# Compile GLMsingle's pairwise repeated-trial cross-validation (calcbadness)
# into fixed per-trial weights and targets.
#
# calcbadness() z-scores betas per session using the reference fit, then for
# every fold sums (z_cand[i] - z_ref[j])^2 over training trials i and held-out
# trials j of the same condition. With W_ij counting those pairs this equals
#   sum_i d_i z_i^2 - 2 sum_i z_i M_i + sum_j w_j z_ref[j]^2,
# with d = rowSums(W), M = W z_ref and w = colSums(W) (glms_cv_compile, C++).
# Only trials with d > 0 are kept ("used").
.glms_cv_compile <- function(geom, ref, zero_sd = "zero") {
  if (is.null(geom$cond_trials)) {
    geom$cond_trials <- split(seq_len(geom$n_trials) - 1L, geom$stimorder)
  }
  test_runs <- vapply(geom$xval_scheme, function(f) seq_len(length(geom$validcolumns)) %in% f,
                      logical(length(geom$validcolumns)))
  cv <- glms_cv_compile(ref, match(geom$session, unique(geom$session))[geom$trial_run] - 1L,
                        geom$trial_run - 1L, unname(geom$cond_trials),
                        matrix(as.integer(test_runs), nrow = length(geom$validcolumns)))
  cv$used <- as.integer(cv$used)
  cv$const <- drop(cv$const)
  cv$d <- drop(cv$d)
  cv$session_used <- as.integer(cv$session_used)
  # candidates: zero-SD voxels get z = 0 ("zero"), or are divided by 1
  # ("python", reproducing GLMsingle's in-place divisor update)
  cv$isd <- if (identical(zero_sd, "python")) ifelse(cv$sd == 0, 1, cv$isd_ref) else cv$isd_ref
  cv
}

# Cross-validation loss (one value per voxel) of candidate betas given only
# on the used trial rows.
.glms_cv_loss <- function(cv, cand_used, ref = FALSE) {
  drop(glms_cv_loss_cpp(cand_used, cv$mu, if (ref) cv$isd_ref else cv$isd,
                        cv$session_used, cv$d, cv$M, cv$const, .glms_nt()))
}

# Loss of the reference fit itself. In GLMsingle the reference's own z-scores
# are always zeroed where the session SD is zero.
.glms_cv_loss_ref <- function(cv, ref_used) .glms_cv_loss(cv, ref_used, ref = TRUE)

# GLMsingle's select_noise_regressors(): first PC count whose improvement
# over 0 PCs is within a factor pcstop of the best improvement.
.glms_select_pcs <- function(xvaltrend, pcstop) {
  curve <- xvaltrend - xvaltrend[1L]
  chosen <- 0L
  best <- -Inf
  for (p in seq_along(curve)) {
    if (curve[p] > best) {
      chosen <- p - 1L
      best <- curve[p]
      if (best * pcstop >= max(curve)) break
    }
  }
  chosen
}

# GLMsingle's robustrange().
.glms_robustrange <- function(m) {
  repeat {
    absmn <- min(m)
    absmx <- max(m)
    vals <- .glms_percentile(m, c(.1, 10, 50, 90, 99.9))
    pmn <- vals[3] - 5 * (vals[3] - vals[2])
    pmx <- vals[3] + 5 * (vals[4] - vals[3])
    rerun <- FALSE
    if (vals[5] <= pmx) {
      finalmx <- if (absmx <= vals[3] + 1.1 * (vals[5] - vals[3])) absmx else vals[5]
    } else {
      rerun <- TRUE
      m <- m[!(m > pmx)]
    }
    if (vals[1] >= pmn) {
      finalmn <- if (absmn >= vals[3] - 1.1 * (vals[3] - vals[1])) absmn else vals[1]
    } else {
      rerun <- TRUE
      m <- m[!(m < pmn)]
    }
    if (!rerun) return(c(finalmn, finalmx))
  }
}

# Two-component Gaussian mixture by EM with deterministic restarts and an
# eps variance floor (GLMsingle MATLAB 91e5b7e).
.glms_gmm2 <- function(x, restarts = 3L, max_iter = 1000L, tol = 1e-10) {
  qs <- list(c(.25, .75), c(.10, .90), c(.50, .95))[seq_len(restarts)]
  best <- NULL
  for (q in qs) {
    mu <- .glms_percentile(x, 100 * q)
    s2 <- rep(stats::var(x), 2)
    w <- c(.5, .5)
    ll_old <- -Inf
    for (it in seq_len(max_iter)) {
      d1 <- w[1] * stats::dnorm(x, mu[1], sqrt(s2[1]))
      d2 <- w[2] * stats::dnorm(x, mu[2], sqrt(s2[2]))
      tot <- d1 + d2
      tot[tot == 0] <- .Machine$double.xmin
      ll <- sum(log(tot))
      r1 <- d1 / tot
      r2 <- 1 - r1
      n1 <- sum(r1); n2 <- sum(r2)
      if (n1 < 1e-8 || n2 < 1e-8) break
      w <- c(n1, n2) / length(x)
      mu <- c(sum(r1 * x) / n1, sum(r2 * x) / n2)
      s2 <- c(sum(r1 * (x - mu[1])^2) / n1, sum(r2 * (x - mu[2])^2) / n2) +
        .Machine$double.eps
      if (abs(ll - ll_old) < tol * max(1, abs(ll))) break
      ll_old <- ll
    }
    if (is.null(best) || ll > best$ll) best <- list(w = w, mu = mu, s2 = s2, ll = ll)
  }
  best
}

# GLMsingle's findtailthreshold(): the point where the posterior of the
# right-hand Gaussian falls to 0.5.
.glms_tail_threshold <- function(v, maxsz = 1e6) {
  v <- v[is.finite(v)]
  # A mixture is unidentified without variation. Use the common R^2 (zero
  # if none is finite); the strict pool/selection inequalities still apply.
  if (!length(v)) return(0)
  if (all(v == v[1L])) return(v[1L])
  if (length(v) > maxsz) v <- v[round(seq(1, length(v), length.out = maxsz))]
  fit <- .glms_gmm2(v)
  rng <- .glms_robustrange(v)
  rng[1] <- min(rng[1], fit$mu)
  rng[2] <- max(rng[2], fit$mu)
  allvals <- seq(rng[1], rng[2], length.out = 500)
  d1 <- fit$w[1] * stats::dnorm(allvals, fit$mu[1], sqrt(fit$s2[1]))
  d2 <- fit$w[2] * stats::dnorm(allvals, fit$mu[2], sqrt(fit$s2[2]))
  post <- cbind(d1, d2) / (d1 + d2)
  post[!is.finite(post)] <- 0.5
  wh <- if (post[500, 1] > 0.5) 1L else 2L
  ix <- 500L
  for (i in 500:1) {
    ix <- i
    if (post[i, wh] <= 0.5) break
  }
  allvals[ix]
}
