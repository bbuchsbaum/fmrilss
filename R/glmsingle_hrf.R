# HRFs used by glmsingle(): the GLMsingle canonical HRF and HRF library,
# reproduced from GLMsingle's getcanonicalhrf() / getcanonicalhrflibrary().

# Monotone piecewise-cubic Hermite interpolation matching
# scipy.interpolate.PchipInterpolator (with extrapolation from the end pieces).
.glms_pchip <- function(x, y, xout) {
  n <- length(x)
  h <- diff(x)
  m <- diff(y) / h
  d <- numeric(n)
  if (n == 2L) {
    d[] <- m[1L]
  } else {
    m0 <- m[-length(m)]
    m1 <- m[-1L]
    h0 <- h[-length(h)]
    h1 <- h[-1L]
    flat <- sign(m1) != sign(m0) | m1 == 0 | m0 == 0
    w1 <- 2 * h1 + h0
    w2 <- h1 + 2 * h0
    whmean <- (w1 / m0 + w2 / m1) / (w1 + w2)
    d[2:(n - 1L)] <- ifelse(flat, 0, 1 / whmean)
    edge <- function(h0, h1, m0, m1) {
      dd <- ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
      if (sign(dd) != sign(m0)) {
        dd <- 0
      } else if (sign(m0) != sign(m1) && abs(dd) > 3 * abs(m0)) {
        dd <- 3 * m0
      }
      dd
    }
    d[1L] <- edge(h[1L], h[2L], m[1L], m[2L])
    d[n] <- edge(h[n - 1L], h[n - 2L], m[n - 1L], m[n - 2L])
  }
  k <- findInterval(xout, x, rightmost.closed = TRUE, all.inside = TRUE)
  t <- xout - x[k]
  hk <- h[k]
  # Hermite form written as PPoly coefficients (as scipy does).
  c0 <- y[k]
  c1 <- d[k]
  c2 <- (3 * m[k] - 2 * d[k] - d[k + 1L]) / hk
  c3 <- (d[k] + d[k + 1L] - 2 * m[k]) / hk^2
  c0 + t * (c1 + t * (c2 + t * c3))
}

.glms_extdata <- function(file) {
  path <- system.file("extdata", file, package = "fmrilss")
  if (!nzchar(path)) stop("Missing package data file: ", file, call. = FALSE)
  path
}

# Convolve a 0.1 s HRF with a stimulus boxcar and resample to the TR.
.glms_hrf_resample <- function(hrf01, duration, tr, sampler) {
  if (duration == 0) duration <- 0.1
  box <- rep(1, max(1L, .glms_alt_round(duration / 0.1)))
  conv <- stats::convolve(hrf01, rev(box), type = "open")
  .glms_pchip((seq_along(conv) - 1) * 0.1, conv, sampler(length(conv)))
}

#' GLMsingle canonical HRF and HRF library
#'
#' Reproduces GLMsingle's canonical HRF (`getcanonicalhrf`) and its library
#' of 20 HRFs (`getcanonicalhrflibrary`) for a stimulus of duration `stimdur`
#' sampled at `tr`. Each HRF is peak-normalised to 1, and the first sample is
#' coincident with stimulus onset.
#'
#' @param stimdur Stimulus duration in seconds (rounded to 0.1 s).
#' @param tr Repetition time in seconds.
#' @return `glmsingle_hrf()` returns a numeric vector; `glmsingle_hrf_library()`
#'   returns a time x 20 matrix.
#' @references Prince, J. S., et al. (2022). Improving the accuracy of
#'   single-trial fMRI response estimates using GLMsingle. eLife, 11, e77599.
#' @export
#' @examples
#' lib <- glmsingle_hrf_library(stimdur = 3, tr = 1)
#' dim(lib)
glmsingle_hrf_library <- function(stimdur, tr) {
  stimdur <- .as_nonnegative_scalar(stimdur, "stimdur")
  tr <- .glms_positive_scalar(tr, "tr")
  lib <- as.matrix(utils::read.table(.glms_extdata("glmsingle_hrflibrary.tsv")))
  sampler <- function(len) seq(0, by = tr, length.out = ceiling(ceiling(len * 0.1) / tr))
  out <- do.call(cbind, lapply(seq_len(ncol(lib)), function(j) {
    .glms_hrf_resample(lib[, j], stimdur, tr, sampler)
  }))
  out <- out / max(out)
  out <- sweep(out, 2L, apply(out, 2L, max), "/")
  dimnames(out) <- list(NULL, paste0("hrf", seq_len(ncol(out))))
  out
}

#' @rdname glmsingle_hrf_library
#' @export
glmsingle_hrf <- function(stimdur, tr) {
  stimdur <- .as_nonnegative_scalar(stimdur, "stimdur")
  tr <- .glms_positive_scalar(tr, "tr")
  basic <- scan(.glms_extdata("glmsingle_basichrf.txt"), quiet = TRUE)
  sampler <- function(len) {
    end <- (len - 1) * 0.1
    seq(0, by = tr, length.out = ceiling(end / tr))
  }
  h <- .glms_hrf_resample(basic, stimdur, tr, sampler)
  h / max(h)
}

.glms_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) {
    stop(sprintf("%s must be a single positive number", name), call. = FALSE)
  }
  as.numeric(x)
}

# Trial design for one run: column i is `hrf` placed at onset i (0-based TR
# index) and truncated at the end of the run, as np.convolve(...)[:ntime].
.glms_trial_design <- function(onsets, n_time, hrf) {
  X <- matrix(0, n_time, length(onsets))
  L <- length(hrf)
  for (i in seq_along(onsets)) {
    rows <- onsets[i] + seq_len(L)
    keep <- rows <= n_time
    X[rows[keep], i] <- hrf[keep]
  }
  X
}
