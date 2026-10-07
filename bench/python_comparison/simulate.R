# Simulate a rapid event-related fMRI dataset shared by the fmrilss and Python
# LSS benchmarks.
#
# Usage: Rscript simulate.R <out_dir> [n_voxels] [seed] [iti_min] [iti_max]
#
# Writes little-endian float64 binaries (column-major, as R stores them) plus
# a JSON manifest with dimensions, so Python can read them with numpy.fromfile.

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "sim"
n_vox <- if (length(args) >= 2) as.integer(args[[2]]) else 2000L
seed <- if (length(args) >= 3) as.integer(args[[3]]) else 1L
iti_range <- if (length(args) >= 5) as.numeric(args[4:5]) else c(2, 8)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(seed)

TR <- 2
n_scans <- 400L
frame_times <- (seq_len(n_scans) - 1) * TR

# --- events: two conditions, jittered design (default ITI 2-8 s, rapid) ----
iti <- runif(400, iti_range[1], iti_range[2])
onsets <- cumsum(c(10, iti))
onsets <- onsets[onsets < max(frame_times) - 30]
n_trials <- length(onsets)
condition <- sample(rep(c("A", "B"), length.out = n_trials))

# --- Glover/SPM canonical HRF, convolved at 0.1 s and sampled at TR --------
dt <- 0.1
t_hrf <- seq(0, 32, by = dt)
hrf <- dgamma(t_hrf, 6) - dgamma(t_hrf, 16) / 6
hrf <- hrf / max(hrf)
t_hi <- seq(0, max(frame_times) + 32, by = dt)
X <- vapply(seq_len(n_trials), function(i) {
  u <- numeric(length(t_hi))
  u[which.min(abs(t_hi - onsets[i]))] <- 1
  # Zero-pad so trials early in the run are not truncated by filter()'s
  # leading NAs.
  pad <- length(hrf) - 1L
  conv <- stats::filter(c(numeric(pad), u), hrf, sides = 1)[-seq_len(pad)]
  stats::approx(t_hi, conv, xout = frame_times)$y
}, numeric(n_scans))

# --- confounds: intercept + DCT high-pass (128 s), 6 motion traces -------
n_dct <- floor(2 * n_scans * TR / 128)
dct <- vapply(seq_len(n_dct), function(k) {
  cos(pi * k * (seq_len(n_scans) - 0.5) / n_scans)
}, numeric(n_scans))
Z <- cbind(1, dct)
motion <- apply(matrix(rnorm(n_scans * 6), n_scans, 6), 2, cumsum)
motion <- scale(motion)

# --- trial amplitudes: condition means differ, trial-level variability ----
cond_mean <- c(A = 1.0, B = 0.25)
mu_vox <- matrix(rnorm(n_vox, 1, 0.3), 1, n_vox)
beta_true <- outer(cond_mean[condition], rep(1, n_vox)) * (mu_vox[rep(1, n_trials), ]) +
  matrix(rnorm(n_trials * n_vox, 0, 0.5), n_trials, n_vox)

# --- AR(1) noise with voxel-specific coefficients, drift and motion ------
phi <- runif(n_vox, 0.1, 0.6)
eps <- matrix(rnorm(n_scans * n_vox), n_scans, n_vox)
noise <- eps
for (t in 2:n_scans) noise[t, ] <- phi * noise[t - 1, ] + eps[t, ]
noise <- noise * rep(sqrt(1 - phi^2), each = n_scans)
noise_sd <- 1.0
drift <- dct[, 1:3] %*% matrix(rnorm(3 * n_vox, 0, 2), 3, n_vox)
mot_eff <- motion %*% matrix(rnorm(6 * n_vox, 0, 0.3), 6, n_vox)
Y <- 100 + X %*% beta_true + drift + mot_eff + noise_sd * noise

wr <- function(m, nm) {
  storage.mode(m) <- "double"
  writeBin(as.vector(m), file.path(out_dir, paste0(nm, ".bin")), endian = "little")
}
wr(Y, "Y"); wr(X, "X"); wr(Z, "Z"); wr(motion, "motion"); wr(beta_true, "beta_true")
writeLines(condition, file.path(out_dir, "condition.txt"))
writeLines(sprintf("%.6f", onsets), file.path(out_dir, "onsets.txt"))
writeBin(as.double(phi), file.path(out_dir, "phi.bin"), endian = "little")

manifest <- sprintf(
  '{"n_scans": %d, "n_trials": %d, "n_vox": %d, "n_z": %d, "n_motion": 6, "TR": %g}',
  n_scans, n_trials, n_vox, ncol(Z), TR
)
writeLines(manifest, file.path(out_dir, "manifest.json"))
cat(manifest, "\n")
