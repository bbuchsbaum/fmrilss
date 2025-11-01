## Benchmark SBHM vs voxel-wise HRF estimation (illustrative)
## Run from an R session after installing/loading fmrilss and fmrihrf.

suppressPackageStartupMessages({
  library(fmrilss)
  library(fmrihrf)
})

set.seed(42)

# Simulation settings (adjust for your machine)
TR     <- 1
Tlen   <- 200          # time points
V      <- 1000         # voxels
ntr    <- 12           # trials
span   <- 32
r      <- 6            # SBHM rank

message(sprintf("Sim config: T=%d, V=%d, trials=%d, r=%d", Tlen, V, ntr, r))

# Sampling frame and time grid
sframe <- sampling_frame(blocklens = Tlen, TR = TR)
times  <- samples(sframe, global = TRUE)

# Build a small HRF library (gamma grid) for SBHM
shapes <- seq(6, 10, by = 2)
rates  <- seq(0.8, 1.2, by = 0.2)
pgrid  <- expand.grid(shape = shapes, rate = rates)
gamma_fun <- function(shape, rate) {
  f <- function(t) fmrihrf::hrf_gamma(t, shape = shape, rate = rate)
  fmrihrf::as_hrf(f, span = span)
}

message("Building SBHM basis...")
sbhm <- sbhm_build(
  library_spec = list(fun = gamma_fun, pgrid = pgrid, span = span, precision = 0.1, method = "conv"),
  r = r,
  sframe = sframe,
  baseline = c(0, 0.5),
  normalize = TRUE
)

hrf_B <- sbhm_hrf(sbhm$B, sbhm$tgrid, sbhm$span)

# Events (safely within acquisition window)
onsets <- seq(12, max(times) - 30, length.out = ntr)
design_spec <- list(sframe = sframe, cond = list(onsets = onsets, duration = 0, span = 30))

# Build per-trial regressors in SBHM basis for simulation
regressors_by_trial <- vector("list", ntr)
for (t in seq_len(ntr)) {
  rr_t <- regressor(onsets = onsets[t], hrf = hrf_B, duration = 0, span = 30, summate = FALSE)
  regressors_by_trial[[t]] <- evaluate(rr_t, grid = times, precision = 0.1, method = "conv")
}

# Simulate data: choose a true library coordinate per voxel; draw trial amplitudes
alpha_idx <- sample(seq_len(ncol(sbhm$A)), V, replace = TRUE)
alpha_mat <- sbhm$A[, alpha_idx, drop = FALSE]     # r×V
amps_true <- matrix(rnorm(ntr * V, mean = 1.5, sd = 0.4), ntr, V)

Y <- matrix(rnorm(Tlen * V, sd = 0.4), Tlen, V)
for (v in seq_len(V)) {
  for (t in seq_len(ntr)) {
    X_t <- regressors_by_trial[[t]]   # T×r
    Y[, v] <- Y[, v] + amps_true[t, v] * as.numeric(X_t %*% alpha_mat[, v])
  }
}

bench <- list()

# 1) SBHM end-to-end
message("Timing SBHM end-to-end...")
bench$sbhm <- system.time({
  res_sbhm <- lss_sbhm(
    Y = Y,
    sbhm = sbhm,
    design_spec = design_spec,
    prepass = list(ridge = list(mode = "fractional", lambda = 0.01)),
    match = list(topK = 1),
    return = "amplitude"
  )
})
cat("SBHM time (sec): ", bench$sbhm[[3]], "\n")

# 2) Voxel-wise HRF (baseline comparator)
# Prefer FIR basis if available in fmrihrf, otherwise fall back to SPMG3.
make_fir <- NULL
if ("make_hrf" %in% getNamespaceExports("fmrihrf")) {
  # Try to construct FIR basis with, say, 20 taps at TR=1
  try({
    make_fir <- fmrihrf::make_hrf("fir", nbasis = 20, span = span)
  }, silent = TRUE)
}
basis_vox <- if (!is.null(make_fir)) make_fir else fmrihrf::HRF_SPMG3

events_df <- data.frame(onset = onsets, duration = 0, condition = "A")

message("Timing voxel-wise HRF estimation + LSS...")
bench$voxhrf <- system.time({
  est <- estimate_voxel_hrf(Y, events = events_df, basis = basis_vox, nuisance_regs = NULL)
  betas_vox <- lss_with_hrf(Y, events = events_df, hrf_estimates = est, engine = "R", verbose = FALSE)
})
cat("Voxel-wise HRF time (sec): ", bench$voxhrf[[3]], "\n")

speedup <- bench$voxhrf[[3]] / bench$sbhm[[3]]
cat(sprintf("Speedup (voxel-wise / SBHM): %.2fx\n", speedup))

## Optional accuracy check (amplitude correlation when bases align)
if (nrow(res_sbhm$amplitude) == nrow(amps_true)) {
  cors <- suppressWarnings(cor(as.vector(res_sbhm$amplitude), as.vector(amps_true)))
  cat(sprintf("Amplitude correlation (SBHM vs truth): %.3f\n", cors))
}

