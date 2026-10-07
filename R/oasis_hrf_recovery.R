#' OASIS HRF Recovery Testing Functions
#'
#' Functions to test OASIS's ability to recover HRF parameters from
#' rapid event-related designs with overlapping HRFs.

#' Generate Rapid Event-Related Design
#'
#' Creates a rapid event-related design with specified ISI range
#' 
#' @param n_events Number of events to generate
#' @param total_time Total time in seconds
#' @param min_isi Minimum inter-stimulus interval in seconds
#' @param max_isi Maximum inter-stimulus interval in seconds
#' @param seed Random seed for reproducibility
#' @return Numeric vector of event onset times
#' @examples
#' generate_rapid_design(n_events = 5, total_time = 40, seed = 1)
#' @importFrom stats runif
#' @export
generate_rapid_design <- function(n_events = 25, 
                                 total_time = 300,
                                 min_isi = 2, 
                                 max_isi = 4,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate random ISIs
  isis <- runif(n_events, min = min_isi, max = max_isi)
  
  # Calculate onset times
  onsets <- cumsum(c(5, isis[-length(isis)]))  # Start at 5s
  
  # Ensure all onsets fit within total time
  valid_onsets <- onsets[onsets < (total_time - 20)]  # Leave 20s at end for HRF tail
  
  return(valid_onsets)
}

#' Generate Synthetic fMRI Data with LWU HRF
#'
#' Creates synthetic fMRI time series using specified LWU HRF parameters
#'
#' @param onsets Vector of event onset times in seconds
#' @param tau LWU tau parameter (time-to-peak)
#' @param sigma LWU sigma parameter (width)
#' @param rho LWU rho parameter (undershoot amplitude)
#' @param TR Repetition time in seconds
#' @param total_time Total scan time in seconds
#' @param n_voxels Number of voxels to simulate
#' @param amplitudes Event amplitudes (scalar or vector)
#' @param noise_sd Standard deviation of noise
#' @param seed Random seed
#' @return List with Y (data matrix), true_hrf, true_betas, and design info
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
#' sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
#' dim(sim$Y)
#' }
#' @importFrom stats rnorm
#' @export
generate_lwu_data <- function(onsets,
                             tau = 6,
                             sigma = 2.5,
                             rho = 0.35,
                             TR = 1.0,
                             total_time = 300,
                             n_voxels = 10,
                             amplitudes = NULL,
                             noise_sd = 0.2,
                             seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!requireNamespace("fmrihrf", quietly = TRUE)) {
    stop("Package 'fmrihrf' is required for HRF simulation")
  }

  # Time grid
  time_points <- seq(0, total_time, by = TR)
  n_time <- length(time_points)
  n_trials <- length(onsets)
  
  # Create HRF
  hrf_times <- seq(0, 30, by = TR)  # 30 second HRF window
  true_hrf <- fmrihrf::hrf_lwu(hrf_times, tau = tau, sigma = sigma, rho = rho, normalize = "height")
  
  # Generate amplitudes if not provided
  if (is.null(amplitudes)) {
    amplitudes <- rnorm(n_trials, mean = 1, sd = 0.1)
  } else if (length(amplitudes) == 1) {
    amplitudes <- rep(amplitudes, n_trials)
  }
  
  # Use the production event builder: fractional onsets are preserved and
  # generator, candidate scoring, and final OASIS fit share precision = 0.1.
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = TR)
  time_points <- fmrihrf::samples(sframe, global = TRUE)
  hrf_obj <- .lwu_recovery_hrf(tau, sigma, rho)
  X <- .oasis_build_X_from_events(list(
    sframe = sframe,
    cond = list(onsets = onsets, hrf = hrf_obj, span = 30,
                amplitude = amplitudes)
  ))$X_trials

  # Generate true betas (trial-specific activations)
  true_betas <- matrix(rnorm(n_trials * n_voxels, mean = 1, sd = 0.3), 
                       nrow = n_trials, ncol = n_voxels)
  
  # Generate data: signal + noise
  signal <- X %*% true_betas
  noise <- matrix(rnorm(n_time * n_voxels, mean = 0, sd = noise_sd), 
                  nrow = n_time, ncol = n_voxels)
  
  # Add AR(1) structure to noise for realism
  ar_coef <- 0.3
  for (v in seq_len(n_voxels)) {
    for (t in 2:n_time) {
      noise[t, v] <- ar_coef * noise[t-1, v] + sqrt(1 - ar_coef^2) * noise[t, v]
    }
  }
  
  Y <- signal + noise
  
  # Create sampling frame for fmrilss
  # blocklens should be number of time points, not total time
  sframe <- fmrihrf::sampling_frame(blocklens = n_time, TR = TR)
  
  return(list(
    Y = Y,
    X = X,
    true_hrf = true_hrf,
    hrf_times = hrf_times,
    true_betas = true_betas,
    onsets = onsets,
    amplitudes = amplitudes,
    sframe = sframe,
    time_points = time_points,
    hrf_params = list(tau = tau, sigma = sigma, rho = rho),
    noise_sd = noise_sd,
    signal = signal,
    noise = noise
  ))
}

#' Create LWU HRF Grid for OASIS Search
#'
#' Generates a grid of LWU HRF models with varying parameters
#'
#' @param tau_range Range of tau values to test
#' @param sigma_range Range of sigma values to test  
#' @param rho_range Range of rho values to test
#' @param n_tau Number of tau values in grid
#' @param n_sigma Number of sigma values in grid
#' @param n_rho Number of rho values in grid
#' @return An object of class `"fmrilss_lwu_grid"`: a list containing the HRF
#'   functions in `hrfs` and their parameter data frame in `parameters`. The
#'   object can be supplied directly as `library_spec$pgrid` to [sbhm_build()].
#' @examples
#' grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
#' nrow(grid$parameters)
#' @export
create_lwu_grid <- function(tau_range = c(4, 8),
                           sigma_range = c(1.5, 3.5),
                           rho_range = c(0.1, 0.6),
                           n_tau = 5,
                           n_sigma = 3,
                           n_rho = 3) {

  # Create parameter grid
  tau_vals <- seq(tau_range[1], tau_range[2], length.out = n_tau)
  sigma_vals <- seq(sigma_range[1], sigma_range[2], length.out = n_sigma)
  rho_vals <- seq(rho_range[1], rho_range[2], length.out = n_rho)
  
  grid <- expand.grid(tau = tau_vals, sigma = sigma_vals, rho = rho_vals)
  
  # Generate HRF models
  hrf_list <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    # Create HRF function with these parameters
    params <- grid[i, ]
    # Create a closure that captures the parameters
    hrf_list[[i]] <- local({
      tau_val <- params$tau
      sigma_val <- params$sigma
      rho_val <- params$rho
      function(t) {
        fmrihrf::hrf_lwu(t, tau = tau_val, sigma = sigma_val, rho = rho_val, normalize = "height")
      }
    })
    # Add attributes for later reference
    attr(hrf_list[[i]], "tau") <- params$tau
    attr(hrf_list[[i]], "sigma") <- params$sigma
    attr(hrf_list[[i]], "rho") <- params$rho
    attr(hrf_list[[i]], "span") <- 30
  }
  
  structure(
    list(hrfs = hrf_list, parameters = grid),
    class = c("fmrilss_lwu_grid", "list")
  )
}

.lwu_recovery_hrf <- function(tau, sigma, rho) {
  force(tau); force(sigma); force(rho)
  structure(function(t) {
    fmrihrf::hrf_lwu(t, tau = tau, sigma = sigma, rho = rho, normalize = "height")
  }, class = c("hrf", "function"), span = 30,
  tau = tau, sigma = sigma, rho = rho)
}

#' Fit OASIS with HRF Grid Search
#'
#' Selects an LWU HRF by joint-model profile fit, then estimates OASIS betas.
#'
#' @details Each candidate is scored by pooled R-squared from unpenalized
#'   joint least squares with a voxel-specific intercept and all trial columns.
#'   The denominator sums squared deviations from each voxel's own mean.
#'   Ridge parameters affect only the final LSS estimates, not HRF selection.
#'   Simulation, scoring, and fitting use the same event construction at
#'   0.1-second precision. This is an in-sample selection criterion; overlapping
#'   events can leave HRF parameters weakly identifiable.
#'
#' @param Y Data matrix (time x voxels)
#' @param onsets Event onset times
#' @param sframe Sampling frame
#' @param hrf_grid List of HRF models to test
#' @param ridge_x Ridge parameter for design matrix
#' @param ridge_b Ridge parameter for aggregator
#' @return List with best HRF index, parameters, and beta estimates
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
#' sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
#' grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
#' fit <- fit_oasis_grid(sim$Y, sim$onsets, sim$sframe, grid)
#' fit$best_params
#' }
#' @export
fit_oasis_grid <- function(Y, onsets, sframe, hrf_grid,
                           ridge_x = 0.01, ridge_b = 0.01) {
  
  if (!is.matrix(Y) || !is.numeric(Y) || !all(is.finite(Y)) ||
      !ncol(Y) || nrow(Y) != length(fmrihrf::samples(sframe, global = TRUE))) {
    stop("Y must be a finite numeric matrix matching the sampling frame", call. = FALSE)
  }
  ss_tot <- sum(sweep(Y, 2L, colMeans(Y))^2)
  if (!is.finite(ss_tot) || ss_tot <= 0) {
    stop("Y must have positive finite within-voxel variation", call. = FALSE)
  }
  n_hrfs <- length(hrf_grid$hrfs)
  if (!n_hrfs || nrow(hrf_grid$parameters) != n_hrfs) {
    stop("hrf_grid must contain matching non-empty HRFs and parameters", call. = FALSE)
  }
  scores <- rep(-Inf, n_hrfs)
  for (i in seq_len(n_hrfs)) {
    scores[i] <- tryCatch({
      params <- hrf_grid$parameters[i, ]
      hrf_obj <- .lwu_recovery_hrf(params$tau, params$sigma, params$rho)
      X_trial <- .oasis_build_X_from_events(list(
        sframe = sframe,
        cond = list(onsets = onsets, hrf = hrf_obj, span = 30)
      ))$X_trials
      # LSS coefficients come from separate trial models. Select HRFs using
      # a joint-model profile residual, then estimate LSS betas for the winner.
      residuals <- stats::lm.fit(cbind(1, X_trial), Y)$residuals
      1 - sum(residuals^2) / ss_tot
    }, error = function(e) -Inf)
  }
  if (!any(is.finite(scores))) {
    stop("No HRF candidate produced a finite profile score", call. = FALSE)
  }
  scores[!is.finite(scores)] <- -Inf
  best_idx <- which.max(scores)
  best_params <- hrf_grid$parameters[best_idx, ]
  best_hrf <- .lwu_recovery_hrf(best_params$tau, best_params$sigma, best_params$rho)

  # Refit with best HRF
  best_beta <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(
          onsets = onsets,
          hrf = best_hrf,
          span = 30
        )
      ),
      ridge_mode = "fractional",
      ridge_x = ridge_x,
      ridge_b = ridge_b
    )
  )
  
  return(list(
    best_idx = best_idx,
    best_params = best_params,
    best_hrf = best_hrf,
    beta = best_beta,
    scores = scores,
    grid = hrf_grid
  ))
}

#' Compare HRF Recovery Methods
#'
#' Compares OASIS, SPMG1, SPMG3, and FIR for HRF recovery
#'
#' @param data Synthetic data from generate_lwu_data
#' @param hrf_grid Optional pre-computed HRF grid for OASIS
#' @return List with results from all methods
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
#' sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
#' grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
#' res <- compare_hrf_recovery(sim, hrf_grid = grid)
#' names(res)
#' }
#' @export
compare_hrf_recovery <- function(data, hrf_grid = NULL) {

  Y <- data$Y
  onsets <- data$onsets
  sframe <- data$sframe
  
  results <- list()
  
  # 1. OASIS with HRF grid search
  if (is.null(hrf_grid)) {
    hrf_grid <- create_lwu_grid()
  }
  
  message("Fitting OASIS with HRF grid search...")
  results$oasis <- fit_oasis_grid(Y, onsets, sframe, hrf_grid)
  
  # 2. Standard SPMG1
  message("Fitting SPMG1...")
  results$spmg1 <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(
          onsets = onsets,
          hrf = fmrihrf::HRF_SPMG1,
          span = 30
        )
      )
    )
  )
  
  # 3. SPMG3 (with derivatives)
  message("Fitting SPMG3...")
  results$spmg3 <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(
          onsets = onsets,
          hrf = fmrihrf::HRF_SPMG3,
          span = 30
        )
      )
    )
  )
  
  # 4. FIR basis
  message("Fitting FIR...")
  fir_hrf <- fmrihrf::hrf_fir_generator(nbasis = 15, span = 30)
  results$fir <- lss(
    Y = Y,
    X = NULL,
    method = "oasis",
    oasis = list(
      design_spec = list(
        sframe = sframe,
        cond = list(
          onsets = onsets,
          hrf = fir_hrf,
          span = 30
        )
      )
    )
  )
  
  # Add ground truth for comparison
  results$true_hrf <- data$true_hrf
  results$true_params <- data$hrf_params
  results$hrf_times <- data$hrf_times %||% data$time_points[
    data$time_points <= 30]
  # Fits use unit-amplitude event columns; include simulated event amplitudes
  # in the coefficient truth on that scale.
  results$true_betas <- data$true_betas * data$amplitudes
  
  return(results)
}

#' Calculate HRF Recovery Metrics
#'
#' Evaluates how well each method recovered the true HRF
#'
#' @param results Output from compare_hrf_recovery
#' @param true_hrf Ground truth HRF sampled on `results$hrf_times`. Shape
#'   comparisons normalize each curve to unit maximum absolute height.
#' @return Data frame with recovery metrics
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
#' sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
#' grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
#' res <- compare_hrf_recovery(sim, hrf_grid = grid)
#' calculate_recovery_metrics(res, sim$true_hrf)
#' }
#' @export
calculate_recovery_metrics <- function(results, true_hrf) {

  hrf_times <- results$hrf_times %||% seq(0, 30, length.out = length(true_hrf))
  if (length(hrf_times) != length(true_hrf)) {
    stop("true_hrf must match the stored HRF time grid", call. = FALSE)
  }
  unit_height <- function(x) as.numeric(x) / max(abs(x))
  true_hrf_eval <- unit_height(true_hrf)

  metrics <- data.frame(
    method = character(),
    mse = numeric(),
    correlation = numeric(),
    peak_time_error = numeric(),
    width_error = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 1. OASIS recovered HRF
  oasis_hrf_eval <- fmrihrf::hrf_lwu(
    hrf_times,
    tau = results$oasis$best_params$tau,
    sigma = results$oasis$best_params$sigma,
    rho = results$oasis$best_params$rho,
    normalize = "height"
  )
  
  oasis_metrics <- list(
    method = "OASIS",
    mse = mean((oasis_hrf_eval - true_hrf_eval)^2),
    correlation = cor(oasis_hrf_eval, true_hrf_eval),
    peak_time_error = abs(hrf_times[which.max(oasis_hrf_eval)] -
                          hrf_times[which.max(true_hrf_eval)]),
    width_error = abs(results$oasis$best_params$sigma - results$true_params$sigma)
  )
  metrics <- rbind(metrics, oasis_metrics)
  
  # 2. SPMG1 HRF
  spmg1_hrf_eval <- unit_height(fmrihrf::evaluate(fmrihrf::HRF_SPMG1, hrf_times))
  spmg1_metrics <- list(
    method = "SPMG1",
    mse = mean((spmg1_hrf_eval - true_hrf_eval)^2),
    correlation = cor(spmg1_hrf_eval, true_hrf_eval),
    peak_time_error = abs(hrf_times[which.max(spmg1_hrf_eval)] - hrf_times[which.max(true_hrf_eval)]),
    width_error = NA  # Can't directly extract width from SPMG1
  )
  metrics <- rbind(metrics, spmg1_metrics)
  
  # 3. SPMG3 HRF (canonical component)
  spmg3_hrf_eval <- unit_height(fmrihrf::evaluate(fmrihrf::HRF_SPMG3, hrf_times)[, 1])  # First basis = canonical
  spmg3_metrics <- list(
    method = "SPMG3",
    mse = mean((spmg3_hrf_eval - true_hrf_eval)^2),
    correlation = cor(spmg3_hrf_eval, true_hrf_eval),
    peak_time_error = abs(hrf_times[which.max(spmg3_hrf_eval)] - hrf_times[which.max(true_hrf_eval)]),
    width_error = NA
  )
  metrics <- rbind(metrics, spmg3_metrics)
  
  # 4. Beta recovery metrics
  metrics$beta_correlation <- NA_real_
  if (!is.null(results$true_betas)) {
    # Calculate beta correlations
    for (method in c("oasis", "spmg1")) {
      if (method == "oasis") {
        beta_cor <- cor(as.vector(results$oasis$beta), 
                        as.vector(results$true_betas))
      } else if (method == "spmg1") {
        beta_cor <- cor(as.vector(results$spmg1), 
                       as.vector(results$true_betas))
      }
      
      idx <- which(metrics$method == toupper(method))
      if (length(idx) > 0) {
        metrics$beta_correlation[idx] <- beta_cor
      }
    }
  }
  
  return(metrics)
}

#' Plot HRF Recovery Comparison
#'
#' Creates visualization comparing true vs recovered HRFs
#'
#' @param results Output from compare_hrf_recovery
#' @param save_path Optional path to save plot
#' @return A `ggplot2` plot object. When `save_path` is supplied, the same plot
#'   is also written to disk.
#' @examplesIf requireNamespace("fmrihrf", quietly = TRUE)
#' \donttest{
#' onsets <- generate_rapid_design(n_events = 4, total_time = 60, seed = 1)
#' sim <- generate_lwu_data(onsets, total_time = 60, n_voxels = 2, seed = 1)
#' grid <- create_lwu_grid(n_tau = 2, n_sigma = 2, n_rho = 2)
#' res <- compare_hrf_recovery(sim, hrf_grid = grid)
#' plot_hrf_comparison(res)
#' }
#' @importFrom ggplot2 ggplot aes geom_line scale_color_manual scale_linetype_manual labs theme_minimal theme ggsave
#' @export
plot_hrf_comparison <- function(results, save_path = NULL) {
  
  time <- hrf <- method <- NULL
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting", call. = FALSE)
  }

  # Time grid for plotting
  hrf_times <- seq(0, 30, by = 0.1)
  
  # True HRF
  true_hrf <- fmrihrf::hrf_lwu(
    hrf_times,
    tau = results$true_params$tau,
    sigma = results$true_params$sigma,
    rho = results$true_params$rho,
    normalize = "height"
  )
  
  # OASIS recovered
  oasis_hrf <- fmrihrf::hrf_lwu(
    hrf_times,
    tau = results$oasis$best_params$tau,
    sigma = results$oasis$best_params$sigma,
    rho = results$oasis$best_params$rho,
    normalize = "height"
  )
  
  # Standard models
  spmg1_hrf <- fmrihrf::evaluate(fmrihrf::HRF_SPMG1, hrf_times)
  spmg3_hrf <- fmrihrf::evaluate(fmrihrf::HRF_SPMG3, hrf_times)[, 1]
  spmg1_hrf <- spmg1_hrf / max(abs(spmg1_hrf))
  spmg3_hrf <- spmg3_hrf / max(abs(spmg3_hrf))
  
  # Create data frame for plotting
  plot_data <- data.frame(
    time = rep(hrf_times, 4),
    hrf = c(true_hrf, oasis_hrf, spmg1_hrf, spmg3_hrf),
    method = rep(c("True", "OASIS", "SPMG1", "SPMG3"), each = length(hrf_times))
  )
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = time, y = hrf, color = method, linetype = method)) +
    ggplot2::geom_line(size = 1.2) +
    ggplot2::scale_color_manual(values = c("True" = "black", 
                                 "OASIS" = "red",
                                 "SPMG1" = "blue", 
                                 "SPMG3" = "green")) +
    ggplot2::scale_linetype_manual(values = c("True" = "solid",
                                    "OASIS" = "dashed",
                                    "SPMG1" = "dotted",
                                    "SPMG3" = "dotdash")) +
    ggplot2::labs(title = "HRF Recovery Comparison",
         subtitle = sprintf("True params: tau=%.1f, sigma=%.1f, rho=%.2f | OASIS: tau=%.1f, sigma=%.1f, rho=%.2f",
                           results$true_params$tau, results$true_params$sigma, results$true_params$rho,
                           results$oasis$best_params$tau, results$oasis$best_params$sigma, results$oasis$best_params$rho),
         x = "Time (seconds)",
         y = "Normalized Response",
         color = "Method",
         linetype = "Method") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  if (!is.null(save_path)) {
    ggplot2::ggsave(save_path, p, width = 10, height = 6)
  }
  
  return(p)
}
