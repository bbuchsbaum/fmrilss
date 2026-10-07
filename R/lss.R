#' Least Squares Separate (LSS) Analysis
#'
#' Computes trial-wise beta estimates using the Least Squares Separate approach
#' of Mumford et al. (2012). This method fits a separate GLM for each trial,
#' with the trial of interest plus a single regressor formed by summing all
#' other trials in a one-basis design. A K-basis OASIS model uses one summed
#' other-trial regressor per basis function.
#'
#' @param Y A numeric matrix of size n × V where n is the number of timepoints
#'   and V is the number of voxels/variables
#' @param X A numeric matrix of size n × T for a one-basis design, with one
#'   column per trial. A raw K-basis OASIS design has n × (T K) columns and
#'   additionally requires `oasis$K`, `oasis$ntrials`, and an explicit
#'   `oasis$trial_basis_map`. When `X` is an unmodified multi-basis
#'   `fmridesign::design_matrix()` with one event term, its column metadata is
#'   used to infer that identity contract and canonicalize rows to trial-major,
#'   basis-within-trial order.
#' @param Z A numeric matrix of size n × F representing common fixed regressors
#'   included in every trial-wise model (e.g., intercept, condition effects,
#'   block effects). Their coefficients are not returned. If NULL, an
#'   intercept-only design is used. Defaults to NULL.
#' @param Nuisance A numeric matrix of size n × N representing nuisance regressors
#'   to be projected out before LSS analysis (e.g., motion parameters, physiological
#'   noise). If NULL, no nuisance projection is performed. Defaults to NULL
#' @param method Character string specifying which implementation to use.
#'   Options are:
#'   \itemize{
#'     \item "r_optimized" - Optimized R implementation (recommended, default)
#'     \item "cpp_optimized" - Optimized C++ implementation with parallel support
#'     \item "r_vectorized" - Standard R vectorized implementation  
#'     \item "cpp" - Standard C++ implementation
#'     \item "naive" - Simple loop-based R implementation (for testing)
#'     \item "oasis" - OASIS method with HRF support and ridge regularization
#'     \item "stglmnet" - overlap-aware elastic-net backend using `glmnet`
#'   }
#' @param block_size An integer specifying the voxel block size for parallel
#'   processing, only applicable when `method = "cpp_optimized"`. Defaults to 96.
#' @param oasis A list of options for the OASIS method (ridge, SE, design
#'   construction, etc.).
#'   See Details and \code{\link{oasis_options}} for the full list.
#'   \strong{Note:} \code{oasis$whiten} is deprecated and ignored with a warning.
#'   Use the \code{prewhiten} parameter instead for all temporal whitening.
#' @param stglmnet A list of options for the `method = "stglmnet"` backend.
#'   See Details and \code{\link{stglmnet_options}} for the common fields.
#' @param prewhiten A list of prewhitening options using the \pkg{fmriAR}
#'   package, or \code{NULL} (no whitening, the default).
#'   See Details and \code{\link{prewhiten_options}} for the full list.
#' @param trial_groups Optional vector (character, factor, or integer) with one
#'   condition label per trial, i.e. per column of `X`. When supplied, each
#'   trial-wise model uses one summed "other trials" regressor per condition
#'   (the LSS-N variant of Turner et al., 2012), with the trial of interest
#'   removed from its own condition's regressor, instead of a single regressor
#'   pooling every other trial. This is the model used by Nilearn's and
#'   NiBetaSeries' LSS beta series and is more accurate when conditions evoke
#'   different responses. Supported by methods `"r_optimized"`,
#'   `"cpp_optimized"`, `"cpp"`, and `"naive"`. Defaults to `NULL` (classic LSS).
#' @param ridge Optional fractional ridge penalty: one number, or two numbers
#'   `c(trial, others)` for the trial-of-interest and the other-trials
#'   coefficients. Each is a fraction of the mean design energy (the
#'   `ridge_mode = "fractional"` convention of OASIS), i.e. the penalty added to
#'   the trial diagonal is `ridge[1] * mean(c_i'c_i)`. Ridge shrinks trial
#'   estimates toward zero and can greatly reduce their variance in rapid
#'   designs where neighbouring trials overlap; it composes with
#'   `trial_groups` and `prewhiten`. Supported by methods `"r_optimized"`,
#'   `"cpp_optimized"`, `"cpp"`, and `"naive"`. Defaults to `NULL` (no ridge).
#'
#' @return Normally, a numeric matrix of trial-wise beta estimates: T × V for
#'   a one-basis design or (T K) × V for OASIS with K basis functions. With
#'   OASIS `return_diag = TRUE` or `return_se = TRUE`, returns
#'   `list(beta, diag?, se?)`; `beta` and `se` have the same row-by-voxel shape.
#'   One-basis diagnostics contain length-T vectors `d`, `alpha`, and `s`;
#'   multi-basis diagnostics contain K × K × T arrays `D`, `C`, and `E`.
#'   Coefficients for common regressors `Z` are not returned. Multi-basis beta
#'   and SE matrices carry the canonical `trial_basis_map` attribute. When
#'   estimated prewhitening is active, the actual fitted `fmriAR_plan` is
#'   available as `attr(result, "whiten_plan")` (and on `result$beta` for
#'   structured returns).
#'
#' @details
#' The LSS approach fits a separate GLM for each trial, where each model includes:
#' \itemize{
#'   \item The trial of interest (from column i of X)
#'   \item For a one-basis design, all other trials combined into one summed
#'     regressor. A K-basis OASIS model uses a K-column summed block.
#'   \item Common fixed regressors (Z matrix), whose coefficients are not returned
#' }
#' 
#' \strong{Computation.} Each LSS estimate is a linear functional of the data,
#' \eqn{\hat\beta_i = w_i^\top y}. The optimized methods build the
#' \eqn{n \times T} weight matrix \eqn{W} from the (small) trial design and
#' compute all trial betas with one matrix product \eqn{W^\top Y}. Because the
#' weights lie in the residual space of the confounds, the data matrix is never
#' residualized or copied, so the cost is a single
#' \eqn{O(nTV)} BLAS call regardless of the number of confounds.
#'
#' If Nuisance regressors are provided, the rank-revealed combined span
#' `cbind(Z, Nuisance)` is projected from both Y and X before fitting. Without a
#' separate Nuisance matrix, Z remains explicitly in every trial-wise model.
#' 
#' When using method="oasis", the following options are available in the oasis list
#' (see also \code{\link{oasis_options}} for a validated constructor):
#' \itemize{
#'   \item \code{design_spec}: A list for building trial-wise designs from event onsets using fmrihrf.
#'     Must contain: \code{sframe} (sampling frame), \code{cond} (list with \code{onsets},
#'     \code{hrf}, and optionally \code{span}), and optionally \code{others} (list of other conditions
#'     to be modeled as nuisances). When provided, X can be NULL and will be constructed automatically.
#'   \item \code{K}: Explicit basis dimension for multi-basis HRF models (e.g., 3 for SPMG3).
#'     A raw multi-basis \code{X} also requires \code{ntrials} and
#'     \code{trial_basis_map}. An unmodified multi-basis
#'     \code{fmridesign::design_matrix()} with one event term is recognized
#'     from its metadata; an ordinary raw \code{X} is otherwise interpreted as K=1.
#'   \item \code{ridge_mode}: Either "fractional" (default) or "absolute". In absolute mode,
#'     ridge_x and ridge_b are used directly as regularization parameters. In fractional mode,
#'     they represent fractions of the mean design energy for adaptive regularization.
#'   \item \code{ridge_x}: Ridge parameter for trial-specific regressors (default 0.05).
#'     Controls
#'     regularization strength for individual trial estimates.
#'   \item \code{ridge_b}: Ridge parameter for the aggregator regressor (default 0.05).
#'     Controls
#'     regularization strength for the sum of all other trials.
#'   \item \code{return_se}: Logical, whether to return model-based standard errors (default FALSE).
#'     This is available only for unpenalized OASIS without estimated prewhitening.
#'   \item \code{return_diag}: Logical, whether to return design diagnostics (default FALSE).
#'     When TRUE, includes diagnostic information about the design matrix structure.
#'   \item \code{block_cols}: Integer, voxel block size for memory-efficient processing (default 4096).
#'     Larger values use more memory but may be faster for systems with sufficient RAM.
#'   \item \code{ntrials}: Required number of trials when a raw K > 1 design is supplied.
#'   \item \code{trial_basis_map}: Required data frame for a raw K > 1 design,
#'     with one row per X column and fields \code{column}, \code{trial}, and \code{basis}.
#'   \item \code{design_spec$hrf_grid}: Candidate HRFs for grid-based selection
#'     within an event-built design. A top-level \code{oasis$hrf_grid} field is
#'     invalid and rejected.
#' }
#'
#' \strong{Prewhitening (temporal autocorrelation correction):}
#'
#' Use the top-level \code{prewhiten} parameter for all temporal whitening.
#' This replaces the old \code{oasis$whiten = "ar1"} syntax, which is now
#' deprecated and ignored.  Do \emph{not} put AR options inside the
#' \code{oasis} list; they belong in \code{prewhiten}.
#'
#' When using \code{method = "stglmnet"}, the backend accepts an additional
#' nested \code{stglmnet=} list for lambda selection, overlap-adaptive
#' penalties, and optional pooled trial parameterizations. The common pattern is
#' \code{stglmnet = stglmnet_options(mode = "cv")} to select lambda by
#' cross-validation, or \code{stglmnet = stglmnet_options(mode = "fixed",
#' lambda = 0.01)} for a fixed elastic-net fit. The backend reuses fmrilss
#' prewhitening and nuisance-projection utilities rather than maintaining a
#' separate whitening path.
#'
#' The \code{prewhiten} list accepts the following fields
#' (see also \code{\link{prewhiten_options}} for a validated constructor):
#' \itemize{
#'   \item \code{method}: Character, \code{"ar"} (default when the list is
#'     non-NULL), \code{"arma"}, or \code{"none"}.
#'     \code{"ar"} fits a pure autoregressive model; \code{"arma"} adds a
#'     moving-average component (requires \code{q > 0}).
#'   \item \code{p}: AR order.
#'     An integer, or \code{"auto"} (default) to select via AIC/BIC up to
#'     \code{p_max}. Use \code{p = 1} for a simple AR(1) model (the most
#'     common choice for fMRI); higher orders are rarely needed but may
#'     help with short TRs or multi-band sequences.
#'   \item \code{q}: Integer MA order for ARMA models (default 0).  Only
#'     relevant when \code{method = "arma"}.
#'   \item \code{p_max}: Integer, maximum AR order when \code{p = "auto"}
#'     (default 6).
#'   \item \code{pooling}: How AR coefficients are estimated across voxels.
#'     One of:
#'     \describe{
#'       \item{\code{"global"}}{(default) A single set of AR coefficients is
#'         estimated from the median autocorrelation across all voxels.
#'         Fast and usually adequate.}
#'       \item{\code{"voxel"}}{Voxel-adaptive noise model. Per-voxel residual
#'         autocorrelations are estimated, voxels are grouped into
#'         \code{voxel_bins} bins (default 50) of similar autocorrelation, and an
#'         AR model is refitted per bin; each bin gets its own filtered design,
#'         as in Nilearn's AR(1) GLM. Supported by methods
#'         \code{"r_optimized"}, \code{"cpp_optimized"} and \code{"cpp"}; other
#'         methods reject it.}
#'       \item{\code{"run"}}{Fit one AR model per run (requires \code{runs}).
#'         Useful when noise structure differs between runs.}
#'       \item{\code{"parcel"}}{Fit one AR model per parcel (requires
#'         \code{parcels}); each parcel is fitted with its own filtered design.
#'         Supported by methods \code{"r_optimized"}, \code{"cpp_optimized"}
#'         and \code{"cpp"}.}
#'     }
#'   \item \code{runs}: Integer vector of length \code{nrow(Y)} giving
#'     run/block labels. Required for \code{pooling = "run"} and recommended
#'     whenever data span multiple runs so that whitening respects run
#'     boundaries.
#'   \item \code{parcels}: Integer vector of length \code{ncol(Y)} giving
#'     parcel labels.  Required for \code{pooling = "parcel"}.
#'   \item \code{exact_first}: Character, \code{"ar1"} (default) or
#'     \code{"none"}.  When \code{"ar1"}, the first observation of each
#'     segment is scaled by \eqn{\sqrt{1 - \phi_1^2}} for the exact
#'     likelihood; \code{"none"} drops the first observation instead.
#'   \item \code{compute_residuals}: Logical (default TRUE). When TRUE,
#'     OLS residuals (see \code{residual_model}) are computed before fitting
#'     the noise model.  Set to FALSE only if Y is already residualized.
#'   \item \code{design}: Optional numeric design matrix whose projection
#'     produced those residuals. Supplying it opts in to fmriAR's correction
#'     for downward bias in residual autocovariance. When
#'     \code{compute_residuals = TRUE}, it must span the same columns as the
#'     full \code{X}/\code{Z}/\code{Nuisance} design, including the intercept
#'     that fmrilss adds when none is already represented.
#'   \item \code{acvf_correction}: Optional correction matrix or list of
#'     matrices from \code{fmriAR::acvf_bias_matrix()}, used instead of
#'     \code{design} when reusing a correction across datasets. The two fields
#'     are mutually exclusive.
#'   \item \code{voxel_bins}: Positive integer number of autocorrelation bins
#'     for \code{pooling = "voxel"} (default 50).
#'   \item \code{residual_model}: \code{"aggregate"} (default) estimates the
#'     noise model from residuals of the confounds plus one summed trial
#'     regressor per \code{trial_groups} level (or per basis function);
#'     \code{"full"} uses residuals of the full trial-wise design, whose many
#'     columns bias the residual autocorrelation downward in rapid designs.
#'     \code{"full"} is implied by \code{design}/\code{acvf_correction}.
#'   \item \code{correction_max_lag}: Positive integer lag budget used when
#'     \code{design} is supplied (default 25). The correction is intended for
#'     high-pass-filtered designs; without high-pass filtering, the required
#'     lag budget can become impractically large. See
#'     \code{fmriAR::fit_noise()} for details.
#' }
#'
#' \strong{Typical prewhiten recipes:}
#' \preformatted{
#'   # Simple AR(1) — good default for most fMRI data
#'   prewhiten = list(method = "ar", p = 1)
#'
#'   # Auto-select AR order (AIC), global pooling
#'   prewhiten = list(method = "ar", p = "auto")
#'
#'   # Per-run AR(1) for multi-run data
#'   prewhiten = list(method = "ar", p = 1, pooling = "run",
#'                    runs = blockids)
#'
#'   # Or use the validated constructor:
#'   prewhiten = prewhiten_options(method = "ar", p = 1, pooling = "run",
#'                                 runs = blockids)
#' }
#'
#' Prewhitening is applied before the LSS analysis to account for temporal
#' autocorrelation in the fMRI time series. Both Y and all design matrices
#' (X, Z, Nuisance) are filtered through the same whitening operator so that
#' OLS on the whitened system is equivalent to GLS on the original data.
#'
#' The OASIS method provides a mathematically equivalent but computationally optimized version
#' of standard LSS. It reformulates the per-trial GLM fitting as a single matrix operation,
#' eliminating redundant computations. This is particularly beneficial for designs with many
#' trials or when processing large datasets. When K > 1 (multi-basis HRFs), the output will
#' have K*ntrials rows, with basis functions for each trial arranged sequentially.
#'
#' @references
#' Mumford, J. A., Turner, B. O., Ashby, F. G., & Poldrack, R. A. (2012).
#' Deconvolving BOLD activation in event-related designs for multivoxel pattern
#' classification analyses. NeuroImage, 59(3), 2636-2643.
#'
#' Turner, B. O., Mumford, J. A., Poldrack, R. A., & Ashby, F. G. (2012).
#' Spatiotemporal activity estimation for multivoxel pattern analysis with
#' rapid event-related designs. NeuroImage, 62(3), 1429-1438.
#'
#' @examples
#' n_timepoints <- 100
#' n_trials <- 10
#' n_voxels <- 50
#'
#' X <- matrix(0, n_timepoints, n_trials)
#' for(i in 1:n_trials) {
#'   start <- (i-1) * 8 + 1
#'   if(start + 5 <= n_timepoints) {
#'     X[start:(start+5), i] <- 1
#'   }
#' }
#'
#' Y <- matrix(rnorm(n_timepoints * n_voxels), n_timepoints, n_voxels)
#' true_betas <- matrix(rnorm(n_trials * n_voxels, 0, 0.5), n_trials, n_voxels)
#' for(i in 1:n_trials) {
#'   Y <- Y + X[, i] %*% matrix(true_betas[i, ], 1, n_voxels)
#' }
#'
#' beta_estimates <- lss(Y, X)
#'
#' Z <- cbind(1, scale(1:n_timepoints))
#' beta_estimates_with_regressors <- lss(Y, X, Z = Z)
#'
#' Nuisance <- matrix(rnorm(n_timepoints * 6), n_timepoints, 6)
#' beta_estimates_clean <- lss(Y, X, Z = Z, Nuisance = Nuisance)
#'
#' \donttest{
#' beta_oasis <- lss(Y, X, method = "oasis",
#'                   oasis = list(ridge_x = 0.1, ridge_b = 0.1,
#'                               ridge_mode = "fractional"))
#'
#' result_with_se <- lss(Y, X, method = "oasis",
#'                      oasis = list(return_se = TRUE,
#'                                   ridge_mode = "absolute",
#'                                   ridge_x = 0, ridge_b = 0))
#' beta_estimates <- result_with_se$beta
#' standard_errors <- result_with_se$se
#'
#'   sframe <- fmrihrf::sampling_frame(blocklens = nrow(Y), TR = 1.0)
#'
#'   beta_auto <- lss(Y, X = NULL, method = "oasis",
#'                    oasis = list(
#'                      design_spec = list(
#'                        sframe = sframe,
#'                        cond = list(
#'                          onsets = c(10, 30, 50, 70),
#'                          hrf = fmrihrf::HRF_SPMG1,
#'                          span = 25
#'                        ),
#'                        others = list(
#'                          list(onsets = c(20, 40, 60, 80))
#'                        )
#'                      )
#'                    ))
#'
#'   beta_multibasis <- lss(Y, X = NULL, method = "oasis",
#'                         oasis = list(
#'                           design_spec = list(
#'                             sframe = sframe,
#'                             cond = list(
#'                               onsets = c(10, 30, 50, 70),
#'                               hrf = fmrihrf::HRF_SPMG3,
#'                               span = 30
#'                             )
#'                           ),
#'                           K = 3
#'                         ))
#' }
#'
#' @export
lss <- function(Y, X, Z = NULL, Nuisance = NULL,
                method = c("r_optimized", "cpp_optimized", "r_vectorized", "cpp", "naive", "oasis", "stglmnet"),
                block_size = 96, oasis = list(), stglmnet = list(), prewhiten = NULL,
                trial_groups = NULL, ridge = NULL) {
  
  method <- match.arg(method)

  if (!is.null(trial_groups) && method %in% c("r_vectorized", "oasis", "stglmnet")) {
    stop(
      "trial_groups (LSS-N) is supported by methods 'r_optimized', ",
      "'cpp_optimized', 'cpp', and 'naive'",
      call. = FALSE
    )
  }
  if (!is.null(ridge) && method %in% c("r_vectorized", "oasis", "stglmnet")) {
    stop(
      "ridge is supported by methods 'r_optimized', 'cpp_optimized', 'cpp', ",
      "and 'naive'; use oasis$ridge_x/ridge_b for method = 'oasis'",
      call. = FALSE
    )
  }
  ridge <- .lss_ridge_arg(ridge)

  if (!is.null(prewhiten)) {
    trusted_plan <- inherits(prewhiten, "fmrilss_internal_prewhiten")
    prewhiten <- .resolve_prewhiten_options(prewhiten, internal = trusted_plan)
  }

  if (method == "cpp_optimized") {
    block_size <- .as_positive_integer(block_size, "block_size")
  }
  
  # Drop legacy oasis$whiten: ignore it and warn if explicitly provided
  if (method == "oasis" && !is.null(oasis$whiten)) {
    warning(
      "oasis$whiten is deprecated and ignored; use the prewhiten parameter instead.",
      call. = FALSE
    )
  }

  # Fast-path to OASIS: it has its own coercion/validation and supports Matrix inputs
  if (method == "oasis") {
    if (is.null(X) && is.null(oasis$design_spec)) {
      stop("For method='oasis', either X or oasis$design_spec must be provided")
    }
    return(.lss_oasis(Y, X, Z, Nuisance, oasis, prewhiten))
  }

  if (method == "stglmnet") {
    return(.lss_stglmnet(Y, X, Z, Nuisance, stglmnet = stglmnet, prewhiten = prewhiten))
  }
  
  # Coerce S4 Matrix/data.frame to base matrices for non-OASIS methods
  to_mat <- function(M) {
    if (is.null(M)) return(NULL)
    if (inherits(M, "Matrix")) return(as.matrix(M))
    if (is.data.frame(M))     return(as.matrix(M))
    M
  }
  Y        <- to_mat(Y)
  X        <- to_mat(X)
  Z        <- to_mat(Z)
  Nuisance <- to_mat(Nuisance)
  
  # Input validation (non-OASIS)
  if (!is.matrix(Y) || !is.numeric(Y)) {
    stop("Y must be a numeric matrix")
  }
  if (nrow(Y) < 1L || ncol(Y) < 1L) {
    stop("Y must have at least one timepoint and one voxel")
  }
  if (!.all_finite(Y)) stop("Y contains non-finite values")
  if (!is.null(Z) && (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != nrow(Y))) {
    stop("Z must be a numeric matrix with the same number of rows as Y")
  }
  if (!is.null(Nuisance) && (!is.matrix(Nuisance) || !is.numeric(Nuisance) || nrow(Nuisance) != nrow(Y))) {
    stop("Nuisance must be a numeric matrix with the same number of rows as Y")
  }
  
  # For non-OASIS methods, X is required
  if (!is.matrix(X) || !is.numeric(X)) {
    stop("X must be a numeric matrix")
  }
  if (nrow(X) < 1L || ncol(X) < 1L) {
    stop("X must have at least one timepoint and one trial regressor")
  }
  if (any(!is.finite(X))) stop("X contains non-finite values")
  if (!is.null(Z) && any(!is.finite(Z))) stop("Z contains non-finite values")
  if (!is.null(Nuisance) && any(!is.finite(Nuisance))) {
    stop("Nuisance contains non-finite values")
  }
  if (nrow(Y) != nrow(X)) {
    stop("Y and X must have the same number of rows (timepoints)")
  }
  trial_names <- .validate_or_default_names(colnames(X), ncol(X), "Trial_", "X column names")
  voxel_names <- .validate_or_default_names(colnames(Y), ncol(Y), "Voxel_", "Y column names")
  colnames(X) <- trial_names
  groups <- .lss_group_codes(trial_groups, ncol(X))
  
  # Set up default experimental regressors (intercept) if not provided
  if (is.null(Z)) {
    Z <- matrix(1, nrow(Y), 1)
    colnames(Z) <- "Intercept"
  }
  
  # Check for zero or near-zero regressors
  .check_zero_regressors(X)
  
  # Step 1: Apply prewhitening if requested
  whiten_plan <- NULL
  prewhiten_active <- !is.null(prewhiten) && !is.null(prewhiten$method) &&
    prewhiten$method != "none"
  if (prewhiten_active && method %in% c("r_optimized", "cpp_optimized", "cpp")) {
    # Weight-matrix path: one filtered design per whitening operator, which
    # also supports voxel- and parcel-specific noise models.
    fit <- .lss_prewhitened(Y, X, Z, Nuisance, prewhiten, groups = groups,
                            method = method, ridge = ridge)
    result <- fit$beta
    rownames(result) <- trial_names
    colnames(result) <- voxel_names
    return(.attach_whiten_plan(result, fit$whiten_plan))
  }
  if (prewhiten_active) {
    if (prewhiten$pooling %in% c("voxel", "parcel")) {
      stop(
        "pooling='", prewhiten$pooling, "' estimates voxel-specific whitening ",
        "operators, which cannot be applied to a shared design matrix with ",
        "method='", method, "'. Use method='r_optimized', 'cpp_optimized' or ",
        "'cpp', or pooling='global'/'run'.",
        call. = FALSE
      )
    }
    whitened <- .prewhiten_data(Y, X, Z, Nuisance, prewhiten,
                                X_noise = .aggregate_trials(X, groups))
    whiten_plan <- whitened$whiten_plan
    Y <- whitened$Y_whitened
    X <- whitened$X_whitened
    if (!is.null(whitened$Z_whitened)) Z <- whitened$Z_whitened
    if (!is.null(whitened$Nuisance_whitened)) Nuisance <- whitened$Nuisance_whitened
  }

  # Step 2: Project out nuisance regressors if provided
  if (!is.null(Nuisance)) {
    # Create full nuisance design matrix
    X_nuisance <- cbind(Z, Nuisance)

    if (method %in% c("naive", "r_vectorized")) {
      # Reference paths operate on explicitly residualized data.
      proj_result <- .project_out_nuisance(Y, X, X_nuisance)
      Y_clean <- proj_result$Y_residual
      X_clean <- proj_result$X_residual
    } else {
      # Weight-matrix paths: the LSS weights of residualized trial regressors
      # lie in the nuisance residual space, so W'Y == W'(QY) and Y never
      # needs to be projected.
      Y_clean <- Y
      X_clean <- .lss_residualize_trials(X, X_nuisance)
    }
  } else {
    Y_clean <- Y
    X_clean <- X
  }

  # Step 3: Run LSS analysis with the chosen method
  result <- switch(method,
    "r_optimized" = .lss_r_optimized(Y_clean, X_clean, Z, groups = groups,
                                     ridge = ridge),
    "cpp_optimized" = .lss_cpp_optimized(Y_clean, X_clean, Z, block_size = block_size,
                                         groups = groups, ridge = ridge),
    "r_vectorized" = .lss_r_vectorized(Y_clean, X_clean, Z),
    "cpp" = .lss_cpp(Y_clean, X_clean, Z, groups = groups, ridge = ridge),
    "naive" = .lss_naive(Y_clean, X_clean, Z, groups = groups, ridge = ridge),
    stop("Unknown method: ", method)
  )
  
  # Add row and column names if available
  rownames(result) <- trial_names
  colnames(result) <- voxel_names

  .attach_whiten_plan(result, whiten_plan)
}

# Helper function to project out nuisance regressors
.project_out_nuisance <- function(Y, X, X_nuisance) {
  # Handle empty or NULL nuisance matrix early
  if (is.null(X_nuisance) || !is.matrix(X_nuisance) || ncol(X_nuisance) == 0L) {
    return(list(Y_residual = Y, X_residual = X))
  }

  # Robust residualization via the estimated-rank QR basis. qr.resid() may
  # retain arbitrary completion directions for rank-deficient designs (for
  # example, duplicate intercept or motion columns) and over-project.
  qrX <- qr(X_nuisance)
  basis <- if (qrX$rank > 0L) {
    qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]
  } else {
    matrix(numeric(), nrow(X_nuisance), 0L)
  }
  Y_residual <- Y - basis %*% crossprod(basis, Y)
  X_residual <- X - basis %*% crossprod(basis, X)

  list(Y_residual = Y_residual, X_residual = X_residual)
}

# Helper function to check for zero or near-zero regressors
.check_zero_regressors <- function(X, eps = 1e-12) {
  trial_names <- if (!is.null(colnames(X))) colnames(X) else paste0("Trial_", seq_len(ncol(X)))
  n <- nrow(X)
  regressor_norm <- sqrt(colSums(X^2))
  regressor_var <- if (n > 1L) {
    colSums(sweep(X, 2L, colMeans(X))^2) / (n - 1)
  } else {
    rep(NA_real_, ncol(X))
  }

  flagged <- which(regressor_norm < eps | (!is.na(regressor_var) & regressor_var < eps))
  for (i in flagged) {
    if (regressor_norm[i] < eps) {
      warning(sprintf("Trial regressor '%s' appears to be zero (norm = %g). This may cause numerical issues or NaN results.",
                     trial_names[i], regressor_norm[i]))
    } else {
      # Non-zero but constant (or nearly constant) regressor
      warning(sprintf("Trial regressor '%s' has very low variance (%g) and may cause numerical instability.",
                     trial_names[i], regressor_var[i]))
    }
  }
}

# Implementation functions (these will call the existing optimized functions)
.lss_r_optimized <- function(Y, X, Z, groups = NULL, ridge = c(0, 0)) {
  .lss_kernel_r(Y, X, Z, groups = groups, ridge = ridge)
}

.lss_cpp_optimized <- function(Y, X, Z, block_size = 96, groups = NULL,
                               ridge = c(0, 0)) {
  block_size <- .as_positive_integer(block_size, "block_size")

  # Pass an orthonormal, rank-revealed confound basis so the C++ projection
  # is exact even for collinear confounds.
  if (is.null(Z)) {
    Zb <- matrix(0, nrow = nrow(Y), ncol = 0)
  } else {
    qrZ <- qr(as.matrix(Z))
    Zb <- qr.Q(qrZ)[, seq_len(qrZ$rank), drop = FALSE]
  }

  # X = trial regressors, Z = confounds, Y = data
  # The C++ function expects: X=confounds, Y=data, C=trials
  lss_fused_optim_cpp(X = Zb, Y = Y, C = X, block_size = block_size,
                      groups = groups, use_omp = !.blas_is_threaded(),
                      ridge_x = ridge[1L], ridge_b = ridge[2L])
}

#' Does the linked BLAS run multithreaded?
#'
#' A single large matrix product is fastest with a multithreaded BLAS, while
#' distributing voxel blocks across OpenMP threads is faster with the
#' reference BLAS. Override with `options(fmrilss.blas_threaded = TRUE/FALSE)`.
#' @keywords internal
#' @noRd
.blas_is_threaded <- function() {
  opt <- getOption("fmrilss.blas_threaded")
  if (!is.null(opt)) return(isTRUE(opt))
  blas <- tryCatch(
    paste(extSoftVersion()[["BLAS"]], La_library()),
    error = function(e) ""
  )
  grepl("openblas|mkl|accelerate|veclib|blis|armpl", blas, ignore.case = TRUE)
}

.lss_r_vectorized <- function(Y, X, Z) {
  # Use existing lss_fast function with use_cpp = FALSE
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_fast(dset = NULL, bdes = bdes, Y = Y, use_cpp = FALSE))
}

.lss_cpp <- function(Y, X, Z, groups = NULL, ridge = c(0, 0)) {
  if (!is.null(groups) || any(ridge > 0)) {
    C_res <- .lss_residualize_trials(X, Z)
    W <- lss_weight_matrix_cpp(C_res, if (is.null(groups)) integer(0) else groups,
                               1e-12, ridge[1L], ridge[2L])
    return(crossprod(W, Y))
  }
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_fast(dset = NULL, bdes = bdes, Y = Y, use_cpp = TRUE))
}

.lss_naive <- function(Y, X, Z, groups = NULL, ridge = c(0, 0)) {
  bdes <- list(
    dmat_base = Z,
    dmat_ran = X,
    dmat_fixed = NULL,
    fixed_ind = NULL
  )
  return(lss_naive(Y, bdes, trial_groups = groups, ridge = ridge))
}

#' Orthogonal Projection Matrix
#'
#' Computes Q = I - X(X'X)^(-1)X' using QR decomposition for numerical stability.
#'
#' @param X Design matrix
#' @return Projection matrix Q
#' @keywords internal
#' @noRd
.Q_project <- function(X) {
  ## returns Q = I − X (XᵀX)⁻¹ Xᵀ   without allocating I
  qrX <- qr(X)
  basis <- if (qrX$rank > 0L) {
    qr.Q(qrX)[, seq_len(qrX$rank), drop = FALSE]
  } else {
    matrix(numeric(), nrow(X), 0L)
  }
  Q   <- diag(nrow(X))                # allocate once
  Q   <- Q - tcrossprod(basis)        # Q = I – QQᵀ for col(X)
  Q
}

#' Vectorized LSS Beta Computation
#'
#' Computes LSS beta estimates without explicit loops using vectorized operations.
#'
#' @param QC Projected trial regressors (n x T)
#' @param Ry Projected data (n x V)
#' @param eps Numerical tolerance
#' @return Beta matrix (T x V)
#' @keywords internal
#' @noRd
.lss_beta_vec <- function(C, Y, eps = 1e-12) {

  T <- ncol(C); V <- ncol(Y)

  # Shared building blocks ----------------------------------------------------
  total  <- rowSums(C)                     # n
  ss_tot <- sum(total^2)                   # scalar

  CtY <- crossprod(C, Y)                   # T × V
  CtC <- colSums(C^2)                      # T
  CtT <- crossprod(C, total)               # T
  total_Y <- drop(crossprod(total, Y))     # 1 × V  (row vector)

  # Per-trial "other" pieces ---------------------------------------------------
  #   b_i^T y  (T × V)
  BtY <- matrix(rep(total_Y, each = T), T, V) - CtY

  #   ||b_i||^2 ,  c_i^T b_i  (length T)
  bt2 <- ss_tot - 2*CtT + CtC
  ctb <- CtT   - CtC

  # Add guard for near-zero other-trial regressors (Fix 3)
  bt2[bt2 < eps] <- Inf

  # Numerator & denominator (memory-efficient version, Fix 1)
  ctb_bt2 <- ctb / bt2                     # T vector
  num <- CtY - sweep(BtY, 1, ctb_bt2, `*`) # Sweep for column-wise multiplication
  den <- CtC - (ctb^2) / bt2               # T vector
  
  # Broadcast den vector for division
  sweep(num, 1, pmax(den, eps), `/`)
}

#' Validate Input Dimensions and Types
#'
#' @param Y Data matrix
#' @param dmat_ran Trial design matrix
#' @keywords internal
#' @noRd
.validate_dims <- function(Y, dmat_ran) {
  if (!is.matrix(Y)) {
    if (is.data.frame(Y)) {
      Y <- as.matrix(Y)
      warning("Converting Y from data.frame to matrix")
    } else {
      stop("Y must be a numeric matrix or data.frame")
    }
  }
  if (nrow(Y) != nrow(dmat_ran)) {
    stop(sprintf("Y has %d timepoints, design has %d.",
                 nrow(Y), nrow(dmat_ran)))
  }
  Y
}

#' @keywords internal
#' @noRd
lss_fast <- function(dset, bdes, Y = NULL, use_cpp = TRUE) {
  # Data preparation and validation
  Y <- if (is.null(Y)) get_data_matrix(dset) else Y
  Y <- .validate_dims(Y, bdes$dmat_ran)

  # Ensure design matrices are matrices
  dmat_base <- as.matrix(bdes$dmat_base)
  dmat_ran <- as.matrix(bdes$dmat_ran)
  dmat_fixed <- if (!is.null(bdes$fixed_ind)) as.matrix(bdes$dmat_fixed) else NULL
  
  # Check for intercept in base design
  if (!any(apply(dmat_base, 2, function(x) all(abs(x - mean(x)) < 1e-10)))) {
    warning("No intercept detected in dmat_base. Consider adding one for proper baseline modeling.")
  }
  
  # Build combined base design matrix
  X_base_fixed <- if (!is.null(dmat_fixed)) {
    cbind(dmat_base, dmat_fixed)
  } else {
    dmat_base
  }
  
  n_events <- ncol(dmat_ran)
  
  # Hot-path early exit for single event
  if (n_events == 1) {
    if (use_cpp) {
      return(lss_compute_cpp(.lss_residualize_trials(dmat_ran, X_base_fixed), Y))
    } else {
      # Simple regression for single event
      Q <- .Q_project(X_base_fixed)
      Ry <- Q %*% Y
      QC <- Q %*% dmat_ran
      beta_matrix <- matrix(crossprod(QC, Ry) / drop(crossprod(QC)), nrow = 1)
      return(beta_matrix)
    }
  }

  if (use_cpp) {
    # Projected trial regressors suffice: the LSS weights lie in the confound
    # residual space, so the data need not be residualized.
    return(lss_compute_cpp(.lss_residualize_trials(dmat_ran, X_base_fixed), Y))
  }

  # --- pure R fall-back ---
  Q  <- .Q_project(X_base_fixed)
  Ry <- Q %*% Y
  QC <- Q %*% dmat_ran
  .lss_beta_vec(QC, Ry)
}

#' Extract Data Matrix from Dataset
#'
#' Internal helper to normalize supported in-memory dataset formats.
#'
#' @param dset Dataset object (format depends on your specific use case)
#' @return A numeric matrix where rows are timepoints and columns are voxels
#' @keywords internal
#' @noRd
get_data_matrix <- function(dset) {
  if (is.matrix(dset)) {
    return(dset)
  } else if (is.data.frame(dset)) {
    return(as.matrix(dset))
  } else {
    stop("Unsupported dataset format. Provide Y as a matrix or data frame.")
  }
}

#' Project Out Confound Variables
#'
#' Computes the orthogonal projection matrix Q = I - X(X'X)^(-1)X' that projects
#' out the space spanned by confound regressors X. This is useful for advanced
#' users who want to cache and reuse projection matrices.
#'
#' @param X Confound design matrix (n x p) where n is number of timepoints
#'   and p is number of confound regressors
#' @return Projection matrix Q (n x n) that projects out the column space of X
#'
#' @details
#' This function uses QR decomposition for numerical stability instead of
#' computing the Moore-Penrose pseudoinverse directly. Only the estimated-rank
#' QR basis is used, so redundant confound columns do not over-project. The
#' resulting matrix Q can be applied to data to remove the influence of
#' confound regressors.
#'
#' @examples
#' \donttest{
#' n <- 100
#' X_confounds <- cbind(1, 1:n)
#' Y_raw <- matrix(rnorm(n * 3), n, 3)
#'
#' Q <- project_confounds(X_confounds)
#'
#' Y_clean <- Q %*% Y_raw
#' }
#'
#' @export
project_confounds <- function(X) {
  .Q_project(as.matrix(X))
}
