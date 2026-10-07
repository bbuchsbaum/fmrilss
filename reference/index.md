# Package index

## Single-trial estimation

- [`LSSBeta`](https://bbuchsbaum.github.io/fmrilss/reference/LSSBeta.md)
  : LSSBeta object

- [`coef(`*`<glmsingle_fit>`*`)`](https://bbuchsbaum.github.io/fmrilss/reference/coef.glmsingle_fit.md)
  : Single-trial betas from a glmsingle fit

- [`glmsingle()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle.md)
  : GLMsingle single-trial response estimation

- [`glmsingle_design()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle_design.md)
  : fmridesign front end for glmsingle()

- [`lsa()`](https://bbuchsbaum.github.io/fmrilss/reference/lsa.md) :
  Least Squares All (LSA) Analysis

- [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) :
  Least Squares Separate (LSS) Analysis

- [`lss_beta_cpp()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_beta_cpp.md)
  : Vectorized LSS Beta Computation Using C++

- [`lss_cpp_optimized()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_cpp_optimized.md)
  : A wrapper for the optimized C++ LSS implementation

- [`lss_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_design.md)
  : LSS Analysis with fmridesign Objects

- [`lss_fit_wrappers`](https://bbuchsbaum.github.io/fmrilss/reference/lss_fit_wrappers.md)
  :

  Convenience wrappers for modern
  [`lss()`](https://bbuchsbaum.github.io/fmrilss/reference/lss.md) usage

- [`lss_fused_optim_cpp()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_fused_optim_cpp.md)
  : Fused Single-Pass LSS Solver (C++)

- [`lss_naive()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_naive.md)
  : Naive Least Squares Separate (LSS) Analysis

- [`lss_naive_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_naive_fit.md)
  : Naive LSS with modern signature

- [`lss_optimized()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_optimized.md)
  : Optimized LSS Analysis (Pure R)

- [`lss_optimized_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_optimized_fit.md)
  : Optimized LSS with modern signature

- [`lss_rank1()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_rank1.md)
  : Rank-1 GLM: joint voxel-wise HRF and trial amplitude estimation

## HRF modeling and recovery

- [`VoxelHRF`](https://bbuchsbaum.github.io/fmrilss/reference/VoxelHRF.md)
  : VoxelHRF object
- [`calculate_recovery_metrics()`](https://bbuchsbaum.github.io/fmrilss/reference/calculate_recovery_metrics.md)
  : Calculate HRF Recovery Metrics
- [`compare_hrf_recovery()`](https://bbuchsbaum.github.io/fmrilss/reference/compare_hrf_recovery.md)
  : Compare HRF Recovery Methods
- [`create_lwu_grid()`](https://bbuchsbaum.github.io/fmrilss/reference/create_lwu_grid.md)
  : Create LWU HRF Grid for OASIS Search
- [`estimate_voxel_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/estimate_voxel_hrf.md)
  : Estimate Voxel-wise HRF Basis Coefficients
- [`generate_lwu_data()`](https://bbuchsbaum.github.io/fmrilss/reference/generate_lwu_data.md)
  : Generate Synthetic fMRI Data with LWU HRF
- [`generate_rapid_design()`](https://bbuchsbaum.github.io/fmrilss/reference/generate_rapid_design.md)
  : OASIS HRF Recovery Testing Functions
- [`glmsingle_hrf_library()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle_hrf_library.md)
  [`glmsingle_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/glmsingle_hrf_library.md)
  : GLMsingle canonical HRF and HRF library
- [`lss_with_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf.md)
  : Perform LSS using Voxel-wise HRFs
- [`lss_with_hrf_pure_r()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_with_hrf_pure_r.md)
  : Least Squares Separate with voxel-wise HRF (basis-weighted)
- [`plot_hrf_comparison()`](https://bbuchsbaum.github.io/fmrilss/reference/plot_hrf_comparison.md)
  : Plot HRF Recovery Comparison

## Shared-basis HRF matching

- [`.sbhm_adaptive_ridge_gls()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_adaptive_ridge_gls.md)
  : Adaptive fractional ridge per voxel from conditioning
- [`.sbhm_auto_gate_threshold()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_auto_gate_threshold.md)
  : Auto-set cond_gate threshold from observed diagnostics
- [`.sbhm_build_trial_regs()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_build_trial_regs.md)
  : Build per-trial basis regressors (each T x r) from SBHM basis +
  design_spec
- [`.sbhm_prewhiten()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_prewhiten.md)
  : Optionally prewhiten Y, regs (stacked), intercept and nuisance
- [`.sbhm_resid()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_resid.md)
  : Residualize columns of M against Z (FWL)
- [`.sbhm_resolve_ridge()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_resolve_ridge.md)
  : Resolve ridge value from spec (absolute or fractional)
- [`.sbhm_solve()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-sbhm_solve.md)
  : Stable linear solve for multiple RHS with tiny ridge
- [`lss_sbhm()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm.md)
  : End-to-End LSS with Shared-Basis HRF Matching (SBHM)
- [`lss_sbhm_design()`](https://bbuchsbaum.github.io/fmrilss/reference/lss_sbhm_design.md)
  : SBHM Pipeline with fmridesign Models
- [`sbhm_amplitude_ls()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_amplitude_ls.md)
  : Single-shape GLM amplitudes given matched coordinates (global LS)
- [`sbhm_amplitude_lss1()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_amplitude_lss1.md)
  : Single-shape LSS amplitudes (2x2 per trial) given matched
  coordinates
- [`sbhm_amplitude_oasis_k1()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_amplitude_oasis_k1.md)
  : OASIS K=1 amplitudes per voxel with matched HRF columns. Returns
  list(beta, se)
- [`sbhm_build()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_build.md)
  : Build a Shared-Basis HRF Library (SBHM)
- [`sbhm_hrf()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_hrf.md)
  : Wrap a Learned Basis as an HRF (SBHM HRF)
- [`sbhm_match()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_match.md)
  : Match Voxels to Library HRFs in Shared Basis (SBHM)
- [`sbhm_prepass()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_prepass.md)
  : SBHM Prepass: Aggregate Fit in a Shared Basis
- [`sbhm_project()`](https://bbuchsbaum.github.io/fmrilss/reference/sbhm_project.md)
  : Project Trial-wise SBHM Coefficients to Scalar Amplitudes

## ITEM workflows

- [`item_build_design()`](https://bbuchsbaum.github.io/fmrilss/reference/item_build_design.md)
  : Build ITEM design metadata
- [`item_compute_u()`](https://bbuchsbaum.github.io/fmrilss/reference/item_compute_u.md)
  : Compute ITEM trial covariance matrix
- [`item_cv()`](https://bbuchsbaum.github.io/fmrilss/reference/item_cv.md)
  : Crossvalidated ITEM decoding
- [`item_fit()`](https://bbuchsbaum.github.io/fmrilss/reference/item_fit.md)
  : Fit ITEM decoder weights
- [`item_from_lsa()`](https://bbuchsbaum.github.io/fmrilss/reference/item_from_lsa.md)
  : Build an ITEM bundle from LS-A estimates
- [`item_predict()`](https://bbuchsbaum.github.io/fmrilss/reference/item_predict.md)
  : Predict targets from ITEM weights
- [`item_slice_fold()`](https://bbuchsbaum.github.io/fmrilss/reference/item_slice_fold.md)
  : Slice an ITEM bundle into train/test fold objects

## Options and supporting functions

- [`.convert_legacy_whiten()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-convert_legacy_whiten.md)
  : Convert old-style whitening options to new format
- [`.mixed_solve_cpp()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-mixed_solve_cpp.md)
  : Mixed Model Solver using C++
- [`.needs_advanced_prewhitening()`](https://bbuchsbaum.github.io/fmrilss/reference/dot-needs_advanced_prewhitening.md)
  : Check if advanced prewhitening is needed
- [`benchmark_mixed_solve()`](https://bbuchsbaum.github.io/fmrilss/reference/benchmark_mixed_solve.md)
  : Benchmark Mixed Model Implementations
- [`fit_oasis_grid()`](https://bbuchsbaum.github.io/fmrilss/reference/fit_oasis_grid.md)
  : Fit OASIS with HRF Grid Search
- [`fmrilss-package`](https://bbuchsbaum.github.io/fmrilss/reference/fmrilss-package.md)
  [`fmrilss`](https://bbuchsbaum.github.io/fmrilss/reference/fmrilss-package.md)
  : fmrilss: Least Squares Separate (LSS) Analysis for fMRI Data
- [`fmrilss_options`](https://bbuchsbaum.github.io/fmrilss/reference/fmrilss_options.md)
  : Option constructors for nested interfaces
- [`mixed_precompute()`](https://bbuchsbaum.github.io/fmrilss/reference/mixed_precompute.md)
  : Precompute Workspace for Optimized Mixed Model
- [`mixed_solve()`](https://bbuchsbaum.github.io/fmrilss/reference/mixed_solve.md)
  [`mixed_solve_cpp()`](https://bbuchsbaum.github.io/fmrilss/reference/mixed_solve.md)
  : Mixed Model Solver
- [`mixed_solve_optimized()`](https://bbuchsbaum.github.io/fmrilss/reference/mixed_solve_optimized.md)
  : Optimized Mixed Model Solver
- [`oasis_options()`](https://bbuchsbaum.github.io/fmrilss/reference/oasis_options.md)
  : Construct OASIS options
- [`prewhiten_options()`](https://bbuchsbaum.github.io/fmrilss/reference/prewhiten_options.md)
  : Construct prewhitening options
- [`project_confounds()`](https://bbuchsbaum.github.io/fmrilss/reference/project_confounds.md)
  : Project Out Confound Variables
- [`project_confounds_cpp()`](https://bbuchsbaum.github.io/fmrilss/reference/project_confounds_cpp.md)
  : Project Out Confounds Using C++
- [`stglmnet_options()`](https://bbuchsbaum.github.io/fmrilss/reference/stglmnet_options.md)
  : Construct stglmnet backend options
