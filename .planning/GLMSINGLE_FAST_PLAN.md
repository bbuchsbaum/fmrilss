# Fast GLMsingle in fmrilss: due diligence and implementation plan

**Status:** plan v2, prepared for third-party review.
**Reference implementation:** GLMsingle Python, `cvnlab/GLMsingle` at
`1ab54a6` (2025-11-09, current HEAD), with `fracridge` 3.0.
**Proposal under review:** "a substantially faster implementation that
preserves GLMsingle's estimator" (external memo, sections cited below as
"proposal §n").

## 0. Goal and decisions already made

Build `glmsingle()` in fmrilss. It must reproduce GLMsingle's single-trial
estimates (types A–D) within floating-point tolerance and run substantially
faster. Changes to the estimator itself (continuous HRFs, SURE, cross-relation
initialisation) form a separate research track, gated behind a separate
accuracy review.

| Decision | Resolution | Rationale |
|---|---|---|
| Where it lives | **fmrilss** | fmrilss is the home for single-trial estimation. It shares the Rcpp/Armadillo/OpenMP toolchain, the fmrihrf/fmridesign integration, and the test and bench infrastructure |
| Parity target | **Python GLMsingle at `1ab54a6`** | It can be scripted and checked without a MATLAB licence. Both ports live in one repo and are actively maintained. Python's latest fix (`19e6617`, 2025-08, pcnum selection + calcbadness session loop) is included. The MATLAB-only commits after it are audited in §2.3 |
| Restrict GLMdenoise CV to PC-selection voxels | **Yes, by default** | Exact: the selected `pcnum` and every estimate are identical (§2.2). It only drops the full-volume `glmbadness` diagnostic, which is available on request. A test asserts that `xvaltrend` and `pcnum` match the unrestricted run (task 3.4) |

## 1. Bottom line

The proposal's core is sound. I checked every source-code claim against the
pinned commit. I also checked the two central exactness claims numerically,
using GLMsingle's own `fracridge` and `calcbadness` code
(`.planning/glmsingle_checks/verify_core_identities.py`; run instructions are in
its header).

| Check | Result |
|---|---|
| Run-blocked fractional ridge vs `fracridge` on GLMsingle's stacked design | coef rel. err 8e-13, alpha rel. err 6e-13 |
| Compiled CV loss `Σ d_i (z_i − m_i)^2 + c` vs `calcbadness` (LORO and grouped folds, 2 sessions) | rel. err 3e-16 |
| `fracridge` grid-interpolated α vs the exact root of proposal §2's equation | median 1.2%, max 4% apart; achieved fraction off by ≤ 0.006 |

### Two corrections to the proposal

1. **`fracridge` does not solve the fraction equation.** It interpolates
   `log(1+α)` on a fixed grid: `10^k`, step 0.2 decades, spanning
   `1e-2·s_min²` to `1e4·s_max²`. Solving the monotone equation exactly, as
   proposal §2 suggests, is a slightly different estimator. For parity we must
   reproduce the grid, built from the *global* `s_min`/`s_max` over all run
   blocks, and the interpolation. With that, the run-blocked version matched
   to 1e-12. The exact root can be an opt-in option.
2. **GLMsingle computes in float32.** It casts data and designs to float32, so
   comparison against Python is limited to about 1e-4 relative. Discrete choices
   (HRF index, `pcnum`, fraction index) flip on near-ties. Exactness is
   therefore proven in two tiers (§4).

### Where the speed actually comes from

The biggest gain is **touching the time series once**. Factorisation savings
are secondary. Details are in §2.2.

## 2. Due diligence

### 2.1 Claim-by-claim

| Proposal § | Claim | Verdict | Notes |
|---|---|---|---|
| 1 | Each run's design carries columns for all trials | **Verified** | `designSINGLE[run]` is `T_r × numtrials` (glmsingle.py:583-612) |
| 1 | Fitter is called once per fraction | **Verified** | The type C/D loop calls `glm_estimatemodel` per frac (≈glmsingle.py:1530). Each call redoes the SVD, the dense projection and the full prediction |
| 1 | Dense `T×T` projectors | **Verified** | `make_projection_matrix` → `combinedmatrix @ data` and `polymatrix @ data` on every call (glm_estimatemodel.py:528-562) |
| 2 | Block-diagonal design; a pooled fraction preserves the estimator | **Verified, with correction** | See §1, correction 1. Also replicate the `tol=1e-10` zeroing of singular values per block |
| 2 | Cost goes from `O((Rn)^3)` to `O(Rn^3)` | Correct, minor | At `Rn≈750` the global SVD costs < 1 s. The waste is elsewhere (§2.2) |
| 3 | Low-rank nuisance application `Y − Q(QᵀY)` | **Verified, high value** | Removes the `T²V` term |
| 3 | Sherman–Morrison PC-prefix updates | Correct (re-derived) but **not needed** | Run blocks are about 60×60, so a fresh solve per PC count costs microseconds. Reusing `QᵀY` prefixes is the real win, and avoiding the update also avoids its cancellation risk |
| 3 | Fit uses poly+extra projection; R² uses poly-only data and the raw design | **Verified** | glm_estimatemodel.py:849-874. `R2 = 100·(1 − Σ‖P_poly(y − X_raw β)‖² / Σ‖P_poly y‖²)`, with no mean subtraction (`calc_cod_stack`) |
| 4 | Compile pairwise repeat CV into fixed targets and weights | **Verified exactly** | Normalisation is fixed from `results[0]` per session. Also holds for the GLMdenoise stage, where the reference is the 0-PC fit |
| 4 | Reconstruct only trials with `d_i > 0` during selection | Correct | Small saving; the fraction norm still uses all trials |
| 5 | Banded Gram plus Woodbury nuisance correction | Correct, **defer** | HRF support is about 25 TRs vs an ISI of about 3 TRs, so the half-bandwidth is about 10 of about 60 trials. Dense BLAS wins |
| 5 | Delay-invariant autocorrelation | Correct, **not applicable** | The GLMsingle library is not a pure delay family |
| 6 | Compressed library with a projector-distance certificate | Correct bound, **probably useless** | R² gaps between library HRFs are often below `δ_h‖y‖²`. Once the data pass below is in place the exact sweep is cheap |
| 7 | VarPro objective and envelope gradient | Correct (re-derived) | Changes the estimator → research track. The HRF-scale normalisation point is important |
| 8 | Polyphase cross-relation initialiser | Mathematically correct | Needs a regular event grid. Polynomial projection breaks commutativity. Research only |
| 9 | SURE / evidence | Correct for fixed λ; the fixed-fraction Jacobian caveat is right | Changes the selection target. Research only |
| 10 | Reject rank-one SVD of an unrestricted basis fit | Agree | |
| 11 | 43.5× kernel timing | Not reproducible (no script) | A hypothesis, and a kernel benchmark rather than end-to-end |

### 2.2 Findings the proposal missed

* **Matched-filter `XᵀY`.** An un-projected trial column has about `L ≈ 25`
  nonzeros. So `X_hᵀY` is a sparse cross-correlation costing `O(n·L·V)` per
  HRF per run, vs `O(n·T·V)` dense. Then `XᵀMY = XᵀY − (XᵀQ)(QᵀY)` exactly.
  The 20-HRF library sweep becomes nearly free relative to reading the data.
* **No time-series predictions.** Every R², R2run, and residual quantity
  GLMsingle reports is a quadratic form in `{β, XᵀPy, XᵀPX, ‖Py‖²}`.
  `glm_predictresponses` is never needed.
* **The GLMdenoise CV only needs the PC-selection voxels.** `glmbadness` is
  computed for every voxel, but `xvaltrend` reads only `ix`: voxels above
  `pcR2cutoff` that also pass the mask, or a top-100 fallback
  (glmsingle.py:1378-1400). Restricting to `ix` is exact for every estimate.
* **Betas only for the winning HRF.** The type-B stage stores betas for all
  `nh` HRFs and discards all but one.
* **All fractions in one grid pass.** In `fracridge`, the per-voxel cost for F
  fractions is one `newlen` evaluation plus one interpolation.
* **`findtailthreshold` is nondeterministic in Python.** It uses sklearn
  `GaussianMixture(n_init=3)` without a seed, so `brainR2`/`pcR2cutoff`, and
  therefore the noise pool, can differ between Python runs. The parity
  harness must record the thresholds Python actually used, and tests inject
  them.

**Cost model:** GLMsingle runs about `1 + nh + (n_pcs+1)·(#HRF groups) + 2·F`
`glm_estimatemodel` passes (≈ 70 with the defaults `nh=20, n_pcs=10, F=20`).
Each pass touches `T²V` for projection and `nTV` for prediction, per run. The
fast path touches the data once for `QᵀY` and `‖Py‖²`, plus 20 sparse
`n·L·V` correlations. Everything after that lives in `n`-space. For
`T≈220, n≈63, L≈25`, a ≥ 10× end-to-end speedup is plausible. That remains
to be measured (task 4.3).

### 2.3 Python vs MATLAB: commits after the Python fix

| MATLAB commit | Change | Effect on the estimator | Plan |
|---|---|---|---|
| `a07f152` (2025-05) | `constructpolynomialmatrix` uses `f*(olsmatrix(f)*v)` instead of a `T×T` projector | None (memory only) | We use an orthonormal basis anyway |
| `91e5b7e` (2025-11) | `findtailthreshold`: GMM `RegularizationValue = eps` | Only in degenerate cases (variance collapsing toward 0) | Adopt it: deterministic EM plus eps variance floor. Document it as the one intentional deviation from Python |

## 3. Architecture: clean, modular, built on existing machinery

### 3.1 Design principles

1. **Layers with one-way dependencies.** Each stage is a pure function from
   plain lists/matrices to plain lists. There is no I/O, no global RNG, and no
   hidden state. The orchestrator is thin. This makes every stage testable in
   isolation, and lets tests inject Python intermediates such as the tail
   thresholds.
2. **Design geometry, voxel statistics, and hyperparameter evaluation are
   separate objects.** Geometry is built once. Voxel statistics are built once
   per voxel tile. Selection only reads statistics.
3. **fmrilss conventions, not GLMsingle's.** Data is `Y` (time × voxels,
   runs concatenated) with run identity from a `sampling_frame` or `run` ids,
   as in `lss_design()`. Betas are trials × voxels, as in `lss()`. Options
   come from an `*_options()` constructor validated with the existing `.as_*`
   helpers. GLMsingle-compatible names are used for every option and output
   field, so users of GLMsingle can map across directly.
4. **One C++ translation unit for the hot path.** It contains only `n`-space
   and data-pass kernels. Orchestration and selection logic stay in R, where
   they are short and readable.
5. **No new hard dependencies.**

### 3.2 Module layout

```
R/glmsingle.R            glmsingle(), glmsingle_design(): public API, orchestration only
R/glmsingle_options.R    glmsingle_options(): defaults mirror GLMsingle (or add to R/options.R)
R/glmsingle_runs.R       run geometry: onset TR indices, trial↔run↔condition↔session maps, folds
R/glmsingle_hrf.R        GLMsingle HRF library (vendored TSV, stimdur conv, pchip→TR, peak-normalise)
R/glmsingle_nuisance.R   per-run orthonormal bases: poly | extra | noise PCs (nested prefixes); PC derivation
R/glmsingle_stats.R      R wrappers around the data pass; quadratic-form R²/R2run
R/glmsingle_ridge.R      run-blocked spectral OLS / fracridge path (fracridge-compatible grid)
R/glmsingle_select.R     compiled CV, tail threshold (GMM), select_noise_regressors, autoscale
R/glmsingle_methods.R    class "glmsingle_fit": print/summary/coef/as.matrix
src/glmsingle_kernels.cpp  data pass (QᵀY, ‖Py‖², matched-filter XᵀY), block spectral path
src/blocked_products.h     shared residualise-and-multiply helper (also used by OASIS; see 3.3)
inst/extdata/glmsingle_hrflibrary.tsv   vendored, BSD-3 attribution in inst/COPYRIGHTS
tools/glmsingle_ref/       pinned Python harness + fixture writer (Rbuildignored)
tests/testthat/helper-glmsingle.R           fixture reader, simulator, dense stage references
tests/testthat/test-glmsingle-{runs,hrf,nuisance,stats,ridge,select,pipeline,parity}.R
```

Dependency order: `options → runs → hrf → nuisance → stats → ridge → select →
glmsingle()`. No module calls upward.

Size budget: about 1,200 lines of R and about 400 lines of C++, excluding
tests. A file approaching 400 lines is a signal to split it.

### 3.3 Reuse map

| Existing machinery | Where | Use in GLMsingle work | Action |
|---|---|---|---|
| `.voxhrf_orthonormal_span()` | R/voxel_hrf.R | Orthonormal nuisance bases with a rank tolerance | **Promote** to `.orthonormal_span()` in R/aaa_utils.R. Update its 3 callers (voxel_hrf, lss_with_hrf). Covered by existing tests |
| `oasis_AtY_SY_blocked`, `oasisk_compute_RY_norm2` | src/oasis_core.cpp | Same pattern as the data pass: blocked `Y − Q(QᵀY)`, products, norms | **Extract** a header-only helper `blocked_products.h`. The GLMsingle data pass uses it. OASIS migrates in its own commit, gated by the OASIS tests |
| `.as_positive_integer`, `.as_scalar_logical`, `.as_nonnegative_scalar`, `.validate_option_names`, `.as_integer_ids` | R/aaa_utils.R | Option and input validation | Reuse as-is |
| `oasis_options()` pattern | R/options.R | Shape of `glmsingle_options()` | Follow it |
| `.set_beta_dimnames`, `.default_trial_names` | R/aaa_utils.R | Naming trials × voxels outputs | Reuse |
| `lss_design()` + `.validate_design_models()` | R/lss_design.R | Template for `glmsingle_design(Y, event_model, baseline_model)`: fmridesign event model → onsets + condition labels; baseline nuisance term → `extra_regressors` | Reuse the validator and mirror the structure. Error if onsets are not on the TR grid, which GLMsingle requires |
| `fmrihrf::sampling_frame`, `blocklens`, `blockids` | fmrihrf | Run structure | Reuse |
| `fmrihrf::evaluate` on user HRF objects | fmrihrf | **Non-parity** option: a user-supplied HRF library as a list of fmrihrf HRFs | Reuse. The default library uses the GLMsingle builder for parity |
| `OPENMP` setup | src/Makevars | Threading over voxel tiles | Reuse |
| bench harness conventions | bench/ | `bench/run_glmsingle_benchmark.R` | Follow them |
| `.item_safe_solve()` (chol→svd→pinv) | R/item_compute_u.R | — | **Do not reuse** in the parity path: its rank rules differ from `olsmatrix2` (drops all-zero columns, else `solve`) and `fracridge` (`tol=1e-10` on singular values). Parity needs those exact rules |
| `generate_rapid_design()` | R/oasis_hrf_recovery.R | — | **Do not reuse**: continuous onsets (GLMsingle needs TR-locked onsets) and it calls global `set.seed`. The simulator lives in the test helper and uses a local RNG |
| Prewhitening (`fmriAR`) | R/prewhitening.R | — | Out of scope: GLMsingle does not prewhiten. The nuisance layer keeps a `whiten` hook point so a later extension does not need restructuring |
| `lss(method = ...)` dispatch | R/lss.R | — | **Not** a new `lss()` method. GLMsingle needs conditions, runs and sessions, and outputs four model types, so it gets a separate entry point that shares the data conventions |

## 4. Definition of "equally accurate"

* **Tier A, parity with Python** (float32 reference): `glmsingle()` vs
  fixtures from Python `1ab54a6`, with Python's tail thresholds injected.
  Continuous outputs must agree to ≤ 1e-4 relative. Discrete choices (HRFindex,
  pcnum, FRACindex) must agree everywhere except voxels whose winning margin is
  below a float32-scaled tie threshold, and those voxels are listed in the test
  output.
* **Tier B, optimisation exactness** (float64): every fast stage vs a compact
  dense reference implementation of that stage. The references are written
  directly from the Python source in `helper-glmsingle.R` and use dense
  projectors, the stacked all-trial design, per-fraction `fracridge` and pairwise
  `calcbadness`. Continuous outputs must agree to ≤ 1e-9 relative, and discrete
  choices must agree exactly except at exact ties. This proves the speedups do
  not change the estimator. The dense references are test-only code and are not
  shipped in `R/`.
* **Tier C, statistical accuracy:** only for research-track changes. Uses a
  ground-truth simulation, as in proposal §12.

## 5. Granular plan

Each task names its deliverable and done-criterion. Phases 0–4 are the
deliverable; Phase 5 is optional research. Each task is roughly one reviewable
commit.

### Phase 0: environment and reference harness

- **0.1 R toolchain.** The container has no R. Install R ≥ 4.3 plus the
  dependencies in DESCRIPTION and devtools, and add a SessionStart hook for
  cloud sessions. *Done when* `devtools::test()` passes on the current branch.
- **0.2 Pinned Python reference.** `tools/glmsingle_ref/requirements.txt`
  pins GLMsingle@1ab54a6, fracridge 3.0, numpy, scipy and scikit-learn. Record
  whether numba JIT was used. Add `tools/` to `.Rbuildignore`. *Done when*
  GLMsingle runs end to end on a synthetic dataset.
- **0.3 Fixture writer** `make_fixtures.py`. It writes inputs and every
  intermediate: onoffR2, meanvol, the tail thresholds used, FitHRFR2(+run),
  HRFindex, pcregressors, the `ix` voxel set, glmbadness, xvaltrend, pcnum,
  rrbadness, FRACvalue, scaleoffset, and betasmd A–D. Use little-endian `.bin`
  files plus a JSON manifest. Scenarios:
  (a) defaults;
  (b) extra regressors;
  (c) 2 sessions + custom grouped `xvalscheme`;
  (d) unequal run lengths;
  (e) a condition with no repeats.
  Keep each under 2 MB in `tests/testthat/fixtures/glmsingle/`. *Done when* the
  fixtures regenerate bit-identically from a fixed seed with thresholds pinned.
- **0.4 Ground-truth simulator** in `helper-glmsingle.R`. It uses TR-locked
  onsets, condition means plus a within-condition behavioural covariate plus
  idiosyncratic trial variation, voxel-specific library HRFs, a structured
  noise pool and drift, and a local RNG. *Done when* it is deterministic and
  documented.
- **0.5 Baseline timing.** Time Python GLMsingle per stage on (a) the simulator
  at 12 runs × 220 TRs × 50k voxels and (b) the GLMsingle example data, if it
  can be downloaded. *Done when* the numbers are in
  `bench/results/glmsingle_baseline.md`.

### Phase 1: shared groundwork (small, independent commits)

- **1.1 Promote `.orthonormal_span()`** to aaa_utils and update its callers.
  *Done when* the existing voxel_hrf and lss_with_hrf tests are green.
- **1.2 Extract `src/blocked_products.h`.** No behaviour change to OASIS yet.
  *Done when* it compiles and the OASIS tests are green.
- **1.3 `glmsingle_options()`** with defaults taken from `glmsingle/defaults.py`: `fracs = seq(1, .05,
  by = -.05)`, `n_pcs = 10`, `pcstop = 1.05`, `brainthresh = c(99, 0.1)`,
  `maxpolydeg` rule, `wantautoscale`, `wantpercentbold`, `chunklen`, and so on.
  Values are validated with the existing helpers. *Done when* the unit tests
  cover every option and its rejection message.

### Phase 2: core modules (each tested at Tier B)

- **2.1 Run geometry** (`glmsingle_runs.R`): from per-run onset TR indices and
  condition labels, plus optional session ids and `xvalscheme`, build an
  immutable object holding trial order (matching designSINGLE's `cnt`
  ordering), `validcolumns`, `stimix`, condition counts, folds, and the
  repeat-availability checks that disable denoise/fracridge.
  *Test:* matches the fixture bookkeeping in every scenario.
- **2.2 HRF library** (`glmsingle_hrf.R`): vendor the TSV, then port
  `getcanonicalhrflibrary` and `getcanonicalhrf` (stimdur boxcar convolution at
  0.1 s, pchip to TR, peak normalisation). This includes a small internal pchip
  that matches scipy. *Test:* agrees with the fixture library to 1e-6
  (float32 source).
- **2.3 Nuisance bases** (`glmsingle_nuisance.R`): per run, `Q = [poly |
  extra | pc₁..pc_k]` built by sequential orthonormalisation so PC prefixes
  nest. Includes the PC derivation: poly-projected noise-pool time series,
  column normalisation, `T×T` SVD, std scaling. *Test:* `I − Q_kQ_kᵀ` equals
  the dense `make_projection_matrix([poly, extra, pcs[:k]])`, including
  rank-deficient extras.
- **2.4 Data pass** (`glmsingle_kernels.cpp` + `glmsingle_stats.R`): per run
  and voxel tile, compute `QᵀY` (all PC columns), `‖y‖²`, and the matched-filter
  `X_hᵀY` for all library HRFs, using `blocked_products.h`. Optional OpenMP
  over tiles. *Test:* agrees with dense products to 1e-12.
- **2.5 Design statistics and quadratic-form fits:** per (run, HRF, prefix
  `k`), `G = XᵀX − (XᵀQ_k)(XᵀQ_k)ᵀ` and `b = XᵀY − (XᵀQ_k)(Q_kᵀY)`, plus the
  poly-only versions for R². OLS follows the `olsmatrix2` zero-column rule.
  R² and R2run come from sufficient statistics. *Test:* betas, R² and R2run
  agree with the dense `glm_estimatemodel` reference.
- **2.6 Spectral ridge path** (`glmsingle_ridge.R` + kernel): eigendecompose
  each run's `G_r` and pool the spectra. Build the `fracridge` grid from the
  global `s_min`/`s_max`. Compute every fraction's α per voxel in one pass.
  Reconstruct coefficients for a requested trial subset. `exact_alpha = TRUE`
  is an opt-in root solve. *Test:* agrees with a direct R port of `fracridge`
  on the stacked design to 1e-10.
- **2.7 Selection** (`glmsingle_select.R`):
  (a) Compiled CV: from the run geometry build `d`, the used-trial mask and
  condition×run sums for the targets `m`, without ever materialising `W`.
  Session normalisation comes from the reference fit, with the
  `zerodiv → 0` rule.
  (b) `select_noise_regressors`.
  (c) A deterministic two-component GMM tail threshold with an eps variance
  floor.
  (d) Autoscale as vectorised 2×2 OLS, including the `h[0] < 0 → (1, 0)` rule.
  (e) Percent-BOLD scaling.
  *Test:* (a) agrees with pairwise `calcbadness` to 1e-12 on random inputs,
  including zero-SD voxels, conditions absent from training, and multi-run
  folds. (c) is within one grid step of Python on the fixtures.

### Phase 3: assembly

- **3.1 `glmsingle()`** orchestrator, about 150 lines: stages A → noise pool →
  PCs → B → PC selection → C → D, chunked over voxel tiles and grouped by HRF
  index within a tile. Returns a `glmsingle_fit` with `typea..typed` fields named
  as in GLMsingle, plus `$hrf_library`, `$options` and `$timing`.
- **3.2 `glmsingle_design(Y, event_model, baseline_model = NULL, ...)`**: the
  fmridesign front-end, mirroring `lss_design()`.
- **3.3 Stage B:** score all HRFs from the statistics, take the argmax, then
  solve once for the winner.
- **3.4 PC selection on `ix` only.** `opt$want_full_glmbadness = TRUE`
  computes all voxels. *Test:* the restricted and full paths give identical
  `xvaltrend`, `pcnum` and type C/D outputs.
- **3.5 Stages C and D:** C is OLS at `pcnum`. D is the full fraction path,
  then compiled CV on used trials, argmin per voxel, full reconstruction for the
  chosen fraction, then autoscale. Also handles the special case where
  frac = 1 is prepended.
- **3.6 Methods:** `print`, `summary` (selected pcnum, fraction histogram, HRF
  index histogram, stage timings), `coef(fit, type = "d")`, and `as.matrix`.
- **3.7 Tests:** `test-glmsingle-pipeline.R` at Tier B, stage by stage.
  `test-glmsingle-parity.R` at Tier A against the fixtures.

### Phase 4: validation, benchmarking, docs

- **4.1 Parity report:** all fixture scenarios plus the GLMsingle example data,
  run locally rather than in CI. Report beta correlations, max relative
  differences, discrete-choice agreement, and the tie-voxel list.
- **4.2 Simulation report:** the same metrics for Python and fmrilss (they
  should be identical up to ties): total and within-condition beta error,
  behavioural-effect recovery, and HRF selection accuracy. This becomes the
  Tier C baseline for Phase 5.
- **4.3 Benchmark** `bench/run_glmsingle_benchmark.R`: per-stage wall time vs
  0.5 at 1 and N threads, plus peak RSS. *Target:* ≥ 5× end to end at 1 thread.
  *Stretch:* ≥ 10×.
- **4.4 Optional:** migrate OASIS onto `blocked_products.h`, gated by the
  OASIS tests.
- **4.5 Docs:** roxygen, a `glmsingle.Rmd` vignette (usage, the fmridesign
  front-end, parity evidence, timing), `_pkgdown.yml` reference section, NEWS,
  and `inst/COPYRIGHTS` (BSD-3, Kendrick Kay) for the vendored library and
  ported algorithms.
- **4.6 `R CMD check`** clean, with the GLMsingle test suite under 60 s.

Out of v1 scope (each errors with a clear message): `hrfmodel = "optimize"`,
`wantlss`, bootstrap/xval resampling modes, figures, and hdf5 output.

### Phase 5: research track (gated on Phase 4; opt-in options only)

Each item needs a Tier C report on the 0.4 simulator, with held-out runs,
showing it is no worse on within-condition and behavioural-effect error.

- 5.1 HRF selection uncertainty: R² margin and profile curvature per voxel;
  optional shrinkage toward neighbourhood HRFs (proposal §7).
- 5.2 Continuous HRF refinement via VarPro with fixed-peak normalisation
  (proposal §7). It could reuse the SBHM library-SVD code in `R/sbhm_build.R`
  for the basis.
- 5.3 Cross-relation initialiser for regular-grid designs (proposal §8).
- 5.4 SURE/evidence fraction selection for designs without repeats (proposal
  §9), including the fixed-fraction Jacobian term.
- 5.5 Compressed library screening with a certificate (proposal §6), only if
  profiling shows the library sweep matters after 2.4.

## 6. Risks

| Risk | Mitigation |
|---|---|
| Discrete-choice flips at Tier A from float32 | Tie-margin reporting. Tier B is the real exactness gate |
| Python tail-threshold nondeterminism | Fixtures record the thresholds and tests inject them. Our GMM is deterministic |
| Rank-deficient custom extras | Rank-tolerant orthonormal span. Dense-reference test with collinear extras. Document behaviour where GLMsingle's `inv` would be ill-posed |
| Data I/O dominates, so speedup falls short | Per-stage baseline first (0.5). The ~70 removed passes over the data are the bulk of compute regardless |
| Upstream GLMsingle changes | Commit pinned in the fixture manifest. Re-audit on upstream releases |
| Scope creep into the research track | Phase 5 is gated on Phase 4 and options-only. Defaults never change without a Tier C report |

## 7. Questions for the third-party reviewer

1. Is the two-tier exactness contract (§4) sufficient evidence of "equally
   accurate" for the default estimator?
2. Is any GLMsingle behaviour missing from §5 that downstream users rely on?
   Examples: FIR diagnostic outputs, `R2run`-based HRF index, `meanvol`, and
   the output file layout.
3. Is adopting MATLAB's eps-regularised GMM (§2.3) acceptable as the single
   intentional deviation from Python?
4. Are the deferral verdicts in §2.1 for proposal §§5, 6 and 3 (Sherman–Morrison)
   correct, given typical NSD-scale designs?
