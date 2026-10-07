GLMsingle acceleration: independent due diligence and reproducibility
Reviewed 2026-10-07

PURPOSE AND SCOPE
This package contains the original ridge microbenchmark, its saved result and a
fresh rerun, the reviewed plan's numerical checker, and new focused checks of
fractional interpolation, finite precision, nuisance rank behavior, HRF support,
score-only algebra, and PC-prefix coefficient computation.

These are source audits and kernel tests. There is no complete glmsingle_fast
implementation here, no end-to-end R/C++ speed measurement, and no measured
real-fMRI accuracy improvement. Numerical differences in this package compare
implementations; they are not errors against physiological ground truth.

REVIEWED REVISIONS
Plan repository: https://github.com/bbuchsbaum/fmrilss
Branch: claude/charming-wozniak-752kte
Plan commit: 40cb4cfbbc5ddb322fc590b56927cf6ccb44f443
Plan: .planning/GLMSINGLE_FAST_PLAN.md
Plan blob: d7d88345161c6ff442ada9ad4f0b51dbbc82d435
Reviewed checker: .planning/glmsingle_checks/verify_core_identities.py
Checker blob: 00d183235366a3049fda489e1dd552417e541a55
GLMsingle commit: 1ab54a65edd3ea41a6133d4b4ecb78a9c7296684
fracridge: the actual PyPI 3.0 wheel, included unchanged.
Wheel SHA256: fc35bf291e6f217735600b5633e6d598186e2d07c2624f499b2fb44cabfaa0bb

HOW TO RUN
Extract the ZIP, then enter the glmsingle_due_diligence directory.
Tested with Python 3.12.14 on Linux x86_64. Install requirements.txt in a fresh
virtual environment if these dependencies are not already available:

    python3 -m pip install -r requirements.txt
    python3 run_all.py

The checks require no network after the Python dependencies are installed.
The source functions are bundled, so no GLMsingle or fracridge installation is
required. Numerical reference calls use jit=False. Timing scripts limit BLAS to
one thread. Timing varies by hardware and concurrent machine activity.
Run serially for meaningful timing; run_all.py does this.

Recorded results are preserved under recorded_results/. Reruns write results
under tmp/glmsingle/ and logs under rerun_logs/. SHA256SUMS.json records every
distributed file other than itself. To verify an untouched extraction:

    python3 verify_manifest.py

WHAT EACH CHECK ESTABLISHES

1. Original benchmark: tmp/glmsingle/check_algebra.py
   This file is byte-for-byte identical to original/check_algebra.py.
   It uses six runs, 220 timepoints and 50 trials per run, 128 voxels, twenty
   fractional candidates and five nuisance regressors. It times repeated
   global Gram-SVD/product/reconstruction versus cached run-block calculations
   for the SAME precomputed voxel-specific penalties. Both reconstruct all
   coefficients. The original timing was 0.363059 s versus 0.008346 s (43.503x).
   The fresh recorded rerun was 0.373214 s versus 0.008186 s (45.593x).
   The timed outputs differ by about 3.19e-15 relative.
   Penalty conversion, design construction, nuisance discovery, HRF selection,
   repeat-CV, diagnostics and file I/O are OUTSIDE the timed section. The
   conversion helper solves the mathematical fraction equation; it does not
   reproduce fracridge's grid interpolation. This limitation does not change
   the equality of the fixed-penalty systems in the timed comparison.
   library_dur3.0.npz is the original retained geometry artifact: the bundled
   library sampled for duration 3 s and TR=4/3 s, peak-normalized, with an SVD
   basis for separate algebra checks. It is not a complete Python32 fixture.

2. Reviewed checker: tmp/glmsingle/review_precision/run_reviewer_checks.py
   It executes the ORIGINAL reviewed test body unchanged. Only dependency
   loading is replaced so that the full original fracridge 3.0 module and real
   pinned calcbadness/zerodiv functions can be imported locally. Result: blocked
   coefficients agree at 8.12e-13 relative; repeated-trial CV agrees at about
   2.7e-16 for both tested fold schemes. This confirms the reviewer's stated
   numbers on that fixture. Its direct block SVD and log1p implementation do
   not test all boundary behavior of upstream Gram SVD and log(1+alpha).
   An additional constant-response probe documents zerodiv's in-place divisor
   mutation; naively zeroing every candidate changes this degenerate case.

3. Precision stress checks: tmp/glmsingle/review_precision/check_precision.py
   The ORIGINAL fracridge 3.0 _do_svd and fracridge function ASTs and constants
   are executed without editing their arithmetic. Independent block statistics
   are optionally injected at the decomposition boundary. This verifies grid,
   interpolation and reconstruction against the original implementation.
   Well-behaved float64 pooled-grid results agree around 4.31e-15 relative.
   A representative exact-root versus interpolated-alpha discrepancy reaches
   3.79%; the exact root is not a literal fracridge replacement.
   Synthetic conditioning stress tests show that a universal 1e-4 Python32
   versus float64 agreement claim, or universal 1e-9 blocked-double claim, is
   invalid. These are intentionally constructed tests, not a prevalence study.
   The constructed near-tie example demonstrates sensitivity, not its frequency.

4. Source-edge probe: tmp/glmsingle/source_edge_probe.py
   Required GLMsingle files are checked against their pinned Git blob hashes,
   then loaded unchanged. The script demonstrates: nuisance projector rank
   differs from a conventional QR cutoff; duplicate nonzero task columns can
   make olsmatrix2 raise; corrected score-only task SSE matches explicit
   predictions with custom nuisances; exact sampled HRF support is not always
   25 TRs; and constant-candidate autoscale cannot use an unconditional 2x2 solve.

5. PC-prefix benchmark: tmp/glmsingle/due_diligence/pc_prefix_benchmark.py
   Compares identical fixed-design coefficient problems with 63 trials and
   eleven nuisance prefixes: cached Cholesky solves, rank-one beta updates,
   and stacked precompiled linear maps. Factor construction, initial data
   projection and repeat-CV are excluded from prefix timings; data projection
   is reported separately. All methods materialize all prefix coefficients.
   At 8,192 voxels, the recorded medians are 97.165, 29.056 and 32.194 ms,
   respectively (about 3.34x and 3.02x faster than repeated cached solves).
   This identifies a useful optimization candidate; it is not a measured
   PC-search or full GLMsingle speedup. The two faster methods' order can vary.

INTERPRETING THE RESULTS
Keep three questions separate: source behavior, mathematical equivalence, and
statistical accuracy. Pooled run spectra plus the same interpolation preserve
the intended fractional estimator in exact arithmetic; finite precision and
rank decisions still need explicit policies. Float64 agreement with a newly
written R reference does not by itself establish parity with pinned Python.
Global adaptive decisions (noise threshold, PC subspace, PC count) can propagate
small numerical differences to many voxels. Test fixed intermediate quantities
and adaptive end-to-end behavior separately.

LICENSES AND ATTRIBUTION
The original GLMsingle source files and sampled HRF data retain their BSD-3
license and Kendrick Kay attribution (licenses/GLMsingle_LICENSE.txt).
fracridge 3.0 retains its BSD-2 license and Kendrick Kay/Ariel Rokem attribution
(licenses/fracridge_LICENSE.txt; also present inside the original wheel).
The reviewed checker is copied from bbuchsbaum/fmrilss at the revision above;
that package declares GPL-3 in DESCRIPTION. The GPL-3 text is in licenses/.
The different items in this archive retain their respective provenance; the
new probes are supplied for evaluating this review, not as upstream code.
