# fmrilss vs Python LSS implementations

Speed and accuracy of `fmrilss::lss()` compared with the LSS approaches used
in the Python fMRI ecosystem, on one shared simulated dataset.

| Python method | What it is |
|---|---|
| `nilearn_ols` | One `nilearn.glm.first_level.run_glm` fit per trial with an LSS-N design (the trial, one "other trials" regressor per condition, confounds), as in Nilearn's beta-series example and NiBetaSeries. |
| `nilearn_ar1` | Same, with Nilearn's default AR(1) noise model (per-voxel AR(1) estimates binned and fitted by GLS per bin). |
| `numpy_loop`, `numpy_loop_n` | One `numpy.linalg.lstsq` per trial: classic LSS / LSS-N. |
| `numpy_closed` | Vectorized closed-form classic LSS in NumPy (an idealized fast Python baseline; we know of no package that ships it). |

| fmrilss call | Model |
|---|---|
| `fmrilss_lss` | `lss(Y, X, Z, motion)`: classic LSS |
| `fmrilss_lss_n` | `+ trial_groups = condition`: LSS-N |
| `fmrilss_lss_n_cpp` | LSS-N with `method = "cpp_optimized"` |
| `fmrilss_lss_n_ar1` | LSS-N + `prewhiten = list(method = "ar", p = 1)` |
| `fmrilss_lss_n_ar1_voxel` | LSS-N + voxel-adaptive AR(1) (`pooling = "voxel"`) |
| `fmrilss_lss_n_ridge` | LSS-N + `ridge = c(1, 0)` |
| `fmrilss_lss_n_ridge_ar1_voxel` | LSS-N + ridge + voxel-adaptive AR(1) |

## Simulation

`simulate.R`: TR 2 s, 400 scans, two conditions (true means 1.0 and 0.25,
trial-level SD 0.5, voxel gain ~N(1, 0.3)), Glover HRF, AR(1) noise with
voxel-specific phi ~ U(0.1, 0.6), DCT drift, and six random-walk motion
traces with voxel-specific effects. The confound model is an intercept, a
128 s DCT high-pass basis and the motion traces.

* `rapid`: ITI ~ U(2, 8) s (150 trials). `slow`: ITI ~ U(8, 14) s (69 trials).
* Accuracy is the mean per-voxel correlation between estimated and true trial
  betas (`r_trial`), the same within condition (`r_within`, isolating
  trial-to-trial variability), and RMSE.

## Results

Linux, Intel Xeon 2.1 GHz (4 cores), R 4.3.3 with OpenBLAS 0.3.26, Python
3.13 with NumPy 2.5.3 and Nilearn 0.14.1, fmriAR 0.3.3. fmrilss times are
the median of 3 runs. Full table: `results/all_results.csv`.

### Speed (seconds)

| Model | Implementation | rapid, 2k vox | slow, 2k vox | rapid, 20k vox |
|---|---|---:|---:|---:|
| LSS-N, OLS | **fmrilss** | **0.017** | **0.009** | **0.054** |
| | Nilearn `run_glm` loop | 1.89 | 1.63 | 11.5 |
| | NumPy `lstsq` loop | 1.87 | 0.86 | 34.1 |
| LSS, OLS | **fmrilss** | **0.010** | **0.007** | **0.054** |
| | NumPy closed form | 0.021 | 0.008 | 0.31 |
| LSS-N, AR(1) | **fmrilss**, global AR | **0.064** | **0.068** | **0.26** |
| | **fmrilss**, voxel-adaptive AR | 0.68 | 0.37 | 1.37 |
| | Nilearn AR(1) loop | 19.8 | 8.3 | 113 |

At 20k voxels fmrilss is ~210x faster than the Nilearn per-trial loop for
OLS LSS-N, and 80-430x faster than Nilearn's AR(1) loop.

### Accuracy (rapid design, 20k voxels)

| Model | r_trial | r_within | RMSE |
|---|---:|---:|---:|
| LSS (fmrilss = NumPy) | 0.530 | 0.432 | 1.066 |
| LSS-N OLS (fmrilss = Nilearn OLS) | 0.528 | 0.442 | 1.032 |
| Nilearn AR(1) | 0.553 | 0.465 | 0.960 |
| fmrilss AR(1), global | 0.549 | 0.461 | 0.969 |
| fmrilss AR(1), voxel-adaptive | **0.553** | **0.465** | 0.956 |
| fmrilss ridge `c(1, 0)` | 0.528 | 0.442 | 0.687 |
| fmrilss ridge + voxel-adaptive AR(1) | **0.553** | **0.465** | **0.660** |

The 2k-voxel rapid and slow designs show the same ordering
(`results/all_results.csv`).

* Algebraically identical estimators agree to < 5e-13 (fmrilss LSS-N vs
  Nilearn OLS and NumPy loops; fmrilss LSS vs NumPy).
* Voxel-adaptive AR(1) matches or slightly exceeds Nilearn's AR(1) accuracy
  in all three scenarios.
* Ridge (shrinking only the trial-of-interest coefficient) lowers RMSE by
  about a third without changing the pattern correlation: it trades a small
  bias for a large variance reduction, which matters for analyses that use
  beta magnitudes. It is not a default; pick the amount by validation on
  your design.
* Before this change, fmrilss estimated the AR model from residuals of the
  full trial-wise design; with 150 trials in 400 scans that gave phi = -0.32
  (truth ~0.35) and AR(1) betas *less* accurate than OLS (r_trial 0.515).
  The default `prewhiten$residual_model = "aggregate"` fixes this.

## Reproduce

```sh
pip install numpy nilearn
R CMD INSTALL .            # from the package root
bash bench/python_comparison/run_all.sh /tmp/lss_bench
```
