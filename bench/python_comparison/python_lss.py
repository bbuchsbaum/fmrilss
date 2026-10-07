"""Python LSS reference implementations for the fmrilss comparison benchmark.

Usage: python python_lss.py <sim_dir> [methods...]

Reads the dataset written by simulate.R and writes, for each method,
``py_<method>.bin`` (n_trials x n_vox float64, column-major) plus a
``py_timings.json`` file. Methods:

* ``nilearn_ols``   - per-trial GLM with nilearn's ``run_glm`` (OLS), LSS-N
                      design (one "other" regressor per condition), as in
                      Nilearn's beta-series example and NiBetaSeries.
* ``nilearn_ar1``   - same, with nilearn's default AR(1) noise model
                      (voxel AR coefficients binned, GLS per bin).
* ``numpy_loop``    - per-trial ``numpy.linalg.lstsq`` on the full design
                      (classic LSS, one pooled "other" regressor).
* ``numpy_loop_n``  - per-trial ``lstsq``, LSS-N design.
* ``numpy_closed``  - vectorized closed-form classic LSS (best-case NumPy).
"""

import json
import os
import sys
import time

import numpy as np


def read_mat(path, nrow, ncol):
    return np.fromfile(path, dtype="<f8").reshape((ncol, nrow)).T.copy()


def load(sim_dir):
    with open(os.path.join(sim_dir, "manifest.json")) as fh:
        m = json.load(fh)
    n, T, V = m["n_scans"], m["n_trials"], m["n_vox"]
    Y = read_mat(os.path.join(sim_dir, "Y.bin"), n, V)
    X = read_mat(os.path.join(sim_dir, "X.bin"), n, T)
    Z = read_mat(os.path.join(sim_dir, "Z.bin"), n, m["n_z"])
    motion = read_mat(os.path.join(sim_dir, "motion.bin"), n, m["n_motion"])
    with open(os.path.join(sim_dir, "condition.txt")) as fh:
        cond = np.array([ln.strip() for ln in fh if ln.strip()])
    return Y, X, np.column_stack([Z, motion]), cond


def trial_design(X, conf, i, cond=None):
    """Design for trial i: [trial, other(s), confounds]."""
    others = np.delete(np.arange(X.shape[1]), i)
    if cond is None:
        other = X[:, others].sum(axis=1, keepdims=True)
    else:
        cols = []
        for c in np.unique(cond):
            idx = others[cond[others] == c]
            if idx.size:
                cols.append(X[:, idx].sum(axis=1))
        other = np.column_stack(cols)
    return np.column_stack([X[:, i], other, conf])


def lss_nilearn(Y, X, conf, cond, noise_model):
    from nilearn.glm.first_level import run_glm

    T, V = X.shape[1], Y.shape[1]
    betas = np.empty((T, V))
    for i in range(T):
        D = trial_design(X, conf, i, cond)
        labels, results = run_glm(Y, D, noise_model=noise_model, n_jobs=1)
        for lab, res in results.items():
            mask = labels == lab
            betas[i, mask] = res.theta[0]
    return betas


def lss_numpy_loop(Y, X, conf, cond=None):
    T = X.shape[1]
    betas = np.empty((T, Y.shape[1]))
    for i in range(T):
        D = trial_design(X, conf, i, cond)
        coef, *_ = np.linalg.lstsq(D, Y, rcond=None)
        betas[i] = coef[0]
    return betas


def lss_numpy_closed(Y, X, conf):
    Q, _ = np.linalg.qr(conf)
    C = X - Q @ (Q.T @ X)
    Yr = Y - Q @ (Q.T @ Y)
    total = C.sum(axis=1)
    CtY = C.T @ Yr
    tY = total @ Yr
    CtC = (C * C).sum(axis=0)
    CtT = C.T @ total
    bt2 = total @ total - 2 * CtT + CtC
    ctb = CtT - CtC
    alpha = ctb / bt2
    den = CtC - ctb ** 2 / bt2
    return ((1 + alpha)[:, None] * CtY - alpha[:, None] * tY[None, :]) / den[:, None]


def main():
    sim_dir = sys.argv[1]
    methods = sys.argv[2:] or [
        "numpy_closed", "numpy_loop", "numpy_loop_n", "nilearn_ols", "nilearn_ar1"
    ]
    Y, X, conf, cond = load(sim_dir)
    runners = {
        "numpy_closed": lambda: lss_numpy_closed(Y, X, conf),
        "numpy_loop": lambda: lss_numpy_loop(Y, X, conf),
        "numpy_loop_n": lambda: lss_numpy_loop(Y, X, conf, cond),
        "nilearn_ols": lambda: lss_nilearn(Y, X, conf, cond, "ols"),
        "nilearn_ar1": lambda: lss_nilearn(Y, X, conf, cond, "ar1"),
    }
    timing_path = os.path.join(sim_dir, "py_timings.json")
    timings = {}
    if os.path.exists(timing_path):
        with open(timing_path) as fh:
            timings = json.load(fh)
    for name in methods:
        t0 = time.perf_counter()
        B = runners[name]()
        timings[name] = time.perf_counter() - t0
        B.T.astype("<f8").tofile(os.path.join(sim_dir, f"py_{name}.bin"))
        print(f"{name:14s} {timings[name]:9.3f}s", flush=True)
    with open(timing_path, "w") as fh:
        json.dump(timings, fh, indent=2)


if __name__ == "__main__":
    main()
