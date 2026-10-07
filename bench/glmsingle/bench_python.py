#!/usr/bin/env python3
"""Time pinned Python GLMsingle on simulated data and save inputs for R.

Usage: OPENBLAS_NUM_THREADS=1 python -I bench_python.py OUTDIR N_RUNS N_TIME N_VOX
Writes OUTDIR/inputs (raw float32 data and designs), OUTDIR/python_result.json
and OUTDIR/python_betas_typed.f32 for comparison with bench/run_glmsingle_benchmark.R.
"""
import json
import os
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "tools", "glmsingle_ref"))
from make_fixtures import simulate  # noqa: E402
from glmsingle.glmsingle import GLM_single  # noqa: E402


def main(outdir, n_runs, n_time, n_vox):
    os.makedirs(os.path.join(outdir, "inputs"), exist_ok=True)
    design, data, _ = simulate(seed=123, n_runs=n_runs, n_time=n_time, n_vox=n_vox, n_cond=40)
    for r, (D, Y) in enumerate(zip(design, data)):
        D.astype("<f4").tofile(os.path.join(outdir, "inputs", f"design_run{r + 1}.f32"))
        Y.astype("<f4").tofile(os.path.join(outdir, "inputs", f"data_run{r + 1}.f32"))
    meta = {"n_runs": n_runs, "n_time": n_time, "n_vox": n_vox, "n_cond": design[0].shape[1],
            "tr": 1.0, "stimdur": 3.0, "brainR2": 5.0, "pcR2cutoff": 5.0}
    params = dict(wantfileoutputs=[0, 0, 0, 0], wantmemoryoutputs=[1, 1, 1, 1],
                  brainR2=meta["brainR2"], pcR2cutoff=meta["pcR2cutoff"])
    with tempfile.TemporaryDirectory() as tmp:
        os.chdir(tmp)
        t0 = time.perf_counter()
        res = GLM_single(params).fit(design, data, meta["stimdur"], meta["tr"],
                                     outputdir=os.path.join(tmp, "out"))
        elapsed = time.perf_counter() - t0
    betas = np.asarray(res["typed"]["betasmd"]).reshape(n_vox, -1)
    betas.astype("<f4").tofile(os.path.join(outdir, "python_betas_typed.f32"))
    meta.update(python_seconds=elapsed, n_trials=betas.shape[1],
                pcnum=int(res["typed"]["pcnum"]))
    with open(os.path.join(outdir, "python_result.json"), "w") as f:
        json.dump(meta, f, indent=1)
    print(json.dumps(meta))


if __name__ == "__main__":
    main(sys.argv[1], *map(int, sys.argv[2:5]))
