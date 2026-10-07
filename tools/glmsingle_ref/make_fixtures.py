#!/usr/bin/env python3
"""Generate glmsingle() test fixtures from pinned Python GLMsingle.

Reference: cvnlab/GLMsingle @ 1ab54a65edd3ea41a6133d4b4ecb78a9c7296684 with
fracridge 3.0 (see requirements.txt). Usage:

    python -I make_fixtures.py OUTDIR

Each scenario writes inputs and GLMsingle outputs as gzipped little-endian
float32 arrays plus manifest.json. Scenarios where fmrilss's default
behaviour differs from GLMsingle are also run with GLMsingle patched to that
behaviour (variant "fmrilss"); the patch is defined in `patched_glmsingle()`:

  * user extra regressors are kept in every GLMdenoise / ridge fit, including
    the 0-PC fit (GLMsingle adds them only when at least one PC is used);
  * zerodiv() no longer modifies its divisor in place, so voxels with zero
    beta SD get z = 0 for every candidate in cross-validation.

Tail thresholds (brainR2, pcR2cutoff) are always passed explicitly so that
GLMsingle's unseeded Gaussian-mixture fit does not enter the fixtures.
"""
import gzip
import hashlib
import importlib
import json
import os
import sys
import tempfile
import types

import numpy as np

import glmsingle
import glmsingle.ssq.calcbadness as cb_mod
from glmsingle.glmsingle import GLM_single
from glmsingle.hrf.gethrf import getcanonicalhrflibrary

PINNED = "1ab54a65edd3ea41a6133d4b4ecb78a9c7296684"


def patched_glmsingle():
    """GLM_single with fmrilss's default extras / zero-SD behaviour."""
    path = os.path.join(os.path.dirname(glmsingle.__file__), "glmsingle.py")
    src = open(path).read()
    for old in ("if n_pc > 0:", "if pcnum > 0:"):
        assert src.count(old) == 1, old
        src = src.replace(old, "if True:")
    mod = types.ModuleType("glmsingle.glmsingle_fmrilss")
    mod.__package__ = "glmsingle"
    mod.__file__ = path
    exec(compile(src, path, "exec"), mod.__dict__)
    return mod.GLM_single


def zerodiv_copy(x, y, val=0, wantcaution=1):
    """zerodiv without the in-place modification of y."""
    from glmsingle.utils import zerodiv as zmod
    return zmod.zerodiv(x, np.array(y, copy=True), val, wantcaution)


def simulate(seed, n_runs=4, n_time=100, n_vox=64, n_cond=12, tr=1.0,
             stimdur=3.0, isi=(3, 5), n_extras=0, collinear_extra=False,
             unrepeated=0, constant_voxel=False, unequal=False):
    rng = np.random.default_rng(seed)
    lib = getcanonicalhrflibrary(stimdur, tr).T.astype(np.float64)
    lib = lib / lib.max(axis=0)
    n_times = [n_time + (13 * r if unequal else 0) for r in range(n_runs)]
    # latent structured noise shared across voxels (what GLMdenoise removes)
    loadings = rng.normal(0, 1, (3, n_vox))
    hrf_idx = rng.integers(0, lib.shape[1], n_vox)
    signal_vox = rng.random(n_vox) < 0.6
    cond_means = rng.normal(1.0, 0.8, (n_cond + unrepeated, n_vox))
    baseline = rng.uniform(500, 1500, n_vox)
    baseline[: n_vox // 10] = rng.uniform(5, 20, n_vox // 10)  # dim voxels
    design, data, extras = [], [], []
    next_unrep = n_cond
    for r in range(n_runs):
        T = n_times[r]
        onsets = []
        t = int(rng.integers(2, 5))
        while t < T - 12:
            onsets.append(t)
            t += int(rng.integers(isi[0], isi[1] + 1))
        conds = rng.integers(0, n_cond, len(onsets))
        if unrepeated and r == 0:
            k = min(unrepeated, len(onsets))
            conds[:k] = np.arange(next_unrep, next_unrep + k)
        D = np.zeros((T, n_cond + unrepeated))
        D[onsets, conds] = 1
        design.append(D)
        Y = np.zeros((n_vox, T))
        for v in range(n_vox):
            stick = np.zeros(T)
            amps = cond_means[conds, v] + rng.normal(0, 0.4, len(onsets))
            stick[onsets] = amps * signal_vox[v]
            Y[v] = np.convolve(stick, lib[:, hrf_idx[v]])[:T]
        latent = np.cumsum(rng.normal(0, 1, (3, T)), axis=1) * 0.3 + rng.normal(0, 1, (3, T))
        tt = np.linspace(-1, 1, T)
        drift = np.outer(rng.normal(0, 2, n_vox), tt) + np.outer(rng.normal(0, 1, n_vox), tt**2)
        noise = rng.normal(0, 1.0, (n_vox, T)) + loadings.T @ latent
        Y = baseline[:, None] * (1 + 0.01 * (Y + noise + drift))
        if n_extras:
            E = rng.normal(0, 1, (T, n_extras)).cumsum(axis=0) * 0.1
            if collinear_extra:
                E = np.c_[E, E[:, 0] + 1e-9 * rng.normal(0, 1, T)]
            Y = Y + baseline[:, None] * 0.01 * (rng.normal(0, 1, (n_vox, E.shape[1])) @ E.T)
            extras.append(E.astype(np.float32).astype(np.float64))
        if constant_voxel:
            Y[-1] = 0.0
        data.append(Y.astype(np.float32))
    return design, data, (extras if n_extras else None)


SCENARIOS = {
    "defaults": dict(sim=dict(seed=1), params={}),
    "extras": dict(sim=dict(seed=2, n_extras=4), params={}, variants=True),
    "extras_pc0": dict(sim=dict(seed=3, n_extras=4), params={"pcstop": 0}, variants=True),
    "sessions": dict(sim=dict(seed=4, n_runs=6),
                     params={"sessionindicator": np.array([1, 1, 1, 2, 2, 2]),
                             "xvalscheme": [np.array([0, 1]), np.array([2, 3]), np.array([4, 5])]}),
    "unequal": dict(sim=dict(seed=5, unequal=True), params={}),
    "unrepeated": dict(sim=dict(seed=6, unrepeated=5), params={}),
    "fast_events": dict(sim=dict(seed=7, isi=(1, 2)), params={}),
    "collinear_extras": dict(sim=dict(seed=8, n_extras=3, collinear_extra=True), params={}, variants=True),
    "zero_voxel": dict(sim=dict(seed=9, constant_voxel=True), params={}, variants=True),
    "single_frac": dict(sim=dict(seed=10), params={"fracs": np.array([0.4])}),
    "no_library": dict(sim=dict(seed=11), params={"wantlibrary": 0}),
}


def write_array(outdir, name, arr, manifest):
    arr = np.asarray(arr)
    if arr.dtype == bool:
        arr = arr.astype(np.float32)
    arr = np.ascontiguousarray(arr.astype("<f4"))
    key = (arr.shape, hashlib.sha256(arr.tobytes()).hexdigest())
    seen = manifest.setdefault("_seen", {})
    if key in seen:  # identical array already written: reference it
        manifest["arrays"][name] = {"file": seen[key], "shape": list(arr.shape)}
        return
    fn = name + ".bin.gz"
    seen[key] = fn
    with gzip.open(os.path.join(outdir, fn), "wb", compresslevel=9) as f:
        f.write(arr.tobytes(order="C"))
    manifest["arrays"][name] = {"file": fn, "shape": list(arr.shape)}


def squeeze_vox(x):
    x = np.asarray(x)
    if x.ndim >= 3 and x.shape[1] == 1 and x.shape[2] == 1:
        x = x.reshape((x.shape[0],) + x.shape[3:])
    return x


def run_glmsingle(cls, design, data, params, tr, stimdur):
    with tempfile.TemporaryDirectory() as tmp:
        cwd = os.getcwd()
        os.chdir(tmp)
        try:
            p = dict(wantfileoutputs=[0, 0, 0, 0], wantmemoryoutputs=[1, 1, 1, 1])
            p.update(params)
            g = cls(p)
            res = g.fit([d.copy() for d in design], [x.copy() for x in data], stimdur, tr,
                        outputdir=os.path.join(tmp, "out"))
        finally:
            os.chdir(cwd)
    return res


def main(outdir):
    os.makedirs(outdir, exist_ok=True)
    tr, stimdur = 1.0, 3.0
    patched_cls = patched_glmsingle()
    original_zerodiv = cb_mod.zerodiv
    index = {"pinned_commit": PINNED, "scenarios": {}}
    for name, spec in SCENARIOS.items():
        design, data, extras = simulate(**spec["sim"])
        params = dict(spec["params"])
        params.setdefault("brainR2", 5.0)
        params.setdefault("pcR2cutoff", 5.0)
        if extras is not None:
            params["extra_regressors"] = extras
        variants = {"upstream": False}
        if spec.get("variants"):
            variants["fmrilss"] = True
        sdir = os.path.join(outdir, name)
        os.makedirs(sdir, exist_ok=True)
        manifest = {"tr": tr, "stimdur": stimdur, "arrays": {}, "variants": list(variants),
                    "params": {k: (v.tolist() if isinstance(v, np.ndarray) else
                                   [x.tolist() for x in v] if isinstance(v, list) and k == "xvalscheme" else v)
                               for k, v in params.items() if k != "extra_regressors"}}
        for r, (D, Y) in enumerate(zip(design, data)):
            write_array(sdir, f"design_run{r + 1}", D, manifest)
            write_array(sdir, f"data_run{r + 1}", Y, manifest)
            if extras is not None:
                write_array(sdir, f"extras_run{r + 1}", extras[r], manifest)
        for variant, patched in variants.items():
            cb_mod.zerodiv = zerodiv_copy if patched else original_zerodiv
            cls = patched_cls if patched else GLM_single
            res = run_glmsingle(cls, design, data, params, tr, stimdur)
            for tname, tres in res.items():
                if patched and tname in ("typea", "typeb"):
                    continue  # identical to upstream by construction
                for key, val in tres.items():
                    if key == "pcregressors" and val is not None:
                        for r, pc in enumerate(val):
                            write_array(sdir, f"{variant}_{tname}_pcregressors_run{r + 1}", pc, manifest)
                        continue
                    if val is None or isinstance(val, (str, dict)):
                        continue
                    arr = np.asarray(val)
                    if arr.dtype == object or arr.size == 0:
                        continue
                    write_array(sdir, f"{variant}_{tname}_{key}", squeeze_vox(arr), manifest)
        cb_mod.zerodiv = original_zerodiv
        manifest.pop("_seen", None)
        with open(os.path.join(sdir, "manifest.json"), "w") as f:
            json.dump(manifest, f, indent=1)
        index["scenarios"][name] = sorted(variants)
        print("wrote", name, flush=True)
    with open(os.path.join(outdir, "index.json"), "w") as f:
        json.dump(index, f, indent=1)


if __name__ == "__main__":
    main(sys.argv[1])
