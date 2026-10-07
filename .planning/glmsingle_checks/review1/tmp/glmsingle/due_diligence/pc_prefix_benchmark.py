#!/usr/bin/env python3
"""Portable PC-prefix coefficient-kernel benchmark, not GLMsingle end to end.

Requirements: Python 3, NumPy, SciPy, threadpoolctl.
Run: python pc_prefix_benchmark.py --output pc_prefix_results.json

Three algebraically equivalent implementations consume the same b=D.T@Y and
u=PC.T@Y, then materialize all 11 sets of coefficients. Geometry/factorization,
noise-pool PCA, CV, HRF selection, I/O, and final output writing are not timed.
The shared data contractions are timed separately and must not be counted as
an incremental-method saving. Prefix-map construction is excluded as well.

The embedded peak-normalized HRF is candidate index 10 (one-based 11) from the
sampled GLMsingle library used in this review, with stimulus duration 3 s and
TR 4/3 s. It removes external file/path dependencies. HRF-library provenance:
cvnlab/GLMsingle, BSD-3-Clause.
"""
import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
import scipy.linalg as la
from scipy.linalg.blas import dger
from threadpoolctl import threadpool_info, threadpool_limits

HRF = np.array([
    0.0, 0.0029637404981514345, 0.12525070852282402, 0.5548706295063759,
    1.0, 0.9810209358115392, 0.6120363274116895, 0.2737876441307163,
    0.08895735880766407, 0.008832374219937003, -0.02617688515055813,
    -0.04587333209700868, -0.06009953799766916, -0.07076013706372968,
    -0.07765446059446478, -0.08052036848990382, -0.07954149962813171,
    -0.07532542108784226, -0.06872961887508551, -0.06067693772207157,
    -0.052013037812258134, -0.04342244674107482, -0.035394813751147836,
    -0.028232453730825852, -0.02207898932525051, -0.016957354660429594,
    -0.012809294689193314, -0.009528985315428003, -0.006989143228762182,
    -0.0050594707924353065, -0.003618170663917071, -0.002558215283183205,
    -0.0017896657602778352, -0.0012396181773880596, -0.0008506501148727532,
    -0.0005786334804985222, -0.0003903633692309231, -0.00026130632612372166,
    -0.00014769079989041987, -5.061655101086586e-05,
], dtype=np.float64)


def run(seed=423, repetitions=7, voxel_counts=(1024, 8192)):
    rng = np.random.default_rng(seed)
    T, n, k = 226, 63, 10
    lag = np.arange(T)[:, None] - 3*np.arange(n)[None, :]
    X = np.where((lag >= 0) & (lag < len(HRF)),
                 HRF[np.clip(lag, 0, len(HRF)-1)], 0)
    poly = la.qr(np.polynomial.legendre.legvander(np.linspace(-1, 1, T), 3),
                 mode='economic')[0]
    D = X-poly@(poly.T@X)
    pcraw = rng.normal(size=(T, k))
    pc = la.qr(pcraw-poly@(poly.T@pcraw), mode='economic')[0]
    C, G = D.T@pc, D.T@D
    chol, vs, ds, maps = [], [], [], []
    for j in range(k+1):
        Gj = G-C[:, :j]@C[:, :j].T
        cf = la.cho_factor(Gj, lower=True)
        chol.append(cf)
        inv = la.cho_solve(cf, np.eye(n))
        mat = np.zeros((n, n+k))
        mat[:, :n] = inv
        mat[:, n:n+j] = -inv@C[:, :j]
        maps.append(mat)
        if j < k:
            v = la.cho_solve(cf, C[:, j])
            vs.append(v)
            ds.append(1-C[:, j]@v)
    stackmap = np.vstack(maps)

    def fresh(b, u):
        bj = np.array(b, order='F', copy=True)
        out = np.empty((k+1, n, b.shape[1]))
        for j in range(k+1):
            if j:
                bj = dger(-1., C[:, j-1], u[j-1, :], a=bj, overwrite_a=1)
            out[j] = la.cho_solve(chol[j], bj, check_finite=False)
        return out

    def incremental(b, u):
        beta = la.cho_solve(chol[0], b, check_finite=False)
        out = np.empty((k+1, n, b.shape[1]))
        out[0] = beta
        for j in range(k):
            e = (C[:, j]@beta-u[j])/ds[j]
            beta = dger(1., vs[j], e, a=beta, overwrite_a=1)
            out[j+1] = beta
        return out

    def batched(b, u):
        return (stackmap@np.vstack([b, u])).reshape(k+1, n, -1)

    def timing(fn, *args):
        fn(*args)
        times = []
        for _ in range(repetitions):
            start = time.perf_counter()
            fn(*args)
            times.append(time.perf_counter()-start)
        return {'median_seconds': float(np.median(times)), 'samples_seconds': times}

    report = {
        'benchmark': 'PC-prefix coefficient kernel; not GLMsingle end to end',
        'seed': seed, 'rng': 'numpy.random.default_rng (PCG64)',
        'dtype': 'float64', 'blas_threads_requested': 1,
        'warmup_calls_per_timed_function': 1, 'timing_repetitions': repetitions,
        'timing_clock': 'time.perf_counter',
        'software': {'python': sys.version, 'numpy': np.__version__,
                     'scipy': scipy.__version__, 'platform': platform.platform(),
                     'machine': platform.machine(), 'threadpools': threadpool_info()},
        'geometry': {
            'T': T, 'n_trials': n, 'n_candidate_PCs': k,
            'n_prefixes_including_zero': k+1, 'event_spacing_TRs': 3,
            'polynomial_degrees': [0, 1, 2, 3],
            'hrf_support_samples': len(HRF),
            'hrf_exact_nonzeros': int(np.count_nonzero(HRF)),
            'hrf_library_candidate_zero_based': 10, 'tr_seconds': 4/3,
            'stimulus_duration_seconds': 3,
            'hrf_float64_little_endian_sha256': hashlib.sha256(
                HRF.astype('<f8').tobytes()).hexdigest(),
            'pc_construction': 'iid N(0,1) columns, polynomial-residualized, QR-orthonormalized',
            'condition_number_G': float(np.linalg.cond(G)),
            'min_schur_denominator': float(min(ds)),
        },
        'timed_operations': {
            'shared_projection': 'Compute D.T@Y and PC.T@Y; reported separately.',
            'fresh_cached_cholesky': 'Increment RHS by rank one, then cached Cholesky solve at each prefix; materialize all prefix betas.',
            'incremental_beta': 'Initial cached solve; then exact rank-one beta updates; materialize all prefix betas.',
            'batched_prefix_maps': 'Stack b and u; one precompiled-map GEMM; materialize all prefix betas.',
        },
        'excluded_operations': [
            'HRF/design construction', 'polynomial/PC basis construction',
            'all Cholesky factorizations', 'all prefix-map construction',
            'noise-pool selection and PCA', 'repeat-CV loss evaluation',
            'HRF selection', 'fractional-ridge selection', 'disk input/output',
            'initial shared projection from prefix timings (timed separately)',
        ],
        'comparison_reference': 'Fresh cached-Cholesky prefix solutions for identical G, b, PCs; not Python GLMsingle.',
        'cases': [],
    }
    methods = [('fresh_cached_cholesky', fresh),
               ('incremental_beta', incremental), ('batched_prefix_maps', batched)]
    for V in voxel_counts:
        Y = np.asarray(rng.normal(size=(T, V)), order='F')
        b, u = D.T@Y, pc.T@Y
        results = {name: fn(b, u) for name, fn in methods}
        baseline = results['fresh_cached_cholesky']
        errors = {name: {
            'relative_frobenius_error': float(np.linalg.norm(res-baseline)/np.linalg.norm(baseline)),
            'maximum_absolute_error': float(np.max(np.abs(res-baseline))),
        } for name, res in results.items()}
        del results, baseline
        timings = {name: timing(fn, b, u) for name, fn in methods}
        projection = timing(lambda y: (D.T@y, pc.T@y), Y)
        t0 = timings['fresh_cached_cholesky']['median_seconds']
        report['cases'].append({
            'Y_shape': [T, V], 'Y_layout': 'Fortran contiguous',
            'output_shape': [k+1, n, V], 'output_bytes_per_method': (k+1)*n*V*8,
            'shared_projection': projection, 'prefix_timings': timings,
            'errors': errors,
            'kernel_speedups_vs_fresh_cached_cholesky': {
                name: t0/info['median_seconds'] for name, info in timings.items()
            },
        })
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path,
                        default=Path(__file__).with_name('pc_prefix_results.json'))
    parser.add_argument('--seed', type=int, default=423)
    parser.add_argument('--repetitions', type=int, default=7)
    parser.add_argument('--voxels', nargs='+', type=int, default=[1024, 8192])
    args = parser.parse_args()
    if args.repetitions < 1 or any(v < 1 for v in args.voxels):
        parser.error('repetitions and voxel counts must be positive')
    with threadpool_limits(limits=1):
        report = run(args.seed, args.repetitions, tuple(args.voxels))
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({
        'results': str(args.output.resolve()),
        'summary': [{
            'V': case['Y_shape'][1],
            'prefix_median_seconds': {name: item['median_seconds'] for name, item in case['prefix_timings'].items()},
            'speedups': case['kernel_speedups_vs_fresh_cached_cholesky'],
            'errors': case['errors'],
        } for case in report['cases']],
    }, indent=2))


if __name__ == '__main__':
    main()
