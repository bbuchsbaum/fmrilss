#!/usr/bin/env python3
"""Reproduce five narrow GLMsingle source/numerical edge observations.

Run next to the unmodified downloaded source files::

    python source_edge_probe.py --output source_edge_results.json

Alternatively, pass --source-dir pointing to the flat source directory. Required
files and their pinned Git blob hashes are listed below. This script performs no
network requests and does not modify the reference kernels. It needs NumPy,
SciPy, and scikit-learn; it does not require installing GLMsingle or fracridge.

These are focused source-edge checks, not an end-to-end estimator benchmark.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import sys
import types

import numpy as np
import scipy
from scipy.linalg import qr
import sklearn


COMMIT = "1ab54a65edd3ea41a6133d4b4ecb78a9c7296684"
REFERENCE_FILES = {
    "olsmatrix.py": (
        "glmsingle/ols/olsmatrix.py",
        "cb22ce1633796b3cb4c731e77cd9124b6ec898cf",
    ),
    "olsmatrix2.py": (
        "glmsingle/ols/olsmatrix2.py",
        "8103b7a71ef7152fcb1e335685662591668d0d6e",
    ),
    "make_poly_matrix.py": (
        "glmsingle/ols/make_poly_matrix.py",
        "5ceeafad5c2da87761c536bf8a66867b07ff75c1",
    ),
    "gethrf.py": (
        "glmsingle/hrf/gethrf.py",
        "4cc08dd63df2f97e212e48515d57cb7321e13f68",
    ),
    "alt_round.py": (
        "glmsingle/utils/alt_round.py",
        "6320592b8b7ea5924f685435d16060cd2ece0c1e",
    ),
    "getcanonicalhrflibrary.tsv": (
        "glmsingle/hrf/getcanonicalhrflibrary.tsv",
        "51247d135662c95d73c7373b5ec97c21e155e1bd",
    ),
}


def source_provenance(source_dir: Path) -> list[dict]:
    records = []
    for filename, (repository_path, expected_blob) in REFERENCE_FILES.items():
        contents = (source_dir / filename).read_bytes()
        blob_header = f"blob {len(contents)}\0".encode("ascii")
        actual_blob = hashlib.sha1(blob_header + contents).hexdigest()
        if actual_blob != expected_blob:
            raise ValueError(f"Reference source hash mismatch: {filename}")
        records.append(
            {
                "local_filename": filename,
                "repository_path": repository_path,
                "pinned_url": (
                    f"https://github.com/cvnlab/GLMsingle/blob/{COMMIT}/"
                    f"{repository_path}"
                ),
                "git_blob_sha1": actual_blob,
                "sha256": hashlib.sha256(contents).hexdigest(),
                "bytes": len(contents),
                "matches_pinned_blob": True,
            }
        )
    return records


def load_reference_kernels(source_dir: Path) -> dict:
    # Provide only import namespaces, then execute the original files unchanged.
    # This avoids running an installed GLMsingle package's __init__ or substituting
    # a reimplementation of any reference kernel used by these checks.
    for name in ("glmsingle", "glmsingle.ols", "glmsingle.utils", "glmsingle.hrf"):
        package = types.ModuleType(name)
        package.__path__ = []
        sys.modules[name] = package

    def load(name: str, filename: str):
        spec = importlib.util.spec_from_file_location(name, source_dir / filename)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"Cannot load {filename}")
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module

    ols = load("glmsingle.ols.olsmatrix", "olsmatrix.py")
    ols2 = load("glmsingle.ols.olsmatrix2", "olsmatrix2.py")
    poly = load("glmsingle.ols.make_poly_matrix", "make_poly_matrix.py")
    load("glmsingle.utils.alt_round", "alt_round.py")
    hrf = load("glmsingle.hrf.gethrf", "gethrf.py")
    return {"ols": ols, "ols2": ols2, "poly": poly, "hrf": hrf}


def relative_norm(error, reference) -> float:
    denominator = max(float(np.linalg.norm(reference)), np.finfo(float).tiny)
    return float(np.linalg.norm(error) / denominator)


def probe_nuisance_rank(kernels: dict) -> dict:
    n = 40
    time = np.linspace(-1, 1, n)
    basis, _ = np.linalg.qr(np.column_stack([np.ones(n), time, time**2]))
    delta = 1e-9
    nuisance = np.column_stack(
        [basis[:, 0], basis[:, 1], basis[:, 1] + delta * basis[:, 2]]
    )
    reference_residual_projector = kernels["poly"].make_projection_matrix(nuisance)
    qr_basis, triangular, pivots = qr(nuisance, mode="economic", pivoting=True)
    qr_cutoff = np.finfo(float).eps * max(nuisance.shape) * abs(triangular[0, 0])
    qr_rank = int(np.sum(abs(np.diag(triangular)) > qr_cutoff))
    qr_projector = np.eye(n) - qr_basis[:, :qr_rank] @ qr_basis[:, :qr_rank].T
    reference_removed_trace = float(np.trace(np.eye(n) - reference_residual_projector))
    return {
        "construction": "Z = [q0, q1, q1 + delta*q2], orthonormal q0:q2",
        "n_timepoints": n,
        "delta": delta,
        "dtype": str(nuisance.dtype),
        "reference_function": "make_projection_matrix -> olsmatrix(mode=0)",
        "reference_policy": (
            "Drop exactly zero columns; unit-normalize; np.linalg.pinv of the "
            "normalized Gram, with no explicit rcond passed by the reference"
        ),
        "reference_removed_subspace_trace": reference_removed_trace,
        "reference_remaining_q2_norm": float(
            np.linalg.norm(reference_residual_projector @ basis[:, 2])
        ),
        "comparison_qr_policy": "eps * max(Z.shape) * abs(R[0,0])",
        "comparison_qr_cutoff": float(qr_cutoff),
        "comparison_qr_rank": qr_rank,
        "comparison_qr_column_pivots": pivots.tolist(),
        "projector_frobenius_difference": float(
            np.linalg.norm(reference_residual_projector - qr_projector)
        ),
        "observed_rank_policy_difference": bool(
            qr_rank == 3 and abs(reference_removed_trace - 2) < 1e-10
        ),
    }


def probe_collinear_task(kernels: dict) -> dict:
    column = np.arange(1.0, 5.0)
    design = np.column_stack([column, column])
    result = {
        "construction": "Two identical nonzero task columns [1,2,3,4]",
        "design": design.tolist(),
        "dtype": str(design.dtype),
        "reference_function": "olsmatrix2(X, lambda_=0)",
        "reference_zero_column_rule": "np.all(X == 0, axis=0)",
    }
    try:
        operator = kernels["ols2"].olsmatrix2(design)
        result.update({"raised": False, "returned_operator": operator.tolist()})
    except np.linalg.LinAlgError as error:
        result.update(
            {"raised": True, "exception_type": type(error).__name__, "message": str(error)}
        )
    result["observed_singular_solve_exception"] = bool(result["raised"])
    return result


def probe_custom_nuisance_score(kernels: dict) -> dict:
    seed = 101
    rng = np.random.default_rng(seed)
    n_time, n_trial, n_voxel, n_extra = 130, 31, 20, 7
    raw_design = rng.normal(size=(n_time, n_trial))
    raw_data = rng.normal(size=(n_time, n_voxel))
    extra = rng.normal(size=(n_time, n_extra))
    polynomial = kernels["poly"].make_polynomial_matrix(n_time, np.arange(3))
    poly_projector = kernels["poly"].make_projection_matrix(polynomial)
    full_projector = kernels["poly"].make_projection_matrix(
        np.column_stack([polynomial, extra])
    )
    # Literal unmodified reference projection and OLS kernels, evaluated in
    # float64 to isolate the algebra from the top-level fitter's float32 casts.
    reference_beta = kernels["ols2"].olsmatrix2(full_projector @ raw_design) @ (
        full_projector @ raw_data
    )
    design = poly_projector @ raw_design
    data = poly_projector @ raw_data
    nuisance_basis, _ = np.linalg.qr(poly_projector @ extra, mode="reduced")
    g = design.T @ nuisance_basis
    d = nuisance_basis.T @ data
    gram = design.T @ design - g @ g.T
    rhs = design.T @ data - g @ d
    lower = np.linalg.cholesky(gram)
    z = np.linalg.solve(lower, rhs)
    a = np.linalg.solve(lower, g).T @ z
    statistic_sse = (
        np.sum(data**2, axis=0)
        - np.sum(z**2, axis=0)
        + np.sum(a**2, axis=0)
        - 2 * np.sum(a * d, axis=0)
    )
    explicit_reference_sse = np.sum((data - design @ reference_beta) ** 2, axis=0)
    naive_sse = np.sum(data**2, axis=0) - np.sum(z**2, axis=0)
    statistic_beta = np.linalg.solve(lower.T, z)
    score_error = relative_norm(statistic_sse - explicit_reference_sse, explicit_reference_sse)
    return {
        "seed": seed,
        "shape": {
            "timepoints": n_time,
            "trials": n_trial,
            "voxels": n_voxel,
            "polynomials": 3,
            "custom_nuisance_columns": n_extra,
        },
        "dtype": str(raw_design.dtype),
        "reference_kernels": ["make_polynomial_matrix", "make_projection_matrix", "olsmatrix2"],
        "scope": (
            "Full-rank OLS algebra using unmodified reference kernels in float64; "
            "not a test of complete GLMsingle float32 execution or autoscaled R2"
        ),
        "statistic_sse_formula": "||e||^2 - ||z||^2 + ||a||^2 - 2*a.T*d",
        "definitions": (
            "A=P_poly*X, e=P_poly*y, U=orth(P_poly*extra), g=A.T*U, d=U.T*e; "
            "G=A.T*A-g*g.T=L*L.T; b=A.T*e-g*d; z=solve(L,b); a=solve(L,g).T*z"
        ),
        "sse_relative_l2_error": score_error,
        "sse_max_absolute_error": float(np.max(abs(statistic_sse - explicit_reference_sse))),
        "beta_relative_frobenius_error": relative_norm(
            statistic_beta - reference_beta, reference_beta
        ),
        "naive_uncorrected_score_relative_l2_error": relative_norm(
            naive_sse - explicit_reference_sse, explicit_reference_sse
        ),
        "corrected_score_agrees": bool(score_error < 1e-10),
    }


def probe_hrf_support(kernels: dict) -> dict:
    rows = []
    for duration in (0.1, 1.0, 3.0):
        for tr in (2.0, 1.333, 1.0, 0.8):
            # Execute the actual library loader, boxcar convolution, PCHIP, and
            # normalization code. Its fpath points to the verified sibling TSV.
            library = kernels["hrf"].getcanonicalhrflibrary(duration, tr)
            nonzero = np.count_nonzero(library, axis=1)
            rows.append(
                {
                    "duration_seconds": duration,
                    "tr_seconds": tr,
                    "n_hrfs": int(library.shape[0]),
                    "n_kernel_samples": int(library.shape[1]),
                    "nonzero_min": int(nonzero.min()),
                    "nonzero_max": int(nonzero.max()),
                    "nonzero_per_hrf": nonzero.tolist(),
                    "last_sample_seconds": float((library.shape[1] - 1) * tr),
                    "dtype": str(library.dtype),
                }
            )
    return {
        "reference_function": "getcanonicalhrflibrary(duration, tr)",
        "notes": (
            "Unmodified actual sampled kernels; exactly nonzero entries counted. "
            "Trial columns near run boundaries may be shorter. No tail threshold "
            "or additional kernel truncation is applied. Later per-HRF peak "
            "normalization does not change these support counts."
        ),
        "rows": rows,
    }


def probe_autoscale_constant(kernels: dict) -> dict:
    candidate = np.full(4, 2, dtype=np.float32)
    target = np.asarray([1, 3, 5, 7], dtype=np.float32)
    design = np.column_stack([candidate, np.ones(candidate.size)]).astype(np.float32)
    reference_parameters = kernels["ols"].olsmatrix(design) @ target
    reference_prediction = design @ reference_parameters
    normal_solve = {"raised": False}
    try:
        normal_solve["parameters"] = np.linalg.solve(design.T @ design, design.T @ target).tolist()
    except np.linalg.LinAlgError as error:
        normal_solve = {
            "raised": True,
            "exception_type": type(error).__name__,
            "message": str(error),
        }
    return {
        "construction": "Constant candidate beta=2 plus intercept; OLS target=[1,3,5,7]",
        "reference_function": "olsmatrix(X) @ target, with X cast to float32 as in autoscale",
        "reference_scale_offset": reference_parameters.tolist(),
        "reference_prediction": reference_prediction.tolist(),
        "reference_scale_is_negative": bool(reference_parameters[0] < 0),
        "simple_2x2_normal_solve": normal_solve,
        "raw_design_pseudoinverse_scale_offset": (np.linalg.pinv(design) @ target).tolist(),
        "observation": (
            "A simple 2x2 solve is singular. The reference's normalized-Gram "
            "pseudoinverse has defined scale/offset; a raw-design pseudoinverse "
            "can give different parameters despite the same fitted constant."
        ),
        "observed_defined_reference_and_singular_normal_solve": bool(
            np.isfinite(reference_parameters).all()
            and normal_solve["raised"]
            and np.allclose(reference_prediction, target.mean(), rtol=1e-6, atol=1e-6)
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("source_edge_results.json"))
    args = parser.parse_args()
    source_dir = args.source_dir.resolve()
    sources = source_provenance(source_dir)
    kernels = load_reference_kernels(source_dir)
    results = {
        "schema_version": 1,
        "scope": "Focused source-edge observations; no end-to-end benchmark or accuracy claim",
        "provenance": {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "reference_repository": "https://github.com/cvnlab/GLMsingle",
            "reference_commit": COMMIT,
            "probe_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scikit_learn": sklearn.__version__,
            "platform": platform.platform(),
            "reference_sources": sources,
            "source_execution": "Unmodified verified files loaded using importlib; package namespaces only are supplied",
        },
        "nuisance_rank_policy": probe_nuisance_rank(kernels),
        "collinear_task_solve": probe_collinear_task(kernels),
        "custom_nuisance_score": probe_custom_nuisance_score(kernels),
        "actual_hrf_support": probe_hrf_support(kernels),
        "constant_candidate_autoscale": probe_autoscale_constant(kernels),
    }
    checks = {
        "nuisance_rank_policy_difference": results["nuisance_rank_policy"]["observed_rank_policy_difference"],
        "collinear_task_raises": results["collinear_task_solve"]["observed_singular_solve_exception"],
        "corrected_score_agrees": results["custom_nuisance_score"]["corrected_score_agrees"],
        "actual_library_has_twenty_hrfs": all(
            row["n_hrfs"] == 20 for row in results["actual_hrf_support"]["rows"]
        ),
        "constant_autoscale_corner_case": results["constant_candidate_autoscale"][
            "observed_defined_reference_and_singular_normal_solve"
        ],
    }
    results["checks"] = checks
    results["all_checks_passed"] = all(checks.values())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(args.output.resolve()), "checks": checks}, indent=2))
    return 0 if results["all_checks_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
