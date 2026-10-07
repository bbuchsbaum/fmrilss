#!/usr/bin/env python3
"""Reproduce numerical behavior of the actual fracridge 3.0 source.

The downloaded wheel and its original module are adjacent to this script.
Only the original _do_svd and fracridge function ASTs plus their constants
are executed, to avoid importing the unrelated sklearn estimator classes.
Every primary comparison runs those original functions with jit=False.
No source arithmetic is edited for the main behavioral reference.
"""
from __future__ import annotations

import ast
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import sys
import warnings
from pathlib import Path

import numpy as np
import scipy
from scipy.linalg import block_diag, svd
from scipy.optimize import brentq
from scipy.stats import gamma
from threadpoolctl import threadpool_info, threadpool_limits

HERE = Path(__file__).resolve().parent
SOURCE = HERE / "fracridge.py"
FRACS = np.arange(0.05, 1.01, 0.05)
FRACS[-1] = 1.0
TOL = 1e-10


def load_original():
    tree = ast.parse(SOURCE.read_text())
    wanted = {"_do_svd", "fracridge"}
    constants = {"BIG_BIAS", "SMALL_BIAS", "BIAS_STEP"}
    body = []
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in wanted:
            body.append(node)
        elif isinstance(node, ast.Assign) and any(
            isinstance(t, ast.Name) and t.id in constants for t in node.targets
        ):
            body.append(node)
    ns = {"np": np, "interp": np.interp, "warnings": warnings, "__name__": "fracridge_precision_ast"}
    exec(compile(ast.Module(body=body, type_ignores=[]), str(SOURCE), "exec"), ns)
    return ns


NS = load_original()
ORIGINAL_SVD = NS["_do_svd"]
ORIGINAL_FRACRIDGE = NS["fracridge"]


def arrays_summary(a):
    a = np.asarray(a)
    return {"dtype": str(a.dtype), "shape": list(a.shape),
            "finite": int(np.isfinite(a).sum()), "size": int(a.size)}


def compare(a, b):
    a, b = np.asarray(a), np.asarray(b)
    good = np.isfinite(a) & np.isfinite(b)
    allgood = bool(good.all())
    out = {"all_finite": allgood,
           "nonfinite_a": int((~np.isfinite(a)).sum()),
           "nonfinite_b": int((~np.isfinite(b)).sum())}
    if allgood:
        delta = a.astype(float) - b.astype(float)
        denom = np.linalg.norm(b.astype(float))
        out.update(relative_l2=float(np.linalg.norm(delta) / max(denom, np.finfo(float).tiny)),
                   max_absolute=float(np.abs(delta).max(initial=0)))
    else:
        out.update(relative_l2=None, max_absolute=None)
    return out


def grid(eig):
    # Intentionally mirror 3.0, including endpoints, raw (unmasked) extrema,
    # np.arange rounding, and log(1+alpha) rather than log1p(alpha).
    val1 = NS["BIG_BIAS"] * eig[0] ** 2
    val2 = NS["SMALL_BIAS"] * eig[-1] ** 2
    val2 = NS["SMALL_BIAS"] if val2 == 0 else val2
    return np.concatenate([np.array([0]), 10 ** np.arange(
        np.floor(np.log10(val2)), np.ceil(np.log10(val1)), NS["BIAS_STEP"])])


def call_actual(X, y, fracs=FRACS, cache=None):
    old = NS["_do_svd"]
    if cache is not None:
        # This injects blockwise sufficient statistics into the UNMODIFIED
        # original grid/interpolation/reconstruction function.
        NS["_do_svd"] = lambda X, y, jit=True: tuple(a.copy() for a in cache)
    try:
        with warnings.catch_warnings(record=True) as ws:
            warnings.simplefilter("always")
            coef, alpha = ORIGINAL_FRACRIDGE(X, y, fracs=fracs, tol=TOL, jit=False)
        return coef, alpha, sorted(set(str(w.message) for w in ws))
    finally:
        NS["_do_svd"] = old


def call_cache_forced_grid(X,y,cache,reference_grid,fracs=FRACS):
    """Control experiment: original interpolation AST with only grid injected.

    This is NOT the unmodified behavioral reference. It isolates interpolation
    grid differences from differences in the computed block spectra.
    """
    ns=load_original()
    tree=ast.parse(SOURCE.read_text())
    fn=next(x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name=="fracridge")
    class InjectGrid(ast.NodeTransformer):
        def visit_Assign(self,node):
            if any(isinstance(t,ast.Name) and t.id=="alphagrid" for t in node.targets):
                node.value=ast.Name(id="INJECTED_GRID",ctx=ast.Load())
            return self.generic_visit(node)
    fn=InjectGrid().visit(fn)
    ns["INJECTED_GRID"]=reference_grid.copy()
    ns["_do_svd"]=lambda X,y,jit=True:tuple(a.copy() for a in cache)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),str(SOURCE),"exec"),ns)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ns["fracridge"](X,y,fracs=fracs,tol=TOL,jit=False)


def cache_actual(X, y):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return ORIGINAL_SVD(X, y, jit=False)


def split_cache(blocks, targets, preserve_global_branch=True, use_global_gram=False):
    X = block_diag(*blocks)
    Y = np.vstack(targets)
    n, p = X.shape
    tall = n > p
    if not tall and sum(min(x.shape) for x in blocks) != min(n, p):
        raise ValueError("mixed-aspect direct-SVD blocks require nullspace padding")
    gram = X.T @ X if use_global_gram and tall else None
    xy = X.T @ Y if use_global_gram and tall else None
    eigs, vs, cs = [], [], []
    col = 0
    with np.errstate(divide="ignore", invalid="ignore"):
        for x, y in zip(blocks, targets):
            if tall:
                g = gram[col:col+x.shape[1], col:col+x.shape[1]] if gram is not None else x.T @ x
                _, ss, vt = svd(g, full_matrices=False)
                eig = np.sqrt(ss)
                n_for_branch = n if preserve_global_branch else x.shape[0]
                if y.shape[-1] >= n_for_branch:
                    yn = np.diag(1.0/eig) @ vt @ x.T @ y
                else:
                    product = xy[col:col+x.shape[1]] if xy is not None else x.T @ y
                    yn = np.diag(1.0/eig) @ vt @ product
            else:
                u, eig, vt = svd(x, full_matrices=False)
                yn = u.T @ y
            coef = (yn.T / eig).T
            eigs.append(eig); vs.append(vt); cs.append(coef)
            col += x.shape[1]
    eig = np.concatenate(eigs)
    vt = block_diag(*vs)
    coef = np.vstack(cs)
    order = np.argsort(-eig, kind="stable")
    return eig[order], vt[order], coef[order]


def exact_fraction_from_same_cache(cache, fracs=FRACS):
    eig, vt, coef = (a.copy() for a in cache)
    coef[eig < TOL] = 0
    d = eig.astype(float) ** 2
    coef = coef.astype(float)
    alpha = np.zeros((len(fracs), coef.shape[1]))
    transformed = np.zeros((len(eig), len(fracs), coef.shape[1]))
    for j in range(coef.shape[1]):
        w = coef[:, j] ** 2
        if not np.isfinite(w).all() or w.sum() == 0 or np.any(d == 0):
            alpha[:, j] = np.nan; transformed[:, :, j] = np.nan
            continue
        for k, f in enumerate(fracs):
            if f == 1:
                a = 0.0
            else:
                fun = lambda a: np.sum((d/(d+a))**2 * w)/w.sum() - f*f
                hi = float(max(d.max(), np.finfo(float).tiny))
                while fun(hi) > 0:
                    hi *= 4
                a = brentq(fun, 0.0, hi, xtol=np.nextafter(0.0, 1.0), rtol=2e-14)
            alpha[k, j] = a
            transformed[:, k, j] = d/(d+a) * coef[:, j]
    beta = (vt.T.astype(float) @ transformed.reshape(len(eig), -1)).reshape(vt.shape[1], len(fracs), coef.shape[1])
    return beta, alpha


def spectrum_blocks(seed, dtype, condition=10, scale=1.0, target_count=13):
    rng = np.random.default_rng(seed)
    blocks, targets = [], []
    for r, (n, p) in enumerate([(82, 17), (91, 19), (76, 16), (84, 18)]):
        u, _ = np.linalg.qr(rng.normal(size=(n, p)))
        v, _ = np.linalg.qr(rng.normal(size=(p, p)))
        vals = scale * [1.0, .3, 2.0, .08][r] * np.geomspace(1, 1/condition, p)
        x = ((u * vals) @ v.T).astype(dtype)
        y = (x.astype(float) @ rng.normal(size=(p, target_count)) + .2*rng.normal(size=(n,target_count))).astype(dtype)
        blocks.append(x); targets.append(y)
    return blocks, targets


def event_blocks(seed, dtype, isi=3.2):
    rng = np.random.default_rng(seed)
    blocks, targets = [], []
    for r, T in enumerate([151, 164, 159, 156]):
        time = np.arange(T, dtype=float)
        ons = np.arange(7., T-32, isi) + rng.uniform(-.15, .15, len(np.arange(7., T-32, isi)))
        tt = time[:,None] - ons[None,:]
        xx = gamma.pdf(tt, 6+.2*r) - gamma.pdf(tt, 16+.3*r)/6
        xx /= gamma.pdf(5+.2*r,6+.2*r)
        raw_y = xx @ rng.normal(size=(len(ons),13)) + .25*rng.normal(size=(T,13))
        nuisance = np.column_stack([np.linspace(-1,1,T)**k for k in range(4)])
        nuisance = np.column_stack([nuisance, rng.normal(size=(T,3))])
        q, _ = np.linalg.qr(nuisance)
        P = (np.eye(T)-q@q.T).astype(dtype)
        blocks.append(P @ xx.astype(dtype))
        targets.append(P @ raw_y.astype(dtype))
    return blocks, targets


def describe_case(name, blocks, targets, root=True):
    X = block_diag(*blocks); Y = np.vstack(targets)
    out = {"name":name, "dtype":str(X.dtype), "shape_X":list(X.shape), "targets":Y.shape[1],
           "condition_X_float64":float(np.linalg.cond(X.astype(float)))}
    b, a, warn = call_actual(X,Y)
    cache = cache_actual(X,Y)
    out["actual_beta"] = arrays_summary(b)
    out["actual_warnings"] = warn
    out["computed_singular_min"] = float(cache[0][-1])
    out["computed_singular_max"] = float(cache[0][0])
    out["masked_modes"] = int(np.sum(cache[0]<TOL))
    out["actual_grid_points"] = int(grid(cache[0]).size)
    out["cache_replay_beta"] = compare(call_actual(X,Y,cache=cache)[0],b)
    try:
        sc = split_cache(blocks,targets)
        bs,as_,ws = call_actual(X,Y,cache=sc)
        out["split_same_global_interp_beta"] = compare(bs,b)
        out["split_same_global_interp_alpha"] = compare(as_,a)
        out["split_grid_points"] = int(grid(sc[0]).size)
        out["split_grid_matches_bitwise"] = bool(np.array_equal(grid(sc[0]),grid(cache[0])))
        out["split_warnings"] = ws
        bf,af=call_cache_forced_grid(X,Y,sc,grid(cache[0]))
        out["split_forced_identical_grid_beta"]=compare(bf,b)
        out["split_forced_identical_grid_alpha"]=compare(af,a)
        if X.shape[0]>X.shape[1]:
            sg=split_cache(blocks,targets,use_global_gram=True)
            bg,ag,_=call_actual(X,Y,cache=sg)
            out["split_from_global_gram_beta"]=compare(bg,b)
            out["split_from_global_gram_alpha"]=compare(ag,a)
    except Exception as exc:
        out["split_exception"] = str(exc)
    if X.dtype == np.float32:
        b64,a64,w64=call_actual(X.astype(float),Y.astype(float))
        out["float64_replay_beta_vs_python32"] = compare(b64,b)
        out["float64_replay_alpha_vs_python32"] = compare(a64,a)
        out["float64_replay_per_fraction"]=[{
            "fraction":float(f),"beta_relative_l2":compare(b64[:,j,:],b[:,j,:])["relative_l2"],
            "reference_beta_norm":float(np.linalg.norm(b[:,j,:]))}
            for j,f in enumerate(FRACS)]
        out["float64_replay_excluding_OLS"]=compare(b64[:,:-1,:],b[:,:-1,:])
        out["float64_replay_up_to_fraction_08"]=compare(b64[:,FRACS<=.8,:],b[:,FRACS<=.8,:])
    if root and np.isfinite(b).all() and np.isfinite(a).all():
        br,ar=exact_fraction_from_same_cache(cache)
        valid=(ar>0)&np.isfinite(ar)&np.isfinite(a)
        rel=np.abs(a[valid]-ar[valid])/ar[valid]
        achieved=np.linalg.norm(b,axis=0)/np.linalg.norm(b[:,-1,:],axis=0)[None,:]
        out["grid_vs_exact_root_alpha_relative"]={"median":float(np.median(rel)),"max":float(np.max(rel))}
        out["grid_fraction_max_absolute_error"]=float(np.max(np.abs(achieved-FRACS[:,None])))
        out["grid_vs_exact_root_beta"]=compare(b,br)
    return out,(b,a,cache,X,Y)


def edge_cases():
    rng=np.random.default_rng(170)
    result=[]
    for dtype in [np.float64,np.float32]:
        for kind in ["exact_zero_column","near_tol_diagonal","duplicated_column","zero_target","tiny_target","tiny_absolute_design"]:
            X=rng.normal(size=(36,8)).astype(dtype)
            Y=rng.normal(size=(36,4)).astype(dtype)
            if kind=="exact_zero_column": X[:,-1]=0
            if kind=="near_tol_diagonal":
                X=np.zeros((36,8),dtype=dtype);X[:8]=np.diag(np.array([1,.1,.01,.001,1e-6,3e-9,2e-10,5e-11],dtype=dtype))
            if kind=="duplicated_column": X[:,-1]=X[:,0]
            if kind=="zero_target":Y[:,0]=0
            if kind=="tiny_target":Y[:,0]*=1e-25
            if kind=="tiny_absolute_design":
                X=np.zeros((36,8),dtype=dtype);X[:8]=np.diag(np.geomspace(1e-8,5e-10,8).astype(dtype))
            row,_=describe_case(kind+"_"+np.dtype(dtype).name,[X],[Y],root=True)
            result.append(row)
    return result


def near_tie_fixture(blocks,targets):
    X=block_diag(*blocks);Y=np.vstack(targets)
    fractions=np.array([.45,.5])
    b32,a32,_=call_actual(X.astype(np.float32),Y.astype(np.float32),fractions)
    b64,a64,_=call_actual(X.astype(np.float32).astype(float),Y.astype(np.float32).astype(float),fractions)
    best=None
    for j in range(Y.shape[1]):
        x0,x1=b64[:,0,j],b64[:,1,j]
        z0,z1=b32[:,0,j],b32[:,1,j]
        mid=(x0+x1)/2; diff=x0-x1
        def delta(u,v,m):return float(np.sum((u-m)**2)-np.sum((v-m)**2))
        delta32=delta(z0,z1,mid)
        if np.dot(diff,diff)==0 or delta32==0:continue
        m=mid+delta32/(4*np.dot(diff,diff))*diff
        l64=np.sum((b64[:,:,j]-m[:,None])**2,axis=0)
        l32=np.sum((b32[:,:,j]-m[:,None])**2,axis=0)
        row={"constructed_not_prevalence_estimate":True,"target":j,"fractions":fractions.tolist(),
             "loss64":l64.tolist(),"loss32":l32.tolist(),"choice64":float(fractions[np.argmin(l64)]),
             "choice32":float(fractions[np.argmin(l32)]),"margin64":float(abs(np.diff(l64)[0])),
             "margin32":float(abs(np.diff(l32)[0])),"relative_beta_change":compare(b64[:,:,j],b32[:,:,j])["relative_l2"]}
        if row["choice64"]!=row["choice32"]:
            best=row
            np.savez_compressed(HERE/"constructed_near_tie.npz",X=X.astype(np.float32),Y=Y.astype(np.float32),target=m,fractions=fractions,beta32=b32,beta64=b64,voxel=j)
            break
    return best


def main():
    result={"metadata":{
        "wheel_sha256":hashlib.sha256((HERE/"fracridge-3.0-py3-none-any.whl").read_bytes()).hexdigest(),
        "source_sha256":hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        "python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"scipy":scipy.__version__,
        "scikit_learn":importlib.metadata.version("scikit-learn"),"numba":"not installed; explicit jit=False",
        "reference":"unmodified AST of fracridge 3.0 _do_svd and fracridge; jit=False",
        "blas_threads":threadpool_info(),"tol":TOL,"fracs":FRACS.tolist(),
        "seeds":{"spectral":117,"events":18,"edge_case_sequence":170,"separate_run_fractions":91},
        "relative_error_definition":"||candidate-reference||_F / ||reference||_F over all finite trials x fractions x targets; reference is original global fracridge at stated dtype",
        "synthetic_stress_fixtures_only":True},"cases":[]}
    for dtype in [np.float64,np.float32]:
        for condition in [10,100,1000,100000]:
            blocks,targets=spectrum_blocks(117,dtype,condition=condition)
            row,_=describe_case("spectrum_condition_"+str(condition)+"_"+np.dtype(dtype).name,blocks,targets)
            row["seed"]=117
            result["cases"].append(row)
        for isi in [3.2,1.2]:
            blocks,targets=event_blocks(18,dtype,isi=isi)
            row,_=describe_case("event_isi_"+str(isi)+"_"+np.dtype(dtype).name,blocks,targets)
            row["seed"]=18
            result["cases"].append(row)
    result["edge_cases"]=edge_cases()
    blocks,targets=spectrum_blocks(117,np.float32,condition=100,target_count=350)
    X=block_diag(*blocks);Y=np.vstack(targets)
    b,a,_=call_actual(X,Y)
    bsmall,asmall,_=call_actual(X,Y[:,:13])
    result["batch_association"]={"total_timepoints":X.shape[0],"large_targets":Y.shape[1],"small_targets":13,
        "first13_beta_difference":compare(b[:,:,:13],bsmall),"first13_alpha_difference":compare(a[:,:13],asmall)}
    targets_middle=[y[:,:100] for y in targets]
    ym=np.vstack(targets_middle)
    bm,am,_=call_actual(X,ym)
    wrongcache=split_cache(blocks,targets_middle,preserve_global_branch=False)
    bw,aw,_=call_actual(X,ym,cache=wrongcache)
    correctcache=split_cache(blocks,targets_middle,preserve_global_branch=True)
    bc,ac,_=call_actual(X,ym,cache=correctcache)
    result["global_vs_local_projection_branch"]={"targets":100,
        "correct_global_branch_beta":compare(bc,bm),"incorrect_local_branch_beta":compare(bw,bm),
        "correct_global_branch_alpha":compare(ac,am),"incorrect_local_branch_alpha":compare(aw,am)}
    blocks,targets=spectrum_blocks(91,np.float64,condition=10)
    X=block_diag(*blocks);Y=np.vstack(targets)
    bg,ag,_=call_actual(X,Y)
    local=[call_actual(x,y) for x,y in zip(blocks,targets)]
    result["wrong_separate_run_fractions"]={"beta_difference":compare(np.concatenate([x[0] for x in local]),bg),
        "note":"This deliberately imposes f separately per run; it is not the same pooled estimator."}
    blocks,targets=spectrum_blocks(117,np.float32,condition=100)
    result["constructed_near_tie"]=near_tie_fixture(blocks,targets)
    # Confirm AST extraction against import of the untouched complete module,
    # with the actual installed sklearn dependency (no dependency stubs).
    spec=importlib.util.spec_from_file_location("actual_fracridge3_full_module",SOURCE)
    full=importlib.util.module_from_spec(spec);spec.loader.exec_module(full)
    X=block_diag(*blocks);Y=np.vstack(targets)
    ab,aa,_=call_actual(X,Y)
    fb,fa=full.fracridge(X,Y,FRACS,tol=TOL,jit=False)
    result["unmodified_full_module_vs_AST"]={"beta":compare(ab,fb),"alpha":compare(aa,fa)}
    def safe_json(x):
        if isinstance(x,dict): return {k:safe_json(v) for k,v in x.items()}
        if isinstance(x,list): return [safe_json(v) for v in x]
        if isinstance(x,float) and not np.isfinite(x): return str(x)
        return x
    result=safe_json(result)
    (HERE/"precision_results.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    for row in result["cases"]:
        print(json.dumps({"name":row["name"],"cond":row["condition_X_float64"],
            "split_beta_rel":row.get("split_same_global_interp_beta",{}).get("relative_l2"),
            "split_alpha_rel":row.get("split_same_global_interp_alpha",{}).get("relative_l2"),
            "replay64_beta_rel":row.get("float64_replay_beta_vs_python32",{}).get("relative_l2"),
            "root_alpha":row.get("grid_vs_exact_root_alpha_relative"),
            "fraction_error":row.get("grid_fraction_max_absolute_error")}))
    for key in ["batch_association","global_vs_local_projection_branch","wrong_separate_run_fractions","constructed_near_tie"]:
        print(json.dumps({key:result[key]}))
    for row in result["edge_cases"]:
        print(json.dumps({"edge":row["name"],"actual_beta":row["actual_beta"],"masked_modes":row["masked_modes"],"singular_min":row["computed_singular_min"],"root_alpha":row.get("grid_vs_exact_root_alpha_relative"),"fraction_error":row.get("grid_fraction_max_absolute_error"),"warnings":row["actual_warnings"]}))


if __name__=="__main__":
    with threadpool_limits(limits=1):
        main()
