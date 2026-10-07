#!/usr/bin/env python3
"""Execute the reviewer's original test body with real pinned dependencies.

Setup only is replaced: import the untouched fracridge 3.0 full module and
the actual calcbadness/zerodiv files at GLMsingle 1ab54a6. No dependency stubs
are registered, and none of the test body's numerical kernels are rewritten.
"""
import ast
import contextlib
import hashlib
import importlib.util
import io
import json
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import brentq
from threadpoolctl import threadpool_limits, threadpool_info

HERE=Path(__file__).resolve().parent
REVIEWER=HERE.parent/"due_diligence/reviewer_verify_core_identities.py"
if len(sys.argv)>1:REVIEWER=Path(sys.argv[1]).resolve()

spec=importlib.util.spec_from_file_location("actual_fracridge3",HERE/"fracridge.py")
fr=importlib.util.module_from_spec(spec);spec.loader.exec_module(fr)
sys.path.insert(0,str(HERE/"pinned"))
from glmsingle.ssq.calcbadness import calcbadness as cb
from glmsingle.utils.zerodiv import zerodiv


def rel(a,b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b),np.finfo(float).tiny))


def run():
    tree=ast.parse(REVIEWER.read_text())
    start=next(i for i,node in enumerate(tree.body) if isinstance(node,ast.Assign) and any(isinstance(t,ast.Name) and t.id=="rng" for t in node.targets))
    body=ast.Module(body=tree.body[start:],type_ignores=[])
    ns={"np":np,"brentq":brentq,"fr":fr,"cb":cb}
    captured=io.StringIO()
    with contextlib.redirect_stdout(captured):
        exec(compile(body,str(REVIEWER),"exec"),ns)
    output=captured.getvalue()
    (HERE/"reviewer_original_stdout.txt").write_text(output)
    coeff,alphas=ns["coef_g"],ns["alpha_g"]
    valid=ns["fracs"]<1
    ar=ns["alpha_exact"]
    metrics={"coefficient_relative_l2":rel(ns["coef_b"],coeff),
        "alpha_relative_l2":rel(ns["alpha_b"],alphas),
        "exactroot_median_abs_logratio":float(np.median(np.abs(np.log(alphas[valid]/ar[valid])))),
        "exactroot_max_abs_logratio":float(np.max(np.abs(np.log(alphas[valid]/ar[valid])))),
        "achieved_fraction_max_error":float(np.max(np.abs(ns["ach"][valid]-ns["fracs"][valid,None]))),
        "note":"Original reviewer uses direct SVD(X) per block and log1p(grid); upstream tall-design code uses SVD(X.T@X) and log(1+grid). This fixture does not test those boundary differences."}
    # Explicitly demonstrate the zero-SD behavior of the actual zerodiv helper.
    columns=[np.arange(4),np.arange(4,8)]
    conditions=[np.arange(4),np.arange(4)]
    results=[np.ones((1,8)),np.full((1,8),.5)]
    actual=cb([np.array(0),np.array(1)],columns,conditions,results,[1,1])
    sd=np.array([0.]);probe=zerodiv(np.zeros((1,8)),sd,val=0,wantcaution=0)
    zero_sd={"upstream_cv_scores":actual.tolist(),"naive_zero_all_candidates_cv_scores":[[0.,0.]],
             "sd_after_actual_zerodiv":sd.tolist(),
             "explanation":"zerodiv aliases tmp=y and replaces zero divisors by1 in-place. The baseline residual is0, but later candidates use divisor1, so they are not all set to0."}
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    out={"metadata":{"seed":1,"source_commit":"GLMsingle 1ab54a6","fracridge":"3.0", "jit":False,
        "python":sys.version,"numpy":np.__version__,"scipy":scipy.__version__,"threads":threadpool_info(),
        "wheel_sha256":sha(HERE/"fracridge-3.0-py3-none-any.whl"),
        "fracridge_source_sha256":sha(HERE/"fracridge.py"),
        "calcbadness_sha256":sha(HERE/"pinned/glmsingle/ssq/calcbadness.py"),
        "zerodiv_sha256":sha(HERE/"pinned/glmsingle/utils/zerodiv.py"),
        "reviewer_script_sha256":sha(REVIEWER),
        "setup_change":"Only dependency-loading setup omitted; original AST from rng initialization onwards executed unchanged; full actual modules and real helpers imported; no stubs"},
        "original_stdout":output,"metrics":metrics,"zero_sd_case":zero_sd}
    (HERE/"reviewer_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(output,end="")
    print(json.dumps({"zero_sd_case":zero_sd}))


if __name__=="__main__":
    with threadpool_limits(limits=1):run()
