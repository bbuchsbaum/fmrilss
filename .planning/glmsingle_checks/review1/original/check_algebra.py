"""Scratch numerical checks for the review; not a GLMsingle benchmark."""
import json
from time import perf_counter
import numpy as np
from scipy import linalg as la
from threadpoolctl import threadpool_limits

rng = np.random.default_rng(20261007)

def rel(a, b):
    return float(la.norm(a-b)/max(la.norm(b), 1e-300))

def local_design(hrf, times, nt):
    x = np.zeros((nt, len(times)))
    for j, t in enumerate(times):
        size = min(len(hrf), nt-t)
        x[t:t+size, j] = hrf[:size]
    return x

def fractions_to_lambdas(d, a, fractions):
    # a are orthonormal coefficient coordinates for the minimum-norm OLS fit.
    weights = a*a
    denominator = weights.sum(axis=0)
    result=[]
    for f in fractions:
        if f == 1:
            result.append(np.zeros(a.shape[1]))
            continue
        lo=np.zeros(a.shape[1])
        hi=np.full(a.shape[1], d.max()*(1/f-1))
        for _ in range(50):
            mid=(lo+hi)/2
            shrink=d[:,None]/(d[:,None]+mid)
            ratio=(weights*shrink**2).sum(axis=0)/denominator
            mask=ratio>f*f
            lo=np.where(mask,mid,lo)
            hi=np.where(mask,hi,mid)
        result.append((lo+hi)/2)
    return np.array(result)

with threadpool_limits(limits=1):
    nr, nt, nv, p = 6, 220, 128, 5
    hrf = np.load('tmp/glmsingle/library_dur3.0.npz')['hrfs'][:,8]
    xs, ys, qs, xr, yr = [], [], [], [], []
    for r in range(nr):
        times=np.arange(3,200,4)+rng.integers(0,2,size=50)
        x=local_design(hrf,times,nt)
        z=np.column_stack([np.ones(nt),np.linspace(-1,1,nt),rng.normal(size=(nt,p-2))])
        q=la.qr(z,mode='economic')[0]
        y=x@rng.normal(size=(len(times),nv))+q@rng.normal(size=(p,nv))+rng.normal(size=(nt,nv))
        xs.append(x);ys.append(y);qs.append(q)
        xr.append(x-q@(q.T@x));yr.append(y-q@(q.T@y))
    xglobal=la.block_diag(*xr)
    yglobal=np.concatenate(yr)
    u,s,vt=la.svd(xglobal,full_matrices=False,check_finite=False)
    ag=(u.T@yglobal)/s[:,None]
    local=[]
    for x,y in zip(xr,yr):
        ur,sr,vtr=la.svd(x,full_matrices=False,check_finite=False)
        local.append((sr,vtr,(ur.T@y)/sr[:,None]))
    ds=np.concatenate([sr*sr for sr,_,_ in local])
    als=np.concatenate([ar for _,_,ar in local])
    fs=np.linspace(.05,1,20)
    lambdas_global=fractions_to_lambdas(s*s,ag,fs)
    lambdas_local=fractions_to_lambdas(ds,als,fs)
    results={'dimensions':{'runs':nr,'T_per_run':nt,'N_per_run':50,'voxels':nv,'fractions':20,'nuisance_rank':p},
             'pooled_fraction_lambda_relative_error':rel(lambdas_local,lambdas_global)}
    maxerr=0
    for lam in lambdas_global:
        direct=vt.T@((s*s)[:,None]/((s*s)[:,None]+lam)*ag)
        split=np.concatenate([vtr.T@((sr*sr)[:,None]/((sr*sr)[:,None]+lam)*ar) for sr,vtr,ar in local])
        maxerr=max(maxerr,rel(split,direct))
    results['pooled_fraction_betas_max_relative_error']=maxerr

    # Timing isolates Gram SVD+data projection+beta reconstruction across an
    # already prepared set of 20 voxel-specific penalties. Both paths produce
    # all betas. Shared fraction root search, design construction, PC discovery,
    # CV, diagnostics, and file I/O are excluded. Not end-to-end GLMsingle.
    # Follow fracridge's tall-design path: factor X.T @ X, not a slower
    # direct tall-matrix SVD. Both implementations use the same strategy.
    def repeated_global():
        out=[]
        for lam in lambdas_global:
            _,dd,vv=la.svd(xglobal.T@xglobal,full_matrices=False,check_finite=False)
            cc=vv@(xglobal.T@yglobal)
            out.append(vv.T@(cc/(dd[:,None]+lam)))
        return out
    def cached_local():
        decomps=[]
        for x,y in zip(xr,yr):
            _,dd,vv=la.svd(x.T@x,full_matrices=False,check_finite=False)
            decomps.append((dd,vv,vv@(x.T@y)))
        return [np.concatenate([vv.T@(cc/(dd[:,None]+lam)) for dd,vv,cc in decomps]) for lam in lambdas_global]
    timed_dense=[];timed_local=[]
    for _ in range(3):
        t=perf_counter();baseline=repeated_global();timed_dense.append(perf_counter()-t)
        t=perf_counter();fast=cached_local();timed_local.append(perf_counter()-t)
    tbaseline=float(np.median(timed_dense));tfast=float(np.median(timed_local))
    results['kernel_timing_seconds']={'repeated_dense_global_gram_svd':tbaseline,'cached_local_gram_svd':tfast,'ratio':tbaseline/tfast,'repetitions':3,'statistic':'median'}
    results['timed_outputs_max_relative_error']=max(rel(a,b) for a,b in zip(baseline,fast))

    # Compact-support event Gram with exact nuisance low-rank downdate.
    x,y,q=xs[0],ys[0],qs[0]
    lam=.6
    g0=x.T@x
    c=x.T@q
    b=x.T@y-c@(q.T@y)
    g=g0-c@c.T
    bandwidth=max(abs(np.nonzero(np.abs(g0)>1e-14)[0]-np.nonzero(np.abs(g0)>1e-14)[1]))
    base=g0+lam*np.eye(g0.shape[0])
    ab=np.zeros((bandwidth+1,len(g0)))
    for k in range(bandwidth+1):
        ab[k,:len(g0)-k]=np.diag(base,k=-k)
    chol=la.cholesky_banded(ab,lower=True,check_finite=False)
    cb=la.cho_solve_banded((chol,True),c,check_finite=False)
    bb=la.cho_solve_banded((chol,True),b,check_finite=False)
    wb=bb+cb@la.solve(np.eye(p)-c.T@cb,c.T@bb,assume_a='pos')
    reference=la.solve(g+lam*np.eye(len(g)),b,assume_a='pos')
    results['banded_nuisance_relative_error']=rel(wb,reference)
    results['event_gram_half_bandwidth']=int(bandwidth)

    # Appending one orthonormal nuisance PC is a rank-one inverse downdate.
    qnext=rng.normal(size=nt);qnext-=q@(q.T@qnext);qnext/=la.norm(qnext)
    uu=x.T@qnext
    d=qnext@y
    aa=g+lam*np.eye(len(g))
    vv=la.solve(aa,uu,assume_a='pos')
    updated=reference+vv[:,None]*(uu@reference-d)[None,:]/(1-uu@vv)
    direct=la.solve(aa-np.outer(uu,uu),b-uu[:,None]*d[None,:],assume_a='pos')
    results['nested_nuisance_relative_error']=rel(updated,direct)

    # Analytic profile gradient from small HRF-basis sufficient statistics.
    hb=np.load('tmp/glmsingle/library_dur3.0.npz')['basis'][:,:3]
    xb=[local_design(hb[:,k],times,nt) for k in range(3)]
    xb=[a-q@(q.T@a) for a in xb]
    yy=yr[0][:,0]
    bv=np.stack([a.T@yy for a in xb])
    gg=np.array([[a.T@b for b in xb] for a in xb])
    theta=np.array([2.0,-.3,.15])
    def profile(theta):
        bb=theta@bv
        aa=np.einsum('i,j,ijab->ab',theta,theta,gg)+lam*np.eye(len(times))
        beta=la.solve(aa,bb,assume_a='pos')
        loss=yy@yy-bb@beta
        grad=[]
        for k in range(len(theta)):
            da=sum(theta[j]*(gg[k,j]+gg[j,k]) for j in range(len(theta)))
            grad.append(-2*bv[k]@beta+beta@da@beta)
        return loss,np.array(grad)
    loss,grad=profile(theta)
    step=1e-5
    fd=np.array([(profile(theta+np.eye(3)[k]*step)[0]-profile(theta-np.eye(3)[k]*step)[0])/(2*step) for k in range(3)])
    results['profile_gradient_relative_error']=rel(grad,fd)
    print(json.dumps(results,indent=2))
    with open('tmp/glmsingle/algebra_results.json','w') as f:
        json.dump(results,f,indent=2)
