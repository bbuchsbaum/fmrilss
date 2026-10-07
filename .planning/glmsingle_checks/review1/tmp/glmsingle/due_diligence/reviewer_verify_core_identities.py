# Numerical checks for GLMSINGLE_FAST_PLAN.md section 1.
# Usage: python3 -I verify_core_identities.py <EXT>
#   <EXT>/GLMsingle  = git clone of cvnlab/GLMsingle at 1ab54a6
#   <EXT>/fr/x       = unzipped fracridge-3.0 wheel
# Needs numpy, scipy, scikit-learn.
import sys, importlib.util, numpy as np
from scipy.optimize import brentq
EXT = sys.argv[1]
def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path); m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m
fr = load("fr", f"{EXT}/fr/x/fracridge/fracridge.py")
# calcbadness needs zerodiv; inject stub module tree
zd = load("zerodiv", f"{EXT}/GLMsingle/glmsingle/utils/zerodiv.py")
import types
pkg = types.ModuleType("glmsingle"); ut = types.ModuleType("glmsingle.utils"); ut.zerodiv = zd
sys.modules.update({"glmsingle": pkg, "glmsingle.utils": ut, "glmsingle.utils.zerodiv": zd})
cb = load("cb", f"{EXT}/GLMsingle/glmsingle/ssq/calcbadness.py").calcbadness

rng = np.random.default_rng(1)
R, n, T, V = 5, 18, 90, 40
fracs = np.linspace(0.05, 1, 20)
t = np.arange(25); hrf = t**5*np.exp(-t)/120; hrf -= 0.1*t**8*np.exp(-t)/40320; hrf/=hrf.max()
def proj(A):
    Q,_ = np.linalg.qr(A); return np.eye(A.shape[0]) - Q@Q.T
Xg, Yg, blocks = [], [], []
for r in range(R):
    on = np.sort(rng.choice(np.arange(2, T-10), n, replace=False))
    S = np.zeros((T, n)); S[on, np.arange(n)] = 1
    Xr = np.apply_along_axis(lambda c: np.convolve(c, hrf)[:T], 0, S)
    P = np.c_[np.ones(T), np.linspace(-1,1,T), np.linspace(-1,1,T)**2, rng.standard_normal((T,3))]
    M = proj(P)
    beta = rng.standard_normal((n, V)) + 1
    Y = Xr@beta + 2*rng.standard_normal((T, V))
    Xm, Ym = M@Xr, M@Y
    full = np.zeros((T, R*n)); full[:, r*n:(r+1)*n] = Xm
    Xg.append(full); Yg.append(Ym); blocks.append((Xm, Ym))
Xg, Yg = np.vstack(Xg), np.vstack(Yg)
coef_g, alpha_g = fr.fracridge(Xg, Yg, fracs, jit=False)   # p x f x b

# block version, replicating fracridge's grid/interp exactly
svs, a_list, V_list = [], [], []
for Xm, Ym in blocks:
    U, s, Vt = np.linalg.svd(Xm, full_matrices=False)
    svs.append(s); a_list.append((U.T@Ym)/s[:,None]); V_list.append(Vt.T)
s_all = np.concatenate(svs); a_all = np.vstack(a_list)
val1 = fr.BIG_BIAS*s_all.max()**2; val2 = fr.SMALL_BIAS*s_all.min()**2
grid = np.concatenate([[0], 10**np.arange(np.floor(np.log10(val2)), np.ceil(np.log10(val1)), fr.BIAS_STEP)])
ssq = s_all**2; sclg_sq = (ssq/(ssq+grid[:,None]))**2
coef_b = np.zeros_like(coef_g); alpha_b = np.zeros_like(alpha_g); alpha_exact = np.zeros_like(alpha_g)
for v in range(V):
    nl = np.sqrt(sclg_sq @ a_all[:,v]**2); nl = nl/nl[0]
    al = np.exp(np.interp(fracs, nl[::-1], np.log1p(grid)[::-1])) - 1
    alpha_b[:, v] = al
    off = 0
    for (s, a, Vm) in zip(svs, a_list, V_list):
        k = len(s); sc = s**2/(s**2+al[:,None])          # f x k
        coef_b[off:off+k, :, v] = Vm @ (sc*a[:,v]).T; off += k
    for i,f in enumerate(fracs):
        if f >= 1: continue
        g = lambda lam: np.sqrt(np.sum(a_all[:,v]**2*(ssq/(ssq+lam))**2)/np.sum(a_all[:,v]**2)) - f
        alpha_exact[i, v] = brentq(g, 0, 1e12)
rel = lambda A,B: np.linalg.norm(A-B)/np.linalg.norm(B)
print("block vs global fracridge coef rel err:", rel(coef_b, coef_g))
print("block vs global alpha rel err:", rel(alpha_b, alpha_g))
m = fracs < 1
print("fracridge grid-interp alpha vs exact root: median |log ratio| = %.3g, max = %.3g" % (
    np.median(np.abs(np.log(alpha_g[m]/alpha_exact[m]))), np.max(np.abs(np.log(alpha_g[m]/alpha_exact[m])))))
# achieved fraction error under fracridge interp
ach = np.array([[np.sqrt(np.sum(a_all[:,v]**2*(ssq/(ssq+alpha_g[i,v]))**2)/np.sum(a_all[:,v]**2)) for v in range(V)] for i in range(len(fracs))])
print("fracridge achieved-vs-requested frac max abs err: %.3g" % np.max(np.abs(ach[m]-fracs[m,None])))

# ---- compiled calcbadness ----
runs = 6; per = 12; ncond = 15
validcolumns, stimix = [], []; c = 0
for r in range(runs):
    validcolumns.append(np.arange(c, c+per)); c += per
    stimix.append(rng.integers(0, ncond, per))
ntr = c; sess = [1,1,1,2,2,2]
results = [rng.standard_normal((V, ntr))*(1+0.1*k) + 0.5 for k in range(8)]
for xvals in ([np.array(i) for i in range(runs)], [np.array([0,1]), np.array([2,3]), np.array([4,5])]):
    ref = cb(xvals, validcolumns, stimix, results, sess)
    # compile: normalization from results[0]
    trial_run = np.concatenate([[r]*per for r in range(runs)]); cond = np.concatenate(stimix)
    mu = np.zeros((V, ntr)); sd = np.ones((V, ntr))
    for s_ in set(sess):
        cols = np.concatenate([validcolumns[r] for r in range(runs) if sess[r]==s_])
        mu[:, cols] = results[0][:, cols].mean(1, keepdims=True); sd[:, cols] = results[0][:, cols].std(1, ddof=1, keepdims=True)
    Z = [(x-mu)/sd for x in results]
    W = np.zeros((ntr, ntr))
    for te in xvals:
        te = np.atleast_1d(te); tr_ = np.setdiff1d(np.arange(runs), te)
        tem = np.isin(trial_run, te); trm = np.isin(trial_run, tr_)
        W += (trm[:,None] & tem[None,:] & (cond[:,None]==cond[None,:]))
    d = W.sum(1); used = d > 0
    mtar = np.where(used[None,:], (Z[0]@W.T)/np.maximum(d,1), 0)
    const = (Z[0]**2)@W.sum(0) - (mtar**2)@d
    comp = np.stack([((Zf[:,used]-mtar[:,used])**2)@d[used] for Zf in Z], 1) + const[:,None]
    print("compiled calcbadness rel err:", rel(comp, ref), " trials used:", used.sum(), "/", ntr)

