"""Pointwise square-root test (PSR) for the faithful proxy b_m = abar(m)^2 R(m).
DR Prop 9.2's D2-step needs, for each non-exceptional k, |sum_m b_m chi_m(k)| << Q^{1/2} * ||b||/sqrt(#m) (sqrt cancellation).
GRH per-mu bound only gives Q^{1/2} * sum_mu |rho_mu|.  We measure max_k |B(k)|^2 / ||b||^2 and the sixth-moment-free tail.
Reuses exp/sieve3.py machinery by exec of its setup part (no k-loop)."""
import sys, os, numpy as np
here = os.path.dirname(os.path.abspath(__file__)); expd = os.path.join(here, '..', 'exp')
sys.path.insert(0, expd)
Q = int(sys.argv[1]); KMAX = int(sys.argv[2]); D = int(sys.argv[3])
src = open(os.path.join(expd, 'sieve3.py')).read().split('vecs = {')[0]
sys.argv = ['x', str(Q), str(KMAX), str(D)]
g = {'__name__': 'sv'}; exec(compile(src, 'sieve3_setup', 'exec'), g)
mods, ks, ka, kb = g['mods'], g['ks'], g['ka'], g['kb']
rho, th, rng = g['rho'], g['th_true'], g['rng']
vecs = {'random': np.exp(2j*np.pi*rng.random(len(mods))), 'TH(ab2*R)': th,
        'ab2*GAUSS': g['alphabar2']*g['g6s'], 'GAUSS': g['g6s']}
B = np.stack(list(vecs.values()), axis=1)
nb = (np.abs(B)**2).sum(axis=0)
mx = np.zeros(B.shape[1]); s2 = np.zeros(B.shape[1]); n = 0
for st in range(0, len(ks), 20000):
    X = g['chi_matrix'](ka[st:st+20000], kb[st:st+20000])
    V = np.abs(X @ B)**2
    mx = np.maximum(mx, V.max(axis=0)); s2 += V.sum(axis=0); n += V.shape[0]
print(f"Q={Q} #m={len(mods)} #k={n} D={D} #mu={len(rho)} sum|rho|={np.abs(rho).sum():.2f} sum|rho|^2={np.sum(np.abs(rho)**2):.2f}")
for (name, mxv, s, b2) in zip(vecs, mx, s2, nb):
    print(f"  {name:12s} mean|B|^2/||b||^2={s/n/b2:.3f}  max|B|^2/||b||^2={mxv/b2:.2f}  (log #k={np.log(n):.1f})")
print(f"  GRH-per-mu pointwise scale (sum|rho|)^2/sum|rho|^2 = {np.abs(rho).sum()**2/np.sum(np.abs(rho)**2):.1f}")
