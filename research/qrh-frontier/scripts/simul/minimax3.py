"""B4 pricing, rigorous class space.  A binding-type row class is labelled by
   delta_N  (zero level of L(psibar_u): bounds the main-slot amplitude q <= delta_N/2, B0.16),
   delta_D  (zero level of L(eta psibar_u): fixes the s-line a = (1+delta_D)/2, B0.53),
   q        (main-slot amplitude, 0 <= q <= delta_N/2).
Count C = min( R_single(delta_N,q), R_single(delta_D,q) [only if q<=delta_D/2], JOINT(delta_N,delta_D,q) ).
JOINT (hypothesis J_lam, e = 1+lam(x-1)_+ on r_i<=1):
  lam > max(delta_N,delta_D):  v = (1-lam + A K)/(1+K),  K = (lam-dN)/(dN-q) + (lam-dD)/(dD-q), A=1-q
      (adversary's best: both witnesses at their MI crossings r_v^i = (A-v)/(d_i-q); regime-checked)
  otherwise: safe bound with t1=t2=1: max over r in [1/2,1]^2 of 1+lam(x-1)-dN r1-dD r2.
F (B0.55 at dyad d) = (1+dD)/2 + h(z0-1/6) - y/2 - l/2 + l q + d(1/4+3C/4-z0) - q min(d,2y-d)/(12k).
theta = min over geometry of max over classes (same outer solver as slack2.solve).
Paper (no joint) is recovered exactly: then max over classes is attained at dN=dD (sanity printed).
"""
import sys, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0, '.'); sys.path.insert(0, '..')
import numpy as np, slack2
from slack2 import DEF

def R_single(d, q, k):
    p = np.where(d > 0, q/d, 0)
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = np.maximum(5/6-d, 0)
    return 1-d+d*C*V/(2*(V+d*C))

def joint(dN, dD, q, lam, L):
    A = 1 - q; sN = dN - q; sD = dD - q
    with np.errstate(all='ignore'):
        K = (lam-dN)/sN + (lam-dD)/sD
        v = (1 - lam + A*K)/(1 + K)
        rN = (A - v)/sN; rD = (A - v)/sD
    ok = (lam > np.maximum(dN, dD)) & (sN > 0) & (sD > 0) & (rN <= 1) & (rD <= 1) & (rN + rD > 1) \
         & ((1-rN)/2 <= L) & ((1-rD)/2 <= L)
    # safe bound (t1=t2=1): affine max over the box [1/2,1]^2
    c1 = lam - dN; c2 = lam - dD
    safe = 1 - lam + np.where(c1 > 0, c1, c1/2) + np.where(c2 > 0, c2, c2/2)
    return np.where(ok, v, safe)

def worst3(beta, x, y, l, c, lam, nd=121, nq=21, L=0.2104):
    k = 2*beta-1; h = 1+l-x
    dD = np.linspace(0.02, k, nd)[:, None, None]
    dN = np.linspace(0.02, k, nd)[None, :, None]
    t = np.linspace(0, 1, nq)[None, None, :]
    q = 0.5*dN*t                                    # q in [0, dN/2]
    C = R_single(dN, q, k)
    C = np.minimum(C, np.where(q <= dD/2 + 1e-15, R_single(dD, q, k), 9))
    if lam is not None:
        C = np.minimum(C, joint(dN, dD, q, lam, L))
    a = (1 + dD)/2
    out = []
    for d in (0.005, h, min(h, y)):
        cap = max(0.0, min(d, 2*y-d))
        Fd = (a + h*(c["z0"]-1/6) - y/2 - l/2 + l*q + d*(0.25 + 0.75*C - c["z0"]) - q*cap/(12*k))
        out.append(Fd.max())
    out.append(0.51 + h*(c["z0"]-1/6) - y/2 - l/2 + l*0.01 + h*(1-c["z0"]))
    return max(out) - beta

def argmax3(beta, x, y, l, lam, nd=241, nq=41, L=0.2104):
    c = DEF; k = 2*beta-1; h = 1+l-x
    dD = np.linspace(0.02, k, nd)[:, None, None]; dN = np.linspace(0.02, k, nd)[None, :, None]
    q = 0.5*dN*np.linspace(0, 1, nq)[None, None, :]
    C = np.minimum(R_single(dN, q, k), np.where(q <= dD/2+1e-15, R_single(dD, q, k), 9))
    if lam is not None: C = np.minimum(C, joint(dN, dD, q, lam, L))
    F = (1+dD)/2 + h*(c["z0"]-1/6) - y/2 - l/2 + l*q + h*(0.25+0.75*C-c["z0"]) - q*(2*y-h)/(12*k)
    i = np.unravel_index(np.argmax(F), F.shape)
    return F[i]-beta, float(dD.ravel()[i[0]]), float(dN.ravel()[i[1]]), float(q.ravel()[i[1]*0] if False else (0.5*dN.ravel()[i[1]]*np.linspace(0,1,nq)[i[2]]))

def solve3(lam, nd=121, nq=21):
    orig = slack2.worst
    slack2.worst = lambda beta, x, y, l, c, nd_=None, npp=None: worst3(beta, x, y, l, c, lam, nd, nq)
    try:
        return slack2.solve()
    finally:
        slack2.worst = orig

if __name__ == "__main__":
    base = slack2.solve()
    print("slack2 baseline:", base[0])
    s0 = solve3(None)
    print("3-param class space, no joint line (must equal baseline up to grid):", s0[0])
    b, x, l = s0[0], s0[1], s0[2]; y = 1-l-x
    print("   argmax (F-b, dD, dN, q):", argmax3(b, x, y, l, None))
    print(f"{'lam':>6} {'theta':>13} {'gain vs 3p-baseline':>20}  argmax(F-b,dD,dN,q)")
    for lam in sys.argv[1:] and [float(a) for a in sys.argv[1:]] or [0.62, 0.60, 0.595, 0.5908, 0.585, 0.58, 0.57, 0.55, 0.53, 0.50]:
        r = solve3(lam)
        b, x, l = r[0], r[1], r[2]; y = 1-l-x
        print(f"{lam:6.4f} {r[0]:13.10f} {s0[0]-r[0]:+20.3e}  {argmax3(b, x, y, l, lam)}", flush=True)
