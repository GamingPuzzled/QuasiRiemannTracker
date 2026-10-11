"""B4 pricing, v2: supply caps on slot capacity (z <= L = l/d), truncations t_i in [tmin, 3/2]
(Lemma 1W's proof is valid for any fixed t0 in (0,3/2]; the paper uses [1,3/2]).
Joint family: e(x) = 1 + lam*(x-1)_+ for x = r1+r2 (x<=1: Poisson range; MI-style slots allowed
there with capacity min((1-x)/2, L)); 'lift' family: e(x) = max(1,(1+5x)/6) (= lam 5/6).
Exact adversary via 1-D level sets (no 2-D grid):
  value(t1,t2) = max v  s.t.  exists r1 in S1(v), r2 in S2(v), J(r1+r2) >= v.
J is V-shaped in x (max of two affine pieces + Poisson-range piece), so for a level v the set
{x : J(x)>=v} is (-inf, xa(v)] U [xb(v), inf); feasibility <=> minS1+minS2 <= xa or maxS1+maxS2 >= xb,
where S_i(v) may be a union of two intervals -> we use exact interval endpoints from a fine grid.
"""
import numpy as np
from functools import lru_cache

def single(r, t, d, q, k, L):
    m = np.maximum(t - r, 0.0)
    zI = np.minimum((1 - r)/2, L)
    AI = 1 - d*r - 2*q*zI
    inv = np.where(r <= 1, AI, (1 + 5*r)/6 - d*r)
    zP = np.minimum((1 - 2*m)/(6*k), L)
    AP = 1 - 2*d*m - 2*q*zP
    return np.minimum(inv, AP)

def Jline(x, d, q, k, L, lam, slots=True):
    if np.isinf(lam):
        J = np.where(x <= 1, 1 - d*x, np.inf)
    else:
        J = 1 + lam*np.maximum(x - 1, 0) - d*x
    if slots:
        zI = np.minimum((1 - x)/2, L)
        J = np.minimum(J, np.where(x <= 1, 1 - d*x - 2*q*zI, np.inf))
    return J

def adversary(t1, t2, d, q, k, L, lam, nr=1201, slots=True):
    r1 = np.linspace(max(0, t1 - .5), t1, nr); c1 = single(r1, t1, d, q, k, L)
    r2 = np.linspace(max(0, t2 - .5), t2, nr); c2 = single(r2, t2, d, q, k, L)
    # candidate levels: bisection on v
    lo, hi = 0.0, min(c1.max(), c2.max())
    xs = np.linspace(0, 3.2, 6401); Jx = Jline(xs, d, q, k, L, lam, slots)
    def feas(v):
        s1 = r1[c1 >= v]; s2 = r2[c2 >= v]
        if len(s1) == 0 or len(s2) == 0: return False
        # sumset of finite grids, coarse: use sorted unique sums of endpoints of runs
        X = (s1[:, None] + s2[None, :]).ravel() if len(s1)*len(s2) < 4e5 else None
        if X is None:
            a = s1[::max(1, len(s1)//600)]; b = s2[::max(1, len(s2)//600)]
            X = (a[:, None] + b[None, :]).ravel()
        return np.any(np.interp(X, xs, np.where(np.isinf(Jx), 9.0, Jx)) >= v)
    for _ in range(40):
        mid = (lo + hi)/2
        lo, hi = (mid, hi) if feas(mid) else (lo, mid)
    return lo

def R_paper(d, p, k):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = max(5/6-d, 0)
    return 1-d+d*C*V/(2*(V+d*C))

def R_new(d, p, k, L, lam, tmin=0.0, nt=31, slots=True, nr=401):
    q = d*p
    ts = np.linspace(max(tmin, 1e-3), 1.5, nt)
    best = (np.inf, None)
    for a in ts:
        for b in ts:
            if b < a: continue
            v = adversary(a, b, d, q, k, L, lam, nr, slots)
            if v < best[0]: best = (v, (a, b))
    h = (1.5 - max(tmin, 1e-3))/(nt - 1)
    for it in range(3):
        a0, b0 = best[1]
        for a in np.linspace(max(tmin, a0 - h), min(1.5, a0 + h), 7):
            for b in np.linspace(max(tmin, b0 - h), min(1.5, b0 + h), 7):
                v = adversary(a, b, d, q, k, L, lam, nr, slots)
                if v < best[0]: best = (v, (a, b))
        h /= 3
    return min(best[0], R_paper(d, p, k)), best

if __name__ == "__main__":
    import sys; sys.path.insert(0, '.')
    from base import record
    P = record(); d0, k0 = float(P['d']), float(P['k']); L0 = float(P['l']/P['h'])
    R0 = R_paper(d0, .5, k0)
    print(f"delta0={d0:.10f} kappa0={k0:.10f} L0={L0:.6f} R_paper={R0:.10f}")
    # sanity: no joint line
    print("sanity (lam=inf, no slots, t in [1,1.5]):", R_new(d0, .5, k0, L0, np.inf, tmin=1.0, slots=False)[1])
    for tmin in (1.0, 0.0):
        print(f"--- truncations t_i >= {tmin}")
        for lam in [np.inf, 1.0, 0.9, 5/6, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45]:
            Rn, (v, tt) = R_new(d0, .5, k0, L0, lam, tmin=tmin)
            print(f"lam={lam:7.4f} R_new={Rn:.6f} dR={Rn-R0:+.4e} adv={v:.6f} t=({tt[0]:.4f},{tt[1]:.4f})", flush=True)
