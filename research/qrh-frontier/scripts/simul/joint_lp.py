"""B4 pricing (exact piecewise-linear version).

Two zero witnesses per binding row (Lemma 1W applied to phi1 = psibar_u and phi2 = eta psibar_u,
truncations t1, t2 in [tmin, 3/2]); witnesses r_i <= t_i, m_i = (t_i - r_i)_+ <= 1/2.
Single-character count lines (paper B0.38/B0.40, with supply cap z <= L = l/d):
   inverse: A_I(r) = 1 - d r - 2q min((1-r)/2, L)    (r <= 1)    long(r) = (1+5r)/6 - d r  (r >= 1)
   plain:   A_P(m) = 1 - 2 d m - 2q min((1-2m)/(6k), L)
Hypothetical joint line  (J_lam):  sum_u |M_u(phi1;U^r1) M_u(phi2;U^r2)|^2 << U^{e(r1,r2)+eps},
   e(r1,r2) = max( E(r1), E(r2), 1 + lam (r1+r2-1)_+ ),   E(r) = max(1,(1+5r)/6)
   (the max with E(r_i) guarantees the hypothesis implies NOTHING new for one factor: no (*) leakage).
   Optional: MI-style slots in the joint Poisson range x=r1+r2<=1 (capacity min((1-x)/2,L)).
value(t1,t2) = max_{r1,r2} min(all lines)  -- solved exactly as a max over linear pieces of LPs.
R_new = min_{t1,t2} value(t1,t2).
"""
import numpy as np
from scipy.optimize import linprog, minimize

def pieces_single(t, d, q, k, L):
    """list of (lo, hi, [affine (a,b): a + b r ]) concave pieces of c(r) on [max(0,t-1/2), t]."""
    lo0 = max(0.0, t - 0.5)
    mcap = (1 - 6*k*L)/2                       # A_P capacity capped for m < mcap
    bps = sorted(set([lo0, t] + [b for b in (1 - 2*L, 1.0, t - mcap) if lo0 < b < t]))
    if len(bps) == 1: bps = [lo0, lo0]          # degenerate interval (t = 0)
    out = []
    for lo, hi in zip(bps[:-1], bps[1:]):
        mid = (lo + hi)/2; aff = []
        if mid <= 1:
            if (1 - mid)/2 <= L: aff.append((1 - q, -(d - q)))          # 1 - d r - q(1-r)
            else:                aff.append((1 - 2*q*L, -d))
        else:
            aff.append((1/6, 5/6 - d))
        m_mid = t - mid                                                  # m = t - r >= 0
        # A_P(t - r)
        if (1 - 2*m_mid)/(6*k) <= L:
            c = q/(3*k); aff.append((1 - c - 2*(d - c)*t, 2*(d - c)))    # 1-c-2(d-c)(t-r)
        else:
            aff.append((1 - 2*q*L - 2*d*t, 2*d))
        out.append((lo, hi, aff))
    return out

def joint_pieces(d, q, L, lam, slots):
    """affine pieces (a, b1, b2) of e(r1,r2) - d(r1+r2), as a MAX of affine functions (each used
    separately: value = max over pieces).  Poisson-range slotted line is a MIN-type improvement:
    returned separately as constraints valid when x<=1."""
    if not np.isfinite(lam): return [None]       # no joint line at all
    return [(1, -d, -d), (1/6, 5/6 - d, -d), (1/6, -d, 5/6 - d), (1 - lam, lam - d, lam - d)]

def lp_max(cons, box):
    """max v s.t. v <= a + b1 r1 + b2 r2 for (a,b1,b2) in cons, r in box. Returns v."""
    A = [[-b1, -b2, 1.0] for (a, b1, b2) in cons]; b = [a for (a, b1, b2) in cons]
    res = linprog([0, 0, -1], A_ub=A, b_ub=b, bounds=[box[0], box[1], (None, None)], method="highs")
    return -res.fun if res.status == 0 else -np.inf

def value(t1, t2, d, q, k, L, lam, slots=True, ret_arg=False):
    best = -np.inf
    S1 = pieces_single(t1, d, q, k, L); S2 = pieces_single(t2, d, q, k, L)
    JP = joint_pieces(d, q, L, lam, slots)
    for (l1, h1, a1) in S1:
        for (l2, h2, a2) in S2:
            base = [(a, b, 0.0) for (a, b) in a1] + [(a, 0.0, b) for (a, b) in a2]
            for jp in JP:
                cons = base + ([jp] if jp is not None else [])
                v = lp_max(cons, [(l1, h1), (l2, h2)])
                if slots and jp is not None and l1 + l2 < 1:
                    # Poisson-range MI-slotted joint line  1 - d x - 2q min((1-x)/2, L)
                    #   = MAX of two affine pieces (convex): treat each piece separately, region x<=1.
                    v_in = -np.inf
                    for extra in [(1 - q, -(d - q), -(d - q)), (1 - 2*q*L, -d, -d)]:
                        cc = cons + [extra]          # exact: J piece AND slotted piece both active
                        A = [[-c[1], -c[2], 1.0] for c in cc] + [[1, 1, 0]]
                        bb = [c[0] for c in cc] + [1.0]
                        r = linprog([0, 0, -1], A_ub=A, b_ub=bb, bounds=[(l1, h1), (l2, h2), (None, None)], method="highs")
                        if r.status == 0: v_in = max(v_in, -r.fun)
                    A = [[-c[1], -c[2], 1.0] for c in cons] + [[-1, -1, 0]]
                    bb = [c[0] for c in cons] + [-1.0]
                    r = linprog([0, 0, -1], A_ub=A, b_ub=bb, bounds=[(l1, h1), (l2, h2), (None, None)], method="highs")
                    v_out = -r.fun if r.status == 0 else -np.inf
                    v = max(v_in, v_out)
                best = max(best, v)
    return best

def R_paper(d, p, k):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = max(5/6-d, 0)
    return 1-d+d*C*V/(2*(V+d*C))

def R_new(d, p, k, L, lam, tmin=0.0, nt=16, slots=True):
    q = d*p
    ts = np.linspace(tmin, 1.5, nt)
    grid = [(value(a, b, d, q, k, L, lam, slots), a, b) for a in ts for b in ts if b >= a]
    grid.sort()
    best = grid[0]
    for (v0, a0, b0) in grid[:4]:
        f = lambda x: value(min(max(x[0], tmin), 1.5), min(max(x[1], tmin), 1.5), d, q, k, L, lam, slots)
        r = minimize(f, [a0, b0], method="Nelder-Mead", options=dict(xatol=1e-7, fatol=1e-10, maxiter=400))
        if r.fun < best[0]: best = (r.fun, *np.clip(r.x, tmin, 1.5))
    return best

if __name__ == "__main__":
    import sys; sys.path.insert(0, '.')
    from base import record
    P = record(); d0, k0 = float(P['d']), float(P['k']); L0 = float(P['l']/P['h'])
    R0 = R_paper(d0, .5, k0); q0 = d0/2
    print(f"delta0={d0:.10f} kappa0={k0:.10f} L0={L0:.6f} R_paper={R0:.12f}")
    # sanity 1: single character value at t* reproduces R_paper
    tst = 1 + d0*(2*(1-.5)*(1-.5/(3*k0))/(3-.5-2*.5/(3*k0)))/(2*((5/6-d0)+d0*(2*(1-.5)*(1-.5/(3*k0))/(3-.5-2*.5/(3*k0)))))
    print("t* =", tst, " value(t*,t*, no joint) =", value(tst, tst, d0, q0, k0, L0, np.inf, slots=False))
    print("min over t (no joint):", R_new(d0, .5, k0, L0, np.inf, tmin=0.0, slots=False))
    import sys as _s
    tmins = [float(a) for a in _s.argv[1:]] or [1.0, 0.0]
    for tmin in tmins:
        print(f"--- t_i >= {tmin}")
        for lam in [5/6, 0.75, 0.65, 0.6, 0.59, 0.58, 0.55, 0.5]:
            v, a, b = R_new(d0, .5, k0, L0, lam, tmin=tmin)
            print(f"lam={lam:.4f}  R_new={v:.9f}  dR={v-R0:+.4e}  t=({a:.4f},{b:.4f})", flush=True)
