"""B4 pricing -> theta.  Closed-form row count with the joint (or fourth-moment) line, then a full
re-solve of the record's binding system.

Closed form (derived in reports/solverI-simul.md, section 2; checked against joint_lp.py):
for lam in (delta, lam_be), with s = delta - q, A = 1 - q  (q = delta p),
    v(lam) = [ (2A - s)(lam - delta) + s(1 - delta) ] / [ s + 2(lam - delta) ],
which is the level where the adversary's largest joint length x = 2 r_v (r_v = (A-v)/s, the MI
crossing) meets the joint line (1-delta) + (lam-delta)(x-1).  R_new = min(R_paper, v) in the regime
  (1 - r_v)/2 <= L (MI capacity not supply-capped), 1/2 < r_v <= 1.
Break-even:  lam_be = delta + (R - (1-delta)) / (2 r_R - 1),  r_R = (A - R)/s.
"""
import numpy as np, sys
sys.path.insert(0, '.')
from mpmath import mp, mpf, findroot, diff
mp.dps = 30

def Rp(d, p, k):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = 5/6 - d
    return 1-d+d*C*V/(2*(V+d*C))

def v_closed(d, p, lam):
    q = d*p; s = d - q; A = 1 - q
    return ((2*A - s)*(lam - d) + s*(1 - d))/(s + 2*(lam - d))

def R_joint(d, p, k, lam, L=0.2104):
    R0 = Rp(d, p, k)
    if lam <= d: return R0      # other regime, not needed for pricing above delta
    v = v_closed(d, p, lam)
    q = d*p; s = d - q; A = 1 - q; rv = (A - v)/s
    ok = (v < R0) and (0.5 < rv <= 1) and ((1 - rv)/2 <= L)
    return v if ok else R0

def lam_be(d, p, k):
    R0 = Rp(d, p, k); q = d*p; s = d - q; A = 1 - q; rR = (A - R0)/s
    return d + (R0 - (1 - d))/(2*rR - 1)

# ---- full re-solve: same binding system as smallgain/model.py, R replaced by R_joint
def F_of(d, p, l, h, lam):
    x = 1 + l - h; y = 1 - l - x; b = x/2 + (y - x)/12 + mpf(5)/6 - x/3 - l/6; k = 2*b - 1
    R = R_joint(d, p, k, lam); q = d*p
    zAh = (2*y - h)/(6*k)
    return (mpf(1)/2 + d/2 - y/2 - l/2 + l*q + h*(mpf(1)/4 + 3*R/4 - mpf(1)/6) - q*zAh/2), b

def solve(lam, guess):
    p = mpf(1)/2
    G = lambda l, d, h: F_of(d, p, l, h, lam)[0] - F_of(d, p, l, h, lam)[1]
    B = lambda l, h: F_of(mpf('0.4'), p, l, h, lam)[1]
    def eqs(l, d, h):
        return [G(l, d, h),
                diff(lambda dd: F_of(dd, p, l, h, lam)[0], d),
                diff(lambda hh: B(l, hh), h)*diff(lambda ll: G(ll, d, h), l)
                - diff(lambda ll: B(ll, h), l)*diff(lambda hh: G(l, d, hh), h)]
    l, d, h = findroot(eqs, guess)
    x = 1 + l - h
    return dict(l=l, d=d, h=h, theta=mpf(11)/12 - l/4, p_check=None)

def worst_check(sol, lam):
    """grid check that the binding point is the global max of F over (delta,p) at d=h."""
    l, h = sol['l'], sol['h']; b = sol['theta']; k = 2*b - 1
    best = -1
    for d in np.linspace(0.02, float(k), 300):
        for p in np.linspace(0, .5, 26):
            F, _ = F_of(mpf(d), mpf(p), l, h, lam)
            best = max(best, float(F))
    return best - float(b)

if __name__ == "__main__":
    from base import record
    P = record(); d0, k0 = float(P['d']), float(P['k'])
    print("break-even lam at binding point:", lam_be(d0, .5, k0))
    print("check vs model: lam=inf ->", end=" ")
    g = (P['l'], P['d'], P['h'])
    s = solve(10.0, g); print(mp.nstr(s['theta'], 15), " paper", mp.nstr(P['b'], 15))
    print(f"{'lam':>7} {'R_new(bind)':>12} {'dR':>11} {'theta':>16} {'gain':>11} {'delta*':>9} {'gridcheck':>10}")
    for lam in [0.60, 0.59, 0.585, 0.58, 0.57, 0.56, 0.55, 0.53, 0.5]:
        try:
            s = solve(lam, g); g = (s['l'], s['d'], s['h'])
            k = 2*float(s['theta']) - 1
            Rn = R_joint(float(s['d']), .5, k, lam); R0 = Rp(float(s['d']), .5, k)
            print(f"{lam:7.4f} {Rn:12.6f} {Rn-R0:+11.3e} {mp.nstr(s['theta'],12):>16} {float(P['b']-s['theta']):+11.3e} {float(s['d']):9.5f} {worst_check(s, lam):+10.2e}")
        except Exception as e:
            print(lam, "fail", e)
