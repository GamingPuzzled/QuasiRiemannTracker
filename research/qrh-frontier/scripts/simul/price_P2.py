"""Full minimax theta for the single-character two-length inverse product hypotheses
   MI2(Poisson range)  +  slope-lam bound beyond it   (family P2), and P0/P1 for comparison.
Row count R(delta,p) = min(R_paper, v_P2 [regime-checked], v_P0 [regime-checked]); for lam <= delta
the safe bound min(R_paper, 1-delta) (t1=t2=1).  Outer optimizer = slack2.solve (harness model)."""
import sys, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0, '.'); sys.path.insert(0, '..')
import numpy as np, slack2
from closed_forms import v_P0, v_P1, v_P2, _coef
R_orig = slack2.R_k

def regime(d, p, k, lam, v, L=0.2104):
    r, m, t, xb, _ = _coef(d, p, k, lam)
    rv = r[0]+r[1]*v; mv = m[0]+m[1]*v; tv = t[0]+t[1]*v; xbv = xb[0]+xb[1]*v
    return (lam > d) & (rv > .5) & (rv <= 1) & ((1-rv)/2 <= L) & (mv > 0) & (mv <= .5) & (tv > 1) \
           & ((1-2*mv)/(6*k) <= L) & (xbv > 1)

def make_R(lam, fam="P2"):
    def R(delta, p, k, c):
        R0 = R_orig(delta, p, k, c)
        d = np.broadcast_to(delta, np.broadcast(delta, p).shape).astype(float); pp = np.broadcast_to(p, d.shape).astype(float)
        out = R0.copy()
        with np.errstate(all='ignore'):
            fs = {"P0": [v_P0], "P1": [v_P0, v_P1], "P2": [v_P0, v_P1, v_P2]}[fam]
            for f in fs:
                v = f(d, pp, k, lam)
                ok = regime(d, pp, k, lam, v) & (v < out)
                out = np.where(ok, v, out)
        out = np.where(lam <= d, np.minimum(out, 1-d), out)
        return out
    return R

if __name__ == "__main__":
    fam = sys.argv[1] if len(sys.argv) > 1 else "P2"
    lams = [float(a) for a in sys.argv[2:]] or [5/6, 0.83, 0.82, 0.81, 0.80, 0.78, 0.75, 0.70, 0.65, 0.60]
    base = slack2.solve()[0]
    print(f"family {fam}; baseline (grid) {base:.10f}")
    for lam in lams:
        slack2.R_k = make_R(lam, fam); r = slack2.solve(); slack2.R_k = R_orig
        print(f"lam={lam:.4f} theta={r[0]:.10f} gain={base-r[0]:+.3e} x={r[1]:.5f} l={r[2]:.6f}", flush=True)
