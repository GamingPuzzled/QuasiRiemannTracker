"""Full minimax re-solve (slack2.solve, the harness's exponent model) with the joint/4th-moment line.
R_k is patched by R_joint (closed form; valid regime-checked, falls back to the paper's R).
For delta >= lam we use the safe bound min(R_paper, 1-delta) (t1=t2=1 forces x>=1)."""
import sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, '.'); sys.path.insert(0, '..')
import numpy as np, slack2
from price_theta import v_closed

R_orig = slack2.R_k
def make_R(lam, L=0.2104):
    def R(delta, p, k, c):
        R0 = R_orig(delta, p, k, c)
        delta = np.broadcast_to(delta, np.broadcast(delta, p).shape); p = np.broadcast_to(p, delta.shape)
        q = delta*p; s = delta - q; A = 1 - q
        with np.errstate(all='ignore'):
            v = ((2*A - s)*(lam - delta) + s*(1 - delta))/(s + 2*(lam - delta))
            rv = (A - v)/s
        ok = (lam > delta) & (v < R0) & (rv > .5) & (rv <= 1) & ((1 - rv)/2 <= L) & (s > 0)
        # lam <= delta: t1=t2=1 forces x=r1+r2>=1 with r_i<=1, joint line <= 1-delta (safe bound)
        low = np.where(lam <= delta, np.minimum(R0, 1 - delta), R0)
        return np.where(ok, v, low)
    return R

if __name__ == "__main__":
    base = slack2.solve()
    print("baseline (grid model):", base[0], " exact record 0.87488664704")
    print(f"{'lam':>6} {'theta':>13} {'gain':>11}")
    for lam in [0.60, 0.5908, 0.59, 0.585, 0.58, 0.57, 0.56, 0.55, 0.53, 0.50]:
        slack2.R_k = make_R(lam)
        r = slack2.solve()
        print(f"{lam:6.4f} {r[0]:13.10f} {base[0]-r[0]:+11.3e}   x={r[1]:.5f} l={r[2]:.6f}", flush=True)
    slack2.R_k = R_orig
