"""LP check of the P1 closed form (Poisson-range joint bound WITHOUT slots) at the binding point."""
import sys; sys.path.insert(0, '.')
from joint_lp import R_new, R_paper
from closed_forms import v_P1
from base import record
P = record(); d0, k0 = float(P['d']), float(P['k']); L0 = float(P['l']/P['h'])
for lam in (0.70, 0.674, 0.65, 0.60):
    v, a, b = R_new(d0, .5, k0, L0, lam, tmin=0.0, slots=False)
    print(f"lam={lam}: LP(no slots) {min(v, R_paper(d0,.5,k0)):.6f}  closed v_P1 {min(v_P1(d0,.5,k0,lam), R_paper(d0,.5,k0)):.6f}  t=({a:.4f},{b:.4f})", flush=True)
