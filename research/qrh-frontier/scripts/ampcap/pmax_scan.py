"""Amplitude-cap lever: (1) theta vs pmax, (2) decomposition of dF/dp at the binding point,
(3) Delta p needed for target gains. Run from scripts/ampcap/."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np
from slack2_ampcap import solve, R_k, DEF

base = solve()[0]
print(f"baseline b0 = {base:.10f}")
for pm in (0.51, 0.5, 0.499, 0.495, 0.49, 0.48, 0.45, 0.40, 0.30):
    r = solve({"pmax": pm})
    print(f"pmax={pm:.3f}  b0={r[0]:.10f}  gain={base-r[0]:+.3e}  x={r[1]:.5f} l={r[2]:.6f}")

# binding point (paper B0.45)
l = 0.1671200785068101428882107
k = (5-3*l)/6
d0 = (25-60*l+27*l**2)/(6*(5+4*l-3*l**2))
C = (6*k-1)/(15*k-2); V = 5/6-d0
Rp = -1 + C*(V**2-d0**2*C)/(2*(V+d0*C)**2)
h = -((1+l)/2 + l/(6*k))/(3*Rp/4 - 1/(24*k))
x = 1+l-h; y = h-2*l
def F4(delta, p):
    R = R_k(delta, p, k, DEF)
    return (1+l)/2 + delta/2 + h*(3*R/4-5/12) + delta*p*(l-(h-4*l)/(12*k))
eps = 1e-6
dRdp = (R_k(d0, .5, k, DEF)-R_k(d0, .5-eps, k, DEF))/eps
amp = d0*l; mark = -d0*(h-4*l)/(12*k); cnt = h*0.75*dRdp
print(f"\nbinding: delta0={d0:.6f} kappa0={k:.6f} h0={h:.6f} F4={F4(d0,.5):.10f} b0={11/12-l/4:.10f}")
print(f"dF/dp at p=1/2: amplitude cost l*delta={amp:+.5f}, numerator mark gain={mark:+.5f}, "
      f"row-count gain (3h/4)dR/dp={cnt:+.5f}, total={amp+mark+cnt:+.5f}")
print(f"  (finite diff total = {(F4(d0,.5)-F4(d0,.5-eps))/eps:+.5f})")
tot = amp+mark+cnt
for g in (1e-6, 3e-5, 1e-4):
    print(f"theta gain {g:.0e} needs Delta p ~ {g/0.0075:.2e} (Delta g_i ~ {d0*g/0.0075:.2e} in prime-sum exponent)")
