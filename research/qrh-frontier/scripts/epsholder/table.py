"""Conditional-payoff table for B5 (robust float solver of cond.py, refined; accuracy calibrated against
the exact baseline).  Lines: lossy 2k-th moments  sum|S|^{2k}|Q|^2 << U^{1+(k-2)(delta-eta)n}."""
from mpmath import mp, mpf
from cond import theta_float, solve_smooth
from eta6 import l6, l8
mp.dps = 30
ex = solve_smooth(0)['theta']; rec = mpf(683505193)/781250000
fb = theta_float()
print(f"exact baseline b0 = {mp.nstr(ex, 18)}   record theta = {mp.nstr(rec, 18)}")
print(f"float baseline    = {fb['theta']:.15f}   |float-exact| = {abs(fb['theta']-float(ex)):.1e}  (solver accuracy)")
cal = solve_smooth(0, H1=mpf(1)/4-mpf('1.3e-6'))['theta']
print(f"calibration H1=1/4-1.3e-6 (exact): theta = {mp.nstr(cal, 16)}  gain {mp.nstr(ex-cal, 4)}")
for name, mk in (("6th", l6), ("8th", l8)):
    for eta in [0.0, 0.002, 0.00236, 0.0024, 0.0025, 0.003, 0.005, 0.01, 0.05]:
        r = theta_float(mk(eta))
        print(f"{name} eta={eta:<8} theta={r['theta']:.13f} gain={float(ex)-r['theta']:+.3e} active-line={r['active']} "
              f"delta*={r['delta']:.5f} h={r['h']:.6f} full-grid maxG={r['full']:.1e}")
# Lindelof-on-average sixth/eighth moments (eta = delta, i.e. e = 1)
import numpy as np
L6 = lambda R, n, d, p, k: (1-2*d*p*(1-2*n)/(6*k))/6 + 5*R/6
r = theta_float(L6)
print(f"6th Lindelof-on-average (e=1, slots kept): theta={r['theta']:.13f} gain={float(ex)-r['theta']:+.3e}")
