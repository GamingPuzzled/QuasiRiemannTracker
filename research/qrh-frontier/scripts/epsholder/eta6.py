"""Hypothetical 'lossy sixth moment':  sum_u |S_u(U^n)|^6 |Q_I|^2 << U^{1+(delta-eta) n}, same capacity zA.
eta = 0 is the trivial bound (max_B|S|^2 <= U^{delta n} times the MP fourth moment) -> no gain.
Holder: (1+(delta-eta)n-2q zA)/6 + 5R/6.  Also 8th: (1+2(delta-eta)n-2q zA)/8 + 7R/8."""
import numpy as np
from mpmath import mp, mpf
from cond import theta_float, solve_smooth, solve_kink, geom, feasible
mp.dps = 40
def l6(eta):
    return lambda R, n, d, p, k: (1+(d-eta)*n-2*d*p*(1-2*n)/(6*k))/6 + 5*R/6
def l8(eta):
    return lambda R, n, d, p, k: (1+2*(d-eta)*n-2*d*p*(1-2*n)/(6*k))/8 + 7*R/8
if __name__ == "__main__":
  b0 = solve_smooth(0)['theta']
  print("baseline", mp.nstr(b0, 16))
  for name, mk in (("6th", l6), ("8th", l8)):
   for eta in ['0', '0.002', '0.0024', '0.0025', '0.003', '0.005', '0.01', '0.02', '0.05', '0.1']:
    fl = theta_float(mk(float(eta)))
    l, h = fl['l'], fl['h']; x = 1+l-h; y = 1-l-x
    msg = f"{name} eta={eta:6s} float theta={fl['theta']:.13f} act={fl['active']} gap={fl['gap']:.1e} y-x={y-x:.3e} delta*={fl['delta']:.5f}"
    if fl['active'] == 2 or fl['gap'] < 1e-4:
        g = (mpf(l), mpf(fl['delta']), mpf(h))
        try:
            ex = solve_smooth(2, mk(mpf(eta)), guess=g) if fl['gap'] > 1e-4 else solve_kink(0, 2, mk(mpf(eta)), guess=g)
            ok = feasible(float(ex['l']), float(ex['h']))
            msg += f" | exact theta={mp.nstr(ex['theta'],15)} gain={mp.nstr(b0-ex['theta'],6)} feasible={ok}"
        except Exception as err:
            msg += " | exact fail " + str(err)[:30]
    print(msg)
