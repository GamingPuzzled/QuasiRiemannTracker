"""Scan the hypothetical sixth-moment exponent e6 (and sharp-large-sieve e6=max(1,3n)) through cond.py."""
from mpmath import mp, mpf
from cond import theta_float, solve_smooth, solve_kink, sixth, Fd, geom
import numpy as np
b0 = solve_smooth(0)['theta']
def run(extra_f, extra_m, label):
    fl = theta_float(extra_f)
    if fl['active'] != 2 and fl['gap'] > 1e-4:
        print(f"{label}: Holder line active, no gain (float theta {fl['theta']:.13f})"); return fl
    try:
        g = (mpf(fl['l']), mpf(fl['delta']), mpf(fl['h']))
        ex = solve_smooth(2, extra_m, guess=g) if fl['gap'] > 1e-4 else solve_kink(0, 2, extra_m, guess=g)
        kind = 'smooth' if fl['gap'] > 1e-4 else 'kink'
        print(f"{label}: {kind} theta={mp.nstr(ex['theta'],16)} gain={mp.nstr(b0-ex['theta'],6)} "
              f"(float {fl['theta']:.13f}, full-grid G={fl['full']:.1e}) l={mp.nstr(ex['l'],10)} h={mp.nstr(ex['h'],10)} "
              f"n={mp.nstr(1-(ex['h']-2*ex['l'])/ex['h'],8)}")
    except Exception as err:
        print(label, "exact solve failed:", str(err)[:60], "float", fl)
    return fl
for e in ['1.19', '1.20', '1.21', '1.22', '1.23', '1.24', '1.25', '1.26']:
    run(sixth(float(e)), sixth(mpf(e)), f"e6={e}")
ls_f = lambda R, n, d, p, k: np.maximum(1, 3*n)/6 + 5*R/6
ls_m = lambda R, n, d, p, k: max(1, 3*n)/6 + 5*R/6
run(ls_f, ls_m, "sharp-LS e6=max(1,3n)")
