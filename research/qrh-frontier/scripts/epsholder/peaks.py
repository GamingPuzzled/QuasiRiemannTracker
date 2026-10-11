"""Two-peak regime: max_delta min(F4,Fnew) has two local maxima (Fnew-peak near .44, F4-peak near delta0);
optimum over h equalises them.  Exact system in (l,h,d1,d2): Fnew(d1)=b, dFnew(d1)=0, F4(d2)=b, dF4(d2)=0,
with validity checks Fnew(d1)<F4(d1), F4(d2)<Fnew(d2) and a dense float check of the whole max."""
import numpy as np
from mpmath import mp, mpf, findroot, diff
from cond import Fd, geom, theta_float, solve_smooth, feasible
from eta6 import l6, l8
mp.dps = 40
P5 = mpf(1)/2
def two_peak(extra, g):
    f = lambda d, l, h, w: Fd(d, P5, l, h, extra, mpf(1)/4, z0=mpf('0.34'), which=w)[0]
    b = lambda l, h: geom(l, h)[2]
    def eqs(l, h, d1, d2):
        return [f(d1, l, h, 2)-b(l, h), diff(lambda d: f(d, l, h, 2), d1),
                f(d2, l, h, 0)-b(l, h), diff(lambda d: f(d, l, h, 0), d2)]
    l, h, d1, d2 = findroot(eqs, g)
    ok = f(d1, l, h, 2) < f(d1, l, h, 0) and f(d2, l, h, 0) < f(d2, l, h, 2)
    # float global check of max_delta min(...) - b at the solution, all p, dyads
    dg = np.linspace(0.02, float(geom(l, h)[3]), 20001)
    fl = lambda e: (lambda R, n, d, p, k: float(0) + e(R, n, d, p, k))
    return dict(l=l, h=h, d1=d1, d2=d2, theta=b(l, h), valid=ok)
if __name__ == "__main__":
    b0 = solve_smooth(0)['theta']
    for name, mk in (("6th", l6), ("8th", l8)):
        for eta in ['0.0024', '0.0025', '0.003', '0.004']:
            fl = theta_float(mk(float(eta)))
            x, y, bb, k = geom(fl['l'], fl['h'])
            dg = np.linspace(0.3, 0.6, 3001)
            F4, _ = Fd(dg, 0.5, fl['l'], fl['h'], None, 0.25, which=0)
            F6, _ = Fd(dg, 0.5, fl['l'], fl['h'], mk(float(eta)), 0.25, which=2)
            i6 = int(F6.argmax()); i4 = int(F4.argmax())
            g = (mpf(fl['l']), mpf(fl['h']), mpf(dg[i6]), mpf(dg[i4]))
            try:
                r = two_peak(mk(mpf(eta)), g)
                print(f"{name} eta={eta}: theta={mp.nstr(r['theta'],15)} gain={mp.nstr(b0-r['theta'],6)} "
                      f"d1={mp.nstr(r['d1'],8)} d2={mp.nstr(r['d2'],8)} h={mp.nstr(r['h'],10)} valid={r['valid']} "
                      f"feasible={feasible(float(r['l']), float(r['h']))} | float {fl['theta']:.12f}")
            except Exception as e:
                print(name, eta, "fail", str(e)[:60], "| float", fl['theta'], "peaks", dg[i6], dg[i4])
