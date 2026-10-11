"""Exact (mpmath) kink solve: theta-optimum where max_delta min(F_i, F_j) sits at the crossing delta_c.
Since b = 11/12 - l/4 depends only on l (M + l = 1), the optimum satisfies  G(l,h)=0, d_h G(l,h)=0,
with G(l,h) = F_i(delta_c(l,h); l,h) - b(l)."""
from mpmath import mp, mpf, findroot, diff
from cond import Fd, geom, theta_float, solve_smooth
from eta6 import l6, l8
mp.dps = 40
P5 = mpf(1)/2
def kink(i, j, extra, g):
    f = lambda d, l, h, w: Fd(d, P5, l, h, extra, mpf(1)/4, z0=mpf('0.34'), which=w)[0]
    state = {'d': g[1]}
    def dc(l, h):
        state['d'] = findroot(lambda d: f(d, l, h, i)-f(d, l, h, j), state['d'], tol=mpf(10)**-35)
        return state['d']
    G = lambda l, h: f(dc(l, h), l, h, i) - geom(l, h)[2]
    l, h = findroot(lambda l, h: [G(l, h), diff(lambda hh: G(l, hh), h)], (g[0], g[2]))
    d = dc(l, h)
    # sanity: crossing is a max of the min (F_i increasing, F_j decreasing, or vice versa)
    si = diff(lambda dd: f(dd, l, h, i), d); sj = diff(lambda dd: f(dd, l, h, j), d)
    return dict(l=l, h=h, delta=d, theta=geom(l, h)[2], slopes=(si, sj))
if __name__ == "__main__":
    b0 = solve_smooth(0)['theta']
    for name, mk in (("6th", l6), ("8th", l8)):
        for eta in ['0.00237', '0.0024', '0.0025', '0.003', '0.004']:
            fl = theta_float(mk(float(eta)))
            try:
                r = kink(0, 2, mk(mpf(eta)), (mpf(fl['l']), mpf(fl['delta']), mpf(fl['h'])))
                print(f"{name} eta={eta}: exact theta={mp.nstr(r['theta'],16)} gain={mp.nstr(b0-r['theta'],8)} "
                      f"delta_c={mp.nstr(r['delta'],10)} slopes(F4,Fnew)=({mp.nstr(r['slopes'][0],4)},{mp.nstr(r['slopes'][1],4)})"
                      f" | float {fl['theta']:.13f}")
            except Exception as e:
                print(name, eta, "fail", str(e)[:50], fl['theta'])
