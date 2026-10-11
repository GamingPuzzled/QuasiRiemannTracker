"""solverJ / B5: exponent model of the record with a modified numerator exponent.

Copy of smallgain/model.py (binding system B0.27, B0.38-B0.41, B0.46-B0.47 with M + l = 1), but the
numerator exponent (U-units) at the binding dyad d = h is   Num = min(H4, PW, extra)   where
    H4 = H1 + (1-H1) R - q zA/2,  H1 = 1/4   (B0.41)
    PW = R + delta n/2                       (B0.15, n = 1 - y/h)
    extra = hypothetical 2k-th moment line e/(2k) + (1-1/(2k)) R   (no slots)
Two stages:
  1. ROBUST float search (numpy): theta = min b(l,h) s.t. max_{delta, p, d-dyads} F - b <= 0. Detects
     which line is active at the optimum and whether the delta-max is smooth or a kink.
  2. EXACT mpmath solve (dps 40) of the corresponding 3x3 system (smooth or kink), 12+ digits.
"""
import numpy as np
from mpmath import mp, mpf, findroot, diff
mp.dps = 40

# ---------------- shared formulas (work for float and mpf) ----------------
def geom(l, h):
    x = 1+l-h; y = 1-l-x
    b = x/2 + (y-x)/12 + 5/6.0 - x/3 - l/6 if isinstance(l, float) else \
        x/2 + (y-x)/12 + mpf(5)/6 - x/3 - l/6
    return x, y, b, 2*b-1

def Rk(d, p, k):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = 5/6.0-d if isinstance(d, (float, np.ndarray)) else mpf(5)/6-d
    return 1-d+d*C*V/(2*(V+d*C))

def lines(d, p, k, n, extra, H1):
    R = Rk(d, p, k); q = d*p; zA = (1-2*n)/(6*k)
    L = [H1 + (1-H1)*R - q*zA/2, R + d*n/2]
    if extra is not None: L.append(extra(R, n, d, p, k))
    return L

def Fd(d, p, l, h, extra=None, H1=0.25, dy=None, z0=0.34, which=None):
    x, y, b, k = geom(l, h)
    dy = h if dy is None else dy
    n = 1-y/dy if y < dy else 0*dy
    L = lines(d, p, k, n, extra, H1)
    if which is None:
        Num = L[0]
        for t in L[1:]: Num = np.minimum(Num, t) if isinstance(Num, np.ndarray) else min(Num, t)
    else:
        Num = L[which]
    return (1+d)/2 + h*(z0-1/6.0) - y/2 - l/2 + l*d*p + dy*(Num - z0), b

# ---------------- stage 1: robust float ----------------
DG = np.linspace(0.02, 0.749, 4000)
def G(l, h, extra, H1, full=False):
    x, y, b, k = geom(l, h)
    dg = DG[DG <= k]
    if not full:
        v, _ = Fd(dg, 0.5, l, h, extra, H1)
        i = int(v.argmax())
        a_, c_ = dg[max(i-1, 0)], dg[min(i+1, len(dg)-1)]
        f = lambda t: -Fd(np.array([t]), 0.5, l, h, extra, H1)[0][0]
        t = gold(f, a_, c_, 50)
        return max(v.max(), -f(t)) - b
    worst = -9
    for p in np.linspace(0, 0.5, 26):
        for dy in list(np.linspace(0.005, h, 40)):
            v, _ = Fd(dg, p, l, h, extra, H1, dy=dy)
            worst = max(worst, v.max() - b)
    return worst

def gold(f, a, c, it=80):
    g = (5**0.5-1)/2
    for _ in range(it):
        x1 = c-g*(c-a); x2 = a+g*(c-a)
        if f(x1) < f(x2): c = x2
        else: a = x1
    return (a+c)/2

def feasible(l, h):
    """paper geometry (slack2.feasible with M+l=1): l<x<=y, y-2x+l<0, l/h>4/21, y<4l."""
    x = 1+l-h; y = 1-l-x
    return l < x <= y and y-2*x+l < 0 and l/h > 4/21 and y < 4*l

def theta_float(extra=None, H1=0.25):
    def Gl(l):
        hs = [hh for hh in np.linspace(0.70, 0.90, 201) if feasible(l, hh)]
        if not hs: return 9.0, None
        vals = [G(l, hh, extra, H1) for hh in hs]
        i = int(np.argmin(vals))
        a_, c_ = hs[max(i-1, 0)], hs[min(i+1, len(hs)-1)]
        hh = gold(lambda t: G(l, t, extra, H1) if feasible(l, t) else 9.0, a_, c_, 60)
        return G(l, hh, extra, H1), hh
    lo, hi = 0.15, 0.20
    for _ in range(45):
        mid = (lo+hi)/2
        (lo, hi) = (mid, hi) if Gl(mid)[0] <= 0 else (lo, mid)
    g, hh = Gl(lo)
    x, y, b, k = geom(lo, hh)
    dg = DG[DG <= k]
    Ls = lines(dg, 0.5, k, 1-y/hh, extra, H1)
    v, _ = Fd(dg, 0.5, lo, hh, extra, H1)
    i = int(v.argmax()); act = int(np.argmin([Lj[i] if np.ndim(Lj) else Lj for Lj in Ls]))
    gaps = sorted([float((Lj[i] if np.ndim(Lj) else Lj)) for Lj in Ls])
    return dict(l=lo, h=hh, theta=b, delta=dg[i], active=act, gap=gaps[1]-gaps[0],
                full=G(lo, hh, extra, H1, full=True))

# ---------------- stage 2: exact mpmath systems ----------------
P5 = mpf(1)/2
def solve_smooth(which, extra=None, H1=mpf(1)/4, guess=None):
    f = lambda d, l, h: Fd(d, P5, l, h, extra, H1, z0=mpf('0.34'), which=which)
    B = lambda l, h: geom(l, h)[2]
    def eqs(l, d, h):
        Fv, b = f(d, l, h)
        Gf = lambda ll, hh: f(d, ll, hh)[0] - B(ll, hh)
        lag = (diff(lambda hh: B(l, hh), h)*diff(lambda ll: Gf(ll, h), l)
               - diff(lambda ll: B(ll, h), l)*diff(lambda hh: Gf(l, hh), h))
        return [Fv-b, diff(lambda dd: f(dd, l, h)[0], d), lag]
    g = guess or (mpf('0.16712007850681'), mpf('0.4693448586981'), mpf('0.7939475973881'))
    l, d, h = findroot(eqs, g)
    return dict(l=l, delta=d, h=h, theta=geom(l, h)[2])

def solve_kink(i, j, extra=None, H1=mpf(1)/4, guess=None):
    """max over delta of min(line i, line j) attained at their crossing delta_c(l,h)."""
    f = lambda d, l, h, w: Fd(d, P5, l, h, extra, H1, z0=mpf('0.34'), which=w)[0]
    B = lambda l, h: geom(l, h)[2]
    def dc(l, h, d0):
        return findroot(lambda d: f(d, l, h, i)-f(d, l, h, j), d0)
    def eqs(l, d, h):
        Gf = lambda ll, hh: f(dc(ll, hh, d), ll, hh, i) - B(ll, hh)
        lag = (diff(lambda hh: B(l, hh), h)*diff(lambda ll: Gf(ll, h), l)
               - diff(lambda ll: B(ll, h), l)*diff(lambda hh: Gf(l, hh), h))
        return [f(d, l, h, i)-f(d, l, h, j), f(d, l, h, i)-B(l, h), lag]
    l, d, h = findroot(eqs, guess)
    return dict(l=l, delta=d, h=h, theta=geom(l, h)[2])

def sixth(e):        # hypothetical sum_rows |S|^6 << U^e (no slots)
    return lambda R, n, d, p, k: e/6 + 5*R/6

if __name__ == "__main__":
    rec = mpf(683505193)/781250000
    b0 = solve_smooth(0)
    print("exact baseline theta =", mp.nstr(b0['theta'], 20), " record-theta =", mp.nstr(rec-b0['theta'], 5))
    fb = theta_float()
    print("float robust baseline theta =", repr(fb['theta']), "active", fb['active'], "full-grid max G =", fb['full'])
    cal = solve_smooth(0, H1=mpf(1)/4-mpf('1.3e-6'))
    print("H1=1/4-1.3e-6: theta =", mp.nstr(cal['theta'], 16), " gain =", mp.nstr(b0['theta']-cal['theta'], 6))
    for e in ['1.1801', '1.18', '1.17', '1.15', '1.10', '1.00']:
        fl = theta_float(sixth(float(e)))
        tag = "smooth" if fl['active'] == 2 and fl['gap'] > 1e-4 else ("kink" if fl['active'] == 2 or fl['gap'] < 1e-4 else "H4")
        try:
            if fl['active'] != 2 and fl['gap'] > 1e-4:
                ex = b0
            elif tag == "smooth":
                ex = solve_smooth(2, sixth(mpf(e)), guess=(mpf(fl['l']), mpf(fl['delta']), mpf(fl['h'])))
            else:
                ex = solve_kink(0, 2, sixth(mpf(e)), guess=(mpf(fl['l']), mpf(fl['delta']), mpf(fl['h'])))
            th = mp.nstr(ex['theta'], 16); gn = mp.nstr(b0['theta']-ex['theta'], 6)
        except Exception as err:
            th, gn = 'fail', str(err)[:40]
        print(f"e6={e}: float theta={fl['theta']:.13f} active={fl['active']} ({tag}) full-grid G={fl['full']:.2e}"
              f" | exact theta={th} gain={gn}")
