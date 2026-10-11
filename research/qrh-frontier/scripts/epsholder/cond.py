"""solverJ / B5: exponent model of the record with a modified numerator exponent.

Copy of smallgain/model.py (binding system B0.27, B0.38-B0.41, B0.46-B0.47 with M + l = 1), but the
numerator exponent (U-units) is   Num = min(H4, PW, extra lines)   where
    H4 = 1/4 + 3R/4 - q zA/2           (B0.41, fourth moment + Holder)
    PW = R + delta*n/2                 (B0.15 pointwise; n = reflected numerator length)
and an optional hypothetical 2k-th moment line  e/(2k) + (1-1/(2k)) R   (no slots).

theta is obtained by a ROBUST nested optimisation (no assumption on which line is active):
    theta = min b(l,h)  s.t.  max_{delta,p,d} F - b <= 0 ,
using a Lagrange/KKT-free bisection on l and golden search on h. mpmath, dps 30.
"""
from mpmath import mp, mpf, findroot, diff
mp.dps = 30
THIRD = mpf(1)/3

def geom(l, h):
    x = 1+l-h; y = 1-l-x
    b = x/2 + (y-x)/12 + mpf(5)/6 - x/3 - l/6
    return x, y, b, 2*b-1

def Rk(d, p, k, V0=mpf(5)/6):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = max(V0-d, 0)
    return 1-d+d*C*V/(2*(V+d*C))

def num(d, p, k, n, extra=None, H1=mpf(1)/4):
    """numerator exponent in U units at reflected length n (0<=n<=1/2)."""
    R = Rk(d, p, k); q = d*p
    zA = (1-2*n)/(6*k)
    lines = [H1 + (1-H1)*R - q*zA/2, R + d*n/2]
    if extra is not None:
        lines.append(extra(R, n, d, p, k))
    return min(lines)

def F(dl, p, l, h, extra=None, dy=None, H1=mpf(1)/4, z0=mpf('0.34')):
    x, y, b, k = geom(l, h)
    dy = h if dy is None else dy
    n = max(0, 1-y/dy)                         # B0.41: m=y/d, reflected length 1-min(m,1)
    a = (1+dl)/2; q = dl*p
    return (a + h*(z0-mpf(1)/6) - y/2 - l/2 + l*q
            + dy*(num(dl, p, k, n, extra, H1) - z0)), b

def inner_max(l, h, extra=None, H1=mpf(1)/4, p=mpf(1)/2, dy=None, grid=60):
    """max over delta in (0.02, kappa] of F (p fixed); golden refinement around grid max."""
    x, y, b, k = geom(l, h)
    f = lambda dl: F(dl, p, l, h, extra, dy, H1)[0]
    lo, hi = mpf('0.02'), k
    pts = [lo + (hi-lo)*i/grid for i in range(grid+1)]
    vals = [f(t) for t in pts]
    i = max(range(len(vals)), key=lambda j: vals[j])
    a_, c_ = pts[max(i-1, 0)], pts[min(i+1, grid)]
    g = (mp.sqrt(5)-1)/2
    for _ in range(110):                       # golden section (function is unimodal locally)
        x1 = c_ - g*(c_-a_); x2 = a_ + g*(c_-a_)
        if f(x1) > f(x2): c_ = x2
        else: a_ = x1
    t = (a_+c_)/2
    return f(t) - b, t

def best_h(l, extra=None, H1=mpf(1)/4, h0=mpf('0.7939475973881361'), w=mpf('0.004')):
    g = (mp.sqrt(5)-1)/2
    a_, c_ = h0-w, h0+w
    G = lambda hh: inner_max(l, hh, extra, H1)[0]
    for _ in range(90):
        x1 = c_ - g*(c_-a_); x2 = a_ + g*(c_-a_)
        if G(x1) < G(x2): c_ = x2
        else: a_ = x1
    hh = (a_+c_)/2
    return G(hh), hh

def theta(extra=None, H1=mpf(1)/4, lo=mpf('0.1665'), hi=mpf('0.1700'), it=60):
    """largest l with min_h max_delta F - b <= 0  -> theta = b."""
    for _ in range(it):
        mid = (lo+hi)/2
        if best_h(mid, extra, H1)[0] <= 0: lo = mid
        else: hi = mid
    G, hh = best_h(lo, extra, H1)
    x, y, b, k = geom(lo, hh)
    return dict(l=lo, h=hh, theta=b, G=G, delta=inner_max(lo, hh, extra, H1)[1])

def sixth(e):        # hypothetical sum_rows |S|^6 << U^e (no slots)
    return lambda R, n, d, p, k: e/6 + 5*R/6

def sixth_LS():      # conjectural sharp large sieve for S^3: e = max(1, 3n)
    return lambda R, n, d, p, k: max(1, 3*n)/6 + 5*R/6

if __name__ == "__main__":
    import sys
    rec = mpf(683505193)/781250000
    base = theta()
    print("baseline (robust)  theta =", mp.nstr(base['theta'], 16), " record - theta =",
          mp.nstr(rec-base['theta'], 4), " delta* =", mp.nstr(base['delta'], 12))
    cal = theta(H1=mpf(1)/4-mpf('1.3e-6'))
    print("H1 = 1/4 - 1.3e-6  theta =", mp.nstr(cal['theta'], 16), " gain =", mp.nstr(base['theta']-cal['theta'], 6))
    for e in ['1.20', '1.1801', '1.18', '1.17', '1.15', '1.1', '1.0']:
        r = theta(extra=sixth(mpf(e)))
        print(f"sixth moment e={e:7s} theta = {mp.nstr(r['theta'],16)}  gain = {mp.nstr(base['theta']-r['theta'],6)}"
              f"  delta*={mp.nstr(r['delta'],8)}")
    r = theta(extra=sixth_LS())
    print("sixth via sharp LS max(1,3n): theta =", mp.nstr(r['theta'], 16), " gain =", mp.nstr(base['theta']-r['theta'], 6))
