"""Vectorized slack model of cycle25 B0.3-B0.7. Constants perturbable via dict."""
import numpy as np
from scipy.optimize import minimize_scalar

DEF = dict(V0=5/6, holder=0.25, capdiv=6.0, pmax=0.5, z0=0.34, Mslack=0.0,
           probe_div=12.0, supply=4/21, damp=1.0)

def R_k(delta, p, k, c):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k))
    V = np.maximum(c["V0"]-delta, 0)
    return 1 - delta + delta*C*V/(2*(V+delta*C))

def worst(beta, x, y, l, c, nd=400, npp=101):
    k = 2*beta-1; h = 1+l-x
    delta = np.linspace(0.02, k, nd)[:, None]
    p = np.linspace(0, c["pmax"], npp)[None, :]
    q = delta*p; a = (1+delta)/2
    R = R_k(delta, p, k, c)
    out = []
    for d in (0.005, h, min(h, y)):
        cap = max(0.0, min(d, 2*y-d))
        Fd = (a + h*(c["z0"]-1/6) - y/2 - l/2 + l*q
              + d*(c["holder"] + (1-c["holder"])*R - c["z0"])
              - c["damp"]*q/(2*c["capdiv"]*k)*2*cap/2*2/2)  # = q*cap/(12k) at defaults
        out.append(Fd.max())
    floor = 0.51 + h*(c["z0"]-1/6) - y/2 - l/2 + l*0.01 + h*(1-c["z0"])
    out.append(floor)
    return max(out) - beta

def probe(x, y, l, c):
    M = x+y
    return x/2 + (y-x)/c["probe_div"] + max(0, 1+3*l-2*M)/8 + 5/6 - x/3 - l/6

def feasible(x, y, l, c):
    h = 1+l-x
    return (l < x <= y and y-2*x+l < 0 and x+y+l <= 1+c["Mslack"]+1e-12
            and l/h > c["supply"] and y < 4*l)

def best_for(x, s, c):
    """Max l (=> min b0) for given x and M-slack s; returns b0 or inf."""
    def ok(l):
        y = 1-l-x+s
        if not feasible(x, y, l, c): return None
        b = probe(x, y, l, c)
        return b if worst(b+1e-12, x, y, l, c) <= 0 else None
    best = np.inf
    ls = np.linspace(0.05, 0.4, 141)
    goods = [l for l in ls if ok(l) is not None]
    if not goods: return np.inf, None
    lo = max(goods); hi = lo+0.0025
    for _ in range(50):
        mid = (lo+hi)/2
        (lo, hi) = (mid, hi) if ok(mid) is not None else (lo, mid)
    return ok(lo), lo

def solve(c=None, s_list=(0.0,)):
    c = {**DEF, **(c or {})}
    best = (np.inf, None, None, None)
    for s in s_list:
        r = minimize_scalar(lambda x: best_for(x, s, c)[0], bounds=(0.2, 0.6), method="bounded",
                            options=dict(xatol=1e-10))
        b, l = best_for(r.x, s, c)
        if b < best[0]: best = (b, r.x, l, s)
    return best

if __name__ == "__main__":
    import sys
    print("baseline", solve())
