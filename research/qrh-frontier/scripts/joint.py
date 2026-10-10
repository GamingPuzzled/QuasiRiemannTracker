"""Joint large-values numerator bound using mixed MP moment (B3.40, two plain sums)."""
import numpy as np
from slack2 import R_k, DEF

def N_old(delta, p, k, d, y):
    R = R_k(delta, p, k, DEF); q = delta*p
    zA = max(0, 2*min(y/d, 1) - 1)/(6*k)
    return 0.25 + 0.75*R - q*zA/2

def N_new(delta, p, k, d, y, mixed=True, nt=101, nr=201, ns=201):
    q = delta*p
    mp = max(0.0, 1 - y/d)                 # reflected numerator length (U units)
    zA = (1 - 2*mp)/(6*k)                  # numerator 4th-moment capacity
    smax = delta*min(mp, 1-mp)/2           # pointwise (B0.15)
    best = np.inf
    sig = np.linspace(0, smax, ns)[None, :]
    for t in np.linspace(1.0, 1.5, nt):
        r = np.linspace(max(0, t-0.5), t, nr)[:, None]
        m = t - r
        inv = np.where(r <= 1, 1 - q - (delta - q)*r, (1 + 5*r)/6 - delta*r)
        AP = 1 - q/(3*k) - 2*(delta - q/(3*k))*m
        num4 = 1 - 4*sig - 2*q*zA
        cnt = np.minimum(np.minimum(inv, AP), num4)
        if mixed:
            zmix = np.maximum(0, 1 - mp - m)/(6*k)
            cnt = np.minimum(cnt, 1 - 2*sig - delta*m - 2*q*zmix)
        val = (np.maximum(cnt, -9) + sig).max()
        best = min(best, val)
    return best

if __name__ == "__main__":
    x, l = 0.373172481118674, 0.16712007850681
    y = 1-l-x; h = 1+l-x; b = 11/12-l/4; k = 2*b-1
    for delta in (0.3, 0.4, 0.469344858698143, 0.55, 0.7):
        o = N_old(delta, .5, k, h, y)
        n0 = N_new(delta, .5, k, h, y, mixed=False)
        n1 = N_new(delta, .5, k, h, y, mixed=True)
        print(f"δ={delta:.3f} old={o:.6f} joint(no mixed)={n0:.6f} joint+mixed={n1:.6f} gain={o-n1:+.5f}")
