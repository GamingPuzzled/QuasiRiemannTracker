import numpy as np
from slack2 import R_k, DEF
x, l = 0.373172481118674, 0.16712007850681
y = 1-l-x; h = 1+l-x; b = 11/12-l/4; k = 2*b-1
delta, p = 0.469344858698143, 0.5; q = delta*p; d = h
mp = 1-y/d; zA = (1-2*mp)/(6*k); smax = delta*min(mp,1-mp)/2
sig = np.linspace(0, smax, 401)[None, :]
print(f"m'={mp:.4f} smax={smax:.4f}")
for t in np.linspace(1.0, 1.5, 11):
    out = {}
    for name, (lo, hi) in {"short r<=1": (max(0,t-.5), min(1,t)), "long r>1": (1, t)}.items():
        if hi < lo: out[name] = None; continue
        r = np.linspace(lo, hi, 401)[:, None]; m = t - r
        inv = np.where(r <= 1, 1-q-(delta-q)*r, (1+5*r)/6-delta*r)
        AP = 1-q/(3*k)-2*(delta-q/(3*k))*m
        num4 = 1-4*sig-2*q*zA
        zmix = np.maximum(0, 1-mp-m)/(6*k)
        mix = 1-2*sig-delta*m-2*q*zmix
        a = (np.minimum(np.minimum(inv, AP), num4)+sig).max()
        bb = (np.minimum(np.minimum(np.minimum(inv, AP), num4), mix)+sig).max()
        cnt = np.minimum(inv, AP).max()
        out[name] = (round(cnt,4), round(a,5), round(bb,5))
    print(f"t={t:.2f}", out)
