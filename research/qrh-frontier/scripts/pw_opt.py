import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import R_k, solve
base = solve()[0]
def worst(beta, x, y, l, c, nd=400, npp=101):
    k = 2*beta-1; h = 1+l-x
    delta = np.linspace(0.02, k, nd)[:, None]; p = np.linspace(0, c["pmax"], npp)[None, :]
    q = delta*p; a = (1+delta)/2; R = R_k(delta, p, k, c)
    out = []
    for d in (0.005, h, min(h, y)):
        cap = max(0.0, min(d, 2*y-d)); m = y/d
        hold = d*(0.25+0.75*R) - q*cap/(12*k)
        pw = d*(R + delta*max(0, min(1-m, m))/2)
        Fd = a + h*(c["z0"]-1/6) - y/2 - l/2 + l*q - d*c["z0"] + np.minimum(hold, pw)
        out.append(Fd.max())
    out.append(0.51 + h*(c["z0"]-1/6) - y/2 - l/2 + l*0.01 + h*(1-c["z0"]))
    return max(out) - beta
slack2.worst = worst
r = solve()
print(f"with min(holder,pointwise): b0={r[0]:.10f} gain={base-r[0]:+.2e} x={r[1]:.5f}")
