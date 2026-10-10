import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import solve
orig = slack2.R_k
base = solve()[0]
def excess(delta, p, k, c):  # R - (1-delta)
    return orig(delta, p, k, c) - (1-delta)
for lam in (1.0, 0.9, 0.75, 0.5, 0.0):
    slack2.R_k = lambda d, p, k, c, lam=lam: (1-d) + lam*excess(d, p, k, c)
    r = solve()
    print(f"R = 1-δ + {lam:.2f}*excess: θ={r[0]:.9f} gain={base-r[0]:+.3e}")
slack2.R_k = orig
x,l=0.373172481118674,0.16712007850681; b=11/12-l/4; k=2*b-1; d0=0.469344858698143
Rv = orig(d0,.5,k,slack2.DEF); a=(1+d0)/2
print(f"binding: a=σ={a:.4f}  R={Rv:.4f}  density exponent R/(1-σ)={Rv/(1-a):.4f} (DH = 2)")
