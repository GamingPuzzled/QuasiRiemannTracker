import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import solve, DEF

def R_general(delta, p, k, c):
    aI, aP, V0 = c.get("aI", 0.5), c.get("aP", 1.0), c["V0"]
    q = delta*p
    alI, beI = 1-2*q*aI, delta-2*q*aI
    cP = q*aP/(3*k)
    alP, beP = 1-cP, 2*delta-2*cP
    # short(t) = alI - beI*r*(t), r* = (alI-alP+beP t)/(beI+beP)
    s1 = beI*beP/(beI+beP); s0 = alI - beI*(alI-alP)/(beI+beP)
    V = np.maximum(V0-delta, 0)
    # long(t) = 1-delta+V(t-1); balance s0 - s1 t = 1-delta+V(t-1)
    t = (s0 - 1 + delta + V)/(s1 + V)
    t = np.clip(t, 1, 1.5)
    return np.maximum(s0 - s1*t, 1-delta+V*(t-1))

# verify against paper formula
k = 0.7498
for d, p in [(0.469, 0.5), (0.3, 0.2), (0.6, 0.5)]:
    print(d, p, R_general(d, p, k, DEF), slack2.R_k(d, p, k, DEF))
base = solve()[0]
slack2.R_k = R_general
print("base(general R)", solve()[0])
for key, val in [("aI", 0.51), ("aP", 1.02), ("V0", 5/6-0.01)]:
    b = solve({key: val})[0]
    dv = val - {"aI":0.5,"aP":1.0,"V0":5/6}[key]
    print(f"{key}: db0={b-base:+.3e} per-unit={(b-base)/dv:+.4f}")
