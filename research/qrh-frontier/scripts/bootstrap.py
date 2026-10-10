import numpy as np
from slack2 import R_k, DEF
x0, l = 0.373172481118674, 0.16712007850681; b = 11/12-l/4; k = 2*b-1
deltas = np.linspace(0.0, k, 751)
def Rall(d):  # rows in bin, any amplitude class: worst p
    if d < 1e-9: return 1.0
    return max(R_k(d, p, k, DEF) for p in np.linspace(0, .5, 51))
RA = np.array([Rall(d) for d in deltas])
R05 = np.array([R_k(d, .5, k, DEF) if d>0 else 1.0 for d in deltas])
for x in (1.0, 1.1, 1.164, 1.25, 1.35, 1.5):
    E = (1+5*x)/6
    ea = (RA + x*deltas); ep = (R05 + x*deltas)
    print(f"x={x:.3f} E={E:.4f}  E'(all p)={ea.max():.4f}@δ={deltas[ea.argmax()]:.3f}  E'(p=.5)={ep.max():.4f}@δ={deltas[ep.argmax()]:.3f}")
