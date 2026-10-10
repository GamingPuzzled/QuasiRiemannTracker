import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import solve, DEF
base = solve()[0]
# E-slope s replaces 5/6 in V = s - delta  (long envelope 1-δ+(s-δ)(t-1))
for s in (5/6, 0.82, 0.80, 0.75, 0.70, 0.60, 0.5):
    print(f"lambda/s={s:.3f}: θ={solve({'V0': s})[0]:.9f} gain={base-solve({'V0': s})[0]:+.2e}")
