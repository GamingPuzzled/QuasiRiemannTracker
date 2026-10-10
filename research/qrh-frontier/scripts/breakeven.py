import warnings; warnings.filterwarnings("ignore")
import slack2, numpy as np
from slack2 import solve
base = solve()[0]
for cst in (0.0, 0.25, 0.5, 0.75, 0.9):
    def probe(x, y, l, c, cst=cst):
        M = x+y; s = max(0, M+l-1)
        return x/2 + (y-x)/12 + cst*s/2 + 5/6 - x/3 - l/6
    slack2.probe = probe
    r = solve({"Mslack": 0.02}, s_list=list(np.linspace(0, 0.02, 5)))
    print(f"energy cost {cst}*s: b0={r[0]:.9f} gain={base-r[0]:+.2e} s={r[3]:.3f}")
