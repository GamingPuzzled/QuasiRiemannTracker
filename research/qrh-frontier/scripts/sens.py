import warnings; warnings.filterwarnings("ignore")
from slack2 import solve, DEF
base = solve()[0]
print(f"base {base:.9f}")
cases = [
 ("Mslack", 0.01, [0,0.005,0.01]),
 ("Mslack", 0.05, [0,0.01,0.02,0.03,0.04,0.05]),
 ("V0", 5/6+0.01, None),
 ("holder", 0.24, None),
 ("capdiv", 5.9, None),
 ("damp", 1.1, None),
 ("pmax", 0.51, None),
 ("probe_div", 12.5, None),
 ("z0", 0.30, None),
 ("supply", 0.17, None),
]
for k, v, sl in cases:
    r = solve({k: v}, s_list=sl or (0.0,))
    print(f"{k:10s}={v:.4f}  b0={r[0]:.9f}  gain={base-r[0]:+.2e}  x={r[1]:.4f} l={r[2]:.5f} s={r[3]}")
