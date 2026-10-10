import warnings; warnings.filterwarnings("ignore")
import slack2
from slack2 import solve
from Rgen import R_general
slack2.R_k = R_general
base = solve()[0]
for c in (6, 5.5, 5, 4, 3, 2):
    r = solve({"capdiv": c, "aP": 6/c})
    print(f"MP slot cost {c}κ: θ={r[0]:.9f} gain={base-r[0]:+.2e} x={r[1]:.4f} ℓ={r[2]:.5f}")
