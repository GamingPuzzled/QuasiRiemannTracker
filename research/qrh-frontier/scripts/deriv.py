import warnings; warnings.filterwarnings("ignore")
from slack2 import solve, DEF
base = solve()[0]
for k, h in [("V0",-0.005),("holder",-0.005),("capdiv",-0.05),("damp",0.05),("pmax",-0.005),("probe_div",0.25)]:
    b = solve({k: DEF[k]+h})[0]
    print(f"{k:10s} step {h:+.3f}: db0={b-base:+.3e}  per-unit={(b-base)/h:+.4f}")
