"""Global break-even of the joint/4th-moment slope lam (full minimax re-solve, slack2 grid model)."""
import sys, warnings; warnings.filterwarnings("ignore"); sys.path.insert(0, '.'); sys.path.insert(0, '..')
import slack2
from price_theta_minimax import make_R, R_orig
base = slack2.solve()[0]
for lam in [0.6, 0.599, 0.598, 0.597, 0.596, 0.595, 0.594, 0.593, 0.592, 0.591]:
    slack2.R_k = make_R(lam); r = slack2.solve(); slack2.R_k = R_orig
    print(f"lam={lam:.4f} theta={r[0]:.10f} gain={base-r[0]:+.3e} x={r[1]:.5f}", flush=True)
