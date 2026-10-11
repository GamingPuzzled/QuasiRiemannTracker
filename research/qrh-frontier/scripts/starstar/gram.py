"""Copy of ../gram.py (c -> theta), extended for solverK-starstar.

Probe Gram term: log_Z of 1 + P_a^c + P_a^2/Y, optionally with a fixed Z-saving eta
on the exceptional term (gram = max(0, c*Pa - eta, 2Pa - y)).
Run:  python3 gram.py
"""
import os, sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import slack2
from slack2 import solve

RECORD = 683505193/781250000

def mk(cexp, eta=0.0):
    def probe(x, y, l, c):
        M = x+y; Pa = y-x  # P_a = Y/X = Y^2/Q up to q_{b*}
        gram = max(0, cexp*Pa - eta, 2*Pa - y)
        E = M + max(0, 1+3*l-2*M)/4
        return -M/2 + (x + gram)/2 + E/2 + 5/6 - x/3 - l/6
    return probe

def theta(cexp, eta=0.0):
    slack2.probe = mk(cexp, eta)
    r = solve()
    return r[0], r[1], r[2]

if __name__ == "__main__":
    base = theta(1/6)[0]
    print(f"record theta = {RECORD:.11f}; model b0(c=1/6) = {base:.11f}")
    print("--- theta as a function of the P_a-exponent c (geometry re-optimized) ---")
    for name, cexp in [("1/6 (record)", 1/6), ("1/6-1e-3", 1/6-1e-3), ("0.16", .16),
                       ("1/7 (perfect pole + L-avg)", 1/7), ("5/42 (litG NC-th+NC-K unif.)", 5/42),
                       ("1/8", 1/8), ("1/9", 1/9), ("1/12", 1/12), ("0", 0.0)]:
        b, x, l = theta(cexp)
        y = 1-l-x
        print(f"c={cexp:.5f} [{name:28s}] b0={b:.9f} gain={base-b:+.3e} x={x:.4f} y={y:.4f} Pa=Z^{y-x:.4f}")
    print("--- theta with a fixed Z-power saving eta on the exceptional term (c=1/6) ---")
    for eta in (1e-4, 5e-4, 1e-3, 3e-3, 7.2e-3, 1.44e-2):
        b, x, l = theta(1/6, eta)
        y = 1-l-x
        print(f"eta={eta:.4g}: b0={b:.9f} gain={base-b:+.3e} x={x:.4f} y={y:.4f} Pa=Z^{y-x:.4f}")
