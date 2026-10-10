import warnings; warnings.filterwarnings("ignore")
import slack2
from slack2 import solve
base = solve()[0]
def mk(cexp):
    def probe(x, y, l, c):
        M = x+y; Pa = y-x  # P_a = Y/X
        gram = max(0, cexp*Pa, 2*Pa - y)       # log_Z of 1 + P_a^c + P_a^2/Y
        # Q^{-1/2} (Q/Y * gram)^{1/2} Z^{E/2}, E = M (+ (1+3l-2M)_+/4)
        E = M + max(0, 1+3*l-2*M)/4
        return -M/2 + (x + gram)/2 + E/2 + 5/6 - x/3 - l/6
    return probe
for cexp in (1/6, 1/8, 1/12, 1/24, 0.0):
    slack2.probe = mk(cexp)
    r = solve()
    print(f"c={cexp:.4f}: b0={r[0]:.9f} gain={base-r[0]:+.3e} x={r[1]:.4f} l={r[2]:.5f} y={1-r[1]-r[2]:.4f}")
