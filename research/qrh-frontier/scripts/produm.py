"""Row count with an extra product-witness constraint PM(x)=e(x)-δx, e(x)=max(1,(1+5x)/6)."""
import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import solve, DEF
orig = slack2.R_k

def R_new(delta, p, k, c, use_pm=True, nt=61, nr=121, nm=61):
    q = delta*p
    best = np.inf
    for t in np.linspace(c.get("tmin", 1.0), 1.5, nt):
        r = np.linspace(0, t, nr)[:, None]
        m = np.linspace(0, 0.5, nm)[None, :]
        ok = (r + m >= t - 1e-12)
        inv = np.where(r <= 1, 1 - q - (delta - q)*r, (1 + 5*r)/6 - delta*r)
        AP = 1 - q/(3*k) - 2*(delta - q/(3*k))*m
        cnt = np.minimum(inv, AP)
        if use_pm:
            x = r + m
            cnt = np.minimum(cnt, np.maximum(1, (1 + 5*x)/6) - delta*x)
        val = np.where(ok, cnt, -9).max()
        best = min(best, val)
    return best

if __name__ == "__main__":
    x, l = 0.373172481118674, 0.16712007850681; b = 11/12-l/4; k = 2*b-1
    for d in (0.3, 0.469344858698143, 0.6):
        print(f"δ={d:.3f} paper R={orig(d,.5,k,DEF):.5f} grid-noPM={R_new(d,.5,k,DEF,False):.5f} withPM={R_new(d,.5,k,DEF):.5f}")
    base = solve()[0]
    slack2.R_k = lambda d, p, k, c: R_new(d, p, k, c, True, nt=31, nr=61, nm=31)
    # coarse solve (slow); reduce grids in worst()
    import slack2 as s
    r = solve()
    print("θ with product UM:", r[0], "gain", base - r[0])
