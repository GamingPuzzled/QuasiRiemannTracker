"""Cross-check lib.py (vectorized dlog/Gauss) against the slow reference code of ../theta_id.py, and basic laws."""
import numpy as np
from lib import prime_elements, symbol_table, gauss_composite, products, Z6, W3, LAM
from eis import sixth_root_table

pel = prime_elements(3000)
targets = [(5, 0), (2, 0), (7, 3), (-1, 0), (11, 6), (1, 2)] + [P["pi"] for P in pel[:40]]
TA = [t[0] for t in targets]; TB = [t[1] for t in targets]
imax = 120
E, Zr, G = symbol_table(pel, imax, TA, TB, POWS=(1, 2), convs=('lam', 'chk'))
err_sym = 0; err_g = 0
for i in range(imax):
    P = pel[i]; p, r = P["p"], P["r"]; tab = sixth_root_table(p, r); e = (p-1)//6
    for t, (x, y) in enumerate(targets):
        res = (x % p + (y % p)*r) % p
        ref = -1 if res == 0 else tab[pow(res, e, p)]
        got = -1 if Zr[i, t] else E[i, t]
        err_sym += ref != got
    piv = P["pi"][0] + P["pi"][1]*W3
    xs = np.arange(1, p); ch = np.array([tab[pow(int(v), e, p)] for v in xs])
    for cv in ('lam', 'chk'):
        tr = (2*np.real(LAM*piv)/(3*p)) if cv == 'lam' else (2*np.real(np.conj(piv))/p)
        for POW in (1, 2):
            ref = np.sum(Z6[(POW*ch) % 6]*np.exp(2j*np.pi*xs*tr))/np.sqrt(p)
            err_g = max(err_g, abs(ref - G[(POW, cv)][i]))
print("symbol mismatches:", err_sym, " max Gauss-sum error:", err_g)
print("|g| for first primes (should be 1):", np.round(np.abs(G[(1, 'lam')][:6]), 10))
# multiplicativity check: chi_pi(xy) = chi_pi(x) chi_pi(y) for targets (5,0),(2,0) -> 10
E2, Z2, _ = symbol_table(pel, imax, [10], [0], POWS=(1,), convs=('lam',))
print("multiplicativity (5*2):", int(np.sum((E[:imax, 0].astype(int) + E[:imax, 1]) % 6 != E2[:imax, 0])), "mismatches")
# cubic Gauss sums g~_3(pi)^3 should be ~ pi-related: |g3|=1 and g(chi^2)^3 / (p * pi-ish) has finite set of values
print("done")
