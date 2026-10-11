"""Exact first-order sensitivity of theta to every exponent coefficient feeding the quartic."""
from mpmath import mp, mpf, diff
from model import solve, PAPER, F_of
mp.dps = 40
base = solve()
th0 = base['theta']
print("theta0 =", mp.nstr(th0, 20))
h = mpf('1e-12')
rows = []
for k_ in PAPER:
    up = solve({k_: PAPER[k_]+h})['theta']; dn = solve({k_: PAPER[k_]-h})['theta']
    der = (up-dn)/(2*h)
    need = mpf('1e-7')/abs(der) if der != 0 else mp.inf
    rows.append((k_, PAPER[k_], der, need))
print(f"{'param':6} {'paper':>10} {'dtheta/dc':>14} {'|dc| for 1e-7 gain':>20}  direction")
for k_, v, der, need in rows:
    print(f"{k_:6} {mp.nstr(v,6):>10} {mp.nstr(der,8):>14} {mp.nstr(need,4):>20}  {'decrease c' if der>0 else 'increase c'}")
# amplitude cap is active: dF/dp > 0 at binding
l, d, hh = base['l'], base['d'], base['h']
print("dF/dp at binding:", mp.nstr(diff(lambda p: F_of(d, p, l, hh, PAPER)[0], mpf(1)/2), 10))
