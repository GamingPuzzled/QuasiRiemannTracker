"""(c) Exponent-level adversarial profile at the binding point.

A 'virtual bin' of U^R rows, each with
    log_U|M_u(U^r)|^2 = mu(r),   log_U|S_u(U^nu)|^2 = s(nu)  (one presentation, = numerator presentation),
    |Q_{u,I}|^2 = U^{2 q z_I} for every slot subcollection of capacity z_I <= L,
all other rows of the family typical (every normalized polynomial of size U^{o(1)}).
We verify that EVERY inequality the paper makes available for rows (UM/B0.9, MI, MP/MP0 incl. mixed
pairs at all lengths, pointwise B0.15/B0.16, amplitude cap, Lemma 1W witnesses for every t0) holds,
while  sum_{bin} |S_u(U^n)| = U^{R + s(n)/2} = U^{H4}.  Hence no level-set/Holder/interpolation
argument built from these inequalities can lower the numerator exponent below H4 (B0.41) at the
binding point: B5 needs an inequality outside this list.
"""
import numpy as np
from mpmath import mp, mpf
from base import record, R as Rfun, tstar
mp.dps = 30
P = record()
l, k, d, h, y = [float(P[s]) for s in ('l', 'k', 'd', 'h', 'y')]
p = 0.5; q = d*p; L = l/h
R = float(Rfun(P['d'], mpf(1)/2, P['k'])); ts = float(tstar(P['d'], mpf(1)/2, P['k']))
n = 1-y/h
E = lambda r: max(1.0, (1+5*r)/6)
zI = lambda r: max(0.0, min((1-r)/2, (3-2*r)/8, L))            # MI capacity (B0.4) & supply
zP = lambda a, b: max(0.0, min(L, (1-a-b)/(6*k)))              # MP capacity, reflected lengths a,b<=1/2
def mu(r):                                                      # inverse profile
    if r <= 1: return min(d*r, 1-R-2*q*zI(r))
    return min(d*r, E(r)-R)
def s(nu):                                                      # plain profile (reflection-symmetric)
    a = min(nu, 1-nu)
    return min(d*a, (1-R-2*q*zP(a, a))/2)
tol = 1e-12; worst = {}
def chk(name, val):                                             # val must be <= 0
    worst[name] = max(worst.get(name, -9), val)
A = np.linspace(0, 0.5, 501)
# C1  MP / MP0 / mixed MP, every pair of lengths, every admissible capacity
for a in A:
    for b in A:
        chk("C1 MP pairs (R+s1+s2+2qz<=1)", R + s(a) + s(b) + 2*q*zP(a, b) - 1)
# C2  MI with slots, every r<=1 and every admissible capacity (max capacity is the worst)
for r in np.linspace(0, 1, 1001):
    chk("C2 MI (R+mu+2qz<=1)", R + mu(r) + 2*q*zI(r) - 1)
# C3  sixth-power-free inverse moment B0.9 for 1<r<=3/2 (UM for r<1 is C2 with z=0)
for r in np.linspace(1, 1.5, 501):
    chk("C3 B0.9 (R+mu<=E(r))", R + mu(r) - E(r))
# C4  pointwise B0.15 / B0.16 (slot amplitude q per unit capacity <= delta/2)
for r in np.linspace(0, 1.5, 1501): chk("C4a pointwise M", mu(r) - d*r)
for a in np.linspace(0, 1, 1001): chk("C4b pointwise S", s(a) - d*min(a, 1-a))
chk("C4c amplitude cap q<=delta/2", q - d/2)
# C5  Lemma 1W witnesses for every t0 in [1,3/2]: exist r<=t0, m in [t0-r,1/2], mu(r)=d r, s(m)=d m
rs = np.linspace(0, 1.5, 3001); ms = np.linspace(0, 0.5, 1001)
fullM = np.array([abs(mu(r)-d*r) < tol for r in rs]); fullS = np.array([abs(s(m)-d*m) < tol for m in ms])
for t0 in np.linspace(1, 1.5, 501):
    ok = False
    for r, fm in zip(rs, fullM):
        if not fm or r > t0 + 1e-12: continue
        need = max(0.0, t0-r)
        if need <= 0.5 and fullS[ms >= need - 1e-12].any(): ok = True; break
    chk("C5 Lemma 1W witness exists (0 = yes)", 0.0 if ok else 1.0)
# C6  the numerator
sig = (1-R-2*q*zP(n, n))/2; H4 = 0.25 + 0.75*R - q*zP(n, n)/2
print(f"R={R:.15f} n={n:.12f} s(n)={s(n):.15f} sigma*={sig:.15f}")
print(f"virtual sum_bin |S| exponent R+s(n)/2 = {R+s(n)/2:.15f};  H4 (B0.41) = {H4:.15f};  diff = {R+s(n)/2-H4:.1e}")
for kname, v in worst.items(): print(f"{kname:42s} max violation = {v:+.3e}  {'OK' if v <= tol else 'VIOLATED'}")
# witness geometry
rI = (1-q-R)/(d-q); mP = (1-q/(3*k)-R)/(2*(d-q/(3*k)))
print(f"witness thresholds r_I={rI:.10f} m_P={mP:.10f}  r_I+m_P={rI+mP:.10f}  t*={ts:.10f}  n-m_P={n-mP:.3e}")
# what the virtual bin does to a hypothetical sixth moment
print(f"virtual bin sixth moment exponent  R+3s(n)+2q zA = {R+3*s(n)+2*q*zP(n,n):.12f} = 1 + delta'*n with delta'={s(n)/n:.10f}"
      f" (trivial bound 1+delta*n, delta={d:.10f}; required saving eta > {d - s(n)/n:.10f})")
