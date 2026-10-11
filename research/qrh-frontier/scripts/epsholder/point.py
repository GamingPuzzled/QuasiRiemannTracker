"""solverJ / B5: exact binding-point quantities for the numerator bound (B0.41) and the
quantitative tests of candidates (a) higher moments and (c) large-values inequalities.

Units: exponents of U = Z^h0 (bin d = h0).  sigma = log_U |S_u(U^n)|^2  (squared level).
"""
from mpmath import mp, mpf
from base import record, R as Rfun, tstar
mp.dps = 40
P = record()
l, k, d, h, y = P['l'], P['k'], P['d'], P['h'], P['y']
p = mpf(1)/2; q = d*p; c3 = q/(3*k)
R = Rfun(d, p, k); ts = tstar(d, p, k)
n = 1 - y/h                      # reflected numerator length (U-units), t = 0 dyad
zA = (1-2*n)/(6*k)               # MP capacity for the numerator 4th moment
L = l/h                          # slot supply
AI = lambda r: 1-q-(d-q)*r       # B0.38
AP = lambda m: 1-c3-2*(d-c3)*m
rI = (1-q-R)/(d-q)               # A_I(rI) = R
mP = (1-c3-R)/(2*(d-c3))         # A_P(mP) = R
sig = (1-R-2*q*zA)/2             # Holder-extremal squared level  (num4 line meets R)
pw = d*n                         # pointwise squared level (B0.15)
H4 = R + sig/2                   # = 1/4 + 3R/4 - q zA/2
out = lambda s, v: print(f"{s:44s} {mp.nstr(v, 16)}")
print("== binding point (record) ==")
for s, v in [("delta0", d), ("kappa0", k), ("h0", h), ("R", R), ("t*", ts), ("r_I (A_I=R)", rI),
             ("m_P (A_P=R)", mP), ("r_I + m_P - t*", rI+mP-ts), ("n = 1 - y/h (numerator dual length)", n),
             ("n - m_P", n-mP), ("z_A", zA), ("supply L", L), ("q = delta p", q),
             ("sigma* (Holder-extremal |S|^2 level)", sig), ("pointwise level delta*n", pw),
             ("pointwise - sigma* (|S|^2 units)", pw-sig), ("H4 = 1/4+3R/4-q zA/2", H4),
             ("pointwise numerator R + delta n/2", R+pw/2), ("gap pointwise-H4 (U units)", R+pw/2-H4),
             ("delta' = sigma*/n (effective level)", sig/n)]:
    out(s, v)

print("\n== (a) higher moments: break-even exponents e_{2k} for sum_rows |S|^{2k} << U^e ==")
print("   Holder value with 2k-th moment, no slots:  e/(2k) + (1-1/(2k)) R ; must be < H4")
for kk in (3, 4, 5):
    e_be = 2*kk*H4 - (2*kk-1)*R
    triv = 1 + (kk-2)*pw - 2*q*zA             # (max_B|S|^2)^{k-2} * MP moment (keeps slots)
    triv_val = triv/(2*kk) + (1-mpf(1)/(2*kk))*R
    lsieve = max(1, kk*n)                     # conjectural SHARP large sieve (U + K)
    bgl = max(1, kk*n, 2*(1+kk*n)/3)          # Blomer-Goldmakher-Louvel sextic sieve
    print(f" k={kk}: break-even e = {mp.nstr(e_be,12)} | pointwise-interp value - H4 = "
          f"{mp.nstr(triv_val-H4,6)} | sharp-LS e = {mp.nstr(lsieve,6)} (excess {mp.nstr(lsieve-e_be,4)})"
          f" | BGL e = {mp.nstr(bgl,6)}")
# Lindelof-on-average sixth moment (e=1) value
for kk in (3, 4):
    print(f" k={kk}: hypothetical Lindelof-on-average e=1: numerator {mp.nstr(1/(2*kk)+(1-mpf(1)/(2*kk))*R, 12)}"
          f" (gain {mp.nstr(H4-(1/(2*kk)+(1-mpf(1)/(2*kk))*R), 6)} U-units)")

print("\n== (c) large values for S^j (length U^{jn}) at the Holder level (squared, normalized N^{j delta'}) ==")
dp = sig/n
print(" Halasz-Montgomery needs normalized level > sqrt(length) off-diagonal, i.e. delta' > 1/2 even with")
print(f"   square-root off-diagonal bounds; here delta' = {mp.nstr(dp,10)} < 1/2  -> no information for any j.")
print(" Burgess (r=2) off-diagonal for k -> chi_k(u v^5), conductor <= U^2: N^{1/2} U^{3/8}:")
for j in (1, 2, 3):
    lenj = j*n
    off = lenj/2 + mpf(3)/8
    print(f"   j={j}: length {mp.nstr(lenj,5)}, off-diag {mp.nstr(off,5)} vs trivial {mp.nstr(lenj,5)}; level {mp.nstr(j*sig,5)}")
print(" Huxley-type count (heuristic transfer, family size U in place of T): N^{x(1-d')} + U N^{x(1-3d')}:")
for j in (1, 2, 3, 4):
    x = j*n
    cnt = max(x*(1-dp), 1+x*(1-3*dp))
    print(f"   j={j}: x={mp.nstr(x,6)} count exponent {mp.nstr(cnt,8)}  vs R={mp.nstr(R,8)}  ({'BETTER' if cnt<R else 'worse'})")
xs = 1/(2*dp)
print(f"   continuous optimum x={mp.nstr(xs,8)} (j={mp.nstr(xs/n,6)}, not an integer power): {mp.nstr(xs*(1-dp),8)}")
