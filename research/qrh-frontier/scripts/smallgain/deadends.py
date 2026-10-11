"""Exact values for the additional variants examined (all at the record binding point)."""
from mpmath import mp, mpf
from base import record, R, tstar
mp.dps = 30
P = record(); d, k, h, y, l = P['d'], P['k'], P['h'], P['y'], P['l']
p = mpf(1)/2; q = d*p; c3 = q/(3*k); R0 = R(d, p, k); ts = tstar(d, p, k)
r = (1-q-(1-c3)+2*(d-c3)*ts)/((d-q)+2*(d-c3)); ms = ts-r; mpr = 1-y/h
zA = (1-2*mpr)/(6*k); holder = mpf(1)/4+3*R0/4-q*zA/2
# 1. sixth-power lift + MI slots for long rows: count (1+5t+10z)/6 - d t - 2 q z
print("lifted MI, long rows: dz-coefficient 5/3-2q =", mp.nstr(mpf(5)/3-2*q, 8), "(>0: slots hurt)")
# 2. lifted MI for short rows: A_I(r) + lp*(5-6q)
print("lifted MI, short rows: d(log P)-coefficient 5-6q =", mp.nstr(5-6*q, 8))
# 3. MP with n2=0 for the plain witness vs pure MP
print("MP(n1=m*,n2=0) count - A_P(m*) =", mp.nstr(ms*(d-c3), 8))
# 4. numerator by 2nd moment with slots (MP n2=0, cap (1-m')/(6k))
z2 = (1-mpr)/(6*k)
print("numerator 2nd-moment exponent", mp.nstr(R0/2+mpf(1)/2-q*z2, 10), "vs Holder", mp.nstr(holder, 10))
# 5. numerator via witness-weighted Cauchy-Schwarz with S(m*): same value (identity)
zmix = (1-ms-mpr)/(6*k)
print("witness-weighted CS exponent", mp.nstr(-d*ms/2+R0/2+mpf(1)/2-q*zmix, 15), "vs Holder", mp.nstr(holder, 15))
# 6. direct full-width MI for long rows (no lift): count r(1-d)
print("full-width MI long count", mp.nstr(ts*(1-d), 8), "vs R", mp.nstr(R0, 8))
# 7. probe energy slack (B0.26) unusable: column branch energy M+s past M+l=1
print("energy margin (1-5l)/4 =", mp.nstr((1-5*l)/4, 6), "but hypothesis M+2z<=N is the active one")
