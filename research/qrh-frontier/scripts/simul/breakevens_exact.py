"""Exact (mpmath) break-even slopes at the record's binding bin (delta0, p=1/2, kappa0)."""
import sys; sys.path.insert(0, '.')
from mpmath import mp, mpf, nstr
from base import record, R, tstar
mp.dps = 40
P = record(); d, k = P['d'], P['k']; p = mpf(1)/2; q = d*p
R0 = R(d, p, k); ts = tstar(d, p, k)
A = 1-q; s = d-q; B = 1-q/(3*k); w = d-q/(3*k); V = mpf(5)/6-d
rR = (A-R0)/s; mR = (B-R0)/(2*w); exc = R0-(1-d)
print("delta0 =", nstr(d, 20), " kappa0 =", nstr(k, 20))
print("R0 =", nstr(R0, 20), " t* =", nstr(ts, 20), " r_R =", nstr(rR, 20), " m_R =", nstr(mR, 20))
print("check r_R + m_R - t* =", nstr(rR+mR-ts, 5))
print("excess R0-(1-delta) =", nstr(exc, 20))
lam0 = d + exc/(2*rR-1)
lam1 = d + exc/((1-R0)/d + 2*mR - ts + rR - 1)
lam2 = d + exc/(ts-1)
print("lam_be P0 (no Poisson-range joint)      =", nstr(lam0, 20))
print("lam_be P1 (UM2 in Poisson range)        =", nstr(lam1, 20))
print("lam_be P2 (MI2 in Poisson range)        =", nstr(lam2, 20), " (5/6 =", nstr(mpf(5)/6, 20), ")")
print("needed exponent at the P2 equilibrium lengths: x = t* =", nstr(ts, 12),
      "-> e(x) < E(t*) = (1+5t*)/6 =", nstr((1+5*ts)/6, 12))
print("binding lengths: short factor ~ m_R =", nstr(mR, 6), ", long factor ~ r_R =", nstr(rR, 6))
