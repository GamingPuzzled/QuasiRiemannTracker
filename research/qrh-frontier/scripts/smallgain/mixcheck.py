"""(c) joint numerator/row-class bound with the mixed MP moment (B3.40, n1=m, n2=m') at binding."""
from mpmath import mp, mpf
from base import record, R, tstar, Ck
mp.dps = 40
P = record(); d, k, h, y = P['d'], P['k'], P['h'], P['y']; p = mpf(1)/2; q = d*p
Rv = R(d, p, k); t = tstar(d, p, k)
c3 = q/(3*k)
# short crossing r*(t): A_I(r) = A_P(t-r)
AI = lambda r: 1-q-(d-q)*r
AP = lambda m: 1-c3-2*(d-c3)*m
r = (1-q - (1-c3) + 2*(d-c3)*t)/((d-q)+2*(d-c3)); m = t-r
print("t*=",mp.nstr(t,12),"r*=",mp.nstr(r,12),"m*=",mp.nstr(m,12),"A_I=",mp.nstr(AI(r),15),"A_P=",mp.nstr(AP(m),15),"R=",mp.nstr(Rv,15))
mprime = 1 - y/h; zA = (1-2*mprime)/(6*k)
sbal = (1-Rv-2*q*zA)/4; smax = d*min(mprime,1-mprime)/2
print("m'=",mp.nstr(mprime,12),"sigma_bal=",mp.nstr(sbal,12),"sigma_pointwise=",mp.nstr(smax,12))
zmix = (1-m-mprime)/(6*k)
mix = 1 - d*m - 2*sbal - 2*q*zmix
print("mixed count at balance:", mp.nstr(mix,15), " minus R:", mp.nstr(mix-Rv,6))
# also mixed with the *largest* admissible numerator m' (h0 shift t=0 is worst)
