"""Independent numeric audit of cycle25 quartic endpoint claims (T1.1-T1.5, B0.45-B0.49, B0.56)."""
from sympy import Rational as Q, symbols, expand, simplify, nsolve, Poly, N, sqrt
import mpmath as mp
mp.mp.dps = 60
L = symbols('l')
p = 927*L**4 - 3135*L**3 + 2433*L**2 + 275*L - 100
ok = lambda name, c: print(("PASS " if c else "FAIL ") + name)

# T1.1 root
lo, hi = Q(16712007850680, 10**14), Q(16712007850682, 10**14)
ok("p(lo)<0<p(hi)", p.subs(L, lo) < 0 < p.subs(L, hi))
ok("p(1/6)<0<p(1/5)", p.subs(L, Q(1,6)) < 0 < p.subs(L, Q(1,5)))
dp = 3708*L**3 - 9405*L**2 + 4866*L + 275
ok("p' formula", expand(p.diff(L) - dp) == 0)
ok("p' lower bound value -9405/25+4866/6+275 = 3549/5", Q(-9405,25)+Q(4866,6)+275 == Q(3549,5))
# check bound is valid: 3708x^3>=0, -9405x^2>=-9405/25, 4866x>=4866/6 on [1/6,1/5]
mins = min(dp.subs(L, Q(1,6)+ (Q(1,5)-Q(1,6))*Q(i,1000)) for i in range(1001))
ok(f"p' min on grid {float(mins):.2f} > 0", mins > 0)
l0 = mp.findroot(lambda x: 927*x**4-3135*x**3+2433*x**2+275*x-100, 0.167)
print("l0 =", mp.nstr(l0, 25))

# T1.5
b0 = mp.mpf(11)/12 - l0/4
th = mp.mpf(683505193)/781250000
ok(f"0<theta-b0<4e-14 (diff={mp.nstr(th-b0,5)})", 0 < th-b0 < 4e-14)
ok("theta<10499/12000", Q(683505193,781250000) < Q(10499,12000))
ok("10499/12000 < 874957019420098946128604623/10^27", Q(10499,12000) < Q(874957019420098946128604623,10**27))
# rigorous: theta - b0 using root interval
b0_hi = Q(11,12) - lo/4; b0_lo = Q(11,12) - hi/4
ok("theta > b0 via interval", Q(683505193,781250000) > b0_hi)
ok("theta - b0 < 4e-14 via interval", Q(683505193,781250000) - b0_lo < Q(4,10**14))

# B0.45 geometry
k0 = (5-3*l0)/6
d0 = (25-60*l0+27*l0**2)/(6*(5+4*l0-3*l0**2))
C0 = (6*k0-1)/(15*k0-2); V0 = mp.mpf(5)/6 - d0
Rp = -1 + C0*(V0**2 - d0**2*C0)/(2*(V0+d0*C0)**2)
h0 = -((1+l0)/2 + l0/(6*k0))/(3*Rp/4 - 1/(24*k0))
x0 = 1+l0-h0; y0 = h0-2*l0
print("delta0,h0,x0,y0 =", [mp.nstr(v, 15) for v in (d0, h0, x0, y0)])
ok("matches paper (.469344858698143,.793947597388136,.373172481118674,.459707440374516)",
   all(abs(a-b) < 1e-14 for a, b in zip((d0,h0,x0,y0),(0.469344858698143,0.793947597388136,0.373172481118674,0.459707440374516))))
ok("l0<x0<y0<2x0", l0 < x0 < y0 < 2*x0)
ok("y0<1/2<h0<2y0", y0 < 0.5 < h0 < 2*y0)
ok("kappa0 = 2 b0 - 1", abs(k0 - (2*b0-1)) < 1e-40)

# R_kappa and F4 at binding
def C(k, pp): return 2*(1-pp)*(1-pp/(3*k))/(3-pp-2*pp/(3*k))
def R(k, d, pp):
    V = mp.mpf(5)/6 - d; c = C(k, pp)
    return 1 - d + d*c*V/(2*(V+d*c))
def F4(k, d, pp, l=l0, h=h0):
    return (1+l)/2 + d/2 + h*(3*R(k,d,pp)/4 - mp.mpf(5)/12) + d*pp*(l - (h-4*l)/(12*k))
ok(f"C0 equals C_kappa(1/2) ({mp.nstr(C0,10)} vs {mp.nstr(C(k0,.5),10)})", abs(C0-C(k0, mp.mpf(1)/2)) < 1e-40)
f = F4(k0, d0, mp.mpf(1)/2)
ok(f"F4 at binding = b0 (diff {mp.nstr(f-b0,5)})", abs(f-b0) < 1e-30)
# B0.48
e1 = (324*k0**2-3)*d0**2 + (-756*k0**2+54*k0+5)*d0 + 300*k0**2 - 40*k0
e2 = d0 - k0*(5-9*l0)/(6*k0+(6*k0+2)*l0)
ok(f"B0.48 eq1 ({mp.nstr(e1,3)})", abs(e1) < 1e-30)
ok(f"B0.48 eq2 ({mp.nstr(e2,3)})", abs(e2) < 1e-30)
ok("R(delta0,1/2) = 5/9+delta0/(18k0)", abs(R(k0,d0,mp.mpf(1)/2) - (mp.mpf(5)/9 + d0/(18*k0))) < 1e-30)
# global max of F4 over delta in (0,k0], p in [0,1/2]
best = max((F4(k0, mp.mpf(i)/4000*k0, mp.mpf(j)/100), i, j) for i in range(1, 4001) for j in range(0, 51, 5))
ok(f"grid max F4 <= b0 (max-b0={mp.nstr(best[0]-b0,5)} at i={best[1]},j={best[2]})", best[0] <= b0 + 1e-12)

# margins
other = mp.mpf(11)/12 - 49*h0/150 + 11*l0/12
ok(f"other endpoint {mp.nstr(other,6)} < b0-.06 ({mp.nstr(b0-0.06,6)})", other < b0 - 0.06)
floor = 0.51 + h0/3 + 0.51*l0
ok(f"floor {mp.nstr(floor,6)} < b0-.01", floor < b0 - 0.01)
z0 = mp.mpf(17)/50; dmin = mp.mpf(1)/200
b56 = h0*(z0-mp.mpf(1)/6) - y0/2 + mp.mpf(63)/50*dmin
ok(f"B0.56 {mp.nstr(b56,6)} < -1/20", b56 < -0.05)
ok(f"capacity l0/h0={mp.nstr(l0/h0,6)} > 1/5 > 4/21", l0/h0 > 0.2 > 4/21)
ok("y0 < 4 l0", y0 < 4*l0)
ok("0 < h0-4l0 < 13/100 and h0 < 4/5", 0 < h0-4*l0 < 0.13 and h0 < 0.8)
