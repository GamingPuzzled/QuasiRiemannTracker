import mpmath as mp
mp.mp.dps = 60
ok = lambda name, c: print(("PASS " if c else "FAIL ") + name)
l0 = mp.findroot(lambda x: 927*x**4-3135*x**3+2433*x**2+275*x-100, 0.167)
b = mp.mpf(11)/12 - l0/4; k = (5-3*l0)/6
d0 = (25-60*l0+27*l0**2)/(6*(5+4*l0-3*l0**2))
C0 = (6*k-1)/(15*k-2); V0 = mp.mpf(5)/6-d0
Rp = -1 + C0*(V0**2-d0**2*C0)/(2*(V0+d0*C0)**2)
h = -((1+l0)/2 + l0/(6*k))/(3*Rp/4 - 1/(24*k)); l = l0
def Cf(kk, pp): return 2*(1-pp)*(1-pp/(3*kk))/(3-pp-2*pp/(3*kk))
def R(kk, d, pp):
    V = mp.mpf(5)/6-d; c = Cf(kk,pp); return 1-d+d*c*V/(2*(V+d*c))
def F4(d, pp): return (1+l)/2+d/2+h*(3*R(k,d,pp)/4-mp.mpf(5)/12)+d*pp*(l-(h-4*l)/(12*k))
j00=300*k-40; j01=120*k+80; j10=24-216*k; j11=144*k-96
B0=24*k-4; B1=48*k; B2=16
T = l-(h-4*l)/(12*k); k0_=b-(1+l)/2-h/3; k1=-mp.mpf(1)/2+3*h/4-T/2
Ks = 8*j10*k1+18*h*B0
def J(d,u): return j00+j01*u+d*(j10+j11*u+96*u**2)
def A1(d): return 8*j01*k0_+(8*(j00*T+j01*k1+j11*k0_)-15*h*B1)*d+(8*(j10*T+j11*k1)+18*h*B1)*d**2
def A2(d): return (8*(j01*T+96*k0_)-15*h*B2)*d+(8*(j11*T+96*k1)+18*h*B2)*d**2
maxerr_R = maxerr_id = mp.mpf(0); Jmin, Jmax = mp.inf, -mp.inf; a1ok = a2ok = True
for i in range(1, 201):
    d = k*mp.mpf(i)/200
    a1, a2 = A1(d), A2(d)
    a1ok &= a1 > 36-174*d+347*d**2; a2ok &= a2 > 38*d+257*d**2
    for jj in range(0, 51):
        u = mp.mpf(jj)/100; pp = mp.mpf(1)/2-u
        Rsub = 1-d+d*(B0+B1*u+B2*u**2)*(5-6*d)/(2*J(d,u))
        maxerr_R = max(maxerr_R, abs(Rsub-R(k,d,pp)))
        lhs = 8*J(d,u)*(b-F4(d,pp)); rhs = Ks*(d-d0)**2+u*a1+u**2*a2+768*T*u**3*d**2
        maxerr_id = max(maxerr_id, abs(lhs-rhs))
        jv = J(d,u); Jmin = min(Jmin, jv); Jmax = max(Jmax, jv)
ok(f"R substitution exact (max err {mp.nstr(maxerr_R,3)})", maxerr_R < 1e-40)
ok(f"B0.49 identity (max err {mp.nstr(maxerr_id,3)})", maxerr_id < 1e-30)
ok(f"69<J<=270 on grid (min {mp.nstr(Jmin,6)}, max {mp.nstr(Jmax,6)})", Jmin > 69 and Jmax <= 270)
ok("A1 > 36-174d+347d^2 on grid", a1ok); ok("A2 > 38d+257d^2 on grid", a2ok)
ok(f"179<K*<180 ({mp.nstr(Ks,8)})", 179 < Ks < 180)
ok(f"117<768T<119 ({mp.nstr(768*T,8)})", 117 < 768*T < 119)
ok("36-174d+347d^2 > 0 always (disc<0)", 174**2 - 4*347*36 < 0)
# T1.3: dC/dk <= 1/18 on p in [0,1/2], k>=7/10;  dR/dC <= 5/12;  dF/dk < 1/10
mx = max(mp.diff(lambda kk: Cf(kk, mp.mpf(j)/100), kk0) for j in range(0,51) for kk0 in (mp.mpf(7)/10, mp.mpf(3)/4, mp.mpf(2)))
ok(f"max dC/dkappa {mp.nstr(mx,6)} < 1/18", mx < mp.mpf(1)/18)
mxF = max(mp.diff(lambda kk: (1+l)/2+d/2+h*(3*R(kk,d,pp)/4-mp.mpf(5)/12)+d*pp*(l-(h-4*l)/(12*kk)), kk0)
          for d in [k*mp.mpf(i)/40 for i in range(1,41)] for pp in [mp.mpf(j)/20 for j in range(0,11)] for kk0 in (mp.mpf(7)/10, k))
ok(f"max dF4/dkappa {mp.nstr(mxF,6)} < 1/10", mxF < mp.mpf(1)/10)
