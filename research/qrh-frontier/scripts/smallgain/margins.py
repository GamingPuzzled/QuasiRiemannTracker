"""(a) Every side condition of the certificate at the exact record point: active or slack?
A constraint with positive margin has zero shadow price (KKT), so sharpening its constant
(1/18, 1/10, 1/20, 4/21, 69/270, .06, .01, ...) cannot move the root of the quartic."""
from mpmath import mp, mpf
from base import record, R, tstar, Ck
mp.dps = 30
P = record(); l, k, d, h, x, y, b = (P[s] for s in ('l','k','d','h','x','y','b'))
p = mpf(1)/2; q = d*p; c3 = q/(3*k); M = x+y
ts = tstar(d, p, k); r = (1-q-(1-c3)+2*(d-c3)*ts)/((d-q)+2*(d-c3)); ms = ts-r
mpr = 1-y/h
C = Ck(p, k); V = mpf(5)/6-d
dCdk = 6*p*(1-p)**2/(3*k*(3-p)-2*p)**2
dRdC = d*V**2/(2*(V+d*C)**2)
rows = [
 ("M+l<=1  (B0.26 hypothesis)",        1-(M+l)),
 ("p<=1/2  (amplitude cap)",           mpf(0)),
 ("l<x",                               x-l),
 ("x<=y",                              y-x),
 ("y-2x+l<0 (Gram P_a^2/Y term)",      -(y-2*x+l)),
 ("energy (1+3l-2M)_+ argument <0",    -(1+3*l-2*M)),
 ("supply l/h>4/21",                   l/h-mpf(4)/21),
 ("numerator cap y<4l",                4*l-y),
 ("kappa<=3/4 (MP proved range)",      mpf(3)/4-k),
 ("kappa>=7/10",                       k-mpf(7)/10),
 ("delta<=kappa",                      k-d),
 ("T1.3 dC/dk<1/18",                   mpf(1)/18-dCdk),
 ("T1.3 dR/dC<=5/12",                  mpf(5)/12-dRdC),
 ("T1.3 dF/dk<1/10 (needs only <1/2)", None),
 ("h0-4l0<13/100",                     mpf(13)/100-(h-4*l)),
 ("h0<4/5",                            mpf(4)/5-h),
 ("other endpoint < b0-.06",           (b-mpf('0.06'))-(mpf(11)/12-49*h/150+11*l/12)),
 ("floor < b0-.01",                    (b-mpf('0.01'))-(mpf('0.51')+h/3+mpf('0.51')*l)),
 ("B0.56 < -1/20",                     -mpf(1)/20-(h*(mpf(17)/50-mpf(1)/6)-y/2+mpf(63)/50/200)),
 ("t* in [1,3/2]",                     min(ts-1, mpf(3)/2-ts)),
 ("r*>=13/21",                         r-mpf(13)/21),
 ("m* in [1/3,1/2]",                   min(ms-mpf(1)/3, mpf(1)/2-ms)),
 ("MI 2nd cond 3-2r-8z at z=(1-r)/2",  2*r-1),
 ("MI cap (1-r*)/2 <= supply l/h",     l/h-(1-r)/2),
 ("numerator m'<1/2",                  mpf(1)/2-mpr),
 ("pointwise - Holder (U units)",      (R(d,p,k)+d*mpr/2)-(mpf(1)/4+3*R(d,p,k)/4-q*(1-2*mpr)/(6*k)/2)),
]
for name, v in rows:
    print(f"{name:40s} {('ACTIVE' if v is not None and abs(v)<mpf('1e-25') else (mp.nstr(v,6) if v is not None else 'slack (see T1.3)')):>16}")
dFdk = mpf(3)/4*mpf(4)/5*mpf(5)/12*dCdk + mpf(5)/12*(h-4*l)/(12*k**2)
print("actual dF/dk bound ingredients: dC/dk =", mp.nstr(dCdk,6), " dR/dC =", mp.nstr(dRdC,6))
print("theta - b0 (rational rounding) =", mp.nstr(mpf(683505193)/781250000-b, 6))
