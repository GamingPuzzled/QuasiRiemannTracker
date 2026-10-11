"""Exact (mpmath) record point of the cycle25 quartic certificate and model pieces."""
from mpmath import mp, mpf, findroot, sqrt
mp.dps = 50
def record():
    p = lambda X: 927*X**4-3135*X**3+2433*X**2+275*X-100
    l = findroot(p, mpf('0.16712007850681'))
    k = (5-3*l)/6
    d = (25-60*l+27*l**2)/(6*(5+4*l-3*l**2))
    C = (6*k-1)/(15*k-2); V = mpf(5)/6-d
    Rp = -1 + C*(V**2-d**2*C)/(2*(V+d*C)**2)
    h = -((1+l)/2 + l/(6*k))/(3*Rp/4 - 1/(24*k))
    x = 1+l-h; y = h-2*l; b = mpf(11)/12 - l/4
    return dict(l=l, k=k, d=d, h=h, x=x, y=y, b=b)
def Ck(p, k): return 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k))
def R(d, p, k, V0=mpf(5)/6):
    C = Ck(p, k); V = max(V0-d, 0)
    return 1-d+d*C*V/(2*(V+d*C))
def tstar(d, p, k, V0=mpf(5)/6):
    C = Ck(p, k); V = V0-d
    return 1 + d*C/(2*(V+d*C))
def F4(d, p, k, P):
    l, h = P['l'], P['h']
    return (1+l)/2 + d/2 + h*(3*R(d,p,k)/4 - mpf(5)/12) + d*p*(l-(h-4*l)/(12*k))
if __name__ == "__main__":
    P = record()
    for kk, v in P.items(): print(kk, mp.nstr(v, 20))
    print("F4 - b0 =", mp.nstr(F4(P['d'], mpf(1)/2, P['k'], P) - P['b'], 5))
    print("R0 =", mp.nstr(R(P['d'], mpf(1)/2, P['k']), 15), " t* =", mp.nstr(tstar(P['d'], mpf(1)/2, P['k']), 15))
