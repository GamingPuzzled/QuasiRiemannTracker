"""Generalized exact model of the cycle25 binding system (B0.27, B0.38-B0.41, B0.47).

Every exponent coefficient that feeds the quartic is a named parameter, so that the exact
first-order effect of each on theta can be computed by the implicit-function theorem.
Paper values in PAPER reproduce l0, delta0, h0, theta to 50 digits.
"""
from mpmath import mp, mpf, findroot, diff, matrix, lu_solve
mp.dps = 40
PAPER = dict(
    aI=mpf(1)/2,     # MI capacity (1-r)*aI          (B0.4: r+2z<=1)
    cP=mpf(6),       # MP slot cost 6*kappa           (B0.5)
    cN=mpf(6),       # MP slot cost used in numerator z_A (B0.41)
    V0=mpf(5)/6,     # long sixth-power slope         (B0.9)
    H1=mpf(1)/4,     # numerator Holder: U^{H1+(1-H1)R}
    pmax=mpf(1)/2,   # amplitude cap g_i<=pmax*delta  (B0.16, B0.37)
    cg=mpf(1)/12,    # probe Gram exponent c/2, c=1/6 (B0.22)
    ca=mpf(1)/2,     # a=(1+ca*2*delta)/2: coefficient of delta in a is ca
    cy=mpf(1)/2,     # Y^{-cy}: numerator normalization
    cl=mpf(1)/2,     # Z^{-cl*l}: slot normalization
    c6=mpf(1)/6,     # 1/6 = zeta(6z) pole location in h(z0-1/6)
)
def b_of(l, x, P):
    y = 1-l-x
    return x/2 + P['cg']*(y-x) + mpf(5)/6 - x/3 - l/6
def R_of(d, p, k, P):
    q = d*p
    # A_I(r)=1-d r-2q aI(1-r);  A_P(m)=1-2 d m-2q(1-2m)/(cP k);  long(t)=1-d+(V0-d)(t-1)
    aI, cP, V0 = P['aI'], P['cP'], P['V0']
    # A_I(r) = alpha_I - beta_I r ; A_P(t-r) = alpha_P - beta_P (t-r)
    aIc, bI = 1-2*q*aI, d-2*q*aI
    aPc, bP = 1-2*q/(cP*k), 2*d-4*q/(cP*k)
    # crossing: aIc - bI r = aPc - bP t + bP r  -> r = (aIc-aPc+bP t)/(bI+bP)
    # short(t) = aIc - bI r(t) = s0 - s1 t
    s1 = bI*bP/(bI+bP); s0 = aIc - bI*(aIc-aPc)/(bI+bP)
    V = V0-d
    t = (s0-1+d+V)/(s1+V)
    return 1-d+V*(t-1), t
def F_of(d, p, l, h, P):
    x = 1+l-h; y = 1-l-x; b = b_of(l, x, P); k = 2*b-1
    R, _ = R_of(d, p, k, P); q = d*p
    zAh = (2*y-h)/(P['cN']*k)          # d*z_A(d) at d=h (2y-h<h here)
    return (mpf(1)/2 + P['ca']*d - P['cy']*y - P['cl']*l + l*q
            + h*(P['H1'] + (1-P['H1'])*R - P['c6']) - q*zAh/2), b
def system(P):
    p = P['pmax']
    def eqs(l, d, h):
        F, b = F_of(d, p, l, h, P)
        return [F-b,
                diff(lambda dd: F_of(dd, p, l, h, P)[0], d),
                lagr(l, d, h)]
    def G(l, d, h):
        F, b = F_of(d, p, l, h, P); return F-b
    def B(l, h): return b_of(l, 1+l-h, P)
    def lagr(l, d, h):
        # optimum of min b s.t. G<=0 over (l,h):  db/dh*dG/dl - db/dl*dG/dh = 0
        return (diff(lambda hh: B(l, hh), h)*diff(lambda ll: G(ll, d, h), l)
                - diff(lambda ll: B(ll, h), l)*diff(lambda hh: G(l, d, hh), h))
    return eqs
def solve(P=None, guess=(mpf('0.16712007850681'), mpf('0.4693448586981'), mpf('0.7939475973881'))):
    P = {**PAPER, **(P or {})}
    l, d, h = findroot(system(P), guess)
    x = 1+l-h
    return dict(l=l, d=d, h=h, theta=b_of(l, x, P))
if __name__ == "__main__":
    s = solve()
    for k_, v in s.items(): print(k_, mp.nstr(v, 25))
    print("paper theta-b0 check:", mp.nstr(mpf(683505193)/781250000 - s['theta'], 5))
