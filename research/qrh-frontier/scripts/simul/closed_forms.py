"""Closed forms for the row count with a two-witness inverse product line (verified against joint_lp.py).

Notation at a bin (delta, q = delta p, kappa):  A = 1-q, s = delta-q, B = 1-q/(3k), w = delta-q/(3k)
  r_v = (A-v)/s          (MI crossing:  A_I(r) = v)
  m_v = (B-v)/(2w)       (MP crossing:  A_P(m) = v)
  t_v = (v-1/6)/(5/6-delta)   (long crossing: (1+5t)/6 - delta t = v)
  joint line beyond Poisson range:  (1-delta) + (lam-delta)(x-1),  x = r1+r2 >= 1.

Family P0 (truncations >= 1 only, or no Poisson-range joint bound):
   2 r_v - 1 = (v-1+delta)/(lam-delta)                        -> v_P0
Family P2 (MI2 in the Poisson range: sum |M(U^r1) M(U^r2) Q_I|^2 << U for r1+r2+2z<=1;
           numerator/second truncation t1 in (m_v + t_s - t_v, x_b - r_v), t_s = r_v+m_v):
   2 (r_v + m_v) - t_v = 1 + (v-1+delta)/(lam-delta)            -> v_P2
   At v = R (paper) this reads t* = 1 + excess/(lam-delta), and excess/(t*-1) = 5/6-delta identically,
   so the P2 break-even is lam_be = 5/6 EXACTLY.
Family P1 (UM2 only in the Poisson range, no slots): x_a = (1-v)/delta replaces r_v in condition (1):
   (t1 - m_v) + (t_v - m_v) > (1-v)/delta  and  t1 + r_v < x_b
   -> (1-v)/delta + 2 m_v - t_v + r_v = x_b(v)                  -> v_P1
"""
import numpy as np

def _coef(d, p, k, lam):
    q = d*p; A = 1-q; s = d-q; B = 1-q/(3*k); w = d-q/(3*k); V = 5/6-d; g = lam-d
    # each quantity is  alpha + beta*v
    r = (A/s, -1/s); m = (B/(2*w), -1/(2*w)); t = (-1/(6*V), 1/V); xb = (1-(1-d)/g, 1/g)
    xa1 = (1/d, -1/d)
    return r, m, t, xb, xa1

def _solve(lhs, rhs):
    # lhs, rhs affine (alpha, beta): alpha1+beta1 v = alpha2+beta2 v
    return (rhs[0]-lhs[0])/(lhs[1]-rhs[1])

def v_P0(d, p, k, lam):
    r, m, t, xb, _ = _coef(d, p, k, lam)
    return _solve((2*r[0], 2*r[1]), xb)

def v_P2(d, p, k, lam):
    r, m, t, xb, _ = _coef(d, p, k, lam)
    lhs = (2*r[0]+2*m[0]-t[0], 2*r[1]+2*m[1]-t[1])
    return _solve(lhs, xb)

def v_P1(d, p, k, lam):
    r, m, t, xb, xa = _coef(d, p, k, lam)
    lhs = (xa[0]+2*m[0]-t[0]+r[0], xa[1]+2*m[1]-t[1]+r[1])
    return _solve(lhs, xb)

def checks(d, p, k, lam, v, L=0.2104):
    """regime conditions used in the derivation"""
    r, m, t, xb, _ = _coef(d, p, k, lam)
    rv = r[0]+r[1]*v; mv = m[0]+m[1]*v; tv = t[0]+t[1]*v; xbv = xb[0]+xb[1]*v
    return dict(rv=rv, mv=mv, tv=tv, xb=xbv, t1=(mv+rv+mv-tv),
                ok=bool(lam > d and 0.5 < rv <= 1 and (1-rv)/2 <= L and 0 < mv <= 0.5 and tv > 1
                        and (1-2*mv)/(6*k) <= L and xbv > 1))

def Rp(d, p, k):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = 5/6-d
    return 1-d+d*C*V/(2*(V+d*C))

if __name__ == "__main__":
    import sys; sys.path.insert(0, '.')
    from base import record
    P = record(); d0, k0 = float(P['d']), float(P['k']); R0 = Rp(d0, .5, k0)
    print("R_paper =", R0)
    print("P2 break-even check: v_P2 at lam=5/6 - R0 =", v_P2(d0, .5, k0, 5/6) - R0)
    print(" lam   v_P0     v_P1     v_P2   (min with R0; LP tmin=1 / tmin=0 for comparison in out_lp_*.txt)")
    for lam in [5/6, 0.8, 0.75, 0.7, 0.674, 0.65, 0.6, 0.59, 0.58, 0.55, 0.5]:
        a, b, c = (min(R0, f(d0, .5, k0, lam)) for f in (v_P0, v_P1, v_P2))
        print(f"{lam:.4f} {a:.6f} {b:.6f} {c:.6f}  P2 regime {checks(d0,.5,k0,lam,c)}")
