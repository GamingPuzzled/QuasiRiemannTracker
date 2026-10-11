"""B4 pricing: row count with two zero witnesses (numerator psibar_u, denominator eta*psibar_u)
and a hypothetical JOINT inverse moment

  (J_e)  sum_{q_u ~ U, 6-free} |M_u(psibar; U^{r1}) M_u(eta psibar; U^{r2})|^2 << U^{e(r1,r2)+eps}.

Count lines (log_U of #rows), at bin (delta, q = delta p, kappa):
  per character i (Lemma 1W with its own truncation t_i, witnesses r_i+m_i >= t_i, m_i <= 1/2):
     A_I(r) = 1-q-(delta-q) r              (MI, r<=1)
     long(r)= (1+5r)/6 - delta r           (B0.9, r>1)
     A_P(m) = 1-q/(3k)-2(delta-q/(3k)) m   (MP)
  joint:  J(x) = e(x) - delta x,  x = r1+r2, with e(x) = 1 + lam*(x-1)_+       ('free' family)
          optional MI-style slots for x<1:  1-q-(delta-q)x  (slot capacity (1-x)/2)
          lam = inf  means: joint bound only in the Poisson range x<=1.
R_new = min_{t1,t2} max_{adversary (r1,r2)} min(all lines).
"""
import numpy as np

def lines_1d(t, delta, q, k, nr):
    lo = max(0.0, t-0.5)
    r = np.linspace(lo, t, nr)
    m = t - r
    inv = np.where(r <= 1, 1-q-(delta-q)*r, (1+5*r)/6-delta*r)
    AP = 1-q/(3*k)-2*(delta-q/(3*k))*m
    return r, np.minimum(inv, AP)

def joint(x, delta, q, lam, slots=True):
    if np.isinf(lam):
        J = np.where(x <= 1, 1-delta*x, np.inf)
    else:
        J = 1 + lam*np.maximum(x-1, 0) - delta*x
    if slots:
        J = np.minimum(J, np.where(x <= 1, 1-q-(delta-q)*x, np.inf))
    return J

def adversary(t1, t2, delta, q, k, lam, nr=301, slots=True):
    r1, c1 = lines_1d(t1, delta, q, k, nr)
    r2, c2 = lines_1d(t2, delta, q, k, nr)
    X = r1[:, None] + r2[None, :]
    C = np.minimum(np.minimum(c1[:, None], c2[None, :]), joint(X, delta, q, lam, slots))
    i = np.unravel_index(np.argmax(C), C.shape)
    return C[i], (r1[i[0]], r2[i[1]])

def R_paper(delta, p, k):
    C = 2*(1-p)*(1-p/(3*k))/(3-p-2*p/(3*k)); V = max(5/6-delta, 0)
    return 1-delta+delta*C*V/(2*(V+delta*C))

def R_new(delta, p, k, lam, tmin=0.6, nt=46, nr=201, slots=True, sym_only=False, refine=True):
    q = delta*p
    ts = np.linspace(tmin, 1.5, nt)
    best = (np.inf, None)
    pairs = [(a, a) for a in ts] if sym_only else [(a, b) for a in ts for b in ts if b >= a]
    for a, b in pairs:
        v, _ = adversary(a, b, delta, q, k, lam, nr, slots)
        if v < best[0]: best = (v, (a, b))
    if refine:
        (a0, b0) = best[1]; h = (1.5-tmin)/(nt-1)
        for a in np.linspace(max(tmin, a0-h), min(1.5, a0+h), 9):
            for b in np.linspace(max(tmin, b0-h), min(1.5, b0+h), 9):
                if sym_only and b != a: continue
                v, _ = adversary(a, b, delta, q, k, lam, 2*nr, slots)
                if v < best[0]: best = (v, (a, b))
    return min(best[0], R_paper(delta, p, k)), best   # paper's single-witness count always available

if __name__ == "__main__":
    import sys
    sys.path.insert(0, '.')
    from base import record
    P = record(); d0, k0 = float(P['d']), float(P['k'])
    R0 = R_paper(d0, .5, k0)
    print(f"binding delta0={d0:.12f} kappa0={k0:.12f}  R_paper={R0:.12f}")
    # sanity: lam huge => joint line never active => R_new == R_paper (grid)
    for lam in [10.0, np.inf, 1.0, 5/6, 0.75, 0.7, 0.65, 0.6, 0.55, 0.5, 0.469, 0.4, 0.3, 0.2, 0.0]:
        Rn, (v, tt) = R_new(d0, .5, k0, lam)
        print(f"lam={lam:7.4f}  R_new={Rn:.6f}  dR={Rn-R0:+.3e}  (raw adversary {v:.6f} at t={tt[0]:.4f},{tt[1]:.4f})")
