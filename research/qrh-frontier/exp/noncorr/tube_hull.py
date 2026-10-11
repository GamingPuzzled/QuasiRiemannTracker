"""Bochner tube-hull bookkeeping (exact LP) for the pole part of C(h') via the triple series
   Z(s,w,u) = sum_{a,m,mu} gt6(a) chi_a(h') conj chi_a(m) * psi abar^2(m) gt6(m) (mu/m)_2 * gt3(mu) alpha(mu)
              * N(a)^-w N(m)^-s N(mu)^-u .
Pole part of C(h') ~ Y^{2/3} * Res_{w=2/3} Z, evaluated (Mellin) at m ~ Q, N(mu) <= D with mu-weight N(mu)^{-1/2}:
   X := |sum_m T(m) rho6(h'^5 m)|  <<  Q^sigma D^(upsilon - 1/2)     for (sigma, 2/3, upsilon) in the holomorphy hull
(assuming polynomial growth propagates to the hull: Phragmen-Lindelof in several variables, as in the MDS literature).
Requirement (remove the P_a^{1/6} loss at h'~1): |C(h')|^2 << Y P_a^{-1/6} E_a with E_a = Q^2, i.e.
   X << Q Y^{-1/6} P_a^{-1/12}.
'True trivial' bound: X <= ||T||_2 ||rho6||_2 ~ Q^{1/2} Q^{5/12} = Q^{11/12}  (sum_m |rho6(m)|^2 ~ Q^{5/6}: B0.22 Gram
exceptional term; ||R||_2^2 = O(log)).  So only a saving Z^{-0.0145} over the trivial bound is needed -- but any
argument that bounds mu-pieces separately pays sum_mu |rho_mu| ~ D^{1/2} (l1 mass) instead of the l2 mass O(1).

Regions where (w-2/3)(w-1/3)Z is holomorphic, with HEURISTIC convexity growth exponents g(x) = hc*max(0, 1-x, 1-2x):
   gPsi(w)  [hc=1/2] untwisted sextic Kubota Psi(w;n)=sum_a gt6(a)conj chi_a(n)N(a)^-w in its index n   (Patterson/Kubota)
   gK(s)    [hc=1/2] ell=2 twisted sextic Kubota sum_m psi abar^2 gt6(m) conj chi_m(n) N(m)^-s: ENTIRE (U-E4 / DR Prop 5.2 analogue)
   gTh(u)   [hc=cth/2] Theta_m(u)=sum_mu gt3 alpha (mu/m)_2 N(mu)^-u: ell=1 twisted cubic Kubota, quadratic twist of
            conductor ~m: ENTIRE; analytic conductor ~N(m)^cth, cth=2 (GL2-like, default) or 1 (optimistic)
   R-m  : sigma >= 1 + gPsi(omega) + gTh(upsilon)                 (m outermost; a and mu inner)
   R-a  : omega >= 1,  sigma >= 1 + gTh(upsilon)                  (a outermost; K by absolute convergence in m)
   R-mu : omega >= 1 + gK(sigma), upsilon >= 1 + 3 gK(sigma)       (a or mu outermost, K via entireness at index a mu^3)
"""
from fractions import Fraction as F
from sympy import Rational
from sympy.solvers.simplex import linprog

LQ, LD, LY = 0.83288, 0.66576, 0.45971   # log_Z of Q, D = Q^2/Z, Y  (paper geometry, solverD/E)
LPA = 2*LY - LQ
NEED = LQ - LY/6 - LPA/12
TRIV = LQ*11/12


def pieces(hc):  # g(x) = max over pieces of (alpha + beta x)
    return [(0, 0), (hc, -hc), (hc, -2*hc)]


def region_constraints(cth, gs):
    """list of regions; each region = list of (a(3-vector), beta) meaning a.x >= beta, x=(sigma, omega, upsilon)."""
    hP, hK, hT = gs*Rational(1, 2), gs*Rational(1, 2), gs*Rational(cth, 2)
    Rm = [((1, -bp, -bt), 1 + ap + at) for (ap, bp) in pieces(hP) for (at, bt) in pieces(hT)]
    Ra = [((0, 1, 0), 1)] + [((1, 0, -bt), 1 + at) for (at, bt) in pieces(hT)]
    Rmu = [((-bk, 1, 0), 1 + ak) for (ak, bk) in pieces(hK)] + [((-3*bk, 0, 1), 1 + 3*ak) for (ak, bk) in pieces(hK)]
    return [Rm, Ra, Rmu]


def hull_min(regions, cvec, w0=Rational(2, 3)):
    """min c.p over p in conv(union regions) with p_omega <= w0. Variables per region: q (3, free), lam >= 0."""
    R = len(regions); nv = 4*R
    A, b = [], []
    for k, reg in enumerate(regions):
        for (a, beta) in reg:                      # a.q_k - beta*lam_k >= 0   ->  -a.q_k + beta*lam_k <= 0
            row = [0]*nv
            for i in range(3):
                row[4*k+i] = -a[i]
            row[4*k+3] = beta
            A.append(row); b.append(0)
    row = [0]*nv                                   # sum_k q_k[omega] <= w0
    for k in range(R):
        row[4*k+1] = 1
    A.append(row); b.append(w0)
    Aeq = [[1 if (j % 4 == 3) else 0 for j in range(nv)]]; beq = [1]
    c = [0]*nv
    for k in range(R):
        for i in range(3):
            c[4*k+i] = cvec[i]
    bounds = [(None, None) if (j % 4 != 3) else (0, None) for j in range(nv)]
    val, x = linprog(c, A, b, Aeq, beq, bounds=bounds)
    p = [sum(x[4*k+i] for k in range(R)) for i in range(3)]
    return val, p


if __name__ == '__main__':
    print(f"geometry (log_Z): Q={LQ}, D={LD}, Y={LY}, P_a={LPA:.4f}.  Requirement X << Z^{NEED:.4f};  true trivial X << Z^{TRIV:.4f}")
    smax = (NEED - LD/2)/LQ
    print(f"per-mu route (bound each mu-piece by Q^sigma N(mu)^theta, then triangle inequality over mu):"
          f" needs sigma <= {smax:.4f} - {LD/LQ:.3f}*theta;  square-root cancellation per mu is sigma = 5/12 = 0.4167,"
          f" per-mu trivial is 11/12")
    for cth, gs, lab in ((2, 1, 'convexity, Theta conductor m^2'), (1, 1, 'convexity, Theta conductor m'),
                         (2, Rational(1, 2), 'all growth exponents halved (subconvexity-type input)')):
        regs = region_constraints(cth, gs)
        v1, p1 = hull_min(regs, (1, 0, 0))
        cq, cd = Rational(LQ).limit_denominator(10**6), Rational(LD).limit_denominator(10**6)
        v2, p2 = hull_min(regs, (cq, 0, cd))
        X = float(v2) - LD/2
        # with mu treated trivially (absolute values over mu: costs D^{1/2}) but per-mu sigma*
        Xtriv = LQ*float(v1) + LD/2
        print(f"[{lab}]\n   per-mu sigma* (u large) = {float(v1):.4f} at (s,w,u)={[round(float(t), 3) for t in p1]}"
              f"  (vs per-mu trivial 0.9167: {'power saving' if float(v1) < 11/12 else 'no saving'}); even ignoring mu-growth,"
              f" mu summed trivially: X << Z^{Xtriv:.4f}  {'OK' if Xtriv < NEED else 'FAILS'}"
              f"\n   joint (s,u) optimum: X << Z^{X:.4f} at (s,w,u)={[round(float(t), 3) for t in p2]}   {'OK' if X < NEED else 'FAILS'} (need {NEED:.4f})")
