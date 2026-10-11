"""Scaling grid for the non-correlation test.  One symbol table up to 2*Ymax; cells (Y, Q) for nested Y, Q lists.

For each cell:  A_m(h') = sum_{a in [Y,2Y]} gt(a) chi_a(h')^P conj chi_a(m)^P over two random halves A0, A1 of the a's.
Statistics for a vector v (see noncorr.py):
   r(h')        = |sum v A|^2 / sum |v|^2 |A|^2                      null ~ 1 (exponential), aligned ~ #m*f
   #m*rho2(h')  = #m * Re(C0 conj C1)/(||v||^2 Re<A0,A1>)            null mean ~ 1, aligned ~ #m
Families of vectors (each evaluated at h'=1 where the pole is strongest; r also pooled over h' != 1):
   nu        : 36 ray-class characters nu(m) = chi_m(lam)^i chi_m(-1)^j chi_m(omega)^k   (untwisted; positive control family)
   ab2*nu    : same times abar(m)^2                                                    (infinity-type test)
   mu        : (mu/m)_2 for each mu (components of R)                                 untwisted / ab2-twisted
   g*mu      : gamma(m)(mu/m)_2  [and conj gamma] untwisted / ab2-twisted  (ab2-twisted = components of T)
   T*nu      : abar^2 gamma R nu  vs  gamma R nu (untwisted)
   RAND      : random phases, same |v| as T (null family, 200 draws)
Output: JSON per cell + printed table.
"""
import sys, json, time, argparse
import numpy as np
from lib import prime_elements, products, symbol_table, gauss_composite, char_matrix, W3

ap = argparse.ArgumentParser()
ap.add_argument('--n', type=int, default=6)
ap.add_argument('--Ys', default='25000,50000,100000')
ap.add_argument('--Qs', default='500,1000,2000')
ap.add_argument('--conv', default='lam')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--nh', type=int, default=12)
ap.add_argument('--nrand', type=int, default=200)
ap.add_argument('--Dexp', type=float, default=0.8)
ap.add_argument('--out', default='out/grid.json')
args = ap.parse_args()
n = args.n; P = 1 if n == 6 else 2
Ys = [int(x) for x in args.Ys.split(',')]; Qs = [int(x) for x in args.Qs.split(',')]
rng = np.random.default_rng(args.seed)
t0 = time.time()
Ymax, Qmax = max(Ys), max(Qs)
Dmax = int(Qmax**args.Dexp)
pel = prime_elements(2*max(Ymax, Qmax))
Msets = {Q: products(pel, Q, 2*Q) for Q in Qs}
mu_all = products(pel, 2, Dmax)
H = [(1, 0)]
cand = sorted((x*x-x*y+y*y, x, y) for x in range(-20, 21) for y in range(-20, 21)
              if 1 < x*x-x*y+y*y <= 400 and (x*x-x*y+y*y) % 6 in (1, 5) and x % 3 == 2 and y % 3 == 0)
H += [(x, y) for (_, x, y) in cand[:args.nh-1]]
UNITS_T = [(-1, 0), (0, 1), (1, 2)]   # -1, omega, lambda = 1+2 omega
targets = []; off = {}
for Q in Qs:
    off[Q] = len(targets); targets += [m[0] for m in Msets[Q]]
off['H'] = len(targets); targets += H
off['mu'] = len(targets); targets += [u[0] for u in mu_all]
off['U'] = len(targets); targets += UNITS_T
TA = [t[0] for t in targets]; TB = [t[1] for t in targets]
imax = len(pel); dummy = imax
iD = sum(1 for Pp in pel if Pp['p'] <= Dmax)
print(f"n={n} conv={args.conv} Ys={Ys} Qs={Qs} #primes={imax} #targets={len(targets)} #mu={len(mu_all)}", flush=True)
E, Zr, G = symbol_table(pel, imax, TA, TB, specs=[(P, args.conv, imax), (2, 'chk', max(iD, 1))])
print(f"  symbol table {time.time()-t0:.0f}s", flush=True)
Gp = G[(P, args.conv)]; G3c = G[(2, 'chk')]
nH = len(H)


def kubota(Y, Q, KV=None):
    """Also returns, for the vectors KV (#k, #m), the a-side energies E_a = sum_a |F_a|^2, F_a = sum_m v_m conj chi_a(m)^P,
    and the a-side norms sum_a |gt(a) chi_a(h')|^2 (h'-dependent), for the a-side Cauchy-Schwarz ratio."""
    A_set = products(pel, Y, 2*Y)
    ga = np.array([gauss_composite(pel, idx, Gp, P) for (_, idx, _) in A_set])
    half = rng.integers(0, 2, len(A_set))
    colM = np.arange(off[Q], off[Q]+len(Msets[Q])); colH = np.arange(off['H'], off['H']+nH)
    cols = np.concatenate([colM, colH]); nM = len(colM)
    A0 = np.zeros((nH, nM), complex); A1 = np.zeros((nH, nM), complex)
    Ea = np.zeros(0 if KV is None else len(KV)); gn = np.zeros(nH)
    for s in range(0, len(A_set), 400):
        X = char_matrix([c[1] for c in A_set[s:s+400]], E, Zr, cols, P, dummy)
        Xm = np.conj(X[:, :nM]); Xh = X[:, nM:]; g = ga[s:s+400]; hb = half[s:s+400]
        A0 += (Xh*(g*(hb == 0))[:, None]).T @ Xm
        A1 += (Xh*(g*(hb == 1))[:, None]).T @ Xm
        gn += np.sum(np.abs(Xh*g[:, None])**2, axis=0)
        if KV is not None:
            Ea += np.sum(np.abs(KV @ Xm.T)**2, axis=1)
    return A0, A1, len(A_set), Ea, gn


def mside(Q):
    M = Msets[Q]; D = int(Q**args.Dexp)
    mus = [k for k, u in enumerate(mu_all) if u[2] <= D]
    cols = np.concatenate([off['mu'] + np.array(mus), off['U'] + np.arange(3)])
    X = char_matrix([m[1] for m in M], E, Zr, cols, 1, dummy)
    quad = X[:, :len(mus)]**3
    um1, uom, ulam = X[:, -3], X[:, -2], X[:, -1]
    mv = np.array([m[0][0] + m[0][1]*W3 for m in M]); ab2 = (np.conj(mv)/np.abs(mv))**2
    gam = np.array([gauss_composite(pel, idx, Gp, P) for (_, idx, _) in M])
    rho = np.array([gauss_composite(pel, mu_all[k][1], G3c, 2)*((mu_all[k][0][0]+mu_all[k][0][1]*W3)/abs(mu_all[k][0][0]+mu_all[k][0][1]*W3))
                    / np.sqrt(mu_all[k][2]) for k in mus])
    Rm = quad @ rho
    nu = np.array([ulam**i * um1**j * uom**k for i in range(6) for j in range(2) for k in range(3)])  # (36, #m)
    fam = {
        'nu (untw, +ctrl)': nu, 'ab2*nu': ab2*nu, 'a2*nu': np.conj(ab2)*nu,
        'mu (untw)': quad.T, 'ab2*mu': ab2*quad.T,
        'g*mu (untw)': gam*quad.T, 'conj(g)*mu (untw)': np.conj(gam)*quad.T,
        'ab2*g*mu [T comps]': ab2*gam*quad.T, 'ab2*conj(g)*mu': ab2*np.conj(gam)*quad.T,
        'g*nu (untw)': gam*nu, 'ab2*g*nu': ab2*gam*nu,
        'g*R*nu (untw)': gam*Rm*nu, 'T*nu = ab2*g*R*nu': ab2*gam*Rm*nu,
        'conj(g)*R*nu (untw)': np.conj(gam)*Rm*nu, 'ab2*conj(g)*R*nu': ab2*np.conj(gam)*Rm*nu,
    }
    T = ab2*gam*Rm
    fam['RAND (null, |T|)'] = np.abs(T)*np.exp(2j*np.pi*rng.random((args.nrand, len(M))))
    return fam, T, rho, quad, ab2, gam


def famstats(V, A0, A1):
    A = A0 + A1
    C = V @ A.T; C0 = V @ A0.T; C1 = V @ A1.T                  # (#v, #h)
    den = (np.abs(V)**2) @ (np.abs(A)**2).T
    r = np.abs(C)**2/den
    v2 = np.sum(np.abs(V)**2, axis=1)
    cross = np.real(np.sum(A0*np.conj(A1), axis=1))
    rho2 = np.real(C0*np.conj(C1))/(v2[:, None]*cross[None, :])
    return r, rho2, C


results = []
for Q in Qs:
    fam, T, rho, quad, ab2, gam = mside(Q)
    nM = len(Msets[Q])
    KVn = ['T', 'one', 'ab2'] + [f'g*nu{j}' for j in range(36)] + [f'T*nu{j}' for j in range(36)] + [f'RAND{j}' for j in range(20)]
    KV = np.vstack([T, np.ones(len(T)), ab2, fam['g*nu (untw)'], fam['T*nu = ab2*g*R*nu'], fam['RAND (null, |T|)'][:20]])
    for Y in Ys:
        A0, A1, nA, Ea, gn = kubota(Y, Q, KV)
        # a-side Cauchy-Schwarz ratio kappa_a(h') = |C(h')|^2 / (sum_a |gt chi_a(h')|^2 * E_a)  (<= 1; 'Y part' of solverE 3.1)
        Cfull = KV @ (A0 + A1).T
        kap = np.abs(Cfull)**2/(Ea[:, None]*gn[None, :])
        kd = {'T': kap[0], 'one': kap[1], 'ab2': kap[2], 'g*nu best': kap[3:39].max(axis=0), 'T*nu mean': kap[39:75].mean(axis=0),
              'RAND mean': kap[75:95].mean(axis=0)}
        Edl = {'T': Ea[0]/(nA*np.sum(np.abs(T)**2)), 'RAND': float(np.mean(Ea[75:95]/(nA*np.sum(np.abs(KV[75:95])**2, axis=1))))}
        A = A0 + A1
        cross = np.real(np.sum(A0*np.conj(A1), axis=1))
        f = 4*cross/np.sum(np.abs(A)**2, axis=1)
        z = cross/np.sqrt(np.sum(np.abs(A0)**2*np.abs(A1)**2, axis=1)/2)
        cell = dict(Y=Y, Q=Q, nA=nA, nM=nM, D=int(Q**args.Dexp), f=f.tolist(), z=z.tolist(), fam={})
        # data-driven pole control: conj(A_other)
        Cp = np.sum(np.conj(A1)*A0, axis=1); denp = np.sum(np.abs(A1)**2*np.abs(A0)**2, axis=1)
        cell['pole_ctrl_r'] = (np.abs(Cp)**2/denp).tolist()
        Cq = np.sum(ab2*np.conj(A1)*A0, axis=1)
        cell['ab2_pole_ctrl_r'] = (np.abs(Cq)**2/denp).tolist()
        print(f"\n[Y={Y} Q={Q}] #a={nA} #m={nM} D={cell['D']}  f(h'=1)={f[0]:.3f} z={z[0]:.1f} | f(h'!=1) mean {f[1:].mean():.3f} "
              f"| POLE ctrl r(1)={cell['pole_ctrl_r'][0]:.1f}  ab2*POLE r(1)={cell['ab2_pole_ctrl_r'][0]:.2f}  ({time.time()-t0:.0f}s)", flush=True)
        print(f"   {'family':24s} {'#v':>4s} {'<r(1)>':>7s} {'max r(1)':>8s} {'<r(!=1)>':>8s} {'<#m rho2(1)>':>12s} {'sd':>6s} {'max':>7s}")
        for name, V in fam.items():
            r, rho2, C = famstats(V, A0, A1)
            d = dict(nv=len(V), r1_mean=float(r[:, 0].mean()), r1_max=float(r[:, 0].max()), rrest_mean=float(r[:, 1:].mean()),
                     mrho1_mean=float(nM*rho2[:, 0].mean()), mrho1_sd=float(nM*rho2[:, 0].std()), mrho1_max=float(nM*rho2[:, 0].max()),
                     r1_list=r[:, 0].tolist() if len(V) <= 80 else None)
            cell['fam'][name] = d
            print(f"   {name:24s} {len(V):4d} {d['r1_mean']:7.2f} {d['r1_max']:8.2f} {d['rrest_mean']:8.2f} {d['mrho1_mean']:12.2f} {d['mrho1_sd']:6.2f} {d['mrho1_max']:7.1f}")
        # T itself, and coherence across mu: compare |sum_mu rho_mu c_mu|^2 to sum |rho_mu|^2 |c_mu|^2 (incoherent) and (sum|rho||c|)^2
        cmu = (ab2*gam*quad.T) @ A[0]          # per-mu correlations at h'=1
        tot = abs(np.sum(rho*cmu))**2
        inc = np.sum(np.abs(rho)**2*np.abs(cmu)**2); worst = np.sum(np.abs(rho)*np.abs(cmu))**2
        rT, rho2T, _ = famstats(T[None, :], A0, A1)
        cell['T'] = dict(r=rT[0].tolist(), mrho1=float(nM*rho2T[0, 0]), coh_tot=float(tot), coh_incoh=float(inc), coh_worst=float(worst))
        cell['kappa_a'] = {k: v.tolist() for k, v in kd.items()}
        cell['Ea_over_diag'] = Edl
        print(f"   T: r(1)={rT[0,0]:.2f}  <r(!=1)>={rT[0,1:].mean():.2f}  #m*rho2(1)={nM*rho2T[0,0]:.2f}  | mu-coherence |sum rho c|^2/sum|rho c|^2 = {tot/inc:.2f}"
              f" (worst-case CS-in-mu ratio {worst/inc:.1f})", flush=True)
        print("   a-side CS ratio kappa_a(h'=1) [<=1; 1/#m ~ %.1e]: " % (1/nM) + "  ".join(f"{k}={v[0]:.2e}" for k, v in kd.items())
              + f" | E_a/(#a||v||^2): T={Edl['T']:.2f} RAND={Edl['RAND']:.2f}", flush=True)
        results.append(cell)
        with open(args.out, 'w') as fh:
            json.dump(dict(args=vars(args), H=H, cells=results), fh)
print(f"done {time.time()-t0:.0f}s")
