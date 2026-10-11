"""Cross-window pooled estimate of the correlation of probe vectors with the Patterson pole direction.

Disjoint a-windows W_i = [Y_i, 2Y_i] give A^{(i)}_m(h') = c_i * pole_m(h') + noise_i, with the SAME pole direction
(residue rho6(h'^5 m) up to a positive window factor c_i ~ Y_i^{2/3}) and independent noise. For a vector v:
    rho2_pool(v) = sum_{i<j} Re(C_i conj C_j) / ( ||v||^2 * sum_{i<j} Re<A_i, A_j> ),   C_i = sum_m v_m A^{(i)}_m
estimates |<v, conj pole>|^2 / (||v||^2 ||pole||^2)  (null: ~1/#m;  full alignment: 1).  Reported as #m*rho2 (null ~ 1).
Also: the pooled pole energy fraction and a pole-direction estimate  P = sum_i A_i / c_i  for display.
"""
import json, time, argparse
import numpy as np
from lib import prime_elements, products, symbol_table, gauss_composite, char_matrix, W3

ap = argparse.ArgumentParser()
ap.add_argument('--n', type=int, default=6)
ap.add_argument('--Ys', default='12500,25000,50000,100000,200000')
ap.add_argument('--Qs', default='500,1000,2000,4000,8000')
ap.add_argument('--conv', default='lam')
ap.add_argument('--nh', type=int, default=6)
ap.add_argument('--nrand', type=int, default=400)
ap.add_argument('--index', default='md', help="'md': conj chi_a(d) -> theta index m*d (correct, (**)); 'mh5': chi_a(h') -> index m*h'^5 (first runs)")
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--out', default='out/pooled.json')
args = ap.parse_args()
n = args.n; P = 1 if n == 6 else 2
Ys = [int(x) for x in args.Ys.split(',')]; Qs = [int(x) for x in args.Qs.split(',')]
rng = np.random.default_rng(args.seed); t0 = time.time()
Dexp = 0.8
pel = prime_elements(2*max(max(Ys), max(Qs)))
Msets = {Q: products(pel, Q, 2*Q) for Q in Qs}
Dmax = int(max(Qs)**Dexp)
mu_all = products(pel, 2, Dmax)
H = [(1, 0)] + [(x, y) for (_, x, y) in sorted((x*x-x*y+y*y, x, y) for x in range(-20, 21) for y in range(-20, 21)
                if 1 < x*x-x*y+y*y <= 400 and (x*x-x*y+y*y) % 6 in (1, 5) and x % 3 == 2 and y % 3 == 0)[:args.nh-1]]
targets = []; off = {}
for Q in Qs:
    off[Q] = len(targets); targets += [m[0] for m in Msets[Q]]
off['H'] = len(targets); targets += H
off['mu'] = len(targets); targets += [u[0] for u in mu_all]
off['U'] = len(targets); targets += [(-1, 0), (0, 1), (1, 2)]
imax = len(pel); dummy = imax
iD = sum(1 for Pp in pel if Pp['p'] <= Dmax)
print(f"n={n} conv={args.conv} Ys={Ys} Qs={Qs} #primes={imax} #targets={len(targets)}", flush=True)
E, Zr, G = symbol_table(pel, imax, [t[0] for t in targets], [t[1] for t in targets],
                        specs=[(P, args.conv, imax), (2, 'chk', max(iD, 1))])
print(f"  symbol table {time.time()-t0:.0f}s", flush=True)
Gp = G[(P, args.conv)]; G3c = G[(2, 'chk')]; nH = len(H)
A_sets = {Y: products(pel, Y, 2*Y) for Y in Ys}
ga = {Y: np.array([gauss_composite(pel, idx, Gp, P) for (_, idx, _) in A_sets[Y]]) for Y in Ys}

out = dict(args=vars(args), H=H, cells=[])
for Q in Qs:
    M = Msets[Q]; nM = len(M); D = int(Q**Dexp)
    colM = np.arange(off[Q], off[Q]+nM); colH = np.arange(off['H'], off['H']+nH)
    cols = np.concatenate([colM, colH])
    As = []
    for Y in Ys:
        A = np.zeros((nH, nM), complex)
        for s in range(0, len(A_sets[Y]), 400):
            X = char_matrix([c[1] for c in A_sets[Y][s:s+400]], E, Zr, cols, P, dummy)
            Xh = X[:, nM:] if args.index == 'mh5' else np.conj(X[:, nM:])
            A += (Xh*ga[Y][s:s+400, None]).T @ np.conj(X[:, :nM])
        As.append(A)
    As = np.array(As)                                           # (#win, #h, #m)
    # m-side vectors
    mus = [k for k, u in enumerate(mu_all) if u[2] <= D]
    Xm = char_matrix([m[1] for m in M], E, Zr, np.concatenate([off['mu'] + np.array(mus), off['U'] + np.arange(3)]), 1, dummy)
    quad = Xm[:, :len(mus)]**3; um1, uom, ulam = Xm[:, -3], Xm[:, -2], Xm[:, -1]
    mv = np.array([m[0][0] + m[0][1]*W3 for m in M]); ab2 = (np.conj(mv)/np.abs(mv))**2
    gam = np.array([gauss_composite(pel, idx, Gp, P) for (_, idx, _) in M])
    rho = np.array([gauss_composite(pel, mu_all[k][1], G3c, 2)*((mu_all[k][0][0]+mu_all[k][0][1]*W3)/abs(mu_all[k][0][0]+mu_all[k][0][1]*W3))
                    / np.sqrt(mu_all[k][2]) for k in mus])
    Rm = quad @ rho; T = ab2*gam*Rm
    nu = np.array([ulam**i * um1**j * uom**k for i in range(6) for j in range(2) for k in range(3)])
    fam = {'T': T[None, :], 'one (untw)': np.ones((1, nM)), 'ab2': ab2[None, :],
           'nu (untw)': nu, 'ab2*nu': ab2*nu, 'T*nu': T*nu, 'g*R*nu (untw)': gam*Rm*nu,
           'ab2*g*mu [T comps]': ab2*gam*quad.T, 'g*nu (untw)': gam*nu, 'ab2*g*nu': ab2*gam*nu,
           'RAND': np.abs(T)*np.exp(2j*np.pi*rng.random((args.nrand, nM)))}
    # pooled cross-window statistics
    W = len(Ys); pairs = [(i, j) for i in range(W) for j in range(i+1, W)]
    crossA = np.array([sum(np.real(np.sum(As[i, h]*np.conj(As[j, h]))) for (i, j) in pairs) for h in range(nH)])
    totA = np.array([sum(np.sum(np.abs(As[i, h])**2) for i in range(W)) for h in range(nH)])
    # pole fraction pooled: for A_i = c_i p + e_i, sum_{i<j} c_i c_j ||p||^2 vs sum_i ||A_i||^2 (approx f with c_i ~ equal scale)
    cell = dict(Q=Q, nM=nM, D=D, cross_over_tot=(crossA/totA).tolist(), fam={})
    # z-score of pooled pole (null sd from random sign flips of windows' noise ~ sqrt(sum_{i<j} sum_m |A_i|^2|A_j|^2 / 2))
    sdA = np.array([np.sqrt(sum(np.sum(np.abs(As[i, h])**2*np.abs(As[j, h])**2) for (i, j) in pairs)/2) for h in range(nH)])
    cell['pole_z'] = (crossA/sdA).tolist()
    print(f"\n[Q={Q}] #m={nM} D={D}  pooled pole z(h')={np.round(crossA/sdA, 1).tolist()}  ({time.time()-t0:.0f}s)")
    print(f"   {'family':22s} {'#v':>4s} {'<#m rho2(1)>':>12s} {'se':>6s} {'max':>7s} {'<#m rho2(h!=1)>':>15s}")
    for name, V in fam.items():
        C = np.einsum('vm,whm->vwh', V, As)                    # (#v, #win, #h)
        num = sum(np.real(C[:, i]*np.conj(C[:, j])) for (i, j) in pairs)   # (#v, #h)
        v2 = np.sum(np.abs(V)**2, axis=1)
        r2 = nM*num/(v2[:, None]*crossA[None, :])
        d = dict(nv=len(V), m1=float(r2[:, 0].mean()), se1=float(r2[:, 0].std()/np.sqrt(len(V))), sd1=float(r2[:, 0].std()),
                 max1=float(r2[:, 0].max()), mrest=float(r2[:, 1:].mean()), vals1=r2[:, 0].tolist() if len(V) <= 64 else None)
        cell['fam'][name] = d
        print(f"   {name:22s} {len(V):4d} {d['m1']:12.2f} {d['se1']:6.2f} {d['max1']:7.1f} {d['mrest']:15.2f}")
    out['cells'].append(cell)
    with open(args.out, 'w') as fh:
        json.dump(out, fh)
print(f"done {time.time()-t0:.0f}s")
