"""Direct non-correlation test (split-prime proxy over Z[omega]).

  A_m(h') = sum_{a ~ Y} gt(a) chi_a(h')^P conj(chi_a(m))^P        (untwisted Kubota sum; carries the Patterson pole)
  C(h')   = sum_{m ~ Q} v(m) A_m(h')

n = 6: P = 1, gt = normalized sextic Gauss sum.   n = 3: P = 2, gt = normalized cubic Gauss sum (validation).
Probe vector T(m) = abar(m)^2 * gamma(m) * R(m),  gamma = gt_n(m),  R(m) = sum_{N(mu)<=D} rho_mu (mu/m)_2,
rho_mu = gt_3(mu, 'chk') alpha(mu) N(mu)^{-1/2}  (same proxy as ../sieve3.py).

Statistics (per h', pooled over h'):
  r(h')   = |C(h')|^2 / sum_m |v_m|^2 |A_m(h')|^2      (= 1 in expectation for random-phase v; ~ #m * rho^2 * f for a main term)
  rho^2   = Re(C0 conj C1) / (||v||^2 Re<A0,A1>)        (A0, A1 = Kubota sums over a random half-split of the a's;
                                                        the common part of A0, A1 is the pole; unbiased estimate of
                                                        |<v, conj pole>|^2 / (||v||^2 ||pole||^2), null mean 0)
  pole fraction f = 4 Re<A0,A1> / ||A0+A1||^2.
"""
import sys, json, time, argparse
import numpy as np
from lib import prime_elements, products, symbol_table, gauss_composite, char_matrix, W3, Z6

ap = argparse.ArgumentParser()
ap.add_argument('--n', type=int, default=6)
ap.add_argument('--Y', type=int, default=20000)
ap.add_argument('--Q', type=int, default=2000)
ap.add_argument('--D', type=int, default=0)
ap.add_argument('--conv', default='lam')
ap.add_argument('--nh', type=int, default=24)
ap.add_argument('--HN', type=int, default=200)
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--nrand', type=int, default=200)
ap.add_argument('--out', default='')
args = ap.parse_args()
n, Y, Q = args.n, args.Y, args.Q
D = args.D or int(Q**0.8)
P = 1 if n == 6 else 2
rng = np.random.default_rng(args.seed)
t0 = time.time()

PMAX = 2*max(Y, Q)
pel = prime_elements(PMAX)
A_set = products(pel, Y, 2*Y)
M_set = products(pel, Q, 2*Q)
mu_set = products(pel, 2, D)
# h' list: 1 and primary elements with N(h) coprime to 6, N <= HN (split and inert)
H = [(1, 0)]
R_ = int(np.sqrt(4*args.HN/3)) + 2
cand = []
for x in range(-R_, R_+1):
    for y in range(-R_, R_+1):
        N = x*x - x*y + y*y
        if 1 < N <= args.HN and N % 2 and N % 3 and x % 3 == 2 and y % 3 == 0:
            cand.append((N, x, y))
cand.sort()
for (N, x, y) in cand[:args.nh-1]:
    H.append((x, y))
nM, nH, nMu = len(M_set), len(H), len(mu_set)
TA = [m[0][0] for m in M_set] + [h[0] for h in H] + [u[0][0] for u in mu_set]
TB = [m[0][1] for m in M_set] + [h[1] for h in H] + [u[0][1] for u in mu_set]
colM = np.arange(nM); colH = np.arange(nM, nM+nH); colMu = np.arange(nM+nH, nM+nH+nMu)
imax = len(pel)
print(f"n={n} Y={Y} Q={Q} D={D} conv={args.conv}: #a={len(A_set)} #m={nM} #h={nH} #mu={nMu} #primes={imax}", flush=True)
iD = sum(1 for Pp in pel if Pp['p'] <= D)
E, Zr, G = symbol_table(pel, imax, TA, TB, specs=[(P, args.conv, imax), (2, 'chk', max(iD, 1))])
print(f"  symbol table {time.time()-t0:.0f}s", flush=True)
Gp = G[(P, args.conv)]
dummy = imax

# ---- Kubota sums over two random halves of the a's
ga = np.array([gauss_composite(pel, idx, Gp, P) for (_, idx, _) in A_set])
half = rng.integers(0, 2, len(A_set))
A0 = np.zeros((nH, nM), complex); A1 = np.zeros((nH, nM), complex)
AR = np.zeros((nH, nM), complex)  # random-phase-coefficient a-sum (no pole) for reference
cols = np.concatenate([colM, colH])
CH = 400
for s in range(0, len(A_set), CH):
    chunk = A_set[s:s+CH]
    X = char_matrix([c[1] for c in chunk], E, Zr, cols, P, dummy)
    Xm = np.conj(X[:, :nM]); Xh = X[:, nM:]
    g = ga[s:s+CH]; hb = half[s:s+CH]
    A0 += (Xh*(g*(hb == 0))[:, None]).T @ Xm
    A1 += (Xh*(g*(hb == 1))[:, None]).T @ Xm
    AR += (Xh*(np.abs(g)*np.exp(2j*np.pi*rng.random(len(g))))[:, None]).T @ Xm
A = A0 + A1
print(f"  Kubota sums {time.time()-t0:.0f}s; |gt(a)|^2 mean {np.mean(np.abs(ga)**2):.3f}", flush=True)

# ---- m-side vectors
mv = np.array([m[0][0] + m[0][1]*W3 for m in M_set]); ab2 = (np.conj(mv)/np.abs(mv))**2
Xm_t = char_matrix([m[1] for m in M_set], E, Zr, np.concatenate([colH, colMu]), 1, dummy)  # chi_m(h), chi_m(mu)
chi_m_h = Xm_t[:, :nH]             # (#m, #h)
quad = Xm_t[:, nH:]**3             # (mu/m)_2   (#m, #mu)
gam = np.array([gauss_composite(pel, idx, Gp, P) for (_, idx, _) in M_set])
G3c = G[(2, 'chk')]
rho = np.array([gauss_composite(pel, idx, G3c, 2)*((u[0]+u[1]*W3)/abs(u[0]+u[1]*W3))/np.sqrt(N) for (u, idx, N) in mu_set])
Rm = quad @ rho
vecs = {
    'T=ab2*g*R': ab2*gam*Rm,
    'ab2*conj(g)*R': ab2*np.conj(gam)*Rm,
    'g*R (no ab2)': gam*Rm,
    'conj(g)*R (no ab2)': np.conj(gam)*Rm,
    'g': gam, 'conj(g)': np.conj(gam),
    'ab2*g': ab2*gam, 'ab2*conj(g)': ab2*np.conj(gam),
    'one': np.ones(nM, complex), 'ab2': ab2,
    'R': Rm.astype(complex), 'ab2*R': ab2*Rm,
}


def stats(V, AA0, AA1):
    """V: (#h, #m) per-h vectors. returns r per h, pooled rho^2, pooled r."""
    AA = AA0 + AA1
    C = np.einsum('hm,hm->h', V, AA); C0 = np.einsum('hm,hm->h', V, AA0); C1 = np.einsum('hm,hm->h', V, AA1)
    den = np.einsum('hm,hm->h', np.abs(V)**2, np.abs(AA)**2)
    r = np.abs(C)**2/den
    v2 = np.sum(np.abs(V)**2, axis=1)
    cross = np.real(np.einsum('hm,hm->h', AA0, np.conj(AA1)))
    num = np.real(C0*np.conj(C1)); dd = v2*cross
    # rho2 reported as: h'=1 value (index 0) and pooled over h' != 1
    return r, (num[0]/dd[0], np.sum(num[1:])/np.sum(dd[1:]))



def pack(r, rho):
    return dict(r_h1=float(r[0]), r_rest=float(r[1:].mean()), r_max=float(r.max()), rho2_h1=float(rho[0]), rho2_rest=float(rho[1]))


res = {}
for k, v in vecs.items():
    r, rho = stats(np.broadcast_to(v, (nH, nM)), A0, A1)
    res[k] = pack(r, rho)
# Gauss controls aligned with chi_m(h)^j (Patterson index shift), best j per h (look-elsewhere factor 6)
for nm, base in (('g*chi_m(h)^j [best j]', gam), ('conj(g)*chi_m(h)^j [best j]', np.conj(gam)),
                 ('ab2*g*chi_m(h)^j [best j]', ab2*gam), ('ab2*conj(g)*chi_m(h)^j [best j]', ab2*np.conj(gam))):
    rr = [stats(base[None, :]*(chi_m_h.T**j), A0, A1) for j in range(6)]
    R = np.array([x[0] for x in rr]); jb = R.argmax(axis=0)
    rho = (rr[jb[0]][1][0], float(np.mean([rr[j][1][1] for j in range(6)])))
    res[nm] = pack(R.max(axis=0), rho)
# data-driven pole controls: v = conj(A1) (other half of the a's), correlate with A0 only
for nm, V in (('POLE: conj(A_other)', np.conj(A1)), ('ab2*conj(A_other)', ab2[None, :]*np.conj(A1)),
              ('|A_other|*T-phase', np.abs(A1)*np.exp(1j*np.angle(vecs['T=ab2*g*R']))[None, :])):
    C = np.einsum('hm,hm->h', V, A0); den = np.einsum('hm,hm->h', np.abs(V)**2, np.abs(A0)**2)
    r = np.abs(C)**2/den
    res[nm] = pack(r, (np.nan, np.nan))
# random-phase nulls with |v| = |T|
absT = np.abs(vecs['T=ab2*g*R'])
nul = []
for t in range(args.nrand):
    v = absT*np.exp(2j*np.pi*rng.random(nM))
    nul.append(pack(*stats(np.broadcast_to(v, (nH, nM)), A0, A1)))
null = {k: (float(np.mean([d[k] for d in nul])), float(np.std([d[k] for d in nul]))) for k in nul[0]}
res['RANDOM (null mean)'] = {k: v[0] for k, v in null.items()}
res['RANDOM (null sd)'] = {k: v[1] for k, v in null.items()}
# pole diagnostics
cross = np.real(np.einsum('hm,hm->h', A0, np.conj(A1)))
f = 4*cross/np.sum(np.abs(A)**2, axis=1)
zsd = np.sqrt(np.einsum('hm,hm->h', np.abs(A0)**2, np.abs(A1)**2)/2)
z = cross/zsd
energy_ratio = np.sum(np.abs(A)**2, axis=1)/np.sum(np.abs(AR)**2, axis=1)
summary = dict(n=n, Y=Y, Q=Q, D=D, conv=args.conv, nA=len(A_set), nM=nM, nH=nH, nMu=nMu, H=H,
               pole_frac=f.tolist(), pole_z=z.tolist(), energy_vs_random=energy_ratio.tolist(),
               res=res, secs=time.time()-t0)
print(f"  pole fraction f per h' (first 8): {np.round(f[:8], 3).tolist()}  z: {np.round(z[:8], 1).tolist()}")
print(f"  ||A||^2/||A_randphase-a||^2 (h'=1): {energy_ratio[0]:.2f};  #m*f(h'=1) = {nM*f[0]:.1f}  (r expected for a vector fully aligned with the pole)")
print(f"  {'vector':34s} {'r(h=1)':>8s} {'r(h!=1)':>8s} {'r_max':>8s} {'#m*rho2(1)':>11s} {'#m*rho2(!=1)':>12s}")
for k, d in res.items():
    print(f"  {k:34s} {d['r_h1']:8.2f} {d['r_rest']:8.2f} {d['r_max']:8.2f} {nM*d['rho2_h1']:11.2f} {nM*d['rho2_rest']:12.2f}")
print(f"  ({time.time()-t0:.0f}s)")
if args.out:
    with open(args.out, 'w') as fh:
        json.dump(summary, fh, indent=1)
