"""Identify the Patterson residue structure: A_m = sum_{a~Y} gt(chi_a^p) * conj(chi_a(m))^p, p=2 (cubic) or 1 (sextic).
Measure correlation of A_m with candidate Gauss-sum structures in m (split-prime proxy)."""
import sys, numpy as np
from eis import split_primes, sixth_root_table, mul
Y = int(sys.argv[1]); Q = int(sys.argv[2]); POW = int(sys.argv[3])
rng = np.random.default_rng(0)
pel = sorted(split_primes(2*max(Y, Q)), key=lambda P: P["p"])
Z6 = np.exp(1j*np.pi/3*np.arange(6)); W3 = np.exp(2j*np.pi/3); LAM = 1+2*W3
def products(lo, hi):
    out = []
    def rec(start, cur, nrm, idx):
        if lo <= nrm <= hi: out.append((cur, tuple(idx), nrm))
        for i in range(start, len(pel)):
            if nrm*pel[i]["p"] > hi: break
            rec(i+1, mul(cur, pel[i]["pi"]), nrm*pel[i]["p"], idx+[i])
    rec(0, (1, 0), 1, []); return out
A_set = products(Y, 2*Y); M_set = products(Q, 2*Q)
def expo_prime(i, xs):  # (x/pi_i)_6 exponent for elements xs=(a,b) arrays; -1 if 0
    P = pel[i]; p, r = P["p"], P["r"]; tab = sixth_root_table(p, r); e = (p-1)//6
    res = (xs[0] % p + (xs[1] % p)*r) % p
    return np.array([-1 if v == 0 else tab[pow(int(v), e, p)] for v in res])
Mx = (np.array([m[0][0] for m in M_set]), np.array([m[0][1] for m in M_set]))
used = sorted({i for _, idx, _ in A_set for i in idx} | {i for _, idx, _ in M_set for i in idx})
EM = {i: expo_prime(i, Mx) for i in used}
def chi_mod(idx, power):  # chi_a(m)^power for all m in M_set, a given by prime indices
    s = np.zeros(len(M_set), int); z = np.zeros(len(M_set), bool)
    for i in idx: s += EM[i]; z |= EM[i] < 0
    return np.where(z, 0, Z6[(power*s) % 6])
gcache = {}
def gauss_prime(i, power, add):  # normalized Gauss sum of (x/pi)_6^power, additive e(x/(lam pi)) if add=='lam' else e-check(x/pi)
    key = (i, power, add)
    if key in gcache: return gcache[key]
    P = pel[i]; p, r = P["p"], P["r"]; tab = sixth_root_table(p, r); e = (p-1)//6
    piv = P["pi"][0] + P["pi"][1]*W3
    tr = (2*np.real(LAM*piv)/(3*p)) if add == 'lam' else (2*np.real(np.conj(piv))/p)
    x = np.arange(1, p); ch = np.array([tab[pow(int(v), e, p)] for v in x])
    g = np.sum(Z6[(power*ch) % 6]*np.exp(2j*np.pi*x*tr))/np.sqrt(p)
    gcache[key] = g; return g
def pair(i, j, power):
    Pj = pel[j]["pi"]; e = expo_prime(i, (np.array([Pj[0]]), np.array([Pj[1]])))[0]
    return 0 if e < 0 else Z6[(power*e) % 6]
def gauss(idx, power, add):
    g = 1
    for i in idx: g *= gauss_prime(i, power, add)
    for a in range(len(idx)):
        for b in range(a+1, len(idx)):
            g *= pair(idx[a], idx[b], power)*pair(idx[b], idx[a], power)
    return g
for add in ('lam', 'chk'):
    A = np.zeros(len(M_set), complex)
    for (a, idx, n) in A_set:
        A += gauss(idx, POW, add)*np.conj(chi_mod(idx, POW))
    rand = np.zeros(len(M_set), complex)
    for (a, idx, n) in A_set:
        rand += np.exp(2j*np.pi*rng.random())*np.conj(chi_mod(idx, POW))
    print(f"[{add}] #a={len(A_set)} #m={len(M_set)}  mean|A|^2={np.mean(abs(A)**2):.1f}  (random-phase {np.mean(abs(rand)**2):.1f})")
    mv = np.array([m[0][0]+m[0][1]*W3 for m in M_set]); ang = mv/abs(mv)
    cands = {}
    for pw in range(1, 6):
        g = np.array([gauss(idx, pw, add) for (_, idx, _) in M_set])
        cands[f"g{pw}"] = g; cands[f"conj g{pw}"] = np.conj(g)
    def corr(u, v):
        return abs(np.vdot(v, u))/np.sqrt(np.vdot(u, u).real*np.vdot(v, v).real)
    res = sorted(((corr(A, c), k) for k, c in cands.items()), reverse=True)
    base = np.mean([corr(rand, c) for c in cands.values()])
    print("   top correlations of A_m with Gauss structures:", ", ".join(f"{k}:{v:.3f}" for v, k in res[:4]), f" | random-A baseline ~{base:.3f}, 1/sqrt(#m)={1/np.sqrt(len(M_set)):.3f}")
