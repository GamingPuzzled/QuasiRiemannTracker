"""Sextic large-sieve experiment over Z[omega] (split-prime proxy).
E(K) = sum_{primary k, N(k)<=K, (N(k),6)=1} |sum_m b_m chi_m(k)|^2, compared with (#k + #m) * ||b||^2.
"""
import sys, itertools, time
import numpy as np
from eis import split_primes, sixth_root_table, mul, primary_associate, norm

Q = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 200000
D = int(sys.argv[3]) if len(sys.argv) > 3 else 60
rng = np.random.default_rng(1)
t0 = time.time()
primes = split_primes(2*Q)
# moduli: squarefree products of distinct prime elements (distinct ideals), N in [Q, 2Q]
pel = sorted(primes, key=lambda P: P["p"])
mods = []
def rec(start, cur, nrm, idx):
    if Q <= nrm <= 2*Q:
        mods.append((cur, tuple(idx)))
    for i in range(start, len(pel)):
        p = pel[i]["p"]
        if nrm*p > 2*Q: break
        rec(i+1, mul(cur, pel[i]["pi"]), nrm*p, idx+[i])
rec(0, (1, 0), 1, [])
mods = [(m, idx) for m, idx in mods]  # m value (a,b) as product of primary primes
used = sorted({i for _, idx in mods for i in idx})
# arguments k: primary, norm <= KMAX, norm coprime to 6
ks = []
R = int(np.sqrt(4*KMAX/3)) + 2
for a in range(-R, R+1):
    for b in range(-R, R+1):
        n = norm(a, b)
        if 0 < n <= KMAX and n % 2 and n % 3 and a % 3 == 2 and b % 3 == 0:
            ks.append((n, a, b))
ks.sort()
kn = np.array([k[0] for k in ks]); ka = np.array([k[1] for k in ks], dtype=np.int64); kb = np.array([k[2] for k in ks], dtype=np.int64)
print(f"Q={Q} #m={len(mods)} #k={len(ks)} primes used={len(used)}  ({time.time()-t0:.1f}s)", flush=True)
# exponent arrays per prime element
Z6 = np.exp(1j*np.pi/3*np.arange(6))
def expo(Pr, A, B):
    p, r = Pr["p"], Pr["r"]
    tab = sixth_root_table(p, r); e = (p-1)//6
    res = (A % p + (B % p) * r) % p
    out = np.full(len(A), -1, dtype=np.int8)
    nz = res != 0
    vals = np.array([pow(int(v), e, p) for v in res[nz]])
    out[nz] = np.array([tab[int(v)] for v in vals], dtype=np.int8)
    return out
def chi_matrix(A, B):
    E = {i: expo(pel[i], A, B) for i in used}
    X = np.empty((len(A), len(mods)), dtype=np.complex128)
    for j, (m, idx) in enumerate(mods):
        s = np.zeros(len(A), dtype=np.int16); z = np.zeros(len(A), dtype=bool)
        for i in idx:
            s += E[i]; z |= E[i] < 0
        X[:, j] = np.where(z, 0, Z6[s % 6])
    return X
# vectors b_m
mval = np.array([m[0] + m[1]*np.exp(2j*np.pi/3) for m, _ in mods])
alphabar2 = (np.conj(mval)/np.abs(mval))**2
mus = [(a, b) for (n, a, b) in ks if n <= D]
MA = np.array([u[0] for u in mus]); MB = np.array([u[1] for u in mus])
Xmu = chi_matrix(MA, MB)  # (#mu, #m): chi_m(mu)
quad = (Xmu**3).sum(axis=0)          # sum_mu (mu/m)_2
sext = Xmu.sum(axis=0)

LAM = 1 + 2*np.exp(2j*np.pi/3)
def gauss_prime(i, power):
    Pr = pel[i]; p, r = Pr["p"], Pr["r"]; tab = sixth_root_table(p, r); e = (p-1)//6
    piv = Pr["pi"][0] + Pr["pi"][1]*np.exp(2j*np.pi/3)
    tr = 2*np.real(LAM*piv)/(3*p)   # Tr(x/(lam pi)) = x*tr
    x = np.arange(1, p)
    ch = np.array([tab[pow(int(v), e, p)] for v in x])
    return np.sum(Z6[(power*ch) % 6] * np.exp(2j*np.pi*x*tr))
def chi_pair(i, j, power):  # chi_{pi_i}(pi_j)^power
    Pj = pel[j]["pi"]; e = expo(pel[i], np.array([Pj[0]]), np.array([Pj[1]]))[0]
    return 0 if e < 0 else Z6[(power*e) % 6]
def gauss_vec(power):
    gp = {i: gauss_prime(i, power) for i in used}
    out = []
    for m, idx in mods:
        g = 1
        for i in idx: g *= gp[i]
        for a in range(len(idx)):
            for b in range(a+1, len(idx)):
                g *= chi_pair(idx[a], idx[b], power) * chi_pair(idx[b], idx[a], power)
        out.append(np.conj(g)/np.sqrt(float(np.prod([pel[i]["p"] for i in idx]))))
    return np.array(out)
g6 = gauss_vec(1); g3 = gauss_vec(2)
m1 = chi_matrix(np.array([-1]), np.array([0]))[0]  # chi_m(-1)
g6s = g6*m1  # g(chi_m-bar) normalized, up to conj convention
POWER = int(sys.argv[4]) if len(sys.argv) > 4 else 1
print("gauss |g|^2 mean (should be 1):", np.mean(np.abs(g6)**2).round(4), np.mean(np.abs(g3)**2).round(4))
vecs = {
  "random": np.exp(2j*np.pi*rng.random(len(mods))),
  "theta-proxy ab2*quad": alphabar2*quad,
  "no-ab2 quad (control)": quad.astype(complex),
  "ab2*sextic": alphabar2*sext,
  "GAUSS (pos ctrl)": (g6s if POWER == 1 else g3),
  "GAUSS conj": np.conj(g6s if POWER == 1 else g3),
}
print("#mu =", len(mus), {k: f"{np.linalg.norm(v)**2:.3g}" for k, v in vecs.items()}, flush=True)
B = np.stack(list(vecs.values()), axis=1).astype(np.complex128)
bn2 = (np.abs(B)**2).sum(axis=0)
checkpoints = [K for K in [Q, 2*Q, 4*Q, 8*Q, 16*Q, 32*Q, 64*Q, 128*Q, 256*Q] if K <= KMAX]
cum = np.zeros(B.shape[1]); start = 0; ci = 0
CH = 20000
rows = []
while start < len(ks) and ci < len(checkpoints):
    end = min(start+CH, len(ks))
    X = chi_matrix(ka[start:end], kb[start:end])
    if POWER != 1: X = X**POWER
    vals = np.abs(X @ B)**2
    cs = np.cumsum(vals, axis=0)
    norms = kn[start:end]
    while ci < len(checkpoints) and (end == len(ks) or norms[-1] > checkpoints[ci]):
        K = checkpoints[ci]
        cnt = np.searchsorted(norms, K, side="right")
        tot = cum + (cs[cnt-1] if cnt > 0 else 0)
        nk = start + cnt
        ratio = tot / ((nk + len(mods)) * bn2)
        rows.append((K, nk, ratio)); print(f"K={K:>8} (K/Q={K/Q:>5.0f}) #k={nk:>7}  ratio:", "  ".join(f"{n}={r:.3f}" for n, r in zip(vecs, ratio)), flush=True)
        if end < len(ks) and norms[-1] <= checkpoints[ci]: break
        ci += 1
    cum += cs[-1]; start = end
print(f"done {time.time()-t0:.0f}s")
