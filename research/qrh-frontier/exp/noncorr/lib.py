"""Fast split-prime machinery over Z[omega] for the non-correlation test (vectorized per rational prime).

Conventions follow ../eis.py and ../theta_id.py:
  * prime elements pi are primary (a = 2 mod 3, b = 0 mod 3), omega = r mod pi;
  * chi_pi(x) = (x/pi)_6 = zeta6^j where x^((p-1)/6) = (-r^2)^j mod p, zeta6 = e^{i pi/3};
  * normalized Gauss sum of chi_pi^POW with additive character
        'lam': e(Tr(x/(lambda pi))),  lambda = 1+2 omega  (paper's convention)
        'chk': e(Tr(x/pi))
  * composite moduli: g(ab) = chi_a(b) chi_b(a) g(a) g(b) (twisted multiplicativity, powers POW).
"""
from __future__ import annotations
import os, sys
import numpy as np
from sympy import primitive_root
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from eis import split_primes, sixth_root_table, mul  # noqa: E402

W3 = np.exp(2j*np.pi/3)
LAM = 1 + 2*W3
Z6 = np.exp(1j*np.pi/3*np.arange(6))


def prime_elements(PMAX):
    pel = sorted(split_primes(PMAX), key=lambda P: (P["p"], P["pi"]))
    return pel


def products(pel, lo, hi):
    """squarefree products of distinct prime elements with norm in [lo, hi]: list of (value, idx tuple, norm)."""
    out = []
    norms = [P["p"] for P in pel]

    def rec(start, cur, nrm, idx):
        if lo <= nrm <= hi:
            out.append((cur, tuple(idx), nrm))
        for i in range(start, len(pel)):
            if nrm*norms[i] > hi:
                break
            idx.append(i)
            rec(i+1, mul(cur, pel[i]["pi"]), nrm*norms[i], idx)
            idx.pop()
    rec(0, (1, 0), 1, [])
    return out


def _powers(g, p):
    B = int(np.ceil(np.sqrt(p - 1)))
    small = np.empty(B, dtype=np.int64); v = 1
    for k in range(B):
        small[k] = v; v = v*g % p
    gB = v
    big = np.empty(B, dtype=np.int64); v = 1
    for j in range(B):
        big[j] = v; v = v*gB % p
    pw = (big[:, None]*small[None, :] % p).ravel()[:p-1]
    return pw


def symbol_table(pel, imax, TA, TB, POWS=(1, 2), convs=('lam', 'chk'), specs=None):
    """For prime elements pel[:imax]: exponent table E[i, t] in 0..5 of chi_pi(target t), mask Zr[i,t] (target divisible),
    and normalized Gauss sums G[(POW, conv)][i].  specs: optional list of (POW, conv, ilimit) restricting which
    Gauss sums are computed (for i < ilimit only)."""
    T = len(TA)
    E = np.zeros((imax + 1, T), dtype=np.int8)   # extra dummy row of zeros (padding)
    Zr = np.zeros((imax + 1, T), dtype=bool)
    if specs is None:
        specs = [(pw, cv, imax) for pw in POWS for cv in convs]
    G = {(pw, cv): np.zeros(imax, complex) for (pw, cv, _) in specs}
    TA = np.asarray(TA, dtype=np.int64); TB = np.asarray(TB, dtype=np.int64)
    cache = {}
    for i in range(imax):
        P = pel[i]; p, r = P["p"], P["r"]
        if p not in cache:
            g = primitive_root(p)
            pw = _powers(g, p)
            dlog = np.empty(p, dtype=np.int64); dlog[pw] = np.arange(p-1)
            cache.clear(); cache[p] = (g, pw, dlog)
        g, pw, dlog = cache[p]
        tab = sixth_root_table(p, r); e = (p-1)//6
        c = tab[pow(int(g), e, p)]
        res = (TA % p + (TB % p)*r) % p
        z = res == 0
        E[i] = np.where(z, 0, (c*dlog[np.where(z, 1, res)]) % 6)
        Zr[i] = z
        piv = P["pi"][0] + P["pi"][1]*W3
        k = np.arange(p-1)
        addc = {}
        for (POW, cv, lim) in specs:
            if i >= lim:
                continue
            if cv not in addc:
                tr = (2*np.real(LAM*piv)/(3*p)) if cv == 'lam' else (2*np.real(np.conj(piv))/p)
                addc[cv] = np.mod(pw*tr, 1.0)
            ph = (POW*c*k % 6)/6.0 + addc[cv]
            G[(POW, cv)][i] = np.exp(2j*np.pi*ph).sum()/np.sqrt(p)
    return E, Zr, G


def pair_expo(pel, i, j):
    """exponent of chi_{pi_i}(pi_j) (or -1 if divisible)."""
    P = pel[i]; p, r = P["p"], P["r"]
    x, y = pel[j]["pi"]
    res = (x % p + (y % p)*r) % p
    if res == 0:
        return -1
    return sixth_root_table(p, r)[pow(res, (p-1)//6, p)]


_pair_cache = {}


def gauss_composite(pel, idx, Gp, POW):
    """normalized Gauss sum for squarefree product of pel[idx], given per-prime normalized sums Gp (array)."""
    g = 1+0j
    for i in idx:
        g *= Gp[i]
    for a in range(len(idx)):
        for b in range(a+1, len(idx)):
            i, j = idx[a], idx[b]
            key = (i, j)
            if key not in _pair_cache:
                e1 = pair_expo(pel, i, j); e2 = pair_expo(pel, j, i)
                _pair_cache[key] = None if (e1 < 0 or e2 < 0) else (e1 + e2)
            s = _pair_cache[key]
            if s is None:
                return 0j
            g *= Z6[(POW*s) % 6]
    return g


def char_matrix(idx_list, E, Zr, cols, POW, dummy):
    """chi_x(target)^POW for each modulus x in idx_list (list of prime-index tuples), targets restricted to cols.
    Returns (len(idx_list), len(cols)) complex array."""
    kmax = max(len(t) for t in idx_list)
    I = np.full((len(idx_list), kmax), dummy, dtype=np.int64)
    for n, t in enumerate(idx_list):
        I[n, :len(t)] = t
    Ec = E[:, cols]; Zc = Zr[:, cols]
    s = np.zeros((len(idx_list), len(cols)), dtype=np.int16)
    z = np.zeros((len(idx_list), len(cols)), dtype=bool)
    for k in range(kmax):
        s += Ec[I[:, k]]; z |= Zc[I[:, k]]
    return np.where(z, 0, Z6[(POW*s) % 6])
