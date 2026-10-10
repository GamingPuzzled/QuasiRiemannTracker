"""Eisenstein-integer arithmetic and sextic residue symbols for split primes (proxy experiment)."""
from __future__ import annotations
import numpy as np
from sympy import primerange, sqrt_mod

W = np.exp(2j*np.pi/3)  # omega

def norm(a: int, b: int) -> int:
    return a*a - a*b + b*b

def is_primary(a: int, b: int) -> bool:
    # convention: a ≡ 2 mod 3, b ≡ 0 mod 3  (associate unique up to sign; fix sign below)
    return a % 3 == 2 and b % 3 == 0

UNITS = [(1, 0), (0, 1), (-1, -1), (-1, 0), (0, -1), (1, 1)]  # 1, w, w^2=-1-w, -1, -w, -w^2=1+w

def mul(x, y):
    a, b = x; c, d = y  # (a+bw)(c+dw) = ac + (ad+bc)w + bd w^2, w^2 = -1-w
    return (a*c - b*d, a*d + b*c - b*d)

def primary_associate(a, b):
    for u in UNITS:
        x = mul((a, b), u)
        if is_primary(*x):
            return x
    return None

def split_primes(P):
    out = []
    for p in primerange(7, P + 1):
        if p % 3 != 1:
            continue
        # find a,b with a^2-ab+b^2 = p
        found = None
        for b in range(1, int(2*np.sqrt(p/3)) + 2):
            # a^2 - b a + b^2 - p = 0
            disc = 4*p - 3*b*b
            if disc < 0: break
            s = int(round(np.sqrt(disc)))
            if s*s == disc and (b + s) % 2 == 0:
                found = ((b + s)//2, b); break
        a, b = found
        for pi in (primary_associate(a, b), primary_associate(*mul((a, b), (0, 0)) if False else (a - b, -b))):
            pass
        pi1 = primary_associate(a, b)
        # conjugate of a+bw is a+b w^2 = (a-b) - b w
        pi2 = primary_associate(a - b, -b)
        for pi in (pi1, pi2):
            # residue map: w -> r mod p with pi ≡ 0
            r = None
            for cand in sqrt_mod(-3, p, all_roots=True):
                rr = (cand - 1) * pow(2, -1, p) % p  # roots of x^2+x+1
                if (pi[0] + pi[1]*rr) % p == 0:
                    r = rr; break
            out.append(dict(p=p, pi=pi, r=r))
    return out

def sixth_root_table(p, r):
    """map residues z^( (p-1)/6 ) to exponent j with value (-w^2)^j  (primitive 6th root e^{i pi/3} = -w^2)."""
    z6 = (-(r*r)) % p
    tab = {}
    v = 1
    for j in range(6):
        tab[v] = j; v = v*z6 % p
    return tab
