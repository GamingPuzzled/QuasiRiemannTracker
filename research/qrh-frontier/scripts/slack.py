"""Numeric model of the quartic QRH exponent optimization (cycle25 paper, B0.3-B0.7).

All parameters are exponents base Z. Constants default to the paper's values and can be
perturbed to measure the sensitivity of the final boundary b0.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, minimize

DEF = dict(
    V0=5 / 6,        # long inverse envelope constant, V = V0 - delta   (B0.40)
    holder=0.25,     # numerator exponent 1/4 + 3R/4                    (B0.41)
    capdiv=6.0,      # plain capacity (1-2m)/(capdiv*kappa)             (B0.38, B0.54)
    pmax=0.5,        # amplitude cap g_i <= pmax*delta                  (B0.37)
    z0=0.34,
    Mslack=0.0,      # allow M + l <= 1 + Mslack                        (B0.26/B0.27)
    probe_div=12.0,  # (y-x)/probe_div in (B0.27)
    supply=4 / 21,
)


def C_k(p, k, c):
    return 2 * (1 - p) * (1 - p / (3 * k)) / (3 - p - 2 * p / (3 * k))


def R_k(delta, p, k, c):
    V = c["V0"] - delta
    C = C_k(p, k, c)
    if V <= 0:
        return 1 - delta + delta * C * 0.5  # t = 3/2 endpoint
    return 1 - delta + delta * C * V / (2 * (V + delta * C))


def F(d, delta, p, beta, x, y, l, c):
    """(B0.55) + beta: exponent of retained nonprincipal rows at U = Z^d."""
    k = 2 * beta - 1
    h = 1 + l - x
    a = (1 + delta) / 2
    q = delta * p
    R = R_k(delta, p, k, c)
    cap = max(0.0, min(d, 2 * y - d))
    return (a + h * (c["z0"] - 1 / 6) - y / 2 - l / 2 + l * q
            + d * (c["holder"] + (1 - c["holder"]) * R - c["z0"])
            - q * 2 / (2 * c["capdiv"] * k) * cap)


def probe(x, y, l, c):
    M = x + y
    k = 5 / 6 - x / 3 - l / 6
    return x / 2 + (y - x) / c["probe_div"] + max(0.0, 1 + 3 * l - 2 * M) / 8 + k


DG = np.linspace(0.0, 1.0, 81)
PG = np.linspace(0.0, 1.0, 41)


def worst(beta, x, y, l, c):
    """max over delta, p, d of F - beta (d at endpoints by convexity, plus kink)."""
    k = 2 * beta - 1
    h = 1 + l - x
    best = -9.0
    dmax = 2 * beta - 1
    for delta in np.linspace(0.02, dmax, 120):
        for p in PG * c["pmax"]:
            for d in (0.005, h, min(h, max(0.005, y))):
                best = max(best, F(d, delta, p, beta, x, y, l, c) - beta)
    # floor bin: a=.51, R=1, q<=.01
    best = max(best, 0.51 + h * (c["z0"] - 1 / 6) - y / 2 - l / 2 + l * 0.01
               + h * (1 - c["z0"]) - beta)
    return best


def feasible(x, y, l, c):
    h = 1 + l - x
    return (l < x <= y and y - 2 * x + l < 0 and x + y + l <= 1 + c["Mslack"] + 1e-12
            and l / h > c["supply"] and y < 4 * l)


def bound(params, c, refine=True):
    x, y, l = params
    if not feasible(x, y, l, c):
        return 9.0
    b = probe(x, y, l, c)
    w = worst(b + 1e-9, x, y, l, c)
    return b + max(0.0, w) * 10  # penalty if nonprincipal rows fail


def solve(c=None, start=(0.37317, 0.45971, 0.16712)):
    c = {**DEF, **(c or {})}
    # paper geometry: y = 1 - l - x (+ slack). Search over (x, l, s) with s in [0, Mslack]
    def obj(v):
        x, l, s = v
        y = 1 - l - x + s
        return bound((x, y, l), c)
    best = None
    for s0 in ([0.0] if c["Mslack"] == 0 else [0.0, c["Mslack"] / 2, c["Mslack"]]):
        r = minimize(obj, [start[0], start[2], s0], method="Nelder-Mead",
                     options=dict(xatol=1e-9, fatol=1e-12, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return best


if __name__ == "__main__":
    r = solve()
    print("baseline", r.fun, r.x)
