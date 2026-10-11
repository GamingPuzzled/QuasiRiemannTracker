"""(b),(c): can mixing / splitting beat the triple-point R or the numerator Holder bound?

1. Symbolic: mixed MP (B3.40, n1=m, n2=m') count line == average of the two pure MP lines.
   Any Holder combination of moment bounds gives a convex combination of count lines.
2. Two simultaneous truncations t1<t* <t2 (Lemma 1W applied twice): worst joint class.
3. Numerator split by level sigma (with pointwise cap and mixed line) == Holder value exactly.
"""
import sympy as sp
from mpmath import mp, mpf
from base import record, R, tstar
mp.dps = 30
# ---- 1. symbolic identity
d, q, k, m, mp_, s = sp.symbols('delta q kappa m mprime sigma', positive=True)
AP = 1 - 2*d*m - 2*q*(1-2*m)/(6*k)
num4 = 1 - 4*s - 2*q*(1-2*mp_)/(6*k)
mix = 1 - d*m - 2*s - 2*q*(1-m-mp_)/(6*k)
print("mix - (AP+num4)/2 =", sp.simplify(mix - (AP+num4)/2))
# ---- 2. two truncations
P = record(); dl, kk = P['d'], P['k']; p = mpf(1)/2; qq = dl*p; c3 = qq/(3*kk)
AI = lambda r: 1-qq-(dl-qq)*r
APf = lambda mm: 1-c3-2*(dl-c3)*mm
long_ = lambda t: (1+5*t)/6-dl*t
def short(t):
    r = (1-qq-(1-c3)+2*(dl-c3)*t)/((dl-qq)+2*(dl-c3)); return AI(r)
R0 = R(dl, p, kk); ts = tstar(dl, p, kk)
best = None
import itertools
grid = [ts + mpf(i)/200 for i in range(-30, 31)]
for t1, t2 in itertools.product(grid, grid):
    if t1 > t2: continue
    # classes: (S1|L1) x (S2|L2); joint count <= min of the two class bounds
    w = max(min(a, b) for a in (short(t1), long_(t1)) for b in (short(t2), long_(t2)))
    if best is None or w < best[0]: best = (w, t1, t2)
print("single-t R          =", mp.nstr(R0, 15))
print("two-truncation best =", mp.nstr(best[0], 15), " at t1,t2 =", mp.nstr(best[1], 6), mp.nstr(best[2], 6))
# ---- 3. numerator split
y, h = P['y'], P['h']; mpr = 1-y/h; zA = (1-2*mpr)/(6*kk); smax = dl*mpr/2
r = (1-qq-(1-c3)+2*(dl-c3)*ts)/((dl-qq)+2*(dl-c3)); mstar = ts-r
holder = mpf(1)/4 + 3*R0/4 - qq*zA/2
def val(sg):
    n4 = 1-4*sg-2*qq*zA
    mx = 1-dl*mstar-2*sg-2*qq*(1-mstar-mpr)/(6*kk)
    return min(R0, n4, mx) + sg
# maximize piecewise-linear concave function on [0,smax] by fine search + breakpoints
cands = [smax, (1-R0-2*qq*zA)/4] + [smax*i/20000 for i in range(20001)]
split = max(val(c) for c in cands if 0 <= c <= smax)
print("Holder numerator    =", mp.nstr(holder, 15))
print("split (incl. mixed) =", mp.nstr(split, 15), " diff =", mp.nstr(split-holder, 5))
print("pointwise-only      =", mp.nstr(R0+smax, 15), " (gap to Holder", mp.nstr(R0+smax-holder, 5), "U-units)")
