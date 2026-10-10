import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import solve
def probe(x, y, l, c):
    M = x+y; Pa = y-x
    gram = max(0, Pa/6, 2*Pa-y)
    E = M + max(0, 1+3*l-2*M)/4
    return -M/2 + (x+gram)/2 + E/2 + 5/6 - x/3 - l/6
slack2.probe = probe
def feas(x, y, l, c):
    h = 1+l-x
    return (l < x and x+y+l <= 1+1e-12 and l/h > c["supply"] and y < 4*l and y > 0)
slack2.feasible = feas
print(solve())
# scan x for fixed s=0
for x in np.linspace(0.30, 0.60, 13):
    print(round(x,3), slack2.best_for(x, 0.0, {**slack2.DEF}))
