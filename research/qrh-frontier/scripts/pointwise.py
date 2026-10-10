import warnings; warnings.filterwarnings("ignore")
import numpy as np, slack2
from slack2 import R_k, DEF
x, l = 0.373172481118674, 0.16712007850681
y = 1-l-x; h = 1+l-x; b = 11/12-l/4; k = 2*b-1
for delta in (0.3, 0.4, 0.469344858698143, 0.55, 0.7):
    p = 0.5; q = delta*p; R = R_k(delta, p, k, DEF); d = h; m = y/d
    holder = d*(0.25+0.75*R) - q*min(d,2*y-d)/(12*k)
    point = d*(R + delta*min(1-m, m)/2)
    print(f"delta={delta:.3f} R={R:.5f} holder={holder:.6f} pointwise={point:.6f} diff={point-holder:+.2e}")
