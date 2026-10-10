import warnings; warnings.filterwarnings("ignore")
import slack2
from slack2 import solve

def energy(M, N, z):
    row = M + max(0, N + 2*z - 2*M)/4
    Td = 2*M + 2*z - N
    col = Td + z - min(Td, z, (Td+z)/3) if Td > M else -9
    return max(row, col)

def probe_general(x, y, l, c):
    M = x + y
    # Cauchy-Schwarz: Q^{-1/2} * (X (Y/X)^{1/6})^{1/2} * Z^{E/2}, Q = XY
    return -M/2 + (x + (y-x)/6)/2 + energy(M, 1+l, l)/2 + 5/6 - x/3 - l/6

slack2.probe = probe_general
import numpy as np
print("check base", solve()[0])
for S in (0.01, 0.03):
    print(S, solve({"Mslack": S}, s_list=list(np.linspace(0, S, 4))))
