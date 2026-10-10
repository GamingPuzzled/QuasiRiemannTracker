# Feasibility of slot schedule in the "framework" model:
# width M runs from 1 down; cumulative amplifier gain G(M)=(1-M)/rho (rho=6);
# removed slot length r(M) nondecreasing; allowance L=G-kappa*r >= 0;
# diagonal: Aeff - r - M <= L, Aeff = 2*min(n, M-n) + z (reflection allowed);
# terminal: r=z reached at some M* >= 0.
import numpy as np
def feasible(n,z,kappa,rho=6.0,N=4001):
    Ms=np.linspace(1,0,N); r=0.0
    for M in Ms:
        G=(1-M)/rho
        Aeff=2*max(0.0,min(n,M-n))+z
        # need (1-kappa) r >= Aeff - M - G  and kappa r <= G ; greedily keep r minimal
        need=(Aeff-M-G)/(1-kappa)
        r=max(r,need,0.0); r=min(r,z) if need<=z else r
        if r>z+1e-12: return False
        if kappa*r>G+1e-12:
            return False
        if Aeff-r-M>G-kappa*r+1e-12: return False
        # removing more early never helps except at terminal; try terminal
        if kappa*z<=G+1e-12 and r<=z: return True
    return False
for n in [0.30,0.42,0.48]:
  for kappa in [0.7,0.7498]:
    lo,hi=0,0.5
    for _ in range(50):
        mid=(lo+hi)/2
        (lo,hi)=(mid,hi) if feasible(n,mid,kappa) else (lo,mid)
    print(n,kappa,"zmax=%.5f"%lo,"(1-2n)/(6k)=%.5f"%((1-2*n)/(6*kappa)))
