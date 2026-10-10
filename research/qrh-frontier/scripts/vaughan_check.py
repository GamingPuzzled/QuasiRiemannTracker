# Vaughan/mollifier route: M_k(D) ~ sum_{X<B<=XH} F_k(B) M_k(D/B), F=sum_{Delta<=X} M(Delta)S(Y), B=Delta*Y.
# Logs base H. D=H^x. Tool T (hypothetical, stronger than framework): sum_k |F_k(B)|^2 << H^{max(1, l)},
# l = delta + min(y,1-y) (length after optimal reflection of the plain factor). Pointwise: |M(A)|^2<<A^kappa.
# Also option: UM on M(A) (if a<1) times pointwise sup of F: kappa*(delta+min(y,1-y)).
import numpy as np
def worst(x,xi,k,n=161):
    w=-1e9
    for d in np.linspace(0,xi,n):
        for y in np.linspace(0,1,n):
            b=d+y
            if b<=xi or b>xi+1: continue
            a=x-b
            if a<0: continue
            l=d+min(y,1-y)
            o1=max(1,l)+k*a                       # T on F, pointwise on M(A)
            o2=(1 if a<1 else (1+5*a)/6)+k*l      # moment on M(A), pointwise on F factors
            w=max(w,min(o1,o2))
    return w
for k in (0.70,0.75):
  for x in (1.05,1.1,1.2):
    best=min((worst(x,xi,k),xi) for xi in np.linspace(x-1,1.0,41))
    lam=(best[0]-1)/(x-1)
    print(f"kappa={k} x={x}: best exponent {best[0]:.4f} at xi={best[1]:.3f}; implied lambda={lam:.3f}; (B0.9)-equivalent {(1+5*x)/6:.4f}")
