# exponent bookkeeping (units of log Z) at the paper geometry used by Solver D
y, M = .45971, .83288           # Y = Z^y (a-length), Q = Z^M (m-length)
D = 2*M - 1                      # reflected mu-length Q^2/Z
Pa = 2*y - M; K = 2*M - y
print(f"P_a={Pa:.4f} K={K:.4f} D={D:.4f}  target E_a=Q^2 -> {2*M:.4f}")
print(f"one-step m-Poisson generic: Q^2 P_a^(1/6) -> {2*M+Pa/6:.4f} (c=1/6); side term Q^(4/3)Y -> {4*M/3+y:.4f} (<2M iff y<2M/3={2*M/3:.4f})")
print(f"HB/BGL a-form: Q(QY)^(2/3) -> {M+2*(M+y)/3:.4f} (c=1/3)")
print(f"DR route, D2 via GRH per mu: Q^2 D -> {2*M+D:.4f}")
print(f"twisted-Kubota pointwise + Lindelof: Y Q D Q^0 -> {y+M+D:.4f}")
print(f"D-threshold for GRH dyadic pieces to beat c=1/6: D_j < P_a^(1/6) = Z^{Pa/6:.4f}")
