"""Exponent bookkeeping for solver D (a-twisted energy). log_Z units."""
import sys; sys.path.insert(0, '..'); sys.path.insert(0, '.')
l=0.16712007850681
k0=(5-3*l)/6; d0=(25-60*l+27*l*l)/(6*(5+4*l-3*l*l))
C0=(6*k0-1)/(15*k0-2); V0=5/6-d0
R0=-1+C0*(V0**2-d0**2*C0)/(2*(V0+d0*C0)**2)
h0=-((1+l)/2+l/(6*k0))/(3*R0/4-1/(24*k0))
x=1+l-h0; y=h0-2*l; M=x+y; Pa=y-x
print(f"x={x:.5f} y={y:.5f} l={l:.5f} M={M:.5f} Pa={Pa:.5f} 1+3l-2M={1+3*l-2*M:.4f}")
K=2*M-y; D=2*M-1; Hh=1+y-M   # K=Q^2/Y, D=Q^2/Z reflected length, H'=ZY/Q
print(f"K=Q^2/Y:{K:.4f} D=Q^2/Z:{D:.4f} H'=ZY/Q:{Hh:.4f}")
tgt=2*M
rows={
 "target E_a<<Q^2": tgt,
 "a-diagonal  Y Q^{7/6}": y+7*M/6,
 "exceptional k=r^6 (no Patterson pole), K^{1/6}D^2*Y/Q": K/6+2*D+y-M,
 "full-mass exceptional (Gauss-phase vector) Y K^{1/6} Q": y+K/6+M,
 "HB a-form Q(YQ)^{2/3}": M+2*(y+M)/3,
 "classical a-form Q(Q+Y^2)": M+max(M,2*y),
 "HB k-form Y(KQ)^{2/3}": y+2*(K+M)/3,
 "label route (B2.42, V=y) column branch per f": 2*M+y-1,
 "trivial Patterson main term (no cancel in mu)": y+4*M/3-y/6+D/2,
}
for k_,v in rows.items(): print(f"{k_:55s} {v:.4f}  excess {v-tgt:+.4f}")
for th in (1, 1/3, 1/6, 1/8, 0):
    print(f"E_a<<Q^2 P_a^{th:.3f}: probe exponent x/2+th*Pa/2 = {x/2+th*Pa/2:.4f}")
print("signal-consistency: c=0 bound compatible iff beta <= 5/6+x/6-l/6 =", 5/6+x/6-l/6)
