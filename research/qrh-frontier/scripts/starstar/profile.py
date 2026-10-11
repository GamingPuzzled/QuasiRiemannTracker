"""h-profile bookkeeping for (**) (solverK-starstar section 3).

Units: exponents of P_a.  h = log H'/log P_a in [0,1] indexes the dyadic range d ~ H' of the
sixth-power-free dual frequency in de Faveri's Lemma 6.1.  For each h the contribution to
E_a/(Q ||c||^2) is P_a^{e(h)} and c = max_h e(h), where e(h) is the best available option:

  generic (de Faveri Thm 1.1 on the d-sum):          g(h) = (1-h)/6
  pole, Cauchy per d, amplitude saving sigma:          p(h) = 1/6 + h/3 - sigma      (h <= h0)
  pole, averaged form, saving s over theta-LS:         pa(h) = 1/6 - h/6 - s         (h <= h0)
  off-pole, Cauchy per d + L-avg, saving tau:          o(h) = 1/12 + 5h/12 - tau     (h <= h0)
  e(h) = min(g, max(o, pole))   (pole/off split only allowed for h <= h0)
"""
import numpy as np

H = np.linspace(0, 1, 200001)

def cval(sigma=None, s=None, tau=0.0, h0=1.0):
    g = (1-H)/6
    o = 1/12 + 5*H/12 - tau
    if s is not None:
        pole = 1/6 - H/6 - s
    else:
        pole = 1/6 + H/3 - sigma
    split = np.maximum(o, pole)
    e = np.where(H <= h0, np.minimum(g, split), g)
    i = int(np.argmax(e))
    return e[i], H[i]

if __name__ == "__main__":
    print("no structure (sigma=0, tau=0):              c=%.5f at h=%.4f" % cval(sigma=0))
    print("perfect pole, L-avg only:                   c=%.5f at h=%.4f  (1/7=%.5f)" % (*cval(sigma=9), 1/7))
    print("litG: pole sigma=1/6, off tau=1/12, all d:  c=%.5f at h=%.4f  (5/42=%.5f)" % (*cval(sigma=1/6, tau=1/12), 5/42))
    print("pole sigma=1/6 alone + L-avg:               c=%.5f at h=%.4f" % cval(sigma=1/6))
    for sig in (0.005, 0.01, 0.02, 0.05, 1/14, 0.1):
        c, h = cval(sigma=sig)
        print(f"per-d pole saving sigma={sig:.4f} (h0=1): c={c:.5f} at h={h:.4f};  1/6-sigma/3={1/6-sig/3:.5f}")
    for h0 in (0.0, 1e-3, 0.01, 0.05):
        c, h = cval(sigma=9, h0=h0)
        print(f"perfect pole only for N(d)<=P_a^{h0}: c={c:.5f}  ((1-h0)/6={(1-h0)/6:.5f})")
    for s in (0.005, 0.01, 0.02, 1/42, 0.05):
        c, h = cval(s=s)
        print(f"averaged pole saving s={s:.4f}: c={c:.5f} at h={h:.4f};  max(1/6-s,1/7)={max(1/6-s,1/7):.5f}")
    # what off-pole saving tau is needed (uniform in d) to get below 1/7, given perfect pole
    for tau in (0.01, 0.02, 1/24, 1/12):
        c, h = cval(sigma=9, tau=tau)
        print(f"perfect pole + off saving tau={tau:.4f}: c={c:.5f} at h={h:.4f}")
