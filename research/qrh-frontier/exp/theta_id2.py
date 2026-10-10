import sys, numpy as np, itertools
ADD = sys.argv[4] if len(sys.argv) > 4 else 'chk'
sys.argv = [sys.argv[0]] + sys.argv[1:4]
exec(open('theta_id.py').read().split("for add in ('lam', 'chk'):")[0])
A = np.zeros(len(M_set), complex); R = np.zeros(len(M_set), complex)
for (a, idx, n) in A_set:
    c = np.conj(chi_mod(idx, POW)); A += gauss(idx, POW, ADD)*c; R += np.exp(2j*np.pi*rng.random())*c
G = {pw: np.array([gauss(idx, pw, ADD) for (_, idx, _) in M_set]) for pw in range(1, 6)}
cands = {}
for pw in range(1, 6):
    cands[f"g{pw}"] = G[pw]; cands[f"~g{pw}"] = np.conj(G[pw])
mv = np.array([m[0][0]+m[0][1]*W3 for m in M_set]); ang = mv/abs(mv)
for k in range(-3, 4):
    cands[f'ang^{k}'] = ang**k
base = dict(cands)
for (i, gi), (j, gj) in itertools.combinations_with_replacement([(k,v) for k,v in base.items() if not k.startswith('ang')], 2):
    cands[f"{i}*{j}"] = gi*gj
def corr(u, v): return abs(np.vdot(v, u))/np.sqrt(np.vdot(u, u).real*np.vdot(v, v).real)
ca = sorted(((corr(A, v), k) for k, v in cands.items()), reverse=True)
cr = sorted((corr(R, v) for v in cands.values()), reverse=True)
print(f"POW={POW} add={ADD} #a={len(A_set)} #m={len(M_set)} E|A|^2={np.mean(abs(A)**2):.0f} vs random {np.mean(abs(R)**2):.0f}")
print("  A top:", ", ".join(f"{k}:{v:.3f}" for v, k in ca[:6]))
print("  random-A top (max-of-many baseline):", ", ".join(f"{v:.3f}" for v in cr[:6]))

# project A onto span of top-K candidates; compare explained fraction with random-A
def expl(u, keys):
    B = np.stack([cands[k] for k in keys], 1); coef, *_ = np.linalg.lstsq(B, u, rcond=None)
    return np.linalg.norm(B@coef)**2/np.linalg.norm(u)**2
keysA = [k for _, k in ca[:3]]
print("  explained by top-3 (A):", round(expl(A, keysA), 3), " same keys on random-A:", round(expl(R, keysA), 3))
