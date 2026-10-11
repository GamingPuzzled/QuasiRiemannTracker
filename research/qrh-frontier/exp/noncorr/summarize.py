"""Compact tables from out/grid_*.json and out/pooled_*.json (for the report)."""
import json, sys, os
import numpy as np


def grid_table(path):
    d = json.load(open(path)); cells = d['cells']
    print(f"\n== {path}: per cell (h'=1) ==")
    print(f"{'Y':>7s} {'Q':>5s} {'#m':>5s} {'f':>6s} {'z':>5s} | {'r(T)':>6s} {'<r>T*nu':>8s} {'<r>Tcomp':>8s} {'<r>RAND':>8s} | "
          f"{'+ctrl max r':>11s} {'#m*f':>6s} | {'kap_a T':>9s} {'kap_a RAND':>10s} {'kap_a +ctrl':>11s} {'1/#m':>8s}")
    rows = []
    for c in cells:
        F = c['fam']
        pos = max(F['nu (untw, +ctrl)']['r1_max'], F['g*nu (untw)']['r1_max'])
        kpos = max(c['kappa_a']['g*nu best'][0], c['kappa_a']['one'][0])
        row = (c['Y'], c['Q'], c['nM'], c['f'][0], c['z'][0], c['T']['r'][0], F['T*nu = ab2*g*R*nu']['r1_mean'],
               F['ab2*g*mu [T comps]']['r1_mean'], F['RAND (null, |T|)']['r1_mean'], pos, c['nM']*c['f'][0],
               c['kappa_a']['T'][0], c['kappa_a']['RAND mean'][0], kpos, 1/c['nM'])
        rows.append(row)
        print(f"{row[0]:7d} {row[1]:5d} {row[2]:5d} {row[3]:6.3f} {row[4]:5.1f} | {row[5]:6.2f} {row[6]:8.2f} {row[7]:8.2f} {row[8]:8.2f} | "
              f"{row[9]:11.1f} {row[10]:6.1f} | {row[11]:9.2e} {row[12]:10.2e} {row[13]:11.2e} {row[14]:8.1e}")
    R = np.array([r[:15] for r in rows], dtype=float)
    for Q in sorted(set(R[:, 1])):
        S = R[R[:, 1] == Q]
        print(f"   Q={int(Q):5d}: mean over Y of r(T)={S[:,5].mean():.2f}, <r>T*nu={S[:,6].mean():.2f}, <r>Tcomp={S[:,7].mean():.2f}, "
              f"+ctrl max r={S[:,9].mean():.1f}, kap_a(T)*#m={np.mean(S[:,11]/S[:,14]):.3f}, kap_a(RAND)*#m={np.mean(S[:,12]/S[:,14]):.3f}, "
              f"kap_a(+ctrl)={S[:,13].mean():.2e}")
    allr = R[:, 5]
    print(f"   r(T) over all {len(allr)} cells: mean {allr.mean():.2f} (null Exp(1): mean 1, sd of mean {1/np.sqrt(len(allr)):.2f}); "
          f"max {allr.max():.2f}; P(max of {len(allr)} Exp(1) >= max) = {1-(1-np.exp(-allr.max()))**len(allr):.3f}")


def pooled_table(path):
    d = json.load(open(path)); cells = d['cells']
    print(f"\n== {path}: pooled cross-window #m*rho2 at h'=1 (null ~1; full alignment = #m) ==")
    names = ['T', 'one (untw)', 'ab2', 'nu (untw)', 'ab2*nu', 'T*nu', 'ab2*g*mu [T comps]', 'g*R*nu (untw)', 'g*nu (untw)', 'ab2*g*nu', 'RAND']
    print(f"{'Q':>5s} {'#m':>5s} {'pole z':>7s} " + " ".join(f"{n[:14]:>16s}" for n in names))
    for c in cells:
        F = c['fam']
        def fmt(n):
            x = F[n]
            if x['nv'] == 1:
                return f"{x['m1']:16.2f}"
            return f"{x['m1']:7.2f}±{x['se1']:4.2f}/{x['max1']:4.0f}"
        print(f"{c['Q']:5d} {c['nM']:5d} {c['pole_z'][0]:7.1f} " + " ".join(fmt(n) for n in names))
    print("   (families: mean±se / max)")


if __name__ == '__main__':
    for p in sys.argv[1:]:
        if 'grid' in os.path.basename(p):
            grid_table(p)
        else:
            pooled_table(p)
