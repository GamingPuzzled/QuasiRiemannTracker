# QRH frontier research notes

Working notes from an attempt to lower the quasi-Riemann bound below Hailey Collet's quartic record
θ = 683505193/781250000 (proof `proofs/cycle25-quartic-20261010/`). **No improvement has been proved.**

- `FINDINGS.md`: running log covering the exponent model, sensitivities, dead ends, obstructions, the numeric audit of the record's endpoint, and the large-sieve and theta-residue experiments.
- `reports/`: solver, grader and literature reports. These are model-generated; their claims are not verified unless stated in `FINDINGS.md`.
- `scripts/`: exponent and slack model (`slack2.py`, `gram.py`, `Rgen.py`, …) and the numeric audit (`audit.py`, `audit2.py`).
- `exp/`: Eisenstein-integer sextic large-sieve and theta-residue experiments (`eis.py`, `sieve*.py`, `theta_id*.py`).

Third-party papers referenced in `reports/literatureE.md` are cited by URL only and are not redistributed.
Run scripts from their own directory. Some reports mention the original local scratchpad paths; the equivalent files live here.

## Current frontier (see the end of `FINDINGS.md`)

The best route to a large improvement (θ ≈ 0.865152, conditional) needs a sharp sextic large sieve for an
angular-twisted, theta-structured vector. Equivalently, it needs non-correlation of the cubic theta sum T
with the sextic theta coefficients ρ₆. The numerics support it, but there is no proof. Open sub-question being worked:
the structure of ρ₆ at squarefree indices.
