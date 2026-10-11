# Solver H — small-gain hunt on the quartic record (θ = 683505193/781250000)

**Verdict: no rigorous improvement found, of any size ≥ 1e-7.** The only "free" slack anywhere in the
chain is the rational rounding θ − b0 = 3.587e-14 (b0 = 11/12 − ℓ0/4 is what is actually proved).
Below, every item is marked **[proved]** (exact computation or a short complete argument) or
**[heuristic]**.

Scripts are in `research/qrh-frontier/scripts/smallgain/` and use mpmath at 30–50 digits:
`base.py` (exact record point), `model.py` (a generalized binding system), `sens_all.py`,
`margins.py`, `triple.py` and `deadends.py`.

## 0. Exact reproduction [proved]
`model.py` solves the three binding equations from first principles: F = b, ∂δF = 0, and the
Lagrange condition in (ℓ,h). Inputs are B0.27, B0.38–B0.41 and B0.47, with M + ℓ = 1 and p = 1/2. Each
exponent coefficient is a named parameter. At the paper's values it reproduces ℓ0, δ0 and h0 to 25 digits,
and θ − b0 = 3.5869e-14, which matches T1.5.

## (a) Slack in the certificate's inequalities [proved: none affects θ]
The root of 927X⁴ − 3135X³ + ... depends only on the active set at the optimum. `margins.py` lists
every side condition at the exact record point. Only these are **ACTIVE**:
- M + ℓ ≤ 1, the hypothesis of B0.26/B2.44;
- p ≤ 1/2, the amplitude cap (∂F/∂p = +0.0324 > 0);
- the δ-max and h-min saddle (B0.48);
- the dyad endpoint d = h.

Every constant the task names has a strictly positive margin, so its KKT multiplier is zero. Sharpening it
changes ℓ0 by exactly 0. The margins:

| condition | margin |
|---|---|
| 1/18 (T1.3, ∂κC) | 0.0205 |
| ∂F/∂κ < 1/10 (only < 1/2 is needed) | large |
| 4/21 supply | 0.0200 |
| B0.56 < −1/20 | 0.0359 |
| other endpoint (< b0 − .06) | 0.0044 |
| floor (< b0 − .01) | 0.0050 |
| y − 2x + ℓ < 0 | 0.1195 |
| energy (1 + 3ℓ − 2M)_+ | argument −0.164 |
| MI second condition 2r + 8z ≤ 3 | 0.491 |
| r* ≥ 13/21 | 0.127 |
| m* ∈ [1/3, 1/2] | 0.082 |
| κ ≤ 3/4 | 2.3e-4 (moves the wrong way for improvements) |

69 < J ≤ 270 is non-binding as well (the audit gives J ∈ [81.5, 269.5]).

The energy slack (1 − 5ℓ)/4 ≈ 0.041 in B0.26 cannot be used. The binding item there is the
hypothesis M + 2z ≤ N, not the energy. Past it, (B2.38)'s column branch gives energy exactly T_d = M + s
(re-derived), the same as Dead end 1.

**Full sensitivity map** (`sens_all.py`, exact implicit-function derivatives). It lists the change |Δc|
needed for a 1e-7 gain:

| coefficient (source) | dθ/dc | \|Δc\| for 1e-7 | nature |
|---|---|---|---|
| Hölder 1/4 in 1/4 + 3R/4 (B0.41) | +0.0749 | 1.3e-6 | estimate (needs > 4th moment) |
| probe Gram 1/12 = c/2, c = 1/6 (B0.22) | +0.0669 | 1.5e-6 | estimate (genuine main term per log) |
| MI capacity ½(1−r) (B0.4) | −0.0083 | 1.2e-5 | estimate (sharp in B3.38 at G = 0) |
| long slope 5/6 (B0.9) | +0.0074 | 1.4e-5 | estimate (open problem ★) |
| amplitude cap ½δ (B0.16) | +0.0075 | 1.3e-5 | sharp (explicit formula) |
| MP cost 6κ in A_P / in z_A | +6.4e-5 / +1.3e-4 | 1.6e-3 / 8e-4 | estimate (Solver A: 6 is a floor) |
| 1/2 in a = (1+δ)/2; Y^{-1/2}; Z^{-ℓ/2}; 1/6 pole | 0.11 / −0.11 / −0.038 / −0.18 | ~1e-6 | identities, not estimates |

## (b) Triple point R = A_I(r*) = A_P(m*) = long(t*) [proved: mixing cannot help]
- **Averaging lemma.** Suppose rows satisfy pointwise lower bounds w_j and moment bounds Σ f_j ≪ U^{e_j}.
  Hölder with weights θ_j gives the count bound Σθ_j(e_j − w_j), a convex combination of the individual
  lines, which is never below their minimum. This covers MI/MP Hölder-mixing, with the same or different
  slot sets.
- **Mixed MP (B3.40, n1 = m, n2 = m′).** Symbolically, its count line is *exactly* (A_P + num4)/2
  (`triple.py`: difference = 0). This is why `joint.py` saw a "quadruple point" (mix = R to 1e-41) and no
  gain.
- **Two truncations t1 < t* < t2.** Applying Lemma 1W twice and taking the worst joint class reproduces R
  exactly at t2 = t*. Classes split by r + m are worst at r + m = t0.
- **Variants checked** (`deadends.py`), all worse:
  - Sixth-power lift combined with MI slots: the z-coefficient is 5/3 − 2q = 1.197 > 0.
  - Lift applied to short rows: the log P coefficient is 5 − 6q > 0.
  - MP with n2 = 0: the count is worse by 0.153.
  - Full-width MI for long rows: 0.618 > R.
  - Refining bins per member of 𝒳(u) (denominator level vs slot level): the adversary equalizes the
    levels, so the worst case is unchanged.

## (c) Splitting the bin instead of taking a min [proved: zero gain]
Splitting by numerator level σ, with the pointwise cap and the mixed line included, gives exactly
min(Hölder, pointwise) = Hölder = 0.688627440499405 (difference 0.0). At δ0, pointwise is *inactive* by
4.985e-4 U-units (3.96e-4 Z-units). Splitting by witness class gives g(max class count) = g(R), because
g is monotone.

Changes confined to δ > 0.49 (where pointwise wins) do not move the optimum. G(h) = max_δ F has a
smooth minimum at h0 driven only by the neighbourhood of δ0, since ∂_hF(δ0) = 0. So the gain is
identically 0, and the earlier 7e-14 was numerics.

## What a ≥1e-7 record would minimally need [heuristic guidance]
Any one of the following would be enough:
- a numerator exponent 1/4 − 1.3e-6;
- a Gram exceptional exponent c ≤ 1/6 − 3e-6;
- an MI capacity of (½ + 1.2e-5)(1−r);
- a long slope 5/6 − 1.4e-5;
- MP cost 6 − 5e-4 (both uses together).

Each is a new estimate, not a bookkeeping fix. One genuinely new direction remains unchecked: binding rows
need simultaneous low zeros of L(s, ψ̄_u) and of the denominator L(s, ηψ̄_u). A UM/MI theorem for
1/(L(ψ̄)L(ηψ̄)) would add a non-averaging count line. It is not available in B3.
