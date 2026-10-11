# Solver J: an ε-improvement of the numerator Hölder step (bottleneck B5)

**Verdict: nothing is proved. θ is unchanged at the record 683505193/781250000 = 0.87488664704.**
- The three candidates (a)–(c) all fail, each for a precise reason.
- The main negative result is an *exponent-level adversarial profile* (§4, Claim 4), checked numerically. It is a
  virtual bin of U^R rows that satisfies **every** row inequality the paper makes available, yet has
  Σ_bin|S_u(U^n)| = U^{H4} exactly, where H4 = 1/4 + 3R/4 − qz_A/2 (B0.41).
- So no ε-gain at B5 can come from level sets, Hölder, interpolation, joint witness classes or mixed moments
  built from UM/(B0.9), MI, MP/MP0 (B3.40), (B0.15)/(B0.16) and Lemma 1W.
- The smallest new input that would suffice is a "lossy sixth moment" on bin rows, priced in §2.

Scripts are in `research/qrh-frontier/scripts/epsholder/`:
- `base.py`: a copy of `smallgain/base.py`.
- `point.py`: binding data, and the (a) and (c) quantities.
- `cond.py`: a copy of the `smallgain/model.py` exponent model. The numerator is min(H4, pointwise, an optional
  hypothetical line). It has a robust nested optimiser that checks the geometry constraints of `slack2.feasible`,
  every p ∈ [0, 1/2] and every dyad d, together with exact mpmath smooth solves.
- `eta6.py` and `table.py` (output in `table.out`): conditional payoffs.
- `scan6.py` and `diag6.py`: fixed-e6 diagnostics.
- `adversary.py`: the profile check (Claim 4).

The float solver reproduces the exact b0 to 3.3e-16. It reproduces solverH's calibration exactly: H1 = 1/4 − 1.3e-6
gives θ = 0.8748865496775439, a gain of 9.736e-8.

## 1. Binding data [exact, `point.py`]

All values are in U-units, with U = Z^{h0}, at d = h0, p = 1/2, δ0, κ0. Write σ = log_U|S_u(U^n)|².

| quantity | value |
|---|---|
| R (B0.40) | 0.5903323535263558 |
| t* | 1.163953576509011 |
| r_I (A_I(r_I) = R) | 0.7456999427242902 |
| m_P (A_P(m_P) = R) | 0.4182536337847208 (and r_I + m_P = t* exactly) |
| n = 1 − y0/h0 = **2ℓ0/h0** (reflected numerator length) | 0.4209851608760782 |
| n − m_P | 0.0027315 |
| z_A = (1 − 2n)/(6κ0) | 0.03512832467618695 |
| σ* = (1 − R − 2qz_A)/2 (Hölder-extremal level) | 0.1965901739460984 |
| pointwise level δ0·n | 0.1975872208453980 |
| δ′ = σ*/n | 0.4669764927984408 |
| H4 | 0.688627440499405 |
| pointwise − H4 | 4.985e-4 |

**Structural fact.** The A_P line evaluated at a plain level σ = δm is *identical* to the numerator's
fourth-moment line, because both are 1 − 2σ − 2q(1−2m)/(6κ). The pointwise bound is therefore inactive by
exactly the amount by which the numerator length n exceeds the largest admissible plain-witness length m_P.
That amount is n − m_P = 0.00273. Exactly, δn − σ*(n) = (δ − q/(3κ))(n − m_P) = 0.36500 × 0.0027315 = 9.970e-4,
because δm_P = σ*(m_P) by the definition of m_P.

## 2. Candidate (a): higher moments, lossy or on a sub-range

**(a1) Trivial interpolation gives nothing [proved].** For k ≥ 2 the bound
Σ_B|S|^{2k}|Q|² ≤ (max_B|S|²)^{k−2}·Σ|S|⁴|Q|² is fed into Hölder with weight 1 − 1/(2k) on |B|. The result is
exactly the convex combination (2/k)·H4 + (1 − 2/k)·PW, where PW = R + δn/2. It is therefore ≥ H4 at the
binding point (identity in Claim 2). The numerical excess over H4 is 1.66e-4 for k = 3, 2.49e-4 for k = 4 and
2.99e-4 for k = 5.

**(a2) What a lossy sixth moment must deliver.** Consider a moment with saving η over the trivial bound:

  H6(η): Σ_{u∈bin}|S_u(U^n)|⁶|Q_{u,I}|² ≪ U^{1+(δ−η)n+ε}  for 2n + 6κz ≤ 1, uniformly in bins.

At fixed record geometry it beats H4 iff η > δ0 − δ′ = **0.0023683659**. That means beating the trivial bound by
U^{9.97e-4}. After re-optimising the geometry (`table.py`, robust solver, every constraint and dyad checked,
maxG ≤ 0), the gains are:

| hypothesis | θ | gain |
|---|---|---|
| H6(0.0023) | 0.8748866470400 | 0 (H4 still active) |
| H6(0.00236) | 0.8748864061460 | 2.41e-7 |
| H6(0.0024) | 0.8748859329701 | 7.14e-7 |
| H6(0.0025) | 0.8748847499834 | 1.90e-6 |
| H6(0.003) | 0.8748788340489 | 7.81e-6 |
| H6(0.01) | 0.8747958353817 | 9.08e-5 |
| H6(0.05) | 0.8743151846637 | 5.72e-4 |
| 8th-moment analogue H8(0.00236) | 0.8748859436965 | 7.03e-7 |
| 8th-moment analogue H8(0.003) | 0.8748749711440 | 1.17e-5 |
| 8th-moment analogue H8(0.01) | 0.8747545824198 | 1.32e-4 |
| Lindelöf-on-average sixth moment (e = 1, slots kept) | 0.8692779055918 | 5.61e-3 |

Notes on the table:
- The θ-threshold lies in (0.0023, 0.00236), slightly below 0.0023684. Past it the optimum jumps to a new branch,
  with h ≈ 0.8011 and binding bin δ ≈ 0.439, where the new line is active.
- An unphysical fixed exponent e6 (independent of n) is misleading. The optimiser then runs into the x ≤ y
  boundary of `slack2.feasible`; this is why we priced η instead (`scan6.py` and `diag6.py`).

**Obstructions to proving H6(η) for any η > 0:**
1. *Length.* S³ has length U^{3n} = U^{1.263} > U. Any general-coefficient mean value gives
   e6 ≥ max(1, 3n) = 1.263. That holds even for the conjectural sharp sieve (U + K); the
   Blomer–Goldmakher–Louvel sextic sieve gives 1.509. Both are *worse than the trivial* 1 + δn = 1.1976.
   A sixth moment would have to use the arithmetic of S_u, beyond the large-sieve barrier by 0.066 in the exponent.
2. *Sub-range.* A sixth moment valid only for 3n ≤ 1 is of no use. The identity n = 2ℓ/h holds because
   y = h − 2ℓ when M + ℓ = 1, and the slot-supply condition ℓ/h > 4/21 used in (B0.40) forces n > 8/21 > 1/3
   at every feasible geometry.
3. *Family.* B3.40 controls second moments of products of **two** plain sums with prime slots. Its induction
   (B3.5–B3.10) is over A = n1 + n2 + z. Σ|S|⁶ needs a three-fold product, which is a new theorem.

## 3. Candidate (b): factorising the numerator through MP0

This route does not help, for three reasons [proved].
- **No factorisation exists.** The plain coefficients n ↦ ν(n)χ_n(u)W(q_n/N) are a character times a smooth weight.
  The only Dirichlet factorisations are 1 = 1 (trivial) and 1 = μ ∗ τ. The second gives M·S·S, and Möbius × plain
  mixed moments do not exist (F4, Solver B).
- **Even granted a factorisation, the moments are too weak.** Suppose S ≈ Σ_j S_1^{(j)}S_2^{(j)} with bounded
  coefficient mass. MP0 bounds Σ|S_1S_2|², which is a *second* moment of the product. Hölder from second moments
  gives Σ_B|S| ≤ U^{(1+R)/2} = U^{0.7952}, worse than H4 by (1 − R)/4 = 0.1024. The current H4 is already MP0
  applied to S·S. Beating it needs Σ|S_1S_2|⁴, a four-fold product.
- **Mixed MP is only an average.** The mixed line with a plain witness of length m (B3.40 with n1 = m, n2 = n) is
  exactly (A_P + num4)/2 (solverH identity, re-checked). At zero slots it also follows from Cauchy–Schwarz of the two
  pure fourth moments. By the averaging lemma (P4) it never beats min(A_P, num4).

θ is unchanged.

## 4. Candidate (c): the equality case and large values

**Is saturation compatible with the bin's witness condition? Yes [Claim 4, `adversary.py`].** Define a virtual bin
of U^R rows, all with one presentation, which is also the numerator presentation:
- μ(r) = log_U|M_u(U^r)|²:
  - for r ≤ 1, μ(r) = min(δr, 1 − R − 2q·z_I(r)), where z_I(r) = min((1−r)/2, (3−2r)/8, L);
  - for 1 < r ≤ 3/2, μ(r) = min(δr, E(r) − R).
- s(ν) = log_U|S_u(U^ν)|² = min(δñ, (1 − R − 2q·z_P(ñ,ñ))/2), where ñ = min(ν, 1−ν) and
  z_P(a,b) = min(L, (1−a−b)/(6κ)).
- |Q_I|² = U^{2qz_I} for every slot subcollection with z_I ≤ L.
- Every other row of the family is typical (size U^{o(1)}).

The script checks, at the exact binding point, that all of the following hold, with maximum violation ≤ 2.2e-16:
- **C1:** MP/MP0/mixed MP for every pair of lengths and every admissible capacity, R + s1 + s2 + 2qz ≤ 1.
  Analytically: s(a) + s(b) ≤ 1 − R − q(z_P(a,a) + z_P(b,b)) = 1 − R − 2q·z_P(a,b).
- **C2:** MI with slots, including r = 0 prime-sum moments.
- **C3:** UM and (B0.9).
- **C4:** pointwise (B0.15) and (B0.16), and the amplitude cap.
- **C5:** Lemma 1W witnesses for every t0 ∈ [1, 3/2]. These are short witnesses with r ≤ r_I and m ≤ m_P when
  t0 ≤ t* = r_I + m_P, and the long witness r = t0 when t0 ≥ t*.

Meanwhile s(n) = σ*, so Σ_bin|S_u(U^n)| = U^{R+σ*/2} = U^{H4} exactly (difference 0.0). In this profile the plain
witnesses are saturated (|S|² = U^{δm}) for m ≤ m_P, and the numerator sits at σ* < δn just past m_P. The
fourth-moment equality case is therefore *consistent* with the witness structure. There is no incompatibility to
exploit.

**Large values (Montgomery–Huxley type) [proved / quantified].**
- **Halász–Montgomery** for D = S^j (length X = U^{jn}, normalised level X^{δ′}) gives
  |B|·X^{δ′} ≤ X + |B|·Off. Information requires X^{δ′} > Off. Even square-root off-diagonal bounds
  (GRH strength, Off ≍ X^{1/2}) require δ′ > 1/2, and here δ′ = 0.46698. So HM gives nothing for any j,
  unconditionally or under GRH-type off-diagonal bounds.
- **Burgess (r = 2).** The off-diagonal k ↦ χ_k(uv⁵) has conductor ≤ U², so Off ≤ X^{1/2}U^{3/8}. This exceeds
  the trivial X for j = 1 and j = 2, and for j = 3 it is far above the level.
- **Huxley-type count** (heuristic transfer, with the family size U in place of T): exponent
  max(x(1−δ′), 1 + x(1−3δ′)), where x = jn. This gives 0.8312, 0.6624, 0.6732 and 0.8976 for j = 1–4, all above
  R = 0.5903. Only the non-integer power j = 2.54 would give 0.5707.

So no large-values inequality gives a count below U^R at the level σ* that the bin needs.

## 5. Precise obstruction (summary)

At (δ0, p = 1/2, κ0, d = h0), H4 is the exact value of the virtual configuration in Claim 4. Any ε-improvement of
B0.41 therefore needs an inequality violated by that configuration. The cheapest such inequalities found are:
- **(i)** H6(η) for some η > 0.00236. This is a bin-restricted sixth moment that beats the trivial interpolation by
  U^{≈1e-3}, about 0.066 beyond what any general-coefficient sieve can give. It needs a three-fold plain-product
  extension of B3.40.
- **(ii)** A count line for {u ∈ bin : |S_u(U^n)|² ≥ U^{σ*}} strictly below R at σ*. This sits below the
  Halász–Montgomery square-root barrier (δ′ < 1/2).
- **(iii)** A cross-moment linking inverse and plain profiles, such as a joint M·S moment for r > 1 (B0.9 type).
  This would break the long-class rows (r = t*, m = 0) that carry count R with no plain witness.

None is available.

## CLAIMS

1. **(Exact data; `point.py`, `base.py`.)** At the record binding point the following hold:
   - n = 1 − y0/h0 = 2ℓ0/h0 = 0.42098516087608, and z_A = 0.03512832467619;
   - σ* = 0.19659017394610 and δ0·n = 0.19758722084540;
   - H4 = 0.688627440499405, and R + δ0n/2 − H4 = 4.985e-4;
   - m_P = 0.41825363378472, and r_I + m_P = t* to 1e-16.
   *θ consequence: none (bookkeeping).*
2. **(Proved, elementary.)** For every k ≥ 2, consider Hölder applied to Σ_B|S|^{2k}|Q|² ≤ (max_B|S|²)^{k−2}·Σ|S|⁴|Q|²,
   with max_B|S|² ≤ U^{δn} (B0.15) and MP (B0.5). Its exponent is (2/k)·H4 + (1 − 2/k)(R + δn/2) ≥ H4, with equality
   iff k = 2. *θ consequence: none.*
3. **(Proved.)** Hypothesis: the numerator plain sum is written as a bounded combination of products of two plain sums.
   Then MP0 (B3.40, z = 0) yields at best Σ_B|S| ≪ U^{(1+R)/2+ε}, which exceeds H4 by (1 − R)/4 + qz_A/2 > 0. The
   mixed MP line equals (A_P + num4)/2 identically. *θ consequence: none.*
4. **(Proved at the exponent level; `adversary.py`.)** The virtual bin of §4 satisfies, up to 2.2e-16, every row
   inequality C1–C5: UM, (B0.9), MI (B0.4), MP/MP0/mixed MP (B3.40) at all lengths and capacities, (B0.15), (B0.16),
   the amplitude cap, and Lemma 1W for every t0 ∈ [1, 3/2]. Its numerator sum is exactly U^{H4}. Hence any bound for
   Σ_bin|S_u(U^n)| derived from these inequalities by dyadic level splitting and Hölder (in any combination, any
   number of truncations t0, any joint witness classes) is ≥ U^{H4−o(1)}. *θ consequence: B5 cannot be achieved
   within the paper's toolset; the record is method-optimal at this step.* The claim concerns the method; it does
   not assert that real L-functions realise the profile.
5. **(Proved.)** For D = S_u(U^n)^j at normalised squared level X^{δ′}, with δ′ = σ*/n = 0.46698 < 1/2, the
   Halász–Montgomery inequality gives no bound on the number of such rows. This holds even with off-diagonal bounds
   of square-root strength. *θ consequence: none.*
6. **(Proved.)** In the paper's regime (M + ℓ = 1 and ℓ/h > 4/21), the reflected numerator length is n = 2ℓ/h > 8/21.
   So a sixth moment restricted to 3n ≤ 1 never applies. Moreover max(1, 3n) > 1 + δn at the binding point, so even
   a sharp general-coefficient sieve is worse than the trivial bound.
7. **(Conditional; `table.py`, robust solver accurate to 3e-16 against the exact baseline.)** Assume H6(η) for
   every nonfloor bin and amplitude class. Then:
   - θ(0.00236) = 0.8748864061460, θ(0.0025) = 0.8748847499834 and θ(0.01) = 0.8747958353817;
   - there is no gain for η ≤ 0.0023;
   - a Lindelöf-on-average sixth moment (e = 1) gives θ = 0.8692779055918.
   The hypothesis is unproved. Its proof needs a three-plain-factor version of B3.40.
8. **(Calibration, exact.)** Lowering the Hölder coefficient to 1/4 − 1.3e-6 gives θ = 0.8748865496775439, a gain
   of 9.736e-8, matching solverH's sensitivity +0.0749.
