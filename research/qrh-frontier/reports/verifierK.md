# Verifier K: independent grading of solverK-starstar.md (B6, (★★))

Verifier inputs: KB.md; the report; `scripts/starstar/{gram,profile}.py`; `scripts/gram.py`, `scripts/slack2.py`;
`exp/noncorr/{noncorr,lib,pooled,grid}.py` and `exp/noncorr/out/*`; de Faveri arXiv:2610.04045 (LaTeX source
`main.tex`, read in full for §§1, 4–7); and paper.md (B0.22–B0.27). No repo file other than this one was edited.

**Overall: no claim is wrong, and θ is unchanged, as the report says.** The numbers reproduce exactly. The
h-profile algebra is correct. The reduction matches de Faveri's lemmas once two bookkeeping items are added: a
lemma-label slip, and the Möbius r-terms of Lemma 5.1. Correction C3 is right about `noncorr.py`, but it is stale
about the repo as a whole.

---

## 0. Independent recomputation

**θ values.** I reran `scripts/starstar/gram.py`, which took 3 min. Its probe is identical to `scripts/gram.py`'s
`mk` at η = 0. I checked that at c = 1/6 it reduces to (B0.27)'s exponent x/2 + (y−x)/12 + (1+3ℓ−2M)₊/8: gram/2
with gram = c(y−x), matching (B0.22)'s P_a^{1/6} term.

`minimize_scalar` is a local method, so I also scanned x independently, with my own script calling `slack2.best_for`:
* a coarse grid on [0.30, 0.40];
* a fine grid on [0.325, 0.345] that includes x = 1/3 ± 1e-9.

| c | report | rerun (gram.py) | independent x-scan |
|---|---|---|---|
| 1/6 | 0.874886616 | 0.87488661568 | 0.874886731 (coarse grid, step .0025) |
| 1/6 − 10⁻³ | 0.874852370 | 0.874852370 | — |
| 1/7 | 0.873641042 | 0.873641042 | **0.873641042** at x = 1/3 |
| 5/42 | 0.872046990 | 0.872046990 | **0.872046990** at x = 1/3 |

All values confirmed. The model baseline sits 3.1e-8 below the record (0.87488664704), as stated in U4.

One observation, not an error: for c ≤ 1/7 the optimum moves to the kink x = 1/3, where a different constraint
becomes active (P_a ≈ Z^{.169}). The model's validity there is inherited and was not re-audited against the
paper's constraint list.

**Profile algebra.** I re-derived the four lines of the §3 table by hand from (R), using ‖c‖² ≍ Q, P_a = Y²/Q and
Θ₆(H′,Y) ≍ Y for H′ ≤ P_a < Y^{1/2} (Thm 1.1 with n = 6):

| line | input | e(h) |
|---|---|---|
| g | (Q/Y)(P_a/H′)^{1/6}·Y·E_a | (1−h)/6 |
| o | (Q/Y)(P_a/H′)^{1/6}·‖c‖²H′QY = Q²P_a^{1/6}H′^{5/6}‖c‖² | 1/12 + 5h/12 |
| p | (Q/Y)(P_a/H′)^{1/6}·Y^{4/3}‖c‖²(H′Q)^{5/6}P_a^{−2σ}; uses Y^{1/3}Q^{5/6} = P_a^{1/6}Q | 1/6 + h/3 − σ |
| p̄ | (Q/Y)(P_a/H′)^{1/6}·Y^{4/3}Q^{5/6}H′^{−1/6}P_a^{−2s}‖c‖² | 1/6 − h/6 − s |

All four match. The crossings also check:
* o = g at h = 1/7, value 1/7;
* o − 1/12 = g at h = 2/7, value 5/42;
* p = g at h = 2σ, value 1/6 − σ/3.

I derived the closed form max(1/6 − min(s, h₀/6), 1/7) case by case (h₀ < 1/7 and h₀ ≥ 1/7). `profile.py`
reproduces every value.

---

## 1. Claim-by-claim

### Claim 1: the reduction to (★★). SCORE 8

I checked it against de Faveri's source. In `main.tex` the lemmas are:
* Lemma 4.1 Symmetry;
* Lemma 4.2 Monotonicity;
* **Lemma 4.3 "Completing the 𝔞 sum" (𝓛₁ ≪ 𝓛₂)**;
* **Lemma 4.4 "GCD sorting"**;
* Lemma 5.1 Poisson (𝓛₃ ≪ A + max (A(B₁B₂)^{1/2}/B)·𝓛₂(K₁,L₁)^{1/2}𝓛₂(K₂,L₂)^{1/2}, K_i ≤ Z^ε B²/(AB_i));
* Lemma 6.1 (𝓛₂(A,B) ≪ max_R R·𝓛₁(A/R^n, B)).

What checks out:
* **Step 1.** E_a = Σ_m c_m·conj G(m), where G(m) = Σ_a F_aχ_a(m). Cauchy over all m ~ Q is valid.
* **Identification with Σ₂(A,B,λ).** Σ_m|G|² is exactly Σ₂(Q, Y, F), with 𝔟 ↔ a and the 𝔞-sum over all of I(S).
* **Diagonal term.** Lemma 5.1's leading A gives Q·E_a.
* **Squarefree a.** This forces 𝔟♯ = 1, so B₁ = B₂ = 1. The root number then splits (Remark 3.4 / [2, Lemma 2.3]),
  and the dual length is K ≤ Z^εB²/A = Z^εP_a.
* **Normalization.** The 1/√N(𝔟⋆) in S₁ and S₂ gives exactly the factor (Q/Y)·Σ|B′|².
* **Pointwise in λ (U2).** I confirmed that the proofs of 4.4, 5.1 and 6.1 use the sup over λ only to bound
  Σ|S_i|² by 𝓛₂ of a derived vector. Keeping that vector instead is legitimate.
* **The weight.** In Lemma 6.1, 𝔞 = 𝔡𝔯^n and χ_𝔟(𝔯)^n is only a coprimality indicator. So the R-count is genuine,
  and the claim "(P_a/N d)^{1/6} is precisely the 6.1 loss, not lossy at the level of (★★)" is correct.
* **Kubota index dm.** I re-derived it with de Faveri's conventions:
  * χ_𝔞(𝔟) = (x m_E/𝔟)_n;
  * Poisson gives ε(χ̃_{b₁,b₂})·conj χ̃(𝔞′), with ε ∝ g(χ_{b₁})·conj g(χ_{b₂});
  * so S₁ ∝ Σ_a F_a g̃(a) χ̄_a(a′) = Σ_m c_m Σ_a g̃(a)χ̄_a(a′m);
  * and g(r,a) = conj(r/a)·g(1,a).

  The index is therefore a′m = dm. This confirms C1.
* **(R) logic.** A term K·E_a gives E_a ≤ K‖c‖², and a term X‖c‖² gives E_a ≤ X^{1/2}‖c‖. Correct.

MAJOR_GAPS: none.

MINOR_GAPS:
1. **Lemma labels.** Step 2 calls Lemma 4.3 "GCD sorting". GCD sorting is **Lemma 4.4**. Lemma 4.3 (𝓛₁ ≤ 𝓛₂) is not
   needed, because the chain starts at Σ₂.
2. **Missing Möbius r-terms of Lemma 5.1.** Lemma 5.1's proof removes (𝔟₁⋆,𝔟₂⋆) = 1 by Möbius and writes
   𝔟ᵢ⋆ = 𝔯𝔟ᵢ′ (eq. endgame_Poisson). For a fixed vector, the r-term carries:
   * F_{rb′} = Σ_m c_mχ̄_r(m)χ̄_{b′}(m);
   * by twisted multiplicativity of g̃(rb′), an extra character in b′ of the shape χ_{b′}(r)χ_r(b′), which changes
     the Kubota index by powers of r.

   These terms are not of the "δ-twist" form carried in B′_{t,δ}. Generically they are smaller by about 1/R, so
   only r ≤ P_a^{O(s)} matters. But (Π_s), and the "equivalence" asserted in §5, should quantify over them too.
3. **β_b identification.** That the root number of χ_a is g̃₆(a) is asserted, not traced through Lemma 3.3's β_𝔟.
   The squarefree case makes this standard.
4. **U1 is open.** Non-squarefree columns are not written out.

### Claim 2: Lemma L (Lindelöf on average). SCORE 7

The claim is correctly labelled **[S]**.

What I checked:
* **First bound.** The AFE length is N(r)^{1/2}, and Θ₆(R, R^{1/2}) ≍ R by Thm 1.1:
  R + R^{1/2} + R^{5/6}R^{1/6} + R^{1/3}R^{5/12} ≍ R. So Σ_{r~R, 6-free}|R_Y|² ≪ RY. The restriction to 6-free r is
  exactly what Θ₆ (rather than Ξ₆) allows, which is a good point.
* **Bootstrap.** Σ|𝒦_Y|² ≤ Θ₆(R,R)·R = R^{13/6} gives β ↦ max(5/6, β − 1/6). That is a self-improving inequality at
  fixed R, and it is valid.
* **Cross-check.** The value R^{5/6} = R^{1−1/n} matches the DFDH expectation (cubic analogue R^{2/3}).

MAJOR_GAPS: none for a sketch.

MINOR_GAPS:
* The n = 6 vector-valued Kubota functional equation, with polynomial uniformity in the shift u and in the
  S-unit/class data, is not written out.
* The starting "convexity value β = 4/3" is not justified.
* The residual terms from both poles (2/3 and the dual 1/3) need to be tracked.

This is the item every c < 1/6 statement depends on, so it cannot be scored ≥ 9 until it is written out.

### Claim 3: the normalizations at d ≍ 1. SCORE 9

* Y^{4/3}Q^{5/6}P_a^{1/6} = Y^{5/3}Q^{2/3} = YQP_a^{1/3}. Fed into (R), this gives c = 1/6 exactly.
* P_a^{1/6} = Z^{.086535/6} = Z^{.01442}, and P_a^{1/12} = Z^{.00721}. Both confirmed.
* This agrees with literatureG §5 ((NC-θ), (NC-K)).
* The lower-half equivalence with sharpness of A^{1−1/n}B^{2/n} (A = R, B = Y gives R^{5/6}Y^{4/3}) checks out.

MAJOR_GAPS: none. MINOR_GAPS: the claim is conditional on Claim 2, as stated.

### Claim 4: d ≍ 1 alone gives nothing. SCORE 9

For h > h₀ only the generic line g(h) = (1−h)/6 is available. The pole Cauchy line p(h) ≥ g there, and the
off-pole line is unusable without the pole. So sup_{h>h₀} = (1−h₀)/6, which tends to 1/6 as h₀ → 0.

This is a correct and important correction to literatureG's "ANY partial saving gives c < 1/6". literatureG's own
text says the requirement "shrinks like (P_a/H′)^{1/6} in mean square", which is true relative to the large sieve.
It does not say the saving may be confined to bounded d.

MINOR_GAP: "c ≥ (1−h₀)/6" is a statement about this method's bookkeeping, not a lower bound on the true E_a. The
report's wording ("[P] as bookkeeping") is adequate.

### Claim 5: the partial-saving formula. SCORE 8

The formula c = max(1/6 − min(s, h₀/6), 1/7) is verified analytically (§0) and numerically. So is the per-d form
max(1/6 − min(σ/3, h₀/6), 1/7). The θ examples are verified to 9 digits.

**Uniformity requirement.** The requirement that s hold uniformly for N(d) ≤ P_a^{6s} follows from
1/6 − s ≤ (1−h₀)/6 ⟺ h₀ ≥ 6s. So it is exactly the condition under which the generic tail beyond h₀ does not bind.

MINOR_GAPS:
* It inherits Claim 1's missing r-terms and U1.
* "Any s > 0 beats the record" holds in the model (θ strictly increasing in c). But for s ≲ 10⁻⁹ the gain (≈ 34·s)
  is below the model's own 3e-8 resolution, so a real certificate would have to be recomputed with exact constants.
  This is a cosmetic point.

### Claim 6: the c = 1/7 floor, and 5/42 for (NC-θ) + (NC-K). SCORE 8

* **Floor.** o(h) = 1/12 + 5h/12 meets g at h = 1/7. So with a perfect pole, c = 1/7.
* **(NC-θ) + (NC-K).** σ = 1/6 and τ = 1/12, uniform in d, give c = 5/42 at h = 2/7. The pole line h/3 < 5/42 there,
  so it does not bind.
* **Range.** The range [P_a^{(12c−1)/5}, P_a^{1−6c}] follows from o(h) > c ⟺ h > (12c−1)/5 and g(h) > c ⟺ h < 1−6c.
  Checked.
* **Interpretation.** Reading literatureG's (NC-K) as a fixed saving τ = 1/12 over "Cauchy in m + Lindelöf on
  average" at every d is the natural reading of its displayed inequality. That inequality is per-d Cauchy against
  (Σ_m|ψ̃(dm)|²)^{1/2}.

MINOR_GAPS:
* "Floor" is method-relative. Off-pole terms are bounded only by per-d Cauchy plus Lemma L, or by the generic line.
  A bilinear treatment of R_Y(dm), for example expanding it through its AFE and applying a large sieve jointly in
  (d, m), is not excluded. The claim's qualifier ("no further off-pole input") covers this.
* The prose "off-pole savings at N(d) ≈ P_a^{1/7…2/7}" (§3, §5) is loose. For c = 5/42 the exact range is
  [P_a^{3/35}, P_a^{2/7}]. The numbered claim gives the correct general range.

### Claim 7: handle (i), sphericity gives no cancellation. SCORE 7

The negative conclusion is right. Sphericity of θ₆ at the complex place is a statement about the archimedean
Whittaker function. It says nothing about the arithmetic sequence r ↦ ϱ₆(r) against arg r.

The torus point is also right in substance. On the KP n-fold cover, the commutator of lifts of diag(y,1) and
diag(y′,1) is a nontrivial power of (y,y′)₆. For n = 6, the square (cubic) symbol is nontrivial. So there is no
Hecke integral unfolding to Σϱλ̄N^{−s}, which is the standard reason metaplectic forms lack Euler products.

I re-did the μ-aspect arithmetic:
* 5/12 + 0.8(1/2 + κ) ≤ 11/12 − s(y−x)/M gives κ ≤ 0.125 − 0.130s. The report writes 0.1253; the exact constant
  is 1/8.
* With conductor N(μ)^k, convexity Q^{−1/12}N(μ)^{k/2} meets this at N(μ) = Q^{0.8} iff k ≤ 3/2.
* At k = 3 the gap over convexity is Q^{0.6}, at analytic conductor Q^{2.4}. Lindelöf (Q^{5/12}) suffices.

The "Lindelöf strength" conclusion stands.

MAJOR_GAPS: none.

MINOR_GAPS:
* The label "[P]" overstates what is shown. The report proves that one standard mechanism fails. It does not prove
  that sphericity cannot help.
* It omits the mechanism that does use sphericity: the Heath-Brown–Patterson fact that a Kubota series twisted by an
  angular Grössencharacter in the *Gauss-sum-modulus* variable is entire. In Π(d), m is also a Gauss-sum modulus
  (through g̃₆(m)ᾱ²(m)). The report folds this into the [H] "twisted A₂ double series" paragraph without saying
  that this is where sphericity actually enters.
* The claim that the "expected k = 3" conductor comes from a GL₃-type residual Eisenstein series is heuristic, and is
  labelled as such.

### Claim 8: no-go for local information. SCORE 8

The construction works.
* Prescribe x on squarefree primary m as conj(c_m)·N(m)^{−1/12}·const.
* Extend it by ω to associates and by the K2 relations to non-squarefree arguments.
* K2 never links two distinct squarefree arguments. The Mp ↔ Mp³ loop involves one non-squarefree value.
* Squarefree values elsewhere can be set to size N^{−1/12}, which gives the R^{5/6} mean square at every scale.
* Cauchy–Schwarz then holds with near-equality, up to the 2^{1/12} dyadic variation of N(m)^{−1/12}.

MINOR_GAPS:
* K2 itself is cited from literatureG / Bröker–Hoffstein §1 without a precise reference. In particular, the claim
  that the 1 ↔ 3 loop composite is consistent (not forcing vanishing) and the exact forms of G₁ and G₃ are not
  pinned down. The claim's own hypothesis clause makes it conditional on this.
* (Also: "g(r,p^k), p ∤ r, vanish unless 6 | k" (§2) is inaccurate. For k ≥ 2 they vanish even when 6 | k, as
  Ramanujan sums. This does not affect any claim.)

### Claim 9: handle (iii), averaging over d. SCORE 9

The averaged pole line at s = 0 is p̄(h) = 1/6 − h/6 = g(h), which I verified algebraically. So perfect orthogonality
of the rows d ↦ ϱ(d·) gives exactly the generic exponent. At h = 0 there is a single d, up to S-units. Correct.

MINOR_GAP: none of substance. The "theta large sieve" normalization (H′+Q)(H′Q)^{−1/6} rests on Claim 2.

### Claim 10: correction C3, the noncorr h′ ≠ 1 tests. SCORE 8

**Checked in the code.**
* `noncorr.py` builds `Xh = X[:, nM:]` (not conjugated) and `Xm = conj(X[:, :nM])`.
* `lib.char_matrix` returns χ_x(t) = (t/x)₆ for modulus x.
* `symbol_table` builds the Gauss sum of (·/π)₆ (POW = 1).

So A_m(h′) = Σ_a g̃(a)(h′/a)₆·conj(m/a)₆ = Σ_a g̃(a)·conj(h′⁵m/a)₆ = Σ_a g̃(h′⁵m, a). The **Kubota index is h′⁵m**,
which is not the (★★) index h′m. Confirmed.

solverG's own report states the residue as ρ₆(h′⁵m) (§4.1, and the `pooled.py` docstring). So the h′ ≠ 1 statistics
in `noncorr.py` and the `mh5` runs do not test the (★★) frequencies.

**The data agree with the vanishing explanation for sextic** (`out/pooled_sextic.txt`, index h′⁵m):
* the pole z at N(h′) = 13 and 19 is about 0 (−2.5 … 0.5);
* at N(h′) = 7 it is 3.5 → 7.8. That is consistent with the ≈ 1/7 of rows where h′ | m, where the index is h′⁶m′ and
  τ₆ ≠ 0. The claim's qualifier "h′ ∤ m" covers this.

MINOR_GAPS:
1. **Stale relative to the repo.** `pooled.py` and `grid.py` already have `--index md` (conj χ_a(d), giving index
   md; default, docstring "correct, (★★)"). They were modified at 00:57, before the report was written at 01:15.
   md outputs exist: `out/pooled_sextic_md.txt`, `grid_sextic_md.txt` and `pooled_cubic_md.txt`. They show a strong
   pole at h′ ≠ 1 (z ≈ 13–20 for sextic), as the index-h′m analysis predicts. So "the binding frequencies d = h′ are
   untested" (summary item 5) is no longer true of exp/noncorr/. It is true only of `noncorr.py` and the solverG
   write-up.
2. **The KP vanishing τ₆(Mp⁵) = 0 is cited, not referenced.** A caution: the cubic analogue in the same repo
   (`pooled_cubic.txt`, index h′²m, j = n−1) shows a clearly non-zero pole at prime h′. Its z(h′)/z(1) ≈ 0.39, 0.27,
   0.23 at N = 7, 13, 19, which is ≈ N(h′)^{−1/2}. The naive local pairing j ↔ n−2−j would predict vanishing at
   j = n−1 for n = 3 just as for n = 6. Either the cubic case differs (convention 'chk', or a different local
   relation), or the pairing is not the right description. The sextic data are consistent with the claim. But the
   explanation "vanishes identically by KP" should be backed by a precise reference. Claim 10's main point, that
   the index is h′⁵m, does not depend on it.

---

## 2. Answers to the advisor's priority questions

1. **c ↦ θ.** Reproduced exactly: θ(1/7) = 0.873641042 and θ(5/42) = 0.872046990, both at x = 1/3, confirmed by an
   independent x-scan. The formula c = max(1/6 − min(s, h₀/6), 1/7) is verified analytically.
2. **Reduction against de Faveri.** It is faithful. Lemmas 4.4 → 5.1 → 6.1 are pointwise in λ until the 𝓛 sup. The
   weight (P_a/N d)^{1/6} is the Lemma 6.1 root count, and the index is dm. Two fixes are needed: the label (GCD
   sorting is Lemma 4.4), and the Möbius r-terms from (𝔟₁⋆,𝔟₂⋆) = 1 in Lemma 5.1's proof, which are not carried.
3. **Uniformity N(d) ≤ P_a^{6s}, and the 1/7 floor.** Correct as bookkeeping. 6s comes from (1−h₀)/6 ≤ 1/6 − s, and
   1/7 is the crossing of the Lemma-L off-pole line with the generic line. Both are method-relative, as stated.
4. **C3.** Correct for `noncorr.py` (index h′⁵m). It is stale about `pooled.py`/`grid.py`, which already run index md
   and show the predicted pole at h′ ≠ 1.
5. **Sphericity.** The conclusion and the non-abelian-torus reason are right. "[P]" overstates it (only a mechanism
   is ruled out), and the Heath-Brown–Patterson angular-twist mechanism in the modulus variable is not discussed
   explicitly.

## 3. Scores

| Claim | Topic | SCORE | Major gaps | Main minor gaps |
|---|---|---|---|---|
| 1 | Reduction to (★★), index dm | 8 | — | Lemma 4.3/4.4 label; Lemma 5.1 Möbius r-terms not carried; U1 |
| 2 | Lemma L (Lindelöf on average) | 7 | — | sketch: n = 6 functional equation uniformity, β = 4/3 start |
| 3 | Normalizations Z^{.0144}, Z^{.0072} | 9 | — | conditional on Claim 2 |
| 4 | d ≍ 1 gives nothing | 9 | — | method-relative wording |
| 5 | c = max(1/6 − min(s, h₀/6), 1/7) | 8 | — | inherits Claim 1 gaps; model resolution |
| 6 | 1/7 floor; 5/42 for NC-θ + NC-K | 8 | — | floor method-relative; loose prose range |
| 7 | Sphericity gives no cancellation | 7 | — | "[P]" overstated; HB–Patterson modulus-twist mechanism not addressed |
| 8 | No-go for local information | 8 | — | K2 relations unreferenced |
| 9 | Averaging over d | 9 | — | — |
| 10 | C3: noncorr index h′⁵m | 8 | — | stale (md runs exist); KP vanishing unreferenced, cubic data differ |

No claim scores ≤ 6.
