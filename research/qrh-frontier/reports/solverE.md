# Solver E: the angular-twisted sextic large sieve (ATS) via Dunn–Radziwiłł §9

**Verdict: no unconditional improvement. The record θ is unchanged, and the paper's c = 1/6 stands.**

What is established:

* **Theorem E1 (conditional on GRH).** Sketched rigorously in §2.3. ATS holds sharply for the pure twisted Gauss vector b_m = ψᾱ²g̃₆μ²W. The twist does remove the Heath-Brown/Patterson obstruction, as the numerics suggested.
* **Our structured vector b = g̃₆·ψᾱ²𝓡 is not covered**, even under GRH.
  * The DR mechanism needs pointwise square-root cancellation of β′ = ψᾱ²𝓡 against every sextic character (DR Axiom (3.2) with η ≤ 1/4).
  * 𝓡 is a combination of D = Q^{0.8} quadratic characters with ℓ¹-mass D^{1/2}, so GRH gives only η ≈ 0.40.
  * Replacing the pointwise bound by a mean square is exactly circular under Poisson.
* **Conditional statement.** Under an explicit hypothesis PSR (numerically true; it is stronger than GRH), c = 0 and θ ≈ 0.865152.
* **The c = 1/6 is explained.** It is the one-step Heath-Brown recursion, and it equals the de Faveri–Dunn–Hoffstein conjectured generic norm Ξ₆ at this geometry. Beating it on that route needs non-correlation of T with **sextic** theta coefficients.

Scripts: `solverE/psr.py` (numerical PSR test) and `solverE/expo.py` (exponent bookkeeping).

---

## 1. Exact formulation

Notation follows FINDINGS and solverD:
* Y = Z^{.4597} is the a-length.
* Q = Z^{.8329} is the m-length.
* P_a = Y²/Q = Z^{.0865}.
* K = Q²/Y = Z^{1.206}.
* D = Q²/Z = Z^{.666}.

The target in solverD is E_a ≪ Q²Z^ε. Let T(m) = κν(m)ᾱ(m)²γ₁(m)𝓡(m), where γ₁ = g̃₆ is the normalized sextic Gauss sum and 𝓡(m) = Σ_{N(μ)≲D} ρ_μ (μ/m)₂ with ρ_μ ≍ g̃₃(μ)α(μ)N(μ)^{-1/2} (solverD Claim 3).

**(F1, the a-form = ATS with H = Y)**  Σ_{a~Y} |Σ_{m~Q} g̃₆(m)β′(m) χ̄_a(m)|² ≪ Q²Z^ε,  with β′ := ψᾱ²𝓡.

This is literally the sextic analogue of the left side of DR (1.13)/(9.6), under the dictionary:

| DR | here |
|---|---|
| A | Y |
| B | Q |
| β_b | β′ = ψᾱ²𝓡 |
| cubic g̃ | g̃₆ |

The range is A ∈ [√B, B]. The note on ATS for generic c_m (c = conj(η_∞)g̃₆ cancels the twist) stands, so ATS is only claimed for the class below.

## 2. Dunn–Radziwiłł §9, step by step

### 2.1 Structure of the proof of DR Prop 9.2 (read from litE/2109.07463.txt, pp. 51–59)

After inclusion–exclusion (9.7)–(9.11), the sieve-minus-correction expands into three terms:
* D, the diagonal (9.13);
* C, the cross term (9.14);
* T, the trivial term (9.15).

**C and T exist only because the correction (2π)^{2/3}g̃(a)N(a)^{-1/6}Σβ_bN(b)^{-1/6}/(3Γ(2/3)) is subtracted.** Without the correction, the sieve is D, up to the (a,b) = 1 bookkeeping.

**D comes from Poisson in a** (Cor 4.1, (9.19)). Twisted multiplicativity (7.13) turns g̃(b₁)·conj g̃(b₂) times the Gauss sum of the character (·/b₁)conj(·/b₂) into c̃r_d(k)·(dπ²γk/b₁)(dπ²γk/b₂)‾. **The Gauss sums cancel completely,** and the dual is a sieve for β itself with k ≤ B²/A. Solver D's Claim 4 is the same computation; it gives Off ≍ (Y/Q)Σ_{k≤K}|Σ_m β′(m)χ_m(k)|².

D then splits:
* **D₁** (perfect-power frequencies dπ²γk = cube). Main term ∝ A^{2/3}|Σ_b β_b N(b)^{-1/6}|² (9.37), coming from the radial Poisson in h (9.25)–(9.35). The non-coprime part D₁˛ (d ≠ 1) is bounded by A^{2/3}B^{5/3}/w using w-roughness (9.21).
* **D₂** (non-power frequencies) is bounded **solely** by Axiom 4 of Def 3.1, i.e. (3.2), applied to each Mellin-separated b-sum, with the sum over k ≪ Z = X^ε(1 + B²/(N(d)A)) estimated trivially. This gives (9.40): D₂ ≪ AB^{2η} + B^{2+2η}.

**Key observation.** With ‖β‖² ≍ B, the bound D₂ ≪ B^{2+2η} is B·‖β‖²·B^{2η}. So the dual sieve is handled by **pointwise** square-root cancellation together with a trivial count of the K = B²/A frequencies. This is sharp exactly because (A/B)·K·B = B². All the depth of DR's improvement is in the C-term (Prop 8.1, GRH), which is needed only to subtract the Patterson bias.

### 2.2 Which hypotheses of Prop 9.2 exclude our vector, and whether each matters

| Hypothesis | Status for β′ = ψᾱ²𝓡 | Essential? |
|---|---|---|
| (1) \|β\| ≤ 1 | fails mildly: 𝓡 is pointwise ≤ D^{1/2}, but mean-square O(log) | no: only ℓ² is used in D₁/D₂ |
| (2) squarefree, w-rough support | squarefree yes, rough no | **no for twisted β**: D₁˛ main term ∝ Σ_b β′_{db}N(b)^{-1/6} has the angular factor ᾱ²(db) and cancels for each d (see 2.3(iii)); DR need roughness only because their D₁˛ has a genuine main term |
| (4), clause ℓ ≠ 0 at ℓ = +2 (undoing ᾱ²), k a cube | fails: Σ_b ψ𝓡(b) has the μ = □ terms principal if ψ is principal on them | **no**: in the proof of Prop 9.2, (3.2) is used only in D₂ with ℓ = 0 (V^† is radial, Lemma 4.4; separation is via N(b)^{-s}). The ℓ ≠ 0 clause is not invoked in 9.2 (my reading of pp. 51–59 and (7.20); flagged U-E1) |
| (4), clause ℓ = 0, with η ≤ 1/4 | **fails** | **yes**, see below |

**The essential failure.** Σ_b β′_b χ_k(b)N(b)^{it} = Σ_μ ρ_μ Σ_m ψᾱ²(m)(μ/m)₂χ_m(k)N(m)^{it}. Every inner sum is a non-principal Hecke character sum, non-principal even for k a sixth power because of ᾱ². So GRH gives ≪ Q^{1/2+ε} per μ, and in total Q^{1/2}Σ_μ|ρ_μ| ≍ Q^{1/2}D^{1/2}. That is η_eff = log D^{1/2}/log Q = 0.40 > 1/4.

Without GRH, (3.2) is not available at all for k up to K. The D₂ bound then becomes B^{2+2η} = Q²·D, giving E_a exponent 2.33 against the target 1.666. That is worse than trivial (expo.py).

**Splitting 𝓡 dyadically does not help.** Each dyadic piece μ ~ D_j carries the full ℓ² mass (Σ_{μ~D_j}|ρ_μ|² ≍ 1). The GRH bound Q²D_j beats the generic Q²P_a^{1/6} only for D_j < Z^{.0144}.

### 2.3 What the DR mechanism does prove: Theorem E1 (conditional on GRH; sketch)

Let b_m = ψ(m)ᾱ(m)^{ℓ}μ²(m)g̃₆(m)W(N(m)/Q), with ℓ ≠ 0 (ℓ = 2 in practice) and ψ a fixed ray class character. Then for Q^{1/2} ≤ H ≤ Q,

  Σ_{a~H} |Σ_m b_m χ̄_a(m)|² ≪ Q² Z^ε  (= (H+Q)‖b‖² Z^ε).

Steps:
* **(i) Positivity.** Extend a to all of 𝒪 in a ray class mod 36 (as in solverD 2.4). Then apply Poisson in a.
* **(ii) Frequency k = 0.** χ̄_{m₁}χ_{m₂} must be principal, so m₁ = m₂ up to units (squarefree). This gives H‖b‖².
* **(iii) Coprime m₁, m₂, k ≠ 0.** By the sextic analogue of (9.19), using twisted multiplicativity of g̃₆ (B1.9) and γ₁γ₋₁ = χ_m(−1), this term is ≍ (H/Q)Σ_{0<N(k)≪K}|Σ_m ψ′ᾱ^{ℓ}μ²W·χ_m(k)|² (solverD's normalization of Off). Each inner sum is a smooth sum of μ² times a Hecke character of infinity type ℓ ≠ 0. That character is **never principal**, including for k a sixth power. So under GRH each inner sum is ≪ Q^{1/2}N(k)^ε, uniformly in t and ℓ as in DR §6, and the total is (H/Q)·K·Q^{1+ε} = H·K·Z^ε = Q^{2}Z^ε, because K = Q²/H.
* **(iii′) Common factor d = (m₁,m₂).** Rescale m_i → dm_i. This is the same computation with the cube-frequency term weighted by Σ_b β_{db}N(b)^{-1/6}; that weight is again an angular Hecke sum (ᾱ^{ℓ}(d) factors out), so there is no main term. Its D₂ part is bounded by GRH with the "u | a" clause of (3.2), exactly as in DR (9.38)–(9.40), and summing over d costs Z^ε.

**Uncertain (U-E2).** The full sextic local table for non-coprime pairs (solverD U3). I did not write it out. DR handle the cubic case through the same rescaling.

This is the sextic, twisted form of DR Prop 9.2 **without** a correction term. Twisted vectors are *easier*: C, T and the bias vanish identically. Theorem E1 is consistent with the numerics: ᾱ²·Gauss drops to the random level.

It does **not** cover our 𝓡, by 2.2.

### 2.4 Conditional result for our vector

**Hypothesis PSR.** Uniformly for 0 < N(k) ≪ KZ^ε, all ray classes, and smooth weights W:

  |Σ_m ψ(m)ᾱ(m)²𝓡(m)W(N(m)/Q)χ_m(k)| ≪ Q^{1/2}Z^ε.

Under PSR, together with solverD's U1–U3, step (iii) above gives Off ≪ Q²Z^ε. Hence E_a ≪ Q²Z^ε, c = 0, and θ ≈ 0.865152 (gain ≈ 9.7e-3; solverD/gram.py figures, not recomputed).

**Numerics (solverE/psr.py, split-prime proxy, Q = 1500, K = 30000 = 20Q):** max_k|B(k)|²/‖b‖² against the mean.

| Vector | mean | max | GRH-per-μ scale (Σ\|ρ\|)²/Σ\|ρ\|² |
|---|---|---|---|
| D = 60, θ-proxy | 0.97 | 8.0 | 11.6 |
| D = 400, θ-proxy | 0.94 | 7.5 | 71 |
| random | 0.94 | 8.2 / 11.0 | — |
| GAUSS (untwisted) | 1.04 | 10.6 | — |

The maximum is Gaussian-like (≈ log #k = 9.1) and does not grow with D. So PSR looks true, and the GRH loss D^{1/2} is an artefact of the method.

PSR is **stronger than GRH**. B(k) = Σ_μ ρ_μ S_μ(k) with |S_μ(k)| ≍ Q^{1/2} under GRH, and PSR needs the phases of S_μ(k) (that is, of the L-functions L(ψᾱ²χ_k(·/μ)₂)) to be incoherent across the ≍ D characters μ, uniformly in k.

**PSR in mean square over k is circular.** Σ_k S_{μ₁}(k)S̄_{μ₂}(k) = Σ_{m₁,m₂}(…)Σ_kχ_{m₁}(k)χ̄_{m₂}(k). Poisson in k has dual length Q²/K = Y and returns the Gauss-sum vector. This is exactly F1: Poisson is an involution, and no information is gained.

## 3. The alternative route: Heath-Brown recursion and twisted Kubota series

### 3.1 Where c = 1/6 comes from (sketch)

By duality, E_a ≤ ‖T‖²·𝔅(Q,Y). Apply Poisson in m (modulus a₁a₂ ~ Y², dual length P_a), then group h = j⁶h′ (HB's cube-grouping, done sextically; BGL Thm 1.3 for the inner norm). This gives

  𝔅(Q,Y) ≪ Z^ε[ Q + (Q/Y)·max_{H′≤P_a}(P_a/H′)^{1/6}·𝔅(H′,Y) ]
         ≪ Z^ε[ Q·P_a^{1/6} + Q^{1/3}Y + Q ].

The terms come from:
* Q·P_a^{1/6}: the "Y" (Cauchy) part of 𝔅(H′,Y) at H′ ≈ 1;
* Q^{1/3}Y: the (H′Y)^{2/3} part at H′ = P_a. It is ≤ Q iff Y ≤ Q^{2/3}, which holds here (.460 < .555).

So E_a ≪ Q²P_a^{1/6}, i.e. **c = 1/6, the paper's value**.

Sanity check: Q^{5/6}Y^{1/3} = Q·P_a^{1/6} is exactly the conjectured optimal generic norm Ξ₆ of DFDH (1.8) at (A,B) = (Y,Q). So c = 1/6 is (conjecturally) the best possible **for general vectors**, which matches solverC's Prop. O. The cubic version of the same computation gives (QY)^{2/3}, consistent with DR Thm 1.4.

**Uncertain (U-E3).** I did not write out non-squarefree h′ (HB Lemma 3.4 type) or gcd(a₁,a₂) > 1. These are standard and cost Z^ε.

### 3.2 What c < 1/6 needs on this route

One must beat the "Y" part for H′ ≤ P_a^{1−δ}, i.e. show Σ_a F_a g̃₆(a)χ_a(h′) has cancellation for small h′. Write

  Σ_a F_a g̃₆(a)χ_a(h′) = Σ_m T(m)·𝒦_Y(m,h′),  with 𝒦_Y = Σ_{a~Y} g̃₆(a)χ_a(h′)χ̄_a(m)W.

The a-sum is **untwisted**: a carries no angular character. So 𝒦_Y has its Patterson pole, with main term ≈ Y^{2/3}·(sextic theta coefficient at h′m⁵).

**(a) Pole part.** One needs Σ_m γ₁(m)ᾱ²(m)ν(m)𝓡(m)·τ₆(h′m⁵) ≈ 0.
* *Cubic analogue: works.* Patterson's explicit τ₃(squarefree m) ∝ conj g̃₃(m)N(m)^{-1/6} cancels the Gauss factor and leaves Σ_m ᾱ²(m)(μ/m)₂χ_m(h′)·smooth. That is a smooth Hecke-character sum with conductor ≤ P_aD = Z^{.75} < Q = Z^{.83}, which is **negligible** (functional equation: the dual length is < 1).
* *Sextic: blocked.* τ₆ at squarefree arguments is unknown for n ≥ 4 (DIPP §1.1.3: the residue at s = ½ + 1/r is open for r ≥ 4). The angular twist alone does not kill it, because γ₁τ₆ is not known to be angularly unbiased.

**(b) Non-pole part.** Cauchy over m with Lindelöf-on-average for Kubota (the DdFDS Lemma 6.6 analogue) gives |·|² ≲ Q²Y, the same as the generic Y‖F‖². So no gain without using T's structure inside the Voronoi dual, and that dual is again the long k-form (circular).

### 3.3 Twisted Kubota series for n = 6 (what is needed)

The twisted sextic Kubota series is **not** the bottleneck.
* In F1 → F2 (Poisson in a), only Hecke characters appear. Solver D's U4 (elementary) suffices.
* In F3 (Poisson in m), the series that matters is the untwisted one, in a.
* Using the twisted series pointwise in a (F_a = Σ_μρ_μ𝒦^{(ℓ=2)}_Q(a⁵μ³)) loses Y·D/Q = Z^{.29} even if Lindelöf is granted (expo.py: E_a exponent 1.958).

For the record, "ℓ ≠ 0 ⇒ entire" for n = 6 should follow from the following argument (my inference, U-E4, not in the literature I have):
1. The ℓ-twisted Kubota series is a Fourier coefficient of the Eisenstein series on the 6-fold cover induced from the torus character that is (z/|z|)^ℓ on the compact part M.
2. Its constant-term factor is a ratio L^S(6s−·, (·/|·|)^{6ℓ}) / L^S(6s−·+1, ·), the analogue of DR's ζ(3s−½; 1_r, 3ℓ) normalization in Prop 5.2 and of DDHL Prop 4.3.
3. That ratio is entire for ℓ ≠ 0.

What would have to be checked is metaplectic Langlands constant-term theory for non-spherical M-types: the cocycle on M ⊂ SU(2) is trivial at the complex place, so the character extends.

## 4. Uncertain steps (all flagged)

* **U-E1.** That (3.2) with ℓ ≠ 0 is not used inside Prop 9.2's proof. I read pp. 51–59; the only (3.2) invocation is at (9.40) and is radial/N(b)^{-s}.
* **U-E2.** The sextic non-coprime local table in Theorem E1 (iii′).
* **U-E3.** Non-squarefree h′ and gcd terms in §3.1.
* **U-E4.** n = 6 twisted Kubota series entire (inference only).
* **Normalization of ρ_μ.** I used |ρ_μ| ≍ N(μ)^{-1/2}, from FINDINGS' proxy. If the true theta weights were N^{-1/6}-normalized, ‖𝓡‖ changes and the analysis must be redone (the obstruction only gets worse).
* **solverD's U1–U3** are inherited by every conditional statement here.
* **Numerics:** split primes only, small Q, sharp cutoffs. This is evidence, not proof.

## 5. Exponents and probe c

| Route / input | E_a exponent | c | Status |
|---|---|---|---|
| Target | 1.6658 | 0 | needs PSR (stronger than GRH) |
| Theorem E1 (pure twisted Gauss, GRH) | 1.6658 | 0 | does not apply to 𝓡 |
| One-step HB recursion, generic | 1.6802 | 1/6 | = paper; = DFDH Ξ₆; sketch |
| BGL Thm 1.3, a-form | 1.6946 | 1/3 | proven |
| DR D₂ via GRH per μ | 2.3315 | — | useless |
| Twisted Kubota pointwise + Lindelöf | 1.9583 | — | useless |

**Net: unconditional c = 1/6 (no gain). Conditional on PSR: c = 0, giving θ ≈ 0.865152.**

## 6. The precise obstruction

After the DR/solverD reduction, ATS for our vector is equivalent (Poisson involution) to a sextic sieve with more arguments than moduli: Σ_{k≤Q²/Y}|Σ_m ψᾱ²𝓡(m)χ_m(k)|² ≪ (Q²/Y)Q. This form can be proved only through **pointwise** cancellation (DR's mechanism). For β′ = ψᾱ²𝓡 that pointwise cancellation requires incoherence of ≍ D = Q^{0.8} GRH-sized Hecke sums across the quadratic twists μ, which goes beyond GRH.

On the other side of the duality (Poisson in m), it requires non-correlation of T = γ₁ᾱ²𝓡 with **sextic** theta coefficients τ₆(h′m⁵) for h′ ≤ P_a, and τ₆ is not known explicitly. In the cubic analogue, step (a) of the same argument goes through via Patterson's explicit τ₃, so the obstruction is genuinely sextic.

The angular twist does everything it can: it kills every main term (D₁, the d ≠ 1 terms, the Patterson bias). What remains is an off-diagonal "Type II" problem that no available tool reaches.
