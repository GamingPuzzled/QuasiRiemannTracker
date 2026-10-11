# Solver K: (★★), running de Faveri's Lemma 6.1 on the probe vector (B6)

**Verdict: nothing is proved unconditionally, so θ is unchanged.** The report does three things:
* it carries out the exact reduction;
* it corrects the bookkeeping of literatureG/B6 in three places that change the payoff;
* it pins down what is known and unknown about ρ₆ for the pole part.

The main new findings:

1. **Savings at d ≍ 1 alone give nothing (correction to B6, "ANY partial saving gives c < 1/6").** Any saving has to hold for every dual frequency d with N(d) ≤ P_a^{h₀}, for some fixed h₀ > 0. The exponent obtained is
   c = max(1/6 − min(s, h₀/6), 1/7),
   where s is the pole saving in averaged form (§3).
2. **(NC-θ) + (NC-K) with literatureG's sizes do not give c = 0.** Even uniformly in d they give c = 5/42 ≈ 0.119. A perfect pole saving combined with the off-pole part at the Lindelöf-on-average level gives **c = 1/7** (θ ≈ 0.873641, a gain of 1.2e-3), not c = 0. Going below 1/7 needs off-pole savings at N(d) ≈ P_a^{1/7…2/7}, not at d ≍ 1.
3. **The d ≍ 1 normalizations of literatureG are correct as a consistent pair.**
   * The pole needs Z^{.0144} = P_a^{1/6} in amplitude; the off-pole part needs Z^{.0072} = P_a^{1/12}.
   * The two heuristic inputs (pole at Y^{2/3}; Σ_{r~R}|ϱ|² ≍ R^{5/6}) enter only through their product. The upper half of that product follows from de Faveri Thm 1.1 together with the L-avg lemma (sketched).
4. **Pole part, handle (i).** "ϱ₆ is spherical at ∞" does **not** make m ↦ ϱ₆(dm) of trivial infinity type in any usable sense. The torus of the 6-fold cover is non-abelian (Hilbert-symbol commutators), so the Hecke integral along the torus does not produce Σ ϱ₆(r)λ(r)N(r)^{-s}.
   * The correct object is the residue of an ᾱ²-twisted A₂ (with 𝓡: A₃-shaped) multiple Dirichlet series. It is not available for n = 6.
   * Even granting it, the μ-aspect of 𝓡 needs Lindelöf-strength uniformity, which agrees with solverG §4.2.
   * Kazhdan–Patterson relations are vacuous at the binding frequencies.
   * Averaging over d (handle iii) cannot help at the binding point h = 0.
5. **New correction C3 (to solverG).** noncorr.py uses A_m(h′) = Σ_a g̃(a)χ_a(h′)χ̄_a(m). That is the Kubota index h′⁻¹m ≡ h′⁵m, not h′m.
   * Its h′ ≠ 1 statistics therefore test the pole at d = h′⁵. At h′ prime that pole vanishes identically by KP (τ₆(Mp⁵) = 0), which explains solverG's "f(h′≠1) ≈ 0".
   * The binding frequencies d = h′ (index h′m, squarefree, non-vanishing pole) are untested. Only h′ = 1 is relevant as measured.

Scripts (new): `scripts/starstar/gram.py` (copy of gram.py, extended with c-values and fixed Z-savings η), `scripts/starstar/profile.py` (h-profile bookkeeping).

Status labels: **[P]** proved (modulo items listed in §6), **[S]** sketched, **[H]** heuristic.

---

## 1. Setup and the exact reduction

Notation follows solverD/E:
* Y = Z^y, Q = Z^M, P_a = Y²/Q = Z^{y−x}, with y − x = .086535 at the record;
* F_a = Σ_{m~Q} c_m χ̄_a(m), c_m = ξΩT(m) = κνψᾱ²γ₁𝓡(m) on squarefree good rows (solverD Claim 3), ‖c‖² ≍ Q;
* E_a = Σ_{a~Y}|F_a|²;
* probe constant c means E_a ≪ Z^ε Q·P_a^{c}·‖c‖². The record has c = 1/6, and gram.py maps c ↦ θ.

**Step 1 (duality with the specific vector) [P].** E_a = Σ_m c_m·conj G(m), with G(m) := Σ_a F_aχ_a(m). Hence
E_a² ≤ ‖c‖²·Σ_{m~Q}|G(m)|².

**Step 2 (de Faveri's chain on Σ_m|G|²) [P modulo §6 U1–U2].** This is de Faveri's Σ₂(A,B,λ) with 𝔞 ↔ m, 𝔟 ↔ a, A = Q, B = Y, λ = F. Apply his Lemma 4.3 (GCD sorting), Lemma 5.1 (Poisson in m) and Lemma 6.1 (h = r⁶d) *without* taking the sup over λ. His proofs are pointwise in λ until the final sup.

For squarefree a we have B₁ = B₂ = 1 and no root-number coupling. The root number of χ_a is g̃₆(a), and the dual length is N(a₁a₂)/Q ≍ P_a. This gives

  Σ_m|G|² ≪ Z^ε [ Q·E_a + (Q/Y)·sup_{|t|≤Z^ε, δ} Σ_{d 6-free, N d≤P_a} (P_a/N d)^{1/6} |B′_{t,δ}(d)|² ],
  B′_{t,δ}(d) := Σ_{a~Y} F^{(δ)}_a g̃₆(a) N(a)^{-it} χ̄_a(d).

Here F^{(δ)} is the GCD-sorted vector: c_m twisted by χ̄_δ(m)1_{(m,δ)=1} with a ↦ aδ. The Δ-penalties of Lemma 4.3 only help.

**The weight (P_a/N d)^{1/6}** is exactly the r-count of Lemma 6.1. The terms B′(r⁶d) = B′(d)·(1 + coprimality corrections) are genuinely almost equal, so that count is **not** lossy at the level of (★★). The whole question is the size of B′(d) itself.

Feeding E_a back in:

  (R)  if, for every dyadic H′ ≤ P_a,  (Q/Y)·Σ_{d~H′}(P_a/H′)^{1/6}|B′(d)|² ≪ Z^ε( Q·P_a^{e}·E_a + Q²P_a^{2e}‖c‖² ),
       then E_a ≪ Z^ε Q‖c‖² P_a^{max(0,e)}.

(A term K·E_a inside E_a² ≤ ‖c‖²(…) gives E_a ≤ K‖c‖²; a term X·‖c‖² gives E_a ≤ X^{1/2}‖c‖.) This is the precise sense in which **(★★) ⇒ c**. It is literatureG's (★★) with the t-twist and the δ-twists made explicit.

**Kernel.** B′(d) = Σ_m c_m 𝒦_Y(dm), with 𝒦_Y(r) := Σ_{a~Y} g̃₆(a)χ̄_a(r)N(a)^{-it}W(N a/Y). The index is dm, as in C1, and I re-derived it:
* the CRT Gauss sum of χ_{a₁}χ̄_{a₂} at frequency h is [χ̄_{a₁}(h)g(χ_{a₁})]·conj[χ̄_{a₂}(h)g(χ_{a₂})] times a class-constant reciprocity factor;
* F_a contributes χ̄_a(m);
* so the product is χ̄_a(hm).

## 2. Pole/off-pole split and the normalizations

**Facts used [P, classical].** Let ψ̃(r,w) = Σ_a g̃₆(r,a)N(a)^{-w} (Kubota 1969; Kazhdan–Patterson 1984).
* It is meromorphic, absolutely convergent for Re w > 1, with a functional equation of conductor ≍ N(r).
* Its only singularity in Re w > 1/2 is a simple pole at w = 1/2 + 1/6 = 2/3, with residue ϱ(r) (the sextic theta coefficient up to normalization).
* The squarefree-a restriction changes it by an Euler-type factor that is holomorphic there; Gauss sums g(r,p^k) with p ∤ r vanish unless 6 | k.

Hence
  𝒦_Y(r) = c₀(t)Y^{2/3−it}ϱ(r) + R_Y(r),  R_Y(r) = (2πi)^{-1}∫_{(1/2+ε)} Y^{w}ψ̃(r,w+it)W̃(w)dw.

**Lemma L (Lindelöf on average) [S].** For R ≥ Y and |t| ≤ Z^ε:
  Σ_{N r~R, 6-free} |R_Y(r)|² ≪ R·Y·(RY)^ε,  and  Σ_{N r~R}|ϱ(r)|² ≪ R^{5/6+ε}.

*Sketch.* Use the AFE of ψ̃(r,1/2+ε+iu) from Kubota's functional equation, with length X = N(r)^{1/2} on each side; the dual coefficients are Gauss sums g̃(r′,a) with r′ = r·(fixed S-unit) on fixed classes. Then:
* By de Faveri Thm 1.1 at (A,B) = (R, R^{1/2}), Θ₆ ≪ R^{1+ε}. Hence both Dirichlet polynomials have mean square ≪ R^{1+ε}.
* Crossing the pole in the AFE leaves a residual term ≍ ϱ(r)X^{1/6}.
* Bootstrap. Suppose Σ|ϱ|² ≪ R^β. Then L-avg holds with an extra R^{β+1/6}. Inserting this into Y^{4/3}Σ|ϱ|² ≤ 2Σ|𝒦_Y|² + 2Σ|R_Y|² at Y = R, with Σ|𝒦_Y|² ≤ Θ₆(R,Y)·Y = O(R^{13/6}), gives β_new = max(5/6, β − 1/6). Starting from the convexity value β = 4/3, three steps give β = 5/6.

Not written out:
* the vector-valued functional equation for n = 6 with uniformity in u;
* the S-unit and residue-class bookkeeping.

This is literatureG's "plausibly short" item. I agree that it is standard, but it is a sketch.

**Normalization check of literatureG's targets [P given Lemma L].** The pole and the norm enter only via Y^{4/3}·Σ|ϱ|². By Lemma L that product is ≪ Y^{4/3}R^{5/6}, so Cauchy–Schwarz on the pole at d = 1 gives
P_a^{1/6}·Y^{4/3}Q^{5/6}‖c‖² = Y·Q‖c‖²·P_a^{1/3}.
Fed into (R), this gives E_a ≪ Q‖c‖²P_a^{1/6}: **exactly c = 1/6, the generic value.** So the pair (2/3, R^{5/6}) is self-consistent with Thm 1.1. Its lower half, Σ|ϱ|² ≫ R^{5/6}, is equivalent to the sharpness of Thm 1.1's A^{1−1/n}B^{2/n} term for n = 6, which is conjectural (DFDH) and proved only for n = 3, 4.

At d ≍ 1:
* the pole needs |Σ_m c_mϱ(dm)| ≪ P_a^{-1/6}‖c‖‖ϱ‖ = Z^{-.01442}‖c‖‖ϱ‖ for c = 0 (**Z^{.0144} confirmed**);
* the off-pole part, by Cauchy and Lemma L, contributes P_a^{1/12} (**Z^{.0072} confirmed**).

solverG §4.1's Z^{.749} vs Z^{.7635} is the same statement.

## 3. The h-profile: what partial savings actually give

Let H′ = P_a^h, 0 ≤ h ≤ 1. For each dyadic d ~ H′, (R) offers the following bounds on the exponent e(h), where E_a ≪ Q‖c‖²P_a^{e}:

| option | input | e(h) |
|---|---|---|
| generic | Thm 1.1 on the d-sum: Σ_{d~H′}\|B′\|² ≪ Θ₆(H′,Y)E_a ≈ Y·E_a (H′ ≤ P_a < Y^{1/2}) | g = (1−h)/6 |
| off-pole | Cauchy in m per d + Lemma L summed over d ~ H′ (Σ_{d,m}\|R_Y(dm)\|² ≤ τ·Σ_{r≪H′Q}\|R_Y\|² ≪ H′QY) | o = 1/12 + 5h/12 |
| pole, per-d Cauchy, amplitude saving P_a^{-σ} | uses Σ_{d~H′}Σ_m\|ϱ(dm)\|² ≪ (H′Q)^{5/6} | p = 1/6 + h/3 − σ |
| pole, averaged form, saving P_a^{-s} over the theta large-sieve size (H′+Q)(H′Q)^{-1/6}‖c‖² | | p̄ = 1/6 − h/6 − s |

**Combination.** The pole/off split is only usable where the savings are assumed (h ≤ h₀). There

  e(h) = min(g, max(o, pole)),

and e(h) = g for h > h₀. Then c = max_h e(h).

**Proposition 3 [P as bookkeeping; inputs as hypotheses].** Computed in closed form and checked by `profile.py` on a grid of 2·10⁵ points.
* (a) No structure (σ = s = 0): c = 1/6.
* (b) Savings only for N(d) ≤ P_a^{h₀}: c ≥ (1−h₀)/6. **Bounded d (h₀ → 0) gives no power gain.**
* (c) Averaged pole saving s, uniform on [0,h₀], with Lemma L: **c = max(1/6 − min(s, h₀/6), 1/7).**
* (d) Per-d pole saving σ: c = max(1/6 − min(σ/3, h₀/6), 1/7). Per-d Cauchy loses H′ relative to the averaged form, hence σ/3.
* (e) Perfect pole (s = ∞) with Lemma L only: **c = 1/7.** The binding point is h = 1/7, where o = g.
* (f) Perfect pole plus an additional off-pole saving P_a^{-τ}, uniform in d: o − τ crosses g at h = (1+12τ)/7, so c = (1 − 2τ)/7. Numerically c = 0.1400 (τ = .01), 0.1310 (τ = 1/24), 5/42 (τ = 1/12).
* (g) literatureG's (NC-θ) [σ = 1/6] and (NC-K) [τ = 1/12], assumed uniform in d: **c = 5/42 ≈ 0.11905, not 0.** At d ≍ 1 only, they give c = 1/6 by (b).

**What c = 0 really needs.** (★★) itself, i.e. for every h ∈ [0,1):
Σ_{d~P_a^h}|B′(d)|² ≪ Z^ε Y (P_a^{h}/P_a)^{1/6}E_a.
This is a mean-square saving (P_a/H′)^{1/6} over the large sieve at every scale, including off-pole savings growing to P_a^{7/24} over per-d Cauchy at h ≈ 1/2. literatureG's "required saving shrinks for larger d" holds relative to the large sieve, but not relative to the per-d pole/off analysis it proposes.

**θ consequences** (`scripts/starstar/gram.py`; geometry re-optimized; the model baseline 0.874886616 sits 3e-8 below the record because of grid error, and gains are relative to it):

| c | source | θ (model b0) | gain |
|---|---|---|---|
| 1/6 | record / no structure | 0.874886616 | 0 |
| 1/6 − 10⁻³ | e.g. s = 10⁻³ (h₀ ≥ .006) | 0.874852370 | 3.4e-5 |
| 0.16 | s = 1/150, or σ = .02 | 0.874627678 | 2.6e-4 |
| **1/7** | perfect pole + Lemma L | **0.873641042** | **1.25e-3** |
| 5/42 | NC-θ + NC-K uniform in d | 0.872046990 | 2.8e-3 |
| 0 | full (★★) | 0.865152069 | 9.7e-3 |

For a fixed Z-power saving Z^{-η} on the whole of (★★) (uniform in h), the gain is ≈ 0.386·η: η = 10⁻³ gives 3.9e-4, and η = .0072 gives 2.8e-3.

Every one of these beats the record. Even s = 10⁻⁴ (in P_a units, i.e. Z^{9e-6}) beats it at the 3e-6 level. But each needs the pole hypothesis of (c) or (d) **uniformly for N(d) ≤ P_a^{6s}**.

## 4. The pole part: what is known about ρ₆, and the three handles

The pole term is Π(d) := Σ_m c_m ϱ(dm) = κΣ_μ ρ_μ Σ_m ψν ᾱ²(m)·g̃₆(m)·(μ/m)₂·ϱ₆(dm)·Ω. It has to be smaller, in the averaged sense of §3(c) and uniformly for N(d) ≤ P_a^{h₀}, than its Cauchy–Schwarz size by P_a^{s}.

### 4.1 What is actually known about ϱ₆ (proven inputs only)

* **K1. Existence and pole** (Kubota; Kazhdan–Patterson). ϱ(r) = Res ψ̃(r,·) at 2/3. The sextic theta function is the residue of the Kubota Eisenstein series on the 6-fold cover. It is spherical at the complex place (where the cover splits) and square-integrable (residual spectrum).
* **K2. Hecke/periodicity relations** (KP; Bröker–Hoffstein §1). For p ∤ M:
  * τ₆(Mp⁶) = N(p)^{1/2}τ₆(M), τ₆(Mp⁵) = 0, τ₆(Mp⁴) = G₁(M,p)τ₆(M);
  * τ₆(Mp²) = G₃τ₆(Mp²), which forces vanishing when the quadratic symbol ≠ 1;
  * τ₆(Mp) ↔ τ₆(Mp³) form a consistent loop and are undetermined.

  So **on squarefree arguments ϱ₆ is unconstrained by any local relation**, and it is not multiplicative.
* **K3. Units.** ϱ(εr) = ω(ε)ϱ(r) for a fixed character ω of the units, from Γ-invariance. So over primary m, Σ_m is a sum over ideals.
* **K4. Mean square.** The upper bound Σ_{N r~R}|ϱ|² ≪ R^{5/6+ε} follows from Thm 1.1 and Lemma L (§2, [S]). Rankin–Selberg for n = 6 is not in the literature found (DFDH: "expected" for n ≥ 5). The lower bound is open.
* **Not usable.** Bröker–Hoffstein Conj. 5.7 and CFH (square indices only) are conjectural; the brief excludes them, and so does this report.

### 4.2 Handle (i): the angular twist versus "radial" theta

**Claim 4.1 [P, negative].** Sphericity of θ₆ at ∞ (K1) does not imply that Σ_r ϱ₆(r)λ(r)N(r)^{-s} is entire for an angular Grössencharacter λ, and neither does any other known theorem. In particular it gives no cancellation in Σ_m ᾱ²(m)(…)ϱ₆(dm).

*Reason.* For a GL₂ automorphic form, "entire for λ of nontrivial infinity type" comes from the Hecke integral ∫_{K^×\A^×}φ(diag(y,1))λ(y)|y|^s d^×y. On the n-fold cover the preimage of the torus is **non-abelian**: the commutator of lifts of diag(y,1) and diag(y′,1) is a power of the Hilbert symbol (y,y′)₆. So y ↦ φ(s(diag(y,1))) is not K^×-invariant up to a character, and the integral does not unfold to Σ_r ϱ(r)λ̄(r).

This is the same mechanism that makes theta coefficients non-multiplicative. Sphericity only says that the archimedean Whittaker function depends on |ry|. Rotation by e^{iφ} is not in G(K), so it imposes nothing on the arithmetic sequence r ↦ ϱ(r) against arg r.

**What would be needed, and its status [H].** Our pole term is not even of the form Σ ϱλ. It involves g̃₆(m)·ϱ₆(dm), which mixes m as a Gauss-sum modulus with m as a theta index. The natural home is the residue at the a-polar line of the twisted A₂ double series
Σ_{a,m} g̃₆(a)g̃₆(m)χ̄_a(dm)ψᾱ²(m)N(a)^{-w}N(m)^{-s},
together with the μ-node from 𝓡 (an A₃-shaped chain of mixed orders 6/2/3; solverG §4.2's triple series).

The constant-term heuristic goes as follows. The torus character has nontrivial angular part on the roots α₂ and α₁+α₂, so the only polar divisor is the untwisted a-line. Its residue (our Π(d) series) is then entire in s. Continuation for this angular-twisted series at n = 6 is not in the literature found: BBF/BBFH treat stable, untwisted-at-∞ series, and U-E4 is still an inference.

**Even granting entireness, the μ-aspect fails.** Write the per-μ pole sum as Π_μ(d) = Σ_m ψν ᾱ²g̃₆(m)(μ/m)₂ϱ₆(dm)Ω. Suppose |Π_μ| ≪ Q^{σ₀}N(μ)^{κ}. Then the μ-sum costs Σ_{μ~M}|ρ_μ|·max|Π_μ| ≍ M^{1/2+κ}Q^{σ₀}. Each dyadic M ≤ D = Q^{.8} carries full ℓ² mass (Σ_{μ~M}|ρ_μ|² ≍ 1), so small-μ truncation saves nothing.

Requiring Σ_μ ≤ Q^{11/12}P_a^{-s} at M = D:
  σ₀ + 0.8(1/2 + κ) ≤ 11/12 − s(y−x)/M, i.e. with σ₀ = 5/12 (square root in m): κ ≤ 0.1253 − 0.130s.

Suppose Π_μ(d) had a functional equation of conductor ≍ N(μ)^k. Convexity for a smooth sum of length Q then gives |Π_μ| ≪ Q^{5/12}(N(μ)^k/Q)^{1/2}. At N(μ) = Q^{0.8} this meets the requirement iff k ≤ 3/2.

The expected value is k = 3 [H: a GL₃-type residual Eisenstein series, induced from θ₆ ⊗ a twisted GL₁ character]. That would need a saving of Q^{0.6} over convexity at analytic conductor ≈ Q^{2.4}. **The requirement is Lindelöf-strength in the μ-aspect.**

Cauchy–Schwarz over μ with Heath-Brown's quadratic large sieve returns exactly the trivial bound, because D < Q. The second Poisson in μ is the PSR involution (solverE §2.4). This agrees with solverG §4.2.

**Claim 4.2 [P, no-go for local information].** Let x(r) be any sequence that satisfies K2 and K3 and has Σ_{r~R}|x|² ≍ R^{5/6}, for example x(m) := conj(c_m)N(m)^{-1/12}(Q/‖c‖²)^{1/2} on squarefree primary m, extended to associates by ω and to non-squarefree r through the K2 relations. For this x, Σ_m c_m x(m) = ‖c‖‖x‖ (Cauchy–Schwarz equality).

So no argument using only the proven local/unit information K2–K3 together with mean-square information K4 can give s > 0. What is needed is **global** information: the joint equidistribution of g̃₆(m)ϱ₆(dm) against the Hecke characters ψᾱ²(·)(μ/·)₂, uniformly in conductor up to Q^{0.8}. That is the missing theorem.

### 4.3 Handle (ii): Bröker–Hoffstein

Not used, being unproven. For the record, it constrains only τ₆(π)² at primes. Our binding index at d = 1 is a squarefree composite m, where even conjecturally nothing multiplicative holds.

### 4.4 Handle (iii): averaging over d

**Claim 4.3 [P, given Lemma L].** Write X_{H′} := Σ_{d~H′}|Π(d)|². The bound X_{H′} ≪ (H′+Q)(H′Q)^{-1/6}‖c‖² (the "theta large sieve", which is what perfect orthogonality of the rows d ↦ ϱ(d·) would give) yields exactly e(h) = g(h), the generic exponent.

So averaging over d does not, by itself, beat the generic bound at any h. More importantly, at the binding point h = 0 (d ∈ {1} ∪ S-units) there is nothing to average.

What averaging does buy is shape. In averaged form the requirement is a **uniform** saving s for h ∈ [0, 6s] (Prop 3(c)), not a saving that grows with d as in per-d form (d). It is a structural relaxation, not a source of savings.

KP relations (K2) also give nothing here: for N(d) ≤ P_a^{h₀} the pairs with (d,m) ≠ 1 carry a proportion ≪ Σ_{p|d}N(p)^{-1} of the mass, and at d = 1 they are absent.

## 5. Precise obstruction

After Steps 1–2, c < 1/6 is **equivalent** (up to Lemma L, and for c ≥ 1/7) to the following, for some s > 0:

  (Π_s) for all 1 ≤ H′ ≤ P_a^{6s}, all |t| ≤ Z^ε and all GCD twists δ:
      Σ_{d 6-free, N d~H′} |Σ_{m~Q} c^{(δ)}_m ϱ₆(dm)|² ≪ Z^ε P_a^{-2s} (H′+Q)(H′Q)^{-1/6} ‖c‖².

The binding instance is d = 1: |Σ_m ψνᾱ²g̃₆(m)𝓡(m)ϱ₆(m)Ω(N m/Q)| ≪ P_a^{-s}Q^{11/12}. This is a correlation between our reflected cubic-theta vector and the sextic theta coefficients at **squarefree** indices. Those coefficients are exactly the ones Kazhdan–Patterson leave undetermined (K2).

The only mechanism that plausibly kills it is the ᾱ²-infinity type of the twisted multiple Dirichlet series. That mechanism is:
* (a) unproven for n = 6 (Claim 4.1 explains why sphericity alone is not it);
* (b) qualitative per μ, while 𝓡 needs Lindelöf-strength uniformity in μ ≤ Q^{0.8} (§4.2), because each dyadic μ-block carries full ℓ² mass.

Below c = 1/7 a second, independent input is needed as well: off-pole savings at N(d) ≈ P_a^{1/7…2/7} beyond Cauchy and Lemma L. That is a Type-II estimate for c_m against the oscillatory part of the Kubota series.

The numerics (solverG) support (Π_s) at d = 1 with full square-root cancellation. Per C3, they do not test d ≠ 1.

## 6. Uncertain steps

* **U1.** Non-squarefree columns a (solverD U1). In de Faveri's chain they enter as B₁, B₂ > 1. The generic bound gives Q P_a^{1/6}B₀^{-1/6} times their mass fraction (≪ q_s^{-1}), which I believe is negligible but did not write out.
* **U2.** Step 2 applies de Faveri Lemmas 4.3, 5.1 and 6.1 to a fixed vector. The proofs are pointwise in λ, as I checked in the source (the sup is taken only in the definitions of 𝓛ᵢ). The GCD twists δ and the t-integral are carried explicitly. The class splitting by R_𝔠 costs O(1).
* **U3.** Lemma L is a sketch (§2). It uses the n = 6 Kubota functional equation with polynomial uniformity, and it is used for every c < 1/6 statement. Without it, the off-pole part cannot be separated from the pole, and (★★) must be attacked as a whole.
* **U4.** The θ values come from the slack2/gram model, which approximates b0 to about 3e-8. A c-change of 10⁻³ moves θ by 3.4e-5, well above that resolution.
* **U5.** solverD U2–U3 (the exact shape of c_m) are inherited. Only "finite ray character × ᾱ²γ₁𝓡" is used.

---

## CLAIMS

1. **(Reduction; [P] modulo U1, U2.)** With F_a, c_m and E_a as in §1, and B′_{t,δ}(d) := Σ_{a~Y}F^{(δ)}_a g̃₆(a)N(a)^{-it}χ̄_a(d):
   E_a² ≪ Z^ε‖c‖²[Q·E_a + (Q/Y) sup_{t,δ} Σ_{d 6-free, N d ≤ P_a}(P_a/N d)^{1/6}|B′_{t,δ}(d)|²].
   The weight (P_a/N d)^{1/6} is precisely de Faveri's Lemma 6.1 loss. The Kubota index in B′(d) = Σ_m c_m𝒦_Y(dm) is dm (confirms C1).
   Consequence: if for every dyadic H′ ≤ P_a the d ~ H′ part of the sum is ≪ Z^ε·Y·P_a^{c}·E_a + Z^ε·(Y/Q)·Q²P_a^{2c}‖c‖², then E_a ≪ Z^ε Q P_a^{c}‖c‖², and θ is given by `scripts/starstar/gram.py` (table in §3).

2. **(Lindelöf on average; [S].)** For R ≥ Y ≥ 1 and |t| ≤ Z^ε: Σ_{N r~R}|R_Y(r)|² ≪ (RY)^{ε}RY and Σ_{N r~R}|ϱ₆(r)|² ≪ R^{5/6+ε}, where R_Y = 𝒦_Y − c₀Y^{2/3−it}ϱ₆. The proof combines the Kubota AFE, de Faveri Thm 1.1 at (R, R^{1/2}) and the bootstrap of §2.

3. **(Normalizations; [P] given Claim 2.)** Cauchy–Schwarz on the pole at d = 1 reproduces c = 1/6 exactly. At d ≍ 1, c = 0 requires amplitude savings Z^{(y−x)/6} = Z^{.0144} on the pole and Z^{(y−x)/12} = Z^{.0072} on the off-pole part, as literatureG states.

4. **(Correction: d ≍ 1 is not enough; [P] as bookkeeping.)** Assume Claims 1–2. If the pole and off-pole savings are only available for N(d) ≤ P_a^{h₀}, then c ≥ (1−h₀)/6. Hence no power gain comes from bounded d.

5. **(Partial-saving formula; [P] conditional on Claims 1–2 and hypothesis (Π_s) of §5 for H′ ≤ P_a^{h₀}.)**
   c = max(1/6 − min(s, h₀/6), 1/7).
   * In per-d form (amplitude saving P_a^{-σ} over ‖c‖(Σ_m|ϱ(dm)|²)^{1/2}): c = max(1/6 − min(σ/3, h₀/6), 1/7).
   * θ examples: c = 1/6 − 10⁻³ gives 0.874852 (gain 3.4e-5); c = 1/7 gives 0.873641 (gain 1.25e-3).
   * Any s > 0 with h₀ ≥ 6s beats the record.

6. **(Correction: the c = 1/7 floor; [P] as bookkeeping.)**
   * Under Claims 1–2 with a perfect pole term (Π_∞) and no further off-pole input, c = 1/7.
   * literatureG's (NC-θ) and (NC-K), even if they hold uniformly for all N(d) ≤ P_a, give c = 5/42 ≈ 0.119 (θ ≈ 0.872047), not c = 0.
   * c < 1/7 needs off-pole mean-square savings for N(d) ∈ [P_a^{(12c−1)/5}, P_a^{1−6c}].

7. **(Handle (i) negative; [P] for the stated reason, [H] for what replaces it.)** Sphericity of θ₆ at the complex place does not imply cancellation of Σ_m ᾱ²(m)(…)ϱ₆(dm). The torus of the 6-fold cover is non-abelian, so there is no Hecke integral producing Σϱ₆λN^{-s}. The required statement is a property of the ᾱ²-twisted A₂/A₃-type multiple Dirichlet series and is not available for n = 6. Even granting it qualitatively, 𝓡 requires per-μ bounds Q^{5/12+ε}N(μ)^{κ} with κ ≤ 0.125 − 0.13s, uniformly for N(μ) ≤ Q^{0.8} (Lindelöf strength).

8. **(No-go for local information; [P], given that K2 relates only values at M·p^j with M fixed and imposes no relation among values at distinct squarefree arguments.)** Kazhdan–Patterson relations, unit equivariance and mean-square bounds on ϱ₆ are jointly consistent with a sequence that attains Cauchy–Schwarz equality against c. So they cannot yield s > 0.

9. **(Handle (iii); [P] given Claim 2.)** Perfect orthogonality across d (the theta large sieve) yields exactly the generic exponent (1−h)/6 at each scale. Averaging over d therefore gives no saving, and nothing at the binding point d = 1. It only makes the requirement uniform in h (Claim 5) instead of growing.

10. **(Correction C3 to solverG; [P].)** noncorr.py's A_m(h′) = Σ_a g̃(a)χ_a(h′)χ̄_a(m) has Kubota index h′⁵m. For prime h′ ∤ m its pole vanishes by KP (τ₆(Mp⁵) = 0). The h′ ≠ 1 statistics therefore do not test the binding frequencies d = h′ of (★★), whose index h′m is squarefree. Only the h′ = 1 data are relevant.

**Net: θ unchanged.** The concrete missing theorem is (Π_s) for some s > 0 (§5): an averaged non-correlation of ψᾱ²g̃₆𝓡 with the squarefree sextic theta coefficients ϱ₆(dm), for N(d) ≤ P_a^{6s}. Its payoff is c = 1/6 − s, down to 1/7.
