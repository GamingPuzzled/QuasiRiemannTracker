# Literature G (arXiv): metaplectic/large-sieve inputs for E_a / ATS / PSR / (★)

Search date: 2026-10-11. arXiv was reachable this time; the export API worked over https only. Sources were saved
under `scratchpad/arxiv/<id>/` (PDF, LaTeX, pdftotext) and treated as untrusted data. I read the full LaTeX for
2403.13151, 2604.01316 and 2610.04045, the full text of 1312.0568, and the relevant parts of 2607.07911. For the other
items I read only the abstract, and those entries are marked [abs].

Exponent conventions are those of FINDINGS/solverD/solverE. Z-exponents: y = .45971, q = .83288 (Q = M),
P_a = Y²/Q = Z^{.0865}, K = Q²/Y = Z^{1.206}, D = Q²/Z = Z^{.666}. The target is E_a ≪ Q²Z^ε (c = 0, θ ≈ 0.865152).
Script: `scratchpad/scripts/dfexp.py`.

---

## 0. Bottom line

1. **No paper proves the input we need.** Nothing found proves ATS/PSR, non-correlation with τ₆, (★), or a
   1/(L·L) mean value. Unconditionally, θ is unchanged.
2. **A new and important item, not in literatureE/F: de Faveri, arXiv:2610.04045 (2 Oct 2026), "Optimal large sieve
   for fixed order characters".** It proves the DFDH-conjectured generic norm for every n ≥ 3 (Thm 1.1). Our
   bookkeeping changes as follows:
   * The "best sextic sieve = BGL (M+N+(MN)^{2/3})" statement in FINDINGS (Tick 3; solverD table) is **outdated**.
   * Solver E's "c = 1/6 via one-step HB recursion (sketch, U-E3)" is now a **theorem** (a black-box corollary of Thm 1.1).
   * Because Thm 1.1 is optimal for general vectors (proved for n = 3, 4; conjectured for n = 6), c = 1/6 is the
     generic barrier, confirmed from the other side. Its proof puts the *entire* gap c = 1/6 → 0 into one step, its
     Lemma 6.1 ("trivial estimation of n-th powers"). That step is where a structured-vector argument would have to win (§5).
3. **Correction to solverE §3.2 / FINDINGS "τ₆(h′m⁵)".** Done carefully (§3.3), the Kubota index in the pole part
   is **h′·m** (linear in m), not h′m⁵. This matters for two reasons:
   * Kazhdan–Patterson gives τ₆(M p⁵) = 0 for p ∤ M (Bröker–Hoffstein §1). With the m⁵ index, the pole part would
     wrongly appear to vanish for free.
   * With the correct index, τ₆ is evaluated at *squarefree* arguments. Those are exactly the coefficients that
     Kazhdan–Patterson leave undetermined.
4. **Bröker–Hoffstein Conj. 5.7 gives τ₆(π)² ≈ ζ₆^{k(π)}·g₃(1,ε²,π)/|π|^{1/2}·(small algebraic factor), with k(π)
   "patternless".** So τ₆ at a prime is a square root of a *cubic* Gauss sum with an unexplained unit. This explains why
   the repo's theta_id experiment found "not a simple Gauss expression". It also suggests a sharper test (§6).

---

## 1. Dunn, arXiv:2403.13151 (ANT 19:9, 2025): "Metaplectic cusp forms and the large sieve"

**Setting.** f is a fixed *cuspidal* Maass form on the cubic cover (Γ₂ ⊂ SL₂(Z[ω])) with coefficients ρ_f(ν). The
paper notes that the only non-cuspidal Hecke eigenform is Patterson's ϑ₃.

**Exact statements** (X ≍ AB; W_K smooth; ‖·‖ are ℓ² norms):
- **Thm 1.1 (mainthm).** 𝒫_f(X,v,u;W_K) := Σ_{λ³ν≡u (v)} ρ_f(ν)Λ(λ³ν)W_K(N(ν)/X) ≪ (XKN(v))^ε K⁸N(v)⁴ X^{1−1/34}.
  Cor 1.2 gives X^{1−1/578+ε} for sharp cutoffs. This is the cubic analogue of Duke–Iwaniec.
- **Lemma 1.3 (trivpoint).** The pointwise Type I sum Σ_{ab≡u(v)} ρ_f(λ^{-3}ab)W_K is ≪ K⁴N(v)^{1/2}N(a)^{1/2}, by Voronoi.
  This is acceptable for A ≪ X^{2/3}.
- **Thm 1.4 (type2a), the "bilinear large sieve".**
  ℬ_f = ΣΣ_{ab≡u(v)} μ²(a)α_aβ_bρ_f(λ^{-3}ab)W_K ≪ (XKN(v))^ε K⁸N(v)⁴((AB)^{1/2} + A^{3/2}B^{1/4})‖μ²α‖_∞‖β‖₂.
  The proof uses the Browning–Vishe circle method, Voronoi on both ν-sums, and Ramanujan sums. The diagonal
  (s₁ = s₂ | (a₁,a₂)) is handled only by averaging over a₁, a₂ with Rankin–Selberg, because no Deligne-type bound is
  known for ρ_f. That is why ‖α‖_∞ appears.
- **Prop 1.5 (type1avgbd).** 𝒜_f = Σ μ²(a)α_aρ_f(λ^{-3}ab)W_K ≪ K^{14/3}N(v)^{5/6}(AX)^{1/3}‖μ²α‖₂.
  The proof applies Voronoi in b, which gives the bilinear form (X/A²)ΣΣ μ²(a)ḡ(a)α_aρ_f(ν)(ν/a)₃. Thm 5.1 (HBcubic) then
  applies Heath-Brown's cubic sieve with one variable not squarefree, in the form
  Σ_d|Σ_c μ²(c)Ψ_c(c/d)₃|² ≪ M^{1/3}(M+N)‖Ψ‖².

**Relevance and transfer.**
- **The kernel is a metaplectic coefficient evaluated at a product ab, not a residue symbol χ_a(m).** It is a Type I/II
  machine for Σ ρ(n)·(arithmetic weight). The closest of our targets is the pole part of solverE §3.2, where we need
  Σ_m c_m τ₆(h′m) to be small. Here c_m = g̃₆(m)ᾱ²(m)ψ(m)𝓡(m).
- **Cubic → sextic.** Voronoi holds for any automorphic form on the 6-fold cover, including the residual ϑ₆. The dual
  coefficients are ϑ₆'s coefficients at other cusps; these are unknown, but are controllable in mean square by
  Rankin–Selberg (DFDH 2607.07911 Prop 6.x does this for n = 4 and "expects" it for n ≥ 5). Lemma 1.3 and Prop 1.5
  therefore plausibly transfer to ϑ₆ for *smooth* or *Type I* weights.
  Two differences: ϑ₆ is non-cuspidal, so the Rankin–Selberg diagonal acquires a pole main term; and the sextic large
  sieve now comes from de Faveri (§4), not HB.
- **Our weight does not fit.** c_m contains g̃₆(m), which is twisted-multiplicative and not a factorable α*β. 𝓡 is a
  long quadratic-character polynomial. Splitting m = m₁m₂ produces a residue-symbol coupling χ_{m₁}(m₂)², not Dunn's
  product kernel. Thm 1.4 does not apply.
- **Hypothetical exponent.** Even if c_m factored as α*β with ‖α‖_∞ ≪ 1, Thm 1.4 saves over trivial only when
  A^{3/2}B^{1/4} ≤ (AB)^{1/2}, i.e. A ≤ B^{1/4}. Its diagonal treatment also needs cuspidality.
  **Assessment: c unchanged (1/6). No adaptation I can see.** Its value here is as a template: Voronoi plus a large
  sieve on the dual, which is what §5 would need for ϑ₆.

## 2. Castillo–de Faveri–Dunn, arXiv:2604.01316 (v2, Jul 2026): "Non-vanishing for quartic Hecke L-functions…"

**Exact statements of the twisted quartic Gauss sums section (§5).**
- Notation: ψ_α(r,s,ξ;v) = Σ_{c≡v(4),(c,α)=1} g₄(r,c)ξ(c)N(c)^{-s}. Here ξ is the angular character (m̄/|m|)^ω; the
  series converges absolutely for Re s > 3/2.
- **Lemma 5.1 (coprimality removal, i–iv), Cor 5.2.** Exact identities express ψ_{abcd}(ab²c³r,…) through ψ(·) with
  Euler-type factors Δ_α, Δ*_α. They follow from the local Gauss-sum evaluations at p^j, j ≤ 4.
- **Lemma 5.3 (poles and convexity).** ψ(r,s,ξ;v) is meromorphic on C.
  * If **ω ≠ 0 it is holomorphic for Re s > 1.**
  * If ω = 0 it has at most a simple pole at s = 5/4.
  * Uniformly for 1+ε ≤ Re s ≤ 3/2+ε: ψ ≪ (N(r)^{1/2}(|s|³+|ω|³+1))^{3/2−Re s+ε}.
  * The proof cites **Diaconu, Invent. Math. 157 (2004), Thm 2.1**, which treats r-th order characters in general,
    and DDHL Prop 4.3.
- **Suzuki relations (5.8)–(5.11).** τ₄(m⁴δ) = τ₄(δ); τ₄(m³δ) = 0; τ₄(m²δ) = ±ḡ₄(δ,m)N(m)^{-3/4}τ₄(δ;mv); τ₄(mδ) = 0
  unless a quadratic condition holds.
- **Lemma 5.4.** τ₄(r) ≪ N(r)^{1/8+ε}.
- **Lemma 5.5 (truncated twisted sums).** Σ_m ξ(m)χ̄_m(ηab²c³d⁴)g₄(m)N(m)^{-s}H(N m/M) = 𝒯_M + ℛ_M.
  * The main term 𝒯_M is ∝ **1_{ω=0}·1_{c=1}**·M^{5/4−s}·τ₄(aη;bv).
  * The error satisfies ℛ_M ≪ M^{1−σ+ε}·∫|ψ(…,1+ε+i(t+y))|(1+|y|)^{-A}dy.
- **Moment results.** Prop 1.3 (first moment) and Prop 1.4 (second moment, error X^{1−1/70}) for mollifiers of length
  X^{tiny}. They use only convexity and the BGL quartic sieve (Thm 4.14 there; now superseded by de Faveri).
  The authors say "our method should generalize to higher order characters".

**Relevance and transfer.**
- **U-E4 (twisted sextic Kubota series is entire for ℓ ≠ 0).** The n = 4 statement is Lemma 5.3. The cited source,
  Diaconu 2004 Thm 2.1, is stated for general order. So the sextic case is very likely proved in the literature.
  I could not read Diaconu (not on arXiv), so this is not verified. It was not the bottleneck anyway (solverE §3.3).
- **Lemma 5.1 transfers mechanically to n = 6** using the Kazhdan–Patterson local relations (§3).
- **What it would give us.**
  * This is the sextic analogue of the pointwise twisted sum F_a = Σ_μ ρ_μ Σ_m ψᾱ²g̃₆(m)χ̄_m(a^{±1}μ³), or S(d) in
    §5 route B. With convexity, N(index)^{1/4} with index ≈ aμ³ is far too lossy. Even Lindelöf gives an E_a
    exponent of 1.958 (solverE expo.py).
  * The mollified moments are central-point results with mollifier length X^{1/(100H)}, which says nothing at σ ≈ .73.
    They are no help for the 1/(L(ψ̄)L(ηψ̄)) target.
  * **Assessment: c unchanged; λ unchanged.**

## 3. Bröker–Hoffstein, arXiv:1312.0568 (Math. Comp. 2016), and the CFH conjecture

### 3.1 Exact statements (B–H §1, §5; τ normalized with τ(1) = 1, bad primes suppressed in §1)
- **Hecke operators (Kazhdan–Patterson).** For (m,p) = 1 and 0 ≤ j ≤ n−2:
  τ_n(mp^j) = G_{j+1}(m,p)τ_n(mp^{n−2−j}). Also τ_n(mp^{n−1}) = 0, and periodicity gives τ_n(mpⁿ) = N(p)^{1/2}τ_n(m).
- **n = 6.**
  * τ₆(Mp⁵) = 0.
  * τ₆(Mp⁴) = G₁(M,p)τ₆(M).
  * τ₆(Mp) = G₂(M,p)τ₆(Mp³), and τ₆(Mp³) = G₄(M,p)τ₆(Mp). This pair is a consistent loop, so it is **undetermined**.
  * τ₆(Mp²) = G₃(M,p)τ₆(Mp²). This is a quadratic self-relation that forces vanishing when the quadratic symbol
    is ≠ 1 (Lemma 5.3).
- **Lemma 5.1 (proved).** For |π| ≡ 1 mod 4, with the representative set V of (5.1):
  τ(π⁴,V)/τ(1,V) = |π|^{-1/2}g₆(1,ε,π). Also τ(1,V) ≈ 0.1358547858696091 numerically.
- **Conj 5.2, 5.4–5.6 (numerical).** τ(π²,V)/τ(1,V) equals 2g₃(1,ε²,π)/|π|^{1/2}, or the inert/7 mod 12 variants,
  up to a 6th/12th root of unity whose exponent "has no pattern found". τ(π²) = 0 for |π| = 37, 313, 373, 661, 769.
- **Conj 5.7.** (τ(π,V)/τ(1,V))² = ζ₆^{k(π)}·g₃(1,ε²,π)/|π|^{1/2}·3^{l(π)}·(1−3ζ₆)^{m(π)}(−2+3ζ₆)^{n(π)}. The
  13-norm variant applies for |π| ≡ 3 mod 4. Constraints: l ∈ {−1,0}; m,n ∈ {0,2}, not both 2; k(π) shows no pattern.
- **CFH conjecture** (Chinta–Friedberg–Hoffstein, Festschrift 2012; not on arXiv). From a double Dirichlet series:
  Σ τ₆(m²)N m^{-u} = Σ τ₃(m)N m^{-u}·Σ G₁^{(3)}(1,d)N d^{-u}. This predicts τ₆(m₁²m₂⁴) in terms of cubic Gauss sums,
  i.e. **square indices only**.

### 3.2 What this means for our targets
- Our pole part needs τ₆ at **squarefree** indices h′m (§3.3). By the relations above, these are precisely the
  *undetermined* coefficients. Only their squares at primes are conjectured. The values at composite squarefree m are
  not multiplicative, and nothing (conjectural or proved) pins them down.
- **No proof route via an explicit τ₆.** The cubic mechanism (Patterson's τ₃ ∝ ḡ₃ cancels γ₁) cannot be copied.
- **Heuristic support for PSR.** Conj 5.7 says τ₆(π) ≈ ±√(unit·g̃₃(π)), with the unit and sign showing no detectable
  pattern. If those phases are independent of Hecke characters and of g̃₆(π), the correlation of γ₁ᾱ²𝓡 with τ₆(h′·)
  should have square-root cancellation. That is consistent with the repo's PSR numerics. It is evidence, not proof.

### 3.3 Index correction (affects solverE §3.2 and FINDINGS FINAL STATUS)
Write F_a = Σ_m c_mχ̄_a(m). Then E_a = Σ_m c̄_m G(m), where G(m) := Σ_a F_aχ_a(m); the dual kernel is the conjugate.
Poisson in m (modulus a₁a₂, coprime) gives g(χ_{a₁}χ̄_{a₂},h) ∝ [χ̄_{a₁}(h)g(χ_{a₁})]·conj[χ̄_{a₂}(h)g(χ_{a₂})].
The dual form is therefore Σ_h|Σ_a F_a g(χ_a,1)χ̄_a(h)|². The a-summand per m is

  c_m·g(χ_a,1)χ̄_a(h)χ̄_a(m) = c_m·g(χ_a, hm).

So the Kubota series in a has index **hm**. The pole gives the theta coefficient τ₆(h′m), up to sixth powers and
units. Starting from Σ_a|Σ_m c_mχ_a(m)|² instead gives the conjugate theta at index hm, so again linear in m.

The label "τ₆(h′m⁵)" is a slip. If taken literally, KP's τ₆(Mp⁵) = 0 would make the pole part vanish for all m
coprime to h′, which is false.

## 4. NEW: de Faveri, arXiv:2610.04045 (2 Oct 2026), "Optimal large sieve for fixed order characters"

**Exact statements.** The setting is any number field K ⊇ μ_n, with the FHL family χ_a. Here Σ* means squarefree and
Σ♭ means n-th-power-free.
- **Thm 1.1.** For n ≥ 3:
  Θ_n(A,B) := sup_{‖λ‖=1} Σ♭_{N a≤A}|Σ♭_{N b≤B} λ_bχ_a(b)|² ≪_ε (AB)^ε(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}).
  Since Ξ_n ≤ Θ_n, this proves the upper half of the DFDH conjecture (1.8). It is optimal for n = 3 (DR, DFDH) and
  n = 4 (DFDH Thm 1.2). For **n = 6** the bound is A + B + A^{5/6}B^{1/3} + A^{1/3}B^{5/6}.
- **Lemma 6.1 (the key step).** Σ_{a~A}|Σ♭_b λ_bχ_b(a)|² ≪ max_{R≤(2A)^{1/n}} R·𝓛₁(A/Rⁿ,B)·‖λ‖². It is obtained by writing
  a = d·rⁿ and *bounding the r-sum trivially*. The extremal term A(B²/A)^{1/n} is the contribution of near-pure n-th
  power frequencies.
- **Prop 7.1 (iteration).** (E_α) ⇒ (E_{2−(n−1)/(nα−1)}), with fixed point 1 + 1/n.
- **"Coupling of root numbers" (§1.3).** For non-squarefree b, ε(χ_{b₁}χ̄_{b₂}) does not split. It is handled by
  absorbing b₂^♯ into the dual variable.
- **Cor 1.6 (zero density, n-th-power-free family).** Σ♭_{N a≤N} N(σ,T,a) ≪ N^{1−δ_n(σ)+ε}·T^{…}, with explicit
  δ_n. Also Thm 1.5 (second moment, n ≥ 4), with error exponent σ_n = 1 − (n−1)/(2(n+1)(n−2)); σ₆ = 51/56.

**Consequences at the record geometry** (`dfexp.py`):

| Route / input | E_a exponent | excess over 2M = 1.6658 | c |
|---|---|---|---|
| BGL (M+N+(MN)^{2/3}), a-form | 1.6946 | .0288 = P_a^{1/3} | 1/3 |
| **de Faveri Thm 1.1, a-form (A,B) = (Y,Q)** | **1.6802** | **.0144 = P_a^{1/6}** | **1/6** |
| de Faveri, k-form (A,B) = (K,Q), for the vector ψᾱ²𝓡 | 1.7424 | .0766 = Y^{1/6} | worse |

The binding term at (Y,Q) is A^{2/n}B^{1−1/n} = Y^{1/3}Q^{5/6} = Q·P_a^{1/6}. This is exactly the record's B0.22
exceptional term and Solver E's HB-recursion value. So:
- **c = 1/6 is now an unconditional black-box consequence of a published theorem.** It is not an improvement, because
  the record already has c = 1/6. It does close solverE U-E3 and lets B0.22 be cross-checked independently.
- **c = 1/6 is the generic barrier.** Any c < 1/6 must use the structure of the vector. That agrees with solverC/E.

**Other levers re-priced with Thm 1.1** (rough, flagged):
- **(★).** A general sieve over q_k ≍ U with Möbius length D = U^r, 1 < r ≤ 3/2, gives
  max(1, r, 5/6 + r/3, 1/3 + 5r/6) = (2+5r)/6. The record's UM gives (1+5r)/6, which is better by U^{1/6}.
  **No gain** for (★): consistent with "UM beats any general sieve".
- **Sixth-moment numerator (B0.41)**, via the sieve on S³ with dual length U^{.42}: about U^{1.383}. This beats
  4th-moment Hölder only when δm > .383. At binding δm ≈ .20, so **no gain**.
- **Zero density, Cor 1.6 at n = 6, σ = .7347.** Family exponent 1 − δ₆ = .750, against the record's R = .5903
  (ratio 2.83 vs 2.225). **Much weaker; no import.**

## 5. The single most promising adaptation: structured n-th-power removal in de Faveri's recursion

**Where the loss sits.** Run the proof of Thm 1.1 on E_a without taking the sup over λ. With ‖c‖² ≍ Q:

  E_a² ≤ ‖c‖²Σ_m|G(m)|²,  Σ_m|G|² ≈ Q·E_a + (Q/Y)·Σ_{0<h≲P_a}|B′(h)|²,  B′(h) := Σ_{a~Y}F_a g̃₆(a)χ̄_a(h).

Write h = r⁶d with d sixth-power-free; B′(r⁶d) = B′(d) up to coprimality. Then **E_a ≪ Q²Z^ε (c = 0) follows from**

  **(★★)**  Σ_{0<N d≤P_a, d 6-free} (P_a/N d)^{1/6}·|Σ_{a~Y}F_a g̃₆(a)χ̄_a(d)|² ≪ Z^ε·Y·Σ_{a~Y}|F_a|²,
  where F_a = Σ_{m~Q}c_mχ̄_a(m) and c_m ∝ γ₁ᾱ²ν𝓡 (solverD Claim 3).

The generic bound (Thm 1.1 / large sieve) is Σ_{d~H′}|B′(d)|² ≤ (Y+H′)‖F‖². It loses exactly the weight
(P_a/H′)^{1/6}, i.e. P_a^{1/6} at d ≍ 1. So (★★) asks to beat the large-sieve diagonal by the factor (P_a/H′)^{1/6}.
That factor is tiny: at most Z^{.0144}. Any partial saving gives c ∈ (0,1/6) and therefore a record. For example,
gram.py gives c = 1/8 → Δθ ≈ 2.4e-3.

**Unfolding (★★).** B′(d) = Σ_m c_m·𝒦_Y(dm), with 𝒦_Y(r) := Σ_{a~Y}g̃₆(a)χ̄_a(r)W. This is the *untwisted* sextic Kubota
partial sum at index dm (§3.3). Its pole sits at normalized exponent 2/3. Write 𝒦_Y(r) ≈ Y^{2/3}·ϱ(r) + Y^{1/2}·ψ̃(r,½+it),
where ϱ(r) ∝ the τ₆ residue. Heuristically ‖ϱ‖²_{m~Q} ≍ Q^{5/6}: this needs the n = 6 analogue of DFDH Prop 6.x,
Σ|τ|² ≍ A^{1−1/n}, which is "expected" there. At d ≍ 1 the two inputs needed are then:

- **(NC-θ) pole part.** |Σ_{m~Q}c_m ϱ(dm)| ≪ Z^{-η}·P_a^{-1/6}·‖c‖‖ϱ‖. That is a saving of P_a^{1/6} = Z^{.0144} ≈ Q^{-.017}
  over Cauchy–Schwarz in the correlation of c with sextic theta coefficients at squarefree indices.
- **(NC-K) off-pole part.** |Σ_m c_mψ̃(dm,½+it)| ≪ P_a^{-1/12}·‖c‖(Σ_m|ψ̃|²)^{1/2}. That is a saving of Z^{.0072},
  together with Lindelöf-on-average Σ_{m~Q}|ψ̃(dm,·)|² ≪ Q^{1+ε}. The latter is open for n = 6. DFDH prove it for
  n = 4 by Rankin–Selberg and note that the DDDS Lemma 6.6 large-sieve proof also works. With de Faveri's Thm 1.1,
  that large-sieve proof now plausibly goes through for n = 6, which is checkable.

For larger d ~ H′ the required saving shrinks like (P_a/H′)^{1/6} in mean square over d.

**Why this is the best lead.** It turns the vague "non-correlation with τ₆" into a quantitative condition with
**tiny required savings**, at the one step (Lemma 6.1) of a now-rigorous optimal argument.

**Why it is still hard.** Both (NC-θ) and (NC-K) are correlations of c_m with GL₂-metaplectic data, and c_m itself
contains a Gauss sum and a long quadratic polynomial 𝓡 = Σ_{μ≤D}ρ_μ(μ/m)₂ with ρ_μ ≈ g̃₃(μ)α(μ)N(μ)^{-1/2}:
- Expanding 𝓡 turns the full object into a **triple Dirichlet series** in (a, m, μ). Its Gauss sums are g̃₆(a),
  g̃₆(m), g̃₃(μ), with couplings χ̄_a(m) (sextic) and (μ/m)₂ (quadratic): an A₃-shaped chain of mixed orders.
- Per μ, the (a,m) part is a twisted A₂ sextic Weyl-group multiple Dirichlet series (Brubaker–Bump–Friedberg type,
  with an angular twist on m). Its a-residue is Σ_m c_mτ₆(dm)N(m)^{-w}.
- If that residue series had continuation with convexity κ in the μ-aspect, then summing over μ would need
  D^{1/2+κ} ≤ Q^{1/2−η}, i.e. **κ ≲ 0.11**. Pointwise convexity (κ ≈ ½) fails; Lindelöf-on-average in μ would suffice.
- Sawin–Whitehead (arXiv:2507.08662) classify which such axiomatic series have enough functional equations, but over
  **F_q(t) only** (their own caveat).

**Honest assessment.** No exponent below 1/6 is in reach with published tools. The concrete missing theorems are
(NC-θ) and (NC-K) above. If both hold with the stated (small) savings, then c = 0 and θ ≈ 0.865152; partial savings
give intermediate c. This inherits solverD's U1–U3 and the heuristic normalizations flagged above (the pole exponent
2/3 and ‖ϱ‖² ≍ Q^{5/6}).

## 6. Smaller leads and suggested checks

- **Experiment (cheap).** In exp/theta_id, test Bröker–Hoffstein Conj 5.7 directly: is ρ₆(p)²·conj(g̃₃(p)) concentrated
  on a few units times {1, 3^{-1}, 7, …}? This validates the τ₆ proxy. The numerical (NC-θ) can then be measured as the
  ratio |Σ_m c_mρ₆(dm)|/(‖c‖‖ρ₆‖) for d ≤ 50, which should decay like Q^{-1/2+ε} if PSR-type heuristics hold.
- **Lindelöf-on-average for sextic Kubota series.** Prove Σ_{N r≤R}|ψ₆(r,1+ε+it)|² ≪ R^{1+ε} via de Faveri Thm 1.1
  plus the DDDS24 Lemma 6.6 argument. This is a standard but needed ingredient for (NC-K), and plausibly short.
- **Ma, arXiv:2606.26440** [abs]. A bound for the bias of higher-order Gauss sums over F_q[t], extending Patterson's
  conjecture beyond n = 3, 4. Function field only: heuristic support for n = 6 pole sizes.
- **Hamdar, arXiv:2509.06152** [abs]. Cubic first moment off the line, with a phase transition at 1/3. It is a moment
  of L, not 1/L. Not useful.
- **Patterson, arXiv:1411.7594** [abs]. Theta coefficients over F_q(t) for all n, in closed form. If a sextic model is
  needed, an F_q(t) computation of τ₆ at squarefree indices is the place to test (NC-θ) exactly.
- **Not found** on arXiv for 2024–2026: any sieve with angular-twisted coefficients (ATS); any explicit squarefree τ₆;
  any zero-density estimate for sextic families near σ ≈ .73 that beats 2.225; any mean value of 1/(L(ψ̄)L(ηψ̄)) or of
  long sextic-twisted Möbius polynomials. Chinta–Friedberg–Hoffstein and Diaconu (Invent. 2004) are not on arXiv and
  were not read.

## 7. Table: per-target status after this round

| Target | Best new input | Effect on exponent |
|---|---|---|
| E_a / ATS (c) | de Faveri Thm 1.1 | c = 1/6 rigorous (= record); generic optimum; gap localized to (★★) |
| τ₆ non-correlation / PSR | B–H Conj 5.7, KP relations | heuristic support only; index corrected to h′m (no free KP vanishing) |
| (★) λ < 5/6 | de Faveri (general sieve) | gives (2+5r)/6, worse than the record's (1+5r)/6: none |
| MP / numerator higher moment | de Faveri on S³ | U^{1.383} at dual length .42; needs δm > .383: none |
| Zero density, sextic family | de Faveri Cor 1.6 | ratio 2.83 vs the record's 2.225: none |
| 1/(L(ψ̄)L(ηψ̄)) mean value | CdFD mollified moments, de Faveri Thm 1.5 | central point, tiny mollifiers: none |
