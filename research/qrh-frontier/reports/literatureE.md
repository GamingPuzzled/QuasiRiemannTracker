# Literature search E — large sieve for n-th order symbols with angular / Grössencharacter-twisted coefficients

Target: conjecture **ATS** (FINDINGS.md last section):
for b_m = η_∞(m)·c_m (η_∞ nontrivial angular Hecke character, c_m of Gauss/theta type),
Σ_{N(h)≤H} |Σ_{N(m)~Q} b_m χ_m(h)|² ≪ (H+Q)‖b‖² Z^ε, for Q^{1/2} ≤ H ≤ Q (χ_m sextic, over Z[ω]).

Legend: **[V]** = I read the statement in the source text (arXiv PDF, pdftotext). **[A]** = from an abstract or
secondary citation only. **[I]** = my inference.
Local copies of the PDFs and text: `scratchpad/litE/*.pdf, *.txt`.

**Bottom line:** I found no published or arXiv result that proves a sharp n-th order large sieve (n = 3, 6 or general)
when the coefficients carry a nontrivial infinity-type twist, and no result that implies ATS. The literature does
support the *mechanism* behind ATS. (i) The extra term in the cubic and quartic large sieves comes from the pole of
Kubota's series, and the vectors that attain it are Gauss-sum vectors. (ii) Kubota series twisted by (c/|c|)^ℓ with
ℓ ≠ 0 are **entire** (verified for n = 3 and n = 4). (iii) For a related class of vectors that cancel against angular
and cubic characters (but not the ATS vectors themselves), Dunn–Radziwiłł get a power saving over (AB)^{2/3}, but not
a sharp bound. The DFDH conjecture (1.8) implies that, without a twist, the sextic sieve is **not** sharp in exactly
the ATS range. So any proof of ATS must use the twist in an essential way.

---

## 1. Dunn–Radziwiłł, *Bias in cubic Gauss sums: Patterson's conjecture*
Ann. of Math. 200 (2024), no. 3, 967–1057. arXiv:2109.07463. https://arxiv.org/abs/2109.07463

- **Thm 1.4 [V]** (GRH for Hecke L-functions over Q(ω)). B(A,B) is the operator norm of the cubic large sieve, with
  squarefree a and b, a ~ A and b ~ B, both ≡ 1 mod 3. Then
  (AB)^{2/3} ≪ B(A,B) ≪ (AB)^{2/3+ε} for A ∈ [√B, B²] \ [B^{1−ε}, B^{1+2ε}],
  and (AB)^{2/3−3ε} ≪ B(A,B) on the excluded middle range. The lower bounds are conditional; the upper bounds are
  Heath-Brown's (Remark 1.1).
- **Extremal vectors [V]** (§1.4, proof in §9, eq. (9.3)): β_b = g̃(b)·W(N(b)/B), the normalized cubic Gauss sum times
  a smooth nonnegative bump. The proof (9.4)–(9.5) uses Voronoi summation (Prop 8.1): Σ_b W g̃(ab) has a main term
  ≍ g̃(a)B^{5/6}/A^{1/6}, coming from the pole of the Kubota series. This gives ≫ A^{2/3}B^{5/3} = (AB)^{2/3}‖β‖².
  The paper says: "Any sequence of the form β_b g̃(b) with β_b non-negative and not correlated with cubic symbols
  would provide a counterexample." The result holds for A ≤ (AB)^{1/2−ξ}; the remaining ranges follow by duality
  B(A,B) = B(B,A) (9.2) and monotonicity (Lemma 9.1).
- **Twisted Kubota series, Prop 5.2 [V].** Let R(s;ψ,ℓ) = Σ_{c≡1(3)} g̃(c)(r/c)_3 (c/|c|)^ℓ N(c)^{−s}. Then
  ζ(3s−½; 1_r, 3ℓ)·R(s;ψ,ℓ) "has a unique pole located at s = 5/6, and it is simple" if ℓ = 0, and "If ℓ ≠ 0 the
  Dirichlet series (5.66) defines an entire function." The paper gives an explicit functional equation with gamma
  factors Γ(5/6+|ℓ|/2−s)Γ(7/6+|ℓ|/2−s)/(Γ(s+|ℓ|/2−1/6)Γ(s+|ℓ|/2+1/6)). The Voronoi formula (Prop 8.1, Cor 8.1) has its
  main term multiplied by δ_{ℓ=0} [V]. **This answers item (a) for n = 3: there is no pole for nontrivial infinity
  type.**
- **"Improved cubic large sieve", Props 9.1–9.2 and (1.13) [V].** Let β be in the class C_η(B,w) (Def. 3.1 [V]):
  |β| ≤ 1, supported on squarefree w-rough b ≡ 1 (3) with N(b) ≍ B, and with square-root-type cancellation
  Σ β_a (a/|a|)^ℓ N(a)^{it} (k/a)_3 ≪ (A/N(u))^{1/2+η+ε} for **every ℓ ≠ 0**, and for ℓ = 0 unless k is a cube.
  Then the large sieve for the vector g̃(b)β_b, *minus the correction term*
  (2π)^{2/3} g̃(a) N(a)^{−1/6} Σ_b β_b N(b)^{−1/6} / (3Γ(2/3)), equals
  (A/B)·Σ_{0<|h|≤B²/A, h≠cube} |Σ_b β_b (b/h)_3|² + O((AB)^{2/3}B/w + (AB)^{2/3−δ}B).
  With square-root cancellation and w > (AB)^ε, the paper concludes ≪ (AB)^{o(1)}(AB + B² + (AB)^{2/3−ε}B). The
  authors add: "in order to beat the cubic large sieve, the correction term alone is not enough; we really need to
  know additional information about the sequence β_b."
- **Relevance to ATS.**
  - [I] The correction term is proportional to Σ_b β_b N(b)^{−1/6}, which is ≈ 0 when β carries a nontrivial angular
    twist. So for twisted vectors the Patterson "bias" term disappears, consistent with the numerics in FINDINGS.
  - [I] Prop 9.2 does **not** apply directly to ATS vectors. Its hypothesis (3.2) asks for cancellation against
    (a/|a|)^ℓ for all ℓ ≠ 0. A vector η_∞·(nonnegative) fails this at ℓ = −k, where k is the twist index.
  - [I] Even where it applies, the gain is only (AB)^{2/3−δ} with a tiny unspecified δ, under GRH. It is cubic only,
    with a Gauss-sum-times-orthogonal vector shape. **It helps (as a template) but does not prove ATS.**

## 2. de Faveri–Dunn–Hoffstein, *Non-orthogonality of the cubic and quartic large sieves via Rankin–Selberg*
arXiv:2607.07911 (v2, 13 Jul 2026). https://arxiv.org/abs/2607.07911

- **Thm 1.1 [V]** (unconditional): Ξ_3(A,B) ≫ (AB)^{−ε}(A + B + (AB)^{2/3}). This removes GRH from DR24 Thm 1.4.
- **Thm 1.2 [V]**: Ξ_4(A,B) ≫ (AB)^{−ε}(A + B + A^{3/4}B^{1/2} + A^{1/2}B^{3/4}). **Conj 1.3 [V]**: this is sharp.
- **General-n conjecture (1.8) [V]:** Ξ_n(A,B) = (AB)^{o(1)}(A + B + A^{1−1/n}B^{2/n} + A^{2/n}B^{1−1/n}). The paper
  also says: "Any power-saving improvement on the (AB)^{2/3+ε} term in (1.5) for any n ≥ 4 would be a substantial
  advance."
- **Extremal vector [V]** (proof of Thm 1.1, §2.1): β_b = g̃_3(b)·1_{b≡1(3)}·H(N(b)/B), the same as DR. The pole term
  P_B(a) ≍ g̃(a) B^{5/6} N(a)^{−1/6} comes from DdFDS24 Lemma 6.7. The second moment of ψ^{(3)}(a,s) is
  Lindelöf-on-average (DdFDS24 Lemma 6.6). The authors say "The main obstruction to perfect orthogonality comes from
  the bias exhibited by Gauss sums."
- **Relevance to ATS [I].**
  - Take n = 6 in (1.8) with {A,B} = {H,Q} and H ≤ Q. The conjectured untwisted norm is ≍ H + Q + Q^{5/6}H^{1/3},
    and Q^{5/6}H^{1/3} ≥ Q exactly when H ≥ Q^{1/2}. So in the ATS range the generic sextic sieve is conjecturally
    **not** sharp: the loss is (H²/Q)^{1/6}, at most Q^{1/6} at H = Q. Compare HB/BGL's (HQ)^{2/3}/Q, which is
    Q^{1/3} at H = Q.
  - So ATS can only hold through the twist, killing the pole-driven extremal vectors. This is consistent with the
    paper's mechanism, but nothing in the paper treats twisted vectors.
  - Note also that ATS, as stated, sums over all h (not only squarefree h). The h = k^6 terms contribute
    ≈ H^{1/6}|Σ_m b_m|². This term is harmless only because the angular twist forces cancellation in Σ_m b_m
    (Heath-Brown already warns that a sharp bound is impossible without squarefree restrictions; quoted in DR §1.4 [V]).

## 3. Blomer–Goldmakher–Louvel, *L-functions with n-th order twists*
IMRN 2014, no. 7, 1925–1955. arXiv:1112.1650. https://arxiv.org/abs/1112.1650

- **Thm 1.3 [V]**: for any number field K ⊇ μ_n (n ≥ 3) and arbitrary λ_b, the sum
  Σ*_{Na≤M} |Σ*_{Nb≤N} λ_b χ_a(b)|² is ≪ (MN)^ε(M + N + (MN)^{2/3}) Σ*|λ_b|². Here * means squarefree ideals coprime
  to S, and χ_a are the Fisher–Friedberg n-th order Hecke characters. They remark: "As in [HB2], the term (MN)^{2/3}
  is not optimal and can most likely be removed" (refuted by DR and DFDH for n = 3).
- **This is the source of the sextic bound in arXiv:2201.01885 [V].** Gao–Zhao, *Moments and one level density of
  sextic Hecke L-functions of Q(ω)*, Lemma 2.9, cite "[2, Theorem 1.3]" = BGL14, for squarefree E-primary m and n with
  the sextic symbol (n/m)_6. **Answers item (c).**
- The proof (§3) is Heath-Brown's recursive argument [V]: inequalities (E_α)/(F_α), a fixed point at α = 4/3, and
  Lemma 3.4 for the non-squarefree reduction. The 2/3 is produced by the recursion, not by any input you could
  replace. There is no special treatment of coefficient vectors and no infinity-type refinement [V by grep: no
  occurrence].
- Relevance: gives the HB-type c = 1/3 loss only. **Does not help ATS beyond the baseline.**

## 4. Heath-Brown, *Kummer's conjecture for cubic Gauss sums*
Israel J. Math. 120 (2000), 97–124. (Not freely accessible. Statement taken from the verbatim quotation in DR §1.4 and
§9, Thm 9.1.)

- **Thm 2 [V, via DR]**: the cubic large sieve with (A + B + (AB)^{2/3})(AB)^ε for squarefree variables. HB writes:
  "It seems possible that the term (AB)^{2/3} could be removed … It should be noted that if the variables are not
  restricted to be squarefree, a result as sharp as Theorem 2 would be impossible." The duality and monotonicity
  lemmas (HB Lemmas 4–5) are reproduced in DR Lemma 9.1 [V].
- **Item (b), unbalanced ranges:** I found **no** published improvement of HB's cubic sieve for general vectors when
  M^{1/2} ≤ N ≤ M.
  - DR Thm 1.4 and DFDH Thm 1.1 show that none is possible: the lower bound (AB)^{2/3} holds on all of [√B, B²]
    (up to ε-losses in the middle range).
  - The only "improvements" are for restricted vector classes: DR Props 9.1–9.2 (power saving δ for the class C_η,
    under GRH, with a correction term), and DdFDS24 Prop 4.1 (§6 below), which is a repackaging, not an improvement.

## 5. David–Dunn–Hamieh–Lin, *Quartic Gauss sums over primes and metaplectic theta functions*
arXiv:2306.11875 (v5, Jan 2026). https://arxiv.org/abs/2306.11875

- **Prop 4.3 [V]**: the quartic Kubota series ψ^{(4)}_{i1}(s,ν,ℓ), which carry the angular twist (c/|c|)^ℓ, continue
  meromorphically. "If ℓ ≠ 0 they are entire, if ℓ = 0 they have at most two possible (simple) poles at s = 3/4 and
  s = 5/4." They satisfy the functional equation (4.28), from Diaconu 2004 (2.48), with convexity bound (4.29), which
  is uniform in ℓ through (|s|² + ℓ² + 1)^{…}. **This answers item (a) for n = 4.**
- **Conj 1.2 [V]**: Σ_{N(c)≤X} g̃_4(c)(c/|c|)^ℓ Λ(c) is b_β X^{3/4} + O(X^{1/2+ε}) if ℓ = 0, and O_ℓ(X^{1/2+ε}) if
  ℓ ≠ 0. So the bias, and with it the large-sieve obstruction, is expected to vanish for nontrivial infinity type.
- **Structural trick [V]** (§1, (1.6)–(1.7)): by twisted multiplicativity, g̃_4(ab) = g̃_4(a)g̃_4(b)(a/b)_2. This turns
  quartic-Gauss bilinear forms into **quadratic** large sieve sums, which are sharp (Onodera; Goldmakher–Louvel).
- [I] For n = 6 the analogous identity has (a/b)_6^2 = (a/b)_3, so sextic-Gauss bilinear forms reduce to the
  **cubic** sieve, which is not sharp. The quartic trick does not transfer to the sextic case. But it suggests
  writing (·)_6 = (·)_2·(·)_3^{−1}, splitting off the sharp quadratic part, and isolating the cubic obstruction.

## 6. David–de Faveri–Dunn–Stucky, *Non-vanishing for cubic Hecke L-functions*
arXiv:2410.03048. https://arxiv.org/abs/2410.03048

- **Prop 4.1 [V]** "Cubic large sieve with decaying coefficients": Σ_{q∈F_3(Q1,Q2)} |Σ_{Nb≤R} a(b)χ_q(b)N(b)^{−1/2}|²
  ≪ (RQ1Q2)^ε (Q1Q2 + Q_*R + Q_*^{1/3}(Q1Q2R)^{2/3}). This is obtained from HB by Cauchy–Schwarz; it is not sharp.
- Lemma 6.6 (Lindelöf on average for the second moment of the Kubota series ψ^{(3)}(a,s)) and Lemma 6.7 (the polar
  term) are the inputs DFDH use [V, via DFDH §2].
- The abstract states that "the (unitary) family considered does not satisfy a perfectly orthogonal large sieve
  bound" [V].
- Relevance: these are technical inputs, such as the Lindelöf-on-average bound for twisted Kubota series. **They do
  not give ATS.**

## 7. Diaconu–Ion–Paşol–Popa, *On the second moment and non-vanishing of central values of Hecke L-functions of r-th order characters*
arXiv:2607.27131 (Jul 2026). https://arxiv.org/abs/2607.27131

- **Thm 4.1.4 [V]**, quoting Brubaker–Bump, *On Kubota's Dirichlet series*, J. reine angew. Math. 598 (2006),
  Thm 1: the general-r Kubota series D̂_S(s,a,ψ) continues meromorphically, "analytic except possibly at s = ½ ± 1/r".
  The functional equation mixes ψξ over characters ξ of order dividing r, with ψ unramified outside S.
- [I] In BB06 the twist Ψ lives on F_S^×/(F_S^×)^r. At a complex place, C^× = (C^×)^r, so **angular characters
  (c/|c|)^ℓ are not covered**. The gamma factor G(s) is fixed and independent of ℓ. I found no general-n analogue of
  DR Prop 5.2 or DDHL Prop 4.3 ("ℓ ≠ 0 ⇒ entire") in BB06-style literature. For n = 6 this is unverified.
- §1.1.3 [V]: the paper states that the current bottleneck for higher moments is HB's cubic large sieve and cites
  DFDH. "Explicitly evaluating the residue … at the simple pole s = ½ + 1/r remains an outstanding open problem for
  r ≥ 4" [V].

## 8. Other items checked

- **Gao–Zhao, arXiv:2201.01885 [V]**: Lemma 2.6 (sextic Kubota series h(r,s;χ), possible pole at s = 7/6, trivial
  infinity type) is proved "from the proof of the Lemma on [27, p. 200]". [27] = Patterson, *The distribution of
  general Gauss sums and similar arithmetic functions at prime arguments*, Proc. LMS (1987). Patterson 1987 and
  Heath-Brown–Patterson (J. reine angew. Math. 310, 1979) prove equidistribution of Gauss-sum angles. [I] That
  requires sums of g̃(c) twisted by angular Grössencharacters, so their Kubota-series continuation must cover
  nontrivial infinity type. I could not access these texts, so the n-th order pole structure for ℓ ≠ 0 is
  **[A/I], not verified**.
- **Baier–Young, *Mean values with cubic characters*, JNT 130 (2010), arXiv:0804.2233 [A]**: a cubic large sieve
  for cubic **Dirichlet** characters over Z (Thm 1.4, with terms such as Q^{11/9} + Q^{2/3}M). This is a different
  family and not applicable here.
- **Gao–Zhao, *Large sieve inequalities for quartic character sums*, Q. J. Math. 63 (2012) [A]**: weaker than BGL
  (per BGL [V]).
- **Gao–Zhao, *Mean values of some Hecke characters*, arXiv:1708.08164 [A]**: smoothed double character sums. Not a
  large sieve.
- **Young, *An improved spectral large sieve inequality for SL_3(Z)*, arXiv:2102.02796 [A]**: notes a structural
  analogy with HB's cubic sieve. DR suggest the GL3 sieve may also be optimal [V, DR §1.4].
- **Livné–Patterson (Invent. 2001), Eckhardt–Patterson (PLMS 1992), Kazhdan–Patterson (1984), Hoffstein,
  Friedberg–Hoffstein–Lieman (Math. Ann. 2003), Wright**: I found no large-sieve statement with
  Grössencharacter-twisted coefficients in any accessible source. I did not read the full texts (paywalled or not
  found).
- **Item (e), "large sieve for Gauss sums twisted by Grössencharacters":** I found **no such result** in any of the
  named authors' work that I could access. The closest are DR Props 9.1–9.2 (cubic, GRH, power saving δ) and the
  entire-ness of twisted Kubota series (DR Prop 5.2 for n = 3; DDHL Prop 4.3 for n = 4).

## Assessment for ATS (Q^{1/2} ≤ H ≤ Q)

1. **Nothing proves ATS.** No sharp or improved n-th order sieve with infinity-type-twisted coefficients exists in
   the literature I searched (to Oct 2026).
2. **Supporting evidence** [V + I]:
   - The only known obstruction comes from the Kubota pole (DR, DFDH).
   - The twisted Kubota series are entire for n = 3 (DR Prop 5.2) and n = 4 (DDHL Prop 4.3), so the explicit
     extremal construction collapses under a twist.
   - The DR correction term is a multiple of Σ_b β_b N(b)^{−1/6}, which vanishes for twisted β.
   These match the numerics in FINDINGS (twist ⇒ random-level ratios).
3. **Obstacles:**
   - (i) Even for cubic vectors that are orthogonal in the strong C_η sense, DR get only (AB)^{2/3−δ} under GRH, and
     explicitly say more input is needed.
   - (ii) DFDH (1.8) for n = 6 predicts an untwisted norm of Q^{5/6}H^{1/3} > Q throughout the ATS range. A proof
     therefore has to show that the twisted vectors avoid *every* near-extremal direction, not only the explicit
     Gauss vector. The extremal directions are not classified anywhere.
   - (iii) For n = 6 the relevant pole structure (theta coefficients, twisted series for ℓ ≠ 0) is not explicitly in
     the accessible literature. Patterson 1987 is the likely source; it is unverified.
4. **Most usable route suggested by the literature [I]:**
   - Adapt DR §9 (Voronoi with the δ_{ℓ=0} main term, then Poisson in the outer variable, as in their D/C/T
     decomposition) to the sextic symbol with an angular twist. With ℓ ≠ 0, the main terms D_1, C and T should all
     vanish, and the problem reduces to the "off-diagonal" dual sieve over |h| ≤ B²/A. That dual sieve is exactly the
     reflected sieve in FINDINGS (dual length Y = Q²/K).
   - Combine with the DDHL twisted-multiplicativity idea, (·)_6 = (·)_2(·)_3^{−1}, to peel off a sharp quadratic
     part.
   - Expected outcome: a δ-power saving, not full sharpness, unless a new idea controls the off-diagonal dual sum.
