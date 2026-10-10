# Solver D: the a-twisted reflected energy E_a

**Verdict: failure, with a precise obstruction.** The target estimate

  E_a ≪ Q(Q+Y)Z^ε

is not proved. The best exponent I can prove by the a-side (column-side) Cauchy–Schwarz route is **c = 1/3**, which is worse than the paper's c = 1/6. So for this geometry the method gives **no improvement**: the record θ is unchanged.

What *is* established is listed below. Each item is either proved or carefully reduced, with uncertain points flagged in §3.

1. **Cauchy–Schwarz over a is not intrinsically lossy.** Solver C's Prop. O applies to Cauchy–Schwarz over m. Over a there is no analogue: the bound is attained on the principal (u = 1) piece itself.
2. **The twist cannot be absorbed into the theta source.** It cannot be absorbed as a B2 label, a shift, or via reciprocity. The (B2.42) label route with V = y goes in the wrong direction.
3. **Reflection gives T a definite phase.** On squarefree good rows R, the reflected T carries the phase ᾱ(R)²γ₁(R) (a sextic Gauss sum times a nontrivial angular Grössencharacter) times a short quadratic sum of length D = Q²/Z.
4. **Poisson in a cancels that Gauss phase exactly.** The sixth-power ("exceptional") frequencies then carry no Patterson residue, because of ᾱ², and they are harmless. So the P_a^{1/6} mechanism of (B0.22) is **absent** from E_a.
5. **What remains is a sharp large sieve that nobody has.** It is a sextic large sieve beyond Heath-Brown, in the range K = Q²/Y > Q, for the specific vector b_m = ψ(m)𝓡(m). Equivalently, it is a non-correlation bound between T and the sextic Patterson theta coefficients.
6. **Poisson in m (the "h-form") is circular.** It lands on the reciprocal-Hecke (Möbius) side of B0.4. Controlling it needs zero-density input, the very thing being proved.

Script: `solverD/expo.py` (exponent bookkeeping at the paper's exact T1.2 geometry).

Geometry: x = .37317, y = .45971, ℓ = .16712, M = .83288, y − x = .08653, 1+3ℓ−2M = −.164.

---

## 1. Exact claims

Notation follows B0.3. Ω is the row profile, Q = q_{b*}XY, Y = Z^y, P_a = Y²/Q. T(m) = B_{m,σ}(Z) is the marked completed theta sum of (B0.6), with total compensator mark length ℓ.

Define

  F_a := Σ_m Ω(q_m/Q) ξ(m) (q_m/Q)^{-iv} χ̄_a(m) T(m),  E_a := Σ_{a∈σ} W_1(q_a/Y) |F_a|².

**Claim 1 (reduction; proved modulo U1).** Uniformly in |v| ≤ Z^ε,

  |Q^{-1/2} Σ_m Ω ξ (q_m/Q)^{-iv} A_{m,σ,v} T(m)| ≪ Q^{-1/2} Y^{-1/2} E_a^{1/2} Z^ε.

Consequently E_a ≪ Q² P_a^{θ} Z^ε (marks included) yields the probe exponent x/2 + θ(y−x)/2 + (1+3ℓ−2M)_+/8. In other words c = θ.

**Claim 2 (no intrinsic loss).** E_a ≥ c·QY|I_{σ,v}|² trivially. On the principal rows u = 1, the vector a ↦ F_a^{(u=1)} is proportional to the probe's own column vector, so Cauchy–Schwarz over a is an equality there.

The target bound would force the principal signal Z^{β−5/6+x/3+ℓ/6} ≪ X^{1/2}, i.e. β ≤ 5/6 + x/6 − ℓ/6 = .8677. That is consistent with the known β* ≤ 7/8, so it is not excluded by any lower bound.

**Claim 3 (reflected phase of T; proved for squarefree good rows R prime to the marks, modulo U2).**

  T(R) = κ · ν₁(R) · ᾱ(R)² · γ₁(R) · 𝓡(R),
  𝓡(R) = Σ_branches Σ_{μ≠0} d(μ)α(μ)ϑ(λ⁴μ) q_μ^{-1/2} χ_R(λ⁴μ)^{-3} V^♯(q_μ Z/q_R²).

Here:
* κ is a bounded scalar fixed by the sector.
* ν₁ is a fixed finite ray character with S-supported conductor, which absorbs γ₃(R) = G(R³) and the χ_p(−1), χ_p(λ) supplementary factors.
* 𝓡 has effective dual length D = Q²/Z = Z^{.666}, by (B2.37).

The same holds with active compensator marks, which change only the μ-side columns (B2.31).

**Claim 4 (Poisson in a; proved for coprime squarefree m₁, m₂, modulo U3).**

  E_a = Diag + Off.

* Diag ≪ Y Q^{7/6} Z^ε. Its exponent is 1.4314 < 2M = 1.6658.
* Off ≍ (Y/Q) Σ_{σ'} Σ_{0≠k, q_k ≪ K Z^ε} |B_{σ'}(k)|², with K = Q²/Y = Z^{1.206}.
* B(k) = Σ_m b_m χ_m(k)·(sep. smooth), and b_m := ξ(m)Ω(q_m/Q) T(m) γ_{-1}(m) = ψ(m)ᾱ(m)² 𝓡(m) for a fixed ray character ψ. **The Gauss phase is gone.**

On the exceptional frequencies k = e r⁶ (e S-supported), Σ_r |B(er⁶)|² ≪ K^{1/6} D² Z^ε. That contributes Z^{1.159} to E_a, well below Q².

**Claim 5 (what the method yields).**

* Target E_a ≪ Q²Z^ε ⇔ Off ≪ Q²Z^ε ⇔ (SLS_T): Σ_{q_k ≪ K} |Σ_{q_m≍Q} b_m χ_m(k)|² ≪ K·Q·Z^ε.
* (SLS_T) would give c = 0, with b0 = .865152 after re-optimizing geometry (gram.py), i.e. a gain of 9.7e-3.
* Available inputs give:

| Input | E_a exponent | c | Probe exponent |
|---|---|---|---|
| Heath-Brown-type sextic sieve, a-form (Y+Q+(YQ)^{2/3})·Σ\|T\|² | 1.6946 (excess +.0288 = P_a^{1/3}) | 1/3 | .2010 |
| Same sieve, k-form Y(KQ)^{2/3} | 1.819 | — | worse |
| Classical multiplicative sieve, Q(Q+Y²) | — | 1 | .2299 |
| Paper's own (B0.22), m-side | — | 1/6 | .1938 |

Net: **c = 1/6, no gain.**

## 2. Proofs

### 2.1 Claim 1

In (B0.21), A_{m,σ,v} = Y^{-3/2} Σ_a c_a W_1 (q_a/Y)^{-1+iv} g_{χ_a}(a,−m), with c_a = χ_a(b*)/(τξ(a)), |c_a| = 1.

**Squarefree good a.** Here g_{χ_a}(a,−m) = χ̄_a(−m) g_{χ_a}(a,1), where |g| = q_a^{1/2} and the value is 0 when (a,m) ≠ 1, consistent with the zero extension. χ_a(−1) is a fixed-class sign. Swapping sums:

  Σ_m(...)A_m T(m) = Y^{-3/2} Σ_a c_a χ_a(−1) g(a,1) W_1 (q_a/Y)^{-1+iv} F_a.

Cauchy–Schwarz over a, with Σ_a |g(a,1)|² W_1² ≪ Y², gives Y^{-3/2}·Y·E_a^{1/2}.

**Non-squarefree a = a₀s, s powerful.** By (B1.13), |𝒢(p^j,h)| ≤ Q^{(j−1)/2} with support v_p(m) = j−1, or a Ramanujan sum when 6 | j. For each fixed s these are characters in m times a bounded factor, restricted to s′ | m. The resulting energies are E_a-type sums over the row sublattice s′ | m, with an outside factor q_s^{-1/2+ε}. Summing over s converges (U1).

Multiplying by Q^{-1/2} gives Claim 1. Taking E_a ≪ Q²P_a^θ gives |I| ≪ Q^{1/2}Y^{-1/2}P_a^{θ/2} = q_{b*}^{1/2} X^{1/2} P_a^{θ/2}.

The compensator (B0.25) and the rescaled subsets of (B0.27) are linear in T with the same columns. Each rescaled piece has its own (X′,Y′), the same P_a, and energy exponent from (B0.26). So the paper's derivation of (B0.27) goes through with (y−x)/12 replaced by θ(y−x)/2.

### 2.2 Claim 2

Poisson in m (B0.4) writes the u = 1 rows of the probe as a sum over a of c_a g(a,1)·(Gauss sum of χ̄_a)·(an a-independent factor). Two cancellations make the a-dependence trivial:
* c_a·χ̄_a(b*)ξ(a) = 1/τ (the twist cancels by design).
* g_{χ_a}g_{χ̄_a} = χ_a(−1)q_a.

So F_a^{(u=1)} is proportional to conj(c_a g(a,1)) times a smooth function, and Cauchy–Schwarz over a is an equality on that piece.

Its size is X^{1/2}Σ_r U(r⁶), with Σ_r ≈ (q_c/X)^{1/6}. That reproduces the paper's X^{1/3} principal scale (B0.28/B0.35, 𝓜(1/6)). The β-consistency figure is printed by expo.py.

Contrast with Prop. O. There, the m-side norm ‖A‖ contains the exceptional mass, and the same Riesz functional drives the signal. On the a-side the column norm is exactly Y, with no exceptional, and all the difficulty sits in E_a. No signal-based lower bound exceeds Q².

### 2.3 Claim 3 (reflected phase)

Take a squarefree good row k = R with every prime active and j_p = 1. Then χ_A(R) = R(A,R)Π_p χ_p(A) by (B2.9), so j_p = 1. By (B2.11), B_{p,1}(x) = χ_p(x)^{-3}.

By (B2.13) and (B2.32), the R-dependent part of ζ is

  ᾱ(c)² Π_{p|R} ξ_{p,1} χ_p(R/p)⁴.

Here ξ_{p,1}χ_p(c/p)⁴ = χ_p(σ_p)^{-2}ω_{p,1}, with ω_{p,1} = τ⁻_{p,1}τ⁺_{p,3}χ_p(ε_p)^{-3}.

Expand using ε_p = −λ^{-5}(c/p)^{-2} and σ_p = λ²c/p:

  χ_p(σ_p)^{-2}χ_p(ε_p)^{-3} = χ_p(λ^{-4}(−λ^{-5})^{-3}) · χ_p(c/p)^{-2+6}.

The first factor is a fixed supplementary character of p. Also:
* τ⁻_{p,1} = χ_p(−1)γ₁(p).
* τ⁺_{p,3} = γ₃(p).

Collect the cross factors:
* Π_{p≠q|R} χ_p(q)⁴ = Π_{p<q}(χ_p(q)χ_q(p))⁴.
* By (B1.9), γ₁(R)γ₃(R) = Πγ₁(p)γ₃(p)·Π_{p<q}(χ_pχ_q)^{1+3}.

So the R-part equals ν(R)ᾱ(R)²γ₁(R)γ₃(R). Since γ₃(R) = G(R³) by (B1.4) and multiplicativity, it is a fixed-class function.

Marks: (B2.33) with j_p = 1, j_q = 0 gives exponent 3, and (q/p)₃³ = 1. So active marks add no row phase. Their factor q_q^{-1/2}χ_q(λ⁴μ)^{-2} lives on the μ side.

### 2.4 Claim 4 (Poisson in a)

Fix the ray class of a modulo 36 so that R(a,m) is a constant (B1.3). Write χ̄_a(m₁)χ_a(m₂) = R·χ̄_{m₁}(a)χ_{m₂}(a) by (B1.5). This is a character of a modulo m₁m₂ (times the fixed modulus). Apply Poisson over a ∈ 𝒪 at length Y.

**k = 0 (Diag).** The term is nonzero only if χ̄_{m₁}χ_{m₂} is principal, i.e. m₂ ∼ m₁ modulo sixth powers and S-units. Each m₁ has ≪ (Q/q_e)^{1/6} ≤ Q^{1/6} partners. Apply |γ₁γ₂| ≤ ½(|γ₁|²+|γ₂|²) and use (B0.26), whose exponent at the optimum is M because 1+3ℓ−2M < 0:

  Diag ≪ Y Q^{1/6} Σ|T|² ≪ Y Q^{7/6} Z^ε.

**k ≠ 0, coprime squarefree m_i.** By CRT the complete sum is

  ψ₁(m₂)ψ₂(m₁) g(ψ₁,k) g(ψ₂,k),  with ψ₁ = χ̄_{m₁}, ψ₂ = χ_{m₂},

and ψ₁(m₂)ψ₂(m₁) = R(m₁,m₂). Also:
* g(χ̄_{m₁},k) = χ_{m₁}(k) q_{m₁}^{1/2} γ_{-1}(m₁).
* g(χ_{m₂},k) = χ̄_{m₂}(k) q_{m₂}^{1/2} γ₁(m₂).

Insert Claim 3, using γ₁γ_{-1}(m) = χ_m(−1) and |γ₁|² = 1. With prefactor Y/(q_{m₁}q_{m₂}), each m_i carries b_{m_i}/√q_{m_i}. The coupling Ŵ(k√Y/(m₁m₂)) is separated by one common Mellin/Fourier inversion, as in (B2.22). R(m₁,m₂) is a bicharacter of ray classes, so after splitting into finitely many classes and applying Cauchy–Schwarz,

  Off ≪ (Y/Q) Σ_{σ'} Σ_{k≠0} |B_{σ'}(k)|².

The dual length is q_k ≪ Q²/Y·Z^ε = K. That is the claimed identity.

**Exceptional k = e r⁶.** Then B(er⁶) = Σ_μ ρ_μ Σ_{(m,r)=1} ψ(m)ᾱ(m)²χ_m(eμ³)w(q_m/Q). For each μ, the inner sum is a smooth sum of a Hecke character with **nontrivial infinity type** ᾱ² and conductor ≪ q_eq_μ ≤ q_eD < Q. So it is ≪ (q_eq_μ)^{1/2}Z^ε, with no main term.

This holds even when eμ³ is a sixth power, because the angular character remains. Equivalently, the Kubota series Σ_m ᾱ(m)²𝒢(m,n)q_m^{-s} has no Patterson pole: the intertwining factor is an L-function of ᾱ^{12} ≠ 1. With Σ|ρ_μ| q_μ^{1/2} ≪ D and (K/q_e)^{1/6} values of r, the contribution to E_a is

  (Y/Q) K^{1/6} D² Z^ε = Z^{1.159}.

Even a full-mass exceptional vector, b_m ≡ 1, would give Y K^{1/6} Q = Z^{1.494} < Q². So the exceptional frequencies are harmless.

### 2.5 Why the twist is not a label or a shift

**Label.** The B2 label enters as χ_A(f)⁴ = χ_A(f⁴). That is a column twist, equivalent to moving the row k → kf⁴. Our twist χ̄_a(m) is a function of the row with modulus a. Two rewrites do not help:
* By reciprocity it equals R·χ̄_m(a): modulus m, argument a.
* By multiplicativity in the modulus it equals χ_{a⁵}(m), which shifts the theta index A → Aa⁵. That destroys the theta coefficients. There is no Hecke operator for a⁵ on the sextic cover; only p⁶ is available.

Even if a label of size Y were usable, (B2.42) with V = y has column branch 2M + y − N = 1.126 > M. A label lengthens the reflected conductor q_c ∼ q_Rq_f, so it can only hurt. Converting a per-row energy back to the bilinear form would also cost a further Cauchy–Schwarz over m, i.e. a factor Q.

**Poisson in m (the h-form).** Here β_cγ₁(c) = μ(c)α(c)G(c) by (B1.10), so the dual coefficients are U(h) = Σ_c μ(c)G(c)χ̄_c(h)/q_c over c ∈ [Xq_h, Z]. Their ℓ²-mass needs cancellation in Möbius sums twisted by χ_h; the large sieve alone gives ≈ 1 instead of 1/X. That is the reciprocal-Hecke side of B0.4, and bounding it needs zero-density input. **Circular.**

## 3. Uncertain steps

* **U1.** Non-squarefree columns a. The (B1.13) support gives the claimed q_s^{-1/2} saving heuristically; I did not write out the row-sublattice energies.
* **U2.** Claim 3 is derived only for squarefree good rows. I did not fully track:
  * powerful and S-supported rows: by (B2.38) they cost O/2 and have frozen non-residual primes with j ≠ 1, whose (B2.33) phases are row scalars;
  * the h₀-branch sum;
  * the exact fixed supplementary characters.

  Only "fixed ray character × ᾱ²γ₁γ₃" matters, and that structure is forced by the computation.
* **U3.** Claim 4 is computed only for coprime squarefree m₁, m₂. Common parts C = (m₁,m₂) need the B1.17–B1.19 local table and projector. I expect the same Gauss-phase cancellation, since L_C has |L_C| ≤ q_C and no γ₁ phase, but I did not verify it.
* **U4.** The claim of no Patterson pole for the angular-twisted sextic Kubota series is standard Eisenstein-series theory (nontrivial infinity type, so the constant-term L-factor is entire). In 2.4 I only use the elementary version: smooth sums of a nonprincipal Hecke character.
* **U5.** The Heath-Brown-type sextic large sieve (M+N+(MN)^{2/3}) is quoted from FINDINGS (via arXiv 2201.01885). I did not check it. It is not needed for the negative conclusion: even granting it, c = 1/3.

## 4. Precise obstruction

After Claims 1–4, the target is **exactly equivalent**, up to Z^ε and the harmless Diag and exceptional terms, to

  (SLS_T)  Σ_{k≠0, q_k≪Q²/Y} |Σ_{q_m≍Q} ψ(m)ᾱ(m)²𝓡(m) χ_m(k)|² ≪ (Q²/Y)·Q·Z^ε.

Dually, in a-form, it is

  Σ_{a~Y} |Σ_m T(m)ξ(m)χ̄_a(m)|² ≪ Q·Σ_m|T(m)|²·Z^ε.

This is a sextic large sieve with no (MN)^{2/3} term, at the ratio Y < Q < Y². The best general bound loses P_a^{1/3} = Z^{.029}.

The extremal vectors for a sextic large sieve at these lengths are the Patterson theta coefficients ρ₆(m). The paper's A^{exc} is exactly Y^{-1/3}ρ₆. So (SLS_T) amounts to

  |Σ_m ρ₆(m)·γ₁(m)ᾱ(m)²𝓡(m)| ≪ (Patterson-norm)·Q^{1/2}·P_a^{-1/12}·...

That is a non-correlation statement between sextic theta coefficients, which are not explicitly known for n = 6, and the reflected cubic theta sum. It is the same input solver C isolated, now in a Gauss-phase-free form.

Neither of the two available tools reaches it:
* B2's reflection only controls |T|² pointwise in m. The phase γ₁(m)ᾱ²(m) of T is invisible to it.
* The a-Poisson that removes the phase produces a dual of length K = Q²/Y > Q, where column-Poisson decay (the B0.23/B0.24 mechanism, with Λ = Q/K < 1) gives nothing.

Treating μ in 𝓡 = Σ_μ ρ_μ χ_m(μ)³ trivially gives E_a exponent 1.8265; this is an excess of 0.16 that would have to come from cancellation over μ. Treating it by Cauchy–Schwarz over μ loses D.

**Bottom line.** No exponent c < 1/6 is obtained. The a-route is sound in principle (no intrinsic loss, no exceptional term), but it needs a new sharp sextic large sieve or a theta-correlation estimate.
