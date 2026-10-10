# Solver B: inverse moment beyond the Poisson range (lemma (*))

## Result: FAILED. No λ < 5/6 is proved.

This report gives (1) the claimed statement, which is none; (2) the reductions I did prove, with full arguments; (3) uncertain points; and (4) the precise obstructions.

Notation follows paper.md B0.1. M_k(D) = D^{-1/2} Σ μ(n)ν(n)χ_n(k)^{±1} W(q_n/D). Logs are in base H with D = H^x and 1 < x ≤ 1.2. For 6-free rows, logs are in base U with D = U^r.

---

## 1. Claimed statement

Nothing new is claimed. The best bound I can justify for 1 < r ≤ 3/2 is still (B0.9):
Σ_{q_u≍U}^{(6)} |M_u(U^r)|^2 ≪ U^{(1+5r)/6+ε}. In the language of (*), this is λ = 5/6.

---

## 2. Proved reductions (rigorous given the paper's UM, (B0.7), (B0.8) and (B0.15))

### 2.1 (*) for all rows is equivalent to its 6-free form, with the same λ

(⇐) Assume Σ_{u≍U}^{(6)} |M_u(D)|^2 ≪ U^{1-s} D^s for U ≤ H. Then every row k ≤ H can be written k = u a^6. By (B0.7) and the sup form (B0.8), M_{ua^6} is controlled by M_u at scales ≤ D, with a divisor loss. There are about (H/U)^{1/6} values of a for each dyadic U. Summing over dyadic U gives max_U (H/U)^{1/6} U^{1-s} D^s = H^{1-s} D^s, because 1-s > 1/6 when s < 5/6.

(⇒) Assume Σ_{k≤H'} |M_k|^2 ≪ H'^{1-λ} D^λ. Lifting as in (B0.9) with H' = U P^6 gives U^{1/6} H'^{5/6-λ} D^λ. This is optimal at H' = U, which gives U^{1-λ} D^λ.

**Consequence:** the lifting u → u a^6 creates no information. It only moves λ between the two forms. The 5/6 in (B0.9) is the lift of the bound "Σ_{k≤H} ≪ max(H, D)" taken at H = D. Indeed U^{1/6} H'^{-1/6} D at H' = D equals U^{1/6} D^{5/6}. At H = D every λ gives the same value D. So (*) is exactly the statement "beat the large-sieve value D when H < D".

### 2.2 Sixth-power rows are not the obstruction under the contradiction hypothesis

Take rows k = u_0 a^6 with q_{u_0} = H^α. There are H^{α + (1-α)/6} of them, and pointwise |M_k(D)|^2 ≪ D^{κ+ε} by (B0.15). Their total contribution is H^{α+(1-α)/6+κx}.

For α near 0, this is at most H^{1-λ+λx} whenever 1/6 + κx ≤ 1 - λ + λx. At κ = 3/4 and x ≤ 1.2, that holds for every λ ≥ 1/3. Pointwise bounds handle any family of rows with α ≤ α_0(λ, κ). As α → 1 the family becomes generic and the pointwise cost H^{1+κx} is hopeless.

So the difficulty is entirely in generic 6-free rows. The suggested strategy "treat u_0 a^6 separately with fixed-L input" is valid, but it removes nothing that matters.

### 2.3 Truncated mollifier (Lemma-1W style) identity: new and exact up to O(H^{-A})

Let M_X(s) = Σ μ(d) w(d/X) ψ_k(d) q_d^{-s} with w = 1 on [0,1]. For n ≍ D > 2X,

  μ(n) = Σ_{bm=n} f(b) μ(m),  where f = δ - 1 * (μ w_X).

f is supported on b > X. Each dyadic block of f is F_k(B) = -Σ_Δ M_k(Δ) S_k(B/Δ). Here S_k(Y) is a smooth plain sum whose conductor (including redundant radical, B3.44) is ≪ H. By the reflection (B3.43)–(B3.45), S_k(Y) = O(H^{-A}) for Y ≥ H^{1+ε} when ψ_k is nonprincipal. Hence

  M_k(D) = Σ_{X < B ≤ X H^{1+ε}} F_k(B) M_k(D/B) + O(H^{-A})   (Mellin-separated).

This is the only nontrivial content of any mollifier or Vaughan identity here: plain factors are shorter than the conductor. Without it the identity is tautological, since μ * 1 * μ_X = μ_X vanishes at n ≍ D.

---

## 3. Uncertain / unverified points

- 2.3 needs nonprincipal ψ_k for every row. Principal rows (k a unit times a sixth power, ν trivial) number H^{1/6}. They are covered by 2.2 only under the contradiction hypothesis applied to ζ_K.
- I did not re-verify UM, MI, MP0 or RE themselves. See remark 4.4.

---

## 4. Precise obstructions

### 4.1 Mollifier / Vaughan route (direction "Lemma 1W mollification")

From 2.3, M(D) ≈ Σ M(Δ) S(Y) M(A) with Δ ≤ X, X < ΔY ≤ XH, Y ≤ H, and A = D/(ΔY).

The framework has no mixed Möbius × plain moment:
- MI allows only prime slots.
- MP and MP0 allow only plain and prime factors.

So F_k(B) is controlled only pointwise, or by Cauchy with MP0 and pointwise·UM:
- Σ|M(Δ)S|^2 ≤ H Δ^{κ/2}, or
- H min(Y, H/Y)^κ.

**Self-dual point.** At Y = √H the plain factor is self-dual: reflection maps length √H to H/√H = √H. So the M·S product cannot be shortened, and the term (Δ = 1, Y = √H, A = D/√H) has coefficient length D > H after optimal reflection. This is the original problem with μ*1 in place of μ. The best available cost is H^{1+κ/4} or worse, compared with the target H^{1+λ(x-1)}, and x - 1 ≤ 0.2.

**Even with a hypothetical sharp tool.** Grant a tool T stronger than the framework: Σ_k |F_k(B)|^2 ≪ max(H, Δ·min(Y, H/Y)). Combine it with pointwise or UM on M(A), and choose X optimally. Numerically (scratchpad/vaughan_check.py, κ = 0.70 and 0.75, x = 1.05, 1.1, 1.2) the best exponent is exactly x, i.e. λ = 1, which is the plain large-sieve value D. That is worse than (B0.9)'s lifted value.

The binding triples are (Δ ≈ X, Y ≈ 1, A ≈ D/X) against (Δ ≈ X, Y ≈ √H) once a ≥ 0 is enforced.

### 4.2 Marked lifting (direction "MI marks + sixth-power averaging")

Lifting Σ_k |M_k(D) Q_k|^2 needs H ≥ D Z_I^2 by the MI condition r' + 2z' < 1 (and 2r' + 8z' ≤ 3 needs z' ≤ 1/4). The bound becomes

  Σ_u |M_u Q_u|^2 / min|Q|^2 ≪ U^{(1+5r)/6 + 5z/3 - 2qz}.

This gains only if q > 5/6. But q ≤ δ/2 < 3/8. The marks lengthen the effective Möbius by Z^2 faster than their amplitude pays.

The same failure occurs for the witness pair M(U^r) S(U^m): a mixed lifted moment would give at best the long envelope at length r+m, which the paper already uses at t.

### 4.3 Large values instead of the second moment (Halász–Montgomery)

Duality over columns, with off-diagonals being plain sums of conductor ≤ U^2 under (B0.15), gives #rows ≤ D/V^2 = U^{r(1-δ)}. That is valid when δr ≥ r/2 + κ(2-r)/2, which is borderline at r = 1.2.

It is still worse than (B0.9) by U^{(r-1)/6}. The Cauchy step over n destroys μ. Keeping μ requires bounding Σ μ(n) Σ_j c_j χ_n(u_j), which is (*) again.

### 4.4 Structural remark (important for any reviewer)

UM is asserted unconditionally (paper line 1774). Restrict (B0.8) to the rows k = a^6 with q_a ≤ H^{1/6} and ν trivial. Using (B0.7) with u = 1, there are about H^{1/6} such rows, and each satisfies |M_1(D)| ≤ d(a) sup_E |M_{a^6}(E)|. Hence

  |D^{-1/2} Σ_{(n,S)=1} μ(n) W(q_n/D)|^2 ≪ H^{5/6+ε} = D^{5(1+θ)/6+ε}

for every θ > 0. Therefore 1/ζ_K is holomorphic on Re s > 11/12, i.e. ζ(s) ≠ 0 for Re s > 11/12, unconditionally. The same argument with rows u_0 a^6 gives the same half-plane for every fixed sextic twist.

This explains b_0 = 11/12 - ℓ_0/4: the whole contradiction machinery improves on what UM alone already asserts. Two consequences:
- (a) UM, through RE in Appendix B2 and the canonical estimate, is itself a quasi-Riemann-strength unconditional claim. It is where verification effort should concentrate.
- (b) Every route to (*) inside the framework must either use UM-type Möbius information at H < D, which the paper never supplies (its UM, MI and canonical estimates are all of large-sieve type, giving max(rows, columns)), or use pointwise κ-input. Pointwise input costs D^κ per row and only pays on families of at most H^{α_0} rows (2.2).

### Bottom line

Lemma (*) with λ < 5/6 is equivalent (2.1) to beating the large-sieve value D for generic 6-free rows when H < D. Every tool in the paper is either:
- large-sieve-shaped (Poisson involution, canonical estimate, MI, MP0), so it is capped at max(H, D); or
- pointwise (cost D^κ per row).

Adding the missing mixed Möbius × plain moment in its sharpest large-sieve form still gives λ = 1 (4.1), because of the self-dual plain length √H. A proof needs a genuinely new input, such as a bound of the form Σ_{k≤H} |M_k(D)|^2 = o(D) for H < D. That would be an average Möbius cancellation statement over generic sextic families, below the large-sieve diagonal.
