# Referee report on solverC.md

S = solverC.md line, P = paper.md line. Numbers were rechecked in Python with (l0, x0, y0) from P625.
- M = 1 - l0 = 0.83288
- 1 + 3l - 2M = -0.1644
- E = x/2 + (y-x)/12 = 0.193797 = b0 - k
- P_a^{1/6} = Z^{0.014422}
- P_a^2/Y = Z^{-0.2866}
- X^{1/2} = Z^{0.186586}
- 0.85 - k = 0.168911

All of the solver's arithmetic (S9, S59, S96) is correct.

## Verdict table

| Claim | Verdict |
|---|---|
| Lemma 1 (S17, S41-48) | **Correct, with fixable gaps.** Most of it is trivial, and the argument for the actual ξ is simpler than written. |
| Lemma 2 (S19-27, S50-59): the asymptotic, hence a **lower** bound | **Unproven.** The main term is plausible; the error control needed for a lower bound is missing. |
| Lemma 3 (S29, S61-68): same functional, same Mellin point | **Correct, with a fixable gap.** The m-side Mellin point checks out; the column side is matched by scaling only, as the solver admits in U2. |
| Prop. O (S31-35, S70-83) as stated ("every CS-type bound") | **Unproven / overstated.** It holds at most for radial row weights, S-class × radial column weights, and per-(σ,v) splits. |
| Prop. O restricted to the paper's own split (B0.22 is sharp for these coefficients) | **Correct, with a fixable gap.** This depends on closing the Lemma 2 error terms. |
| Refutation of (iii) (S94) | **Correct** as a non-sequitur objection. |
| Heuristic \|I\| << Z^{0.169} (S96) | **Unproven.** The numbers are correct, but the cited derivation does not apply at β = 0.85 < b0. |
| Dual-route remarks (S106-110) | **Unverified assertions.** None of them is load-bearing. |

## 1. Lemma 1 (|Col(e0)| = |Col'(1)|)

**Tame p | b\* (S42): correct.**
- 1 + pO is pro-p and p ∤ 6, so a character of order dividing 6 factors through (O/p)^×.
- (O/p)^× is cyclic with 6 | q_p - 1, so the character is χ_p(·)^j.
- B1.5 (P937-940) then gives χ_p(n) = R(p,n)χ_n(p).
- R depends on n mod 4 (P923-928), so it is constant on σ mod 36.
- Unstated hypotheses: the primes of b\* must be primary generators, so that B1.5 applies, and must lie in S, so that e0 is one of the e that actually occur in the exceptional sum (P396: "e supported on C and the fixed primes"). Both appear to hold because b\* is "the product of fixed excluded prime generators" (P346), but the solver should say so.

**λ above 3 (S44): correct but superfluous.**
- log: 1 + 3O_λ ≅ 3O_λ ≅ Z_3^2, and (1 + 3O)^3 = 1 + 9O. So every order-3 character of primary n has conductor dividing 9.
- 9 | 36, so the character is already constant on mod-36 classes. No symbol is needed, and the claim that χ_n(ω) and χ_n(λ) "span" is unverified and irrelevant.
- Caveat: the step "primary n ∈ 1 + 3O" depends on the paper's normalization of "primary". The example in P505, n ≡ -2-3ω, is ≡ 1 mod 3, which is consistent, but neither the solver nor the cited section fixes this.

**2 (S46): correct.**
- Order-3 characters factor through F_4^× (conductor 2), so they are constant mod 36 with nothing to do. Since ξ is "cubic at norm four" (P346), this is the only case the paper actually uses.
- Quadratic characters: O_2^×/squares has order 8. The squares in 1 + 4O are 1 + 4·F_2 + 8O, because the Artin–Schreier image is the trace-zero F_2. So the conductor-4 subgroup has index 2, as claimed.
- (2/n)_2 = (2,n)_{K_2} = (2, Nn)_{Q_2} by the product formula; the λ-symbol is tame and trivial. Its conductor is 8 because N(1+4ω) = 13 ≡ 5 mod 8. Checked.

**Role of ψ and b\*.** ψ(n) = χ_n(b\*)ξ̄(n) is itself a character of n mod b\*·(4-data), by B1.5. So e0 = (b\*)^5·e_ξ, reduced mod sixth powers, which is legitimate.

**Exactness.** The identity |Col(e0)| = |Col'(1)| is exact class by class for **any** column weight f. Two qualifications follow.
- (a) const_σ depends on σ. For vectors that mix classes, |Col(e0)| ≠ |Col'(1)| in general; see Prop. O, item (d) below.
- (b) The true r-sum carries the coprimality condition (r, n1n2) = 1 (via χ_n(r^6) = 1_{(r,n)=1}). This puts a density factor φ(n)/q_n into **both** Col and the signal (χ_a(r^6) = 1_{(a,r)=1}, S64). The definition of Col at S26 omits it, and S55 calls it "1 + O(q_r-divisor terms)". That is wrong as stated: φ(n)/q_n is not 1 + o(1). It is harmless because the same factor appears on both sides and factorizes when (n1, n2) = 1.

## 2. Lemma 2: is the exceptional mass bounded BELOW? (focus item 1)

**What is right.**
1. **Hermitian structure.** With e and -e both in E_S (units are included), the C = d = 1 summand pairs to (χ_{n1}(-1) + χ_{n2}(-1))R(n1,n2)χ_{n1}(e)χ̄_{n2}(e). Within one σ this equals 2, using R(σ,σ) = χ_σ(-1) from B1.4 (P932) and bimultiplicativity. So the per-class C = 1 form is PSD, and S53 is correct.
2. **The r-sum (S54).** Ω' is radial and |r|^6 = (r r̄)^3 is a polynomial, so the lattice Poisson sum gives the main term c·|α|^{-1/3}κ(Ω') plus an O(1) r = 0 correction. The main term factorizes as (q_{n1}q_{n2})^{1/6}, and κ(Ω') = π∫Ω'(ρ)ρ^{-1/6}dρ > 0. Correct.
   - The constant is off. In R^2 with e^{-2πi}, the Fourier transform of |w|^{-5/3} is π^{2/3}Γ(1/6)/Γ(5/6)·|x|^{-1/3}, not π^{-1/3}… (S24). This is presumably a measure convention, and the sign is unaffected.
3. **Size of the competing terms.** The diagonal is Q/Y, the non-exceptional off-diagonal is (Q/Y)Z^{-0.287+ε}, the column-sum errors are of relative size Y^{-1/2}, and the r = 0 and Möbius-on-r errors carry relative factor τ(n1n2)P_a^{-1/6}. All of these are power-smaller than (Q/Y)P_a^{1/6} = (Q/Y)Z^{0.0144}. Fine.

**The gap: the C > 1 and d > 1 terms.**
- The solver bounds these with the paper's **upper** bound, "dyadic weight C_0^{-1/6+ε}D_0^{-1+ε}" (P396), reinterpreted as per-modulus decay q_C^{-7/6}, q_d^{-2} (S56).
- That bound carries Z^ε: divisor factors and the root count "≪ Z^ε(RP_a/C_0)^{1/6}" are absorbed into Z^ε (P387, P396).
- For a lower bound the C > 1 error must be smaller than the C = 1 main term (Q/Y)P_a^{1/6}·w_{e0}|Col'(1)|^2 ≍ (Q/Y)P_a^{1/6}. The available bound is (Q/Y)P_a^{1/6}·Z^ε·Σ_{C_0 > s}C_0^{-1/6}. This is not smaller than the main term for fixed S as Z → ∞.
- An aggregate dyadic bound does not yield per-modulus bounds. The C > 1 local values L_C (table at P1032-1037) are signed, so the C > 1 terms are not positive and could cancel the main term.
- A rough local count suggests the relative per-prime contribution could be as large as q_p^{-1} rather than q_p^{-7/6}. With L_p ~ q_p when p | k, the Euler product "1 + O(η_S)" may then not converge absolutely. I did not settle this, but the solver did not either.
- **Fix required:** an explicit main-term evaluation of the exceptional C > 1 terms, giving per-prime local factors with an ε-free decay q_p^{-1-δ}, and the same for the d-projector terms.

The solver's U1 (S87) partly acknowledges this, but the headline at S5 and S19 states the result as a theorem. Until this is done, Lemma 2 gives only the upper bound the paper already has. The lower bound, which is the whole content of Prop. O, is **unproven**.

**Restriction on v.** The lower bound needs |Col'_{σ,v}(1)| ≫ 1. That holds only for bounded v, since the error terms carry (1+|v|)^J (P375) while Col' decays in v. Harmless after integration against Ŵ_0, but it is not stated.

## 3. Lemma 3: "same Mellin point" (focus item 2)

**m side: verified.** κ(Ω) = π∫Ω(ρ)ρ^{-1/6}dρ, and the paper's 𝓜(1/6) ∝ ∫W_0(r)r^{-1/6}dr (P517-519). These are the same point. In the v-separated form κ(Ω_v) = π∫Ω(ρ)ρ^{-1/6-iv}dρ, and |κ(Ω_v)| ≤ κ(Ω') for Ω' ≥ |Ω| follows from the positivity of the kernel. Correct.

**Column side: matched by scaling only.**
- Both the Gram and the signal get t^{1/6} from the same r^6 lattice count against a radial transform, so the column weight W_1(t)t^{-5/6+iv} is the same on both sides by construction.
- The solver did not reconcile this with the paper's own normalization: c_S contains Ŵ_1(1) (P512), and B0.28 has Y^{w-1}, with the power of q_a hidden in the w-parametrization.
- The claim at S29 that the residue is "a constant times κ·Col' and nothing else" leaves out the φ(a)/q_a density (§1(b)) and the T-side factors (H_η, 1/L). The comparison is A-factor against A-factor only, so this is acceptable, but it should be stated.

**Verdict:** correct with a fixable gap. A one-page residue computation would close it.

## 4. Prop. O: does it cover ALL CS splits? (focus item 4)

The proof covers exactly the following: row-CS on the v-separated, single-σ vector A_{σ,v}; radial smooth row weights Ω'(q_m/Q); and column weights of the form (S-class function) × (radial profile). The phrase "every CS-type bound" (S5, S31-33) is not established. Specific holes:

**(a) Non-radial or arithmetic row weights.** κ(Ω') is defined only for radial Ω' (S24). Weights depending on arg m, on m mod a non-fixed modulus, on sixth-power-freeness of m, or on 1_{m ∉ supp of large ρ(m)} change the Poisson structure. They are not addressed. "Any positive CS weights" (S31) is unproved.

**(b) General column weights f(a)**, such as sieve weights, μ^2, primes, or amplifiers. The identity |Col(e0)| = |Col'(1)| survives, but Lemma 2 does not. The paper explicitly says its Gram proof "does not give the same bound for arbitrary column coefficients" (P396), and a lower bound needs the signed non-exceptional terms controlled.

**(c) CS without v-separation.** U3 (S89) argues that cancellation in v only strengthens the obstruction. That addresses the signal side only. For a CS applied to the un-separated vector ∫Ŵ_0(iv)(q_m/Q)^{-iv}A_{m,σ,v}dv, the exceptional mass is a double v-integral with cross kernels κ(Ω'·x^{-i(v1-v2)}). It is not ∫‖A_v‖ and not |∫Sig_v|^2. The triangle-inequality argument at S75 bounds only the per-v split from below.

**(d) Recombined classes (S83).** "Lemma 1 applies" is not enough.
- const_σ varies with σ, and for mixed classes the exceptional matrix is (χ_{σ1}(-1) + χ_{σ2}(-1))R(σ1,σ2)·Σ_e θ_{σ1}(e)θ̄_{σ2}(e).
- This is PSD and diagonal only if the symbols χ_n(e), e ∈ E_S, separate the mod-4 square classes on which R depends. The candidates are (-1/n)_2 and (λ/n)_2, and that they separate those classes was not checked.
- Otherwise the off-diagonal sign pattern of R could make the form indefinite. This is a potential escape route the solver closed too quickly.

**(e) Partial summation or Abel summation in q_m before CS, Hölder with higher moments, and the amplified or a-twisted split.** These are not covered. The solver's own §4(d) (S100-110) identifies the dual CS over columns a as **not** covered. That contradicts "every" at S5 and S33.

**(f) The u_0 twist (S82)** is a heuristic size count, not a proof.

**(g) The hypothesis "T controlled only through (B0.26)"** does all the work. Prop. O is then equivalent to the statement "‖A‖_{Ω'} ≥ c(Q/Y)^{1/2}P_a^{1/12}|Sig_A|/κ under the listed reweightings", i.e. a sharpness statement for B0.22. That is a useful and modest result: it shows the (y-x)/12 term in B0.27 cannot be removed by improving B0.22 alone. It should be stated that way, not as an obstruction to every CS split.

**Status.** The modest statement is correct with a fixable gap, conditional on §2. The universal statement is unproven.

## 5. Refutation of (iii) (S94)

g_{χ_n}(n, -mr^6) = χ̄_n(r^6)g(n,-m) = g(n,-m) for (r,n) = 1. This holds for every column, so invariance under m ↦ mr^6 cannot single out A^{exc}.

Small caveat: the invariance is only of the coefficients. The row weight changes by q_r^6, and columns sharing primes with r break it. That strengthens the solver's point. **Correct.**

## 6. Heuristic |I| << Z^{0.169} (focus item 5)

**Arithmetic.** 0.85 - k = 0.168911 < x/2 = 0.186586. Correct.

**Derivation: invalid.**
- B0.51 is proved under the contradiction hypothesis β\* > b0 (P691; μ = β\* - b0 > 0). The numerator estimates give F^{(4)}_κ - β\* ≤ -4(β\* - b0)/5 (P686).
- At β = 0.85 < b0 = 0.87489 the right side is **+0.0199**. The non-principal numerator error is then of size Z^{b0-k+O(·)}, not Z^{β-k}. So "the paper's own B0.51 machinery" yields only |I| ≪ Z^{b0-k} = Z^{0.1938}, which is the CS bound itself.
- B0.51 bounds J - f, not I. Bounding the u ≠ 1 rows of B0.28 under a zero hypothesis needs a separate argument. The u-sum converges only for Re z ≥ 0.34 (P498-500). A naive contour placement (s = 0.85, z = 0.34) gives X^{1/2-z}Z^{s+z-1} = Z^{0.25} before any w-savings. That is not obviously below Z^{0.169}.

**Verdict:** the conclusion "the loss is purely a method loss" (S96) is a plausible GRH-type expectation, not a consequence of the paper. **Unproven.**

## 7. Minor items

- **S106-110.** The bound "operator norm ≥ QP_a^{1/6}" and the exponents 0.833 and 0.2010 are asserted without computation. A direct sixth-power test vector for the character large sieve gives the weaker Δ ≥ YQ^{1/6}, so the stronger claim needs the Gauss-sum/character identification spelled out.
- **S98.** The Patterson-residue identification A^{exc}_m ≈ Y^{-1/3}ρ_ψ(m) is heuristic. It is not used in any proof.
- **S56.** "C_0^{-1/6}D_0^{-1} after summing over each dyad ⇒ q_C^{-7/6} per modulus" turns an averaged upper bound into a pointwise one, which is not valid in general (see §2).

## Bottom line

The solver has correctly identified why the P_a^{1/6} term appears: the r^6 exceptional frequencies and the principal-signal residue are both Riesz-1/6 lattice counts against the same column functional. Lemma 1 holds, and the refutation of (iii) holds. Proposition O, though, rests on a lower bound for the Gram exceptional mass that the solver has not proved: the C > 1 and d > 1 terms are controlled only by Z^ε-lossy aggregate upper bounds. Even granting that lower bound, Prop. O covers only radial-weight, per-(σ,v) row CS with S-class × radial column weights, not "every" CS split. The Z^{0.169} heuristic misapplies B0.51 below b0.
