# Solver C: can the P_a^{1/6} exceptional term in the probe bound (B0.27) be avoided?

**Verdict: no improvement is proved.** What I did establish:

* **Obstruction theorem (Prop. O).** Consider every bound that goes through a Cauchy–Schwarz (CS) step of the form |<A,T>| <= ||A||_{Omega'} ||T||_{Omega''}. This covers any positive row weights, any column reweighting by S-classes or radial profile, any admissible twist xi, and the v-separation used in B0.3. Every such bound loses a factor >> P_a^{1/12} relative to the principal signal. The reason is that the exceptional Gram mass and the principal signal are controlled by the same linear functional of the column weights. Rigour level: main terms are derived, but a few steps are only sketched (see section 3).
* **Two flaws in the earlier arguments.** Item (iii), "T is aligned with A^{exc}", does not follow from what was given. Item (ii) is correct, and I closed the 2-adic and 3-adic gaps it left open.
* **The one route left open needs a new bilinear input.** That input is a correlation bound between sextic theta coefficients (Patterson residues) and the cubic theta sum T, or equivalently a T-specific sextic large sieve. Nothing in the paper supplies it, and neither do the published large sieves.

Numerical identity used throughout. At the optimum M = x+y = 1-ℓ, which makes (1+3ℓ-2M)_+ = 0 (its value is -0.164). Then b0 = E + k = 5/6 + x/6 + c(y-x)/2 - ℓ/6 with c = 1/6. Python check: E + k = 0.87488664704 = b0. The exponent of P_a is y-x = 0.0865, so P_a^{1/6} = Z^{0.0144}.

---

## 1. Exact claims

Notation follows the paper. A_{m,σ,v} is the column vector from B0.3 restricted to the ray class σ (mod 36 O). c_a = χ_a(b*)ξ̄(a)/τ is the column twist, and ψ(n) := χ_n(b*)ξ̄(n). Poisson in m is done against a row weight Ω' >= 0 (any positive weight used in CS).

**Lemma 1 (every admissible twist is a symbol twist; strengthens (ii)).** Let λ be any character of order dividing 6 on primary ideals n coprime to S, with conductor supported on S. Then there is an e in E_S (S-supported elements with valuations 0..5, units included) such that n ↦ λ(n)χ_n(e) is constant on every class mod 36 O.

**Lemma 2 (the exceptional Gram term is a positive form, dominated by its C = d = 1 part).** Assume P_a >= Z^{δ}, with fixed δ > 0 and P_a^2/Y <= P_a^{1/6-δ}. Then for each σ and v:

Σ_m Ω'(q_m/Q)|A_{m,σ,v}|^2 = (Q/Y)·P_a^{1/6}·κ(Ω')·(1+O(η_S))·Σ_{e∈E_S} w_e |Col_{σ,v}(e)|^2 + O((Q/Y)(1+P_a^2/Y+P_a^{1/6-δ'})Z^ε).

Here:
* κ(Ω') = c_0 ∫_{C} Ω'(|x|^2)|x|^{-1/3} dμ(x) > 0 is the Riesz moment. It is the Fourier transform of |z|^{-5/3} in the plane: π^{-1/3}Γ(1/6)/Γ(5/6) > 0.
* w_e ≍ q_e^{-1/6} > 0.
* Col_{σ,v}(e) := Y^{-1}Σ_{n∈σ} W_1(q_n/Y)(q_n/Y)^{-1+iv}(q_n/Y)^{1/6}ψ(n)χ_n(e).
* η_S = Σ_{p∉S} q_p^{-7/6+ε}. This tends to 0 as the small-prime cutoff of S grows.

**Lemma 3 (the same functional drives the principal signal).** In the v-separated form of I_η, the u = 1, z = 1/6 residue for the piece (σ, v) equals a nonzero constant times κ(Ω_v) times Col_{σ,v}(1)' and nothing else. Here Ω_v(x) = Ω(x)x^{-iv}, and Col' is Col with ψ replaced by 1, because the twist is cancelled exactly by the m-Poisson transform. By Lemma 1 there is an e_0 with ψχ(e_0) ≡ const_σ, so |Col_{σ,v}(e_0)| = |Col'_{σ,v}(1)|.

**Proposition O (obstruction).** Allow all of the following to vary: the column weights f(a) (an S-class function times any smooth radial profile W_1), any admissible ξ (primitive mod b*, with ξ^6 principal), any extra column twist by an order-6 S-character or by χ̄_a(u_0) with u_0 ∉ S, and any positive CS weights Ω', Ω''. Assume P_a >= Z^{δ}. Then every CS-type bound on I^{comp}, with T controlled only through (B0.26), satisfies

bound / |principal signal coefficient| >> P_a^{1/12} × (the paper's normalized X^{1/2} scale).

In words, the effective probe exponent cannot drop below x/2 + (y-x)/12 + (1+3ℓ-2M)_+/8 by any modification of this kind.

**Non-claim.** I do **not** claim that |I^{comp}| really has size X^{1/2}P_a^{1/12}. Section 4(b) argues the opposite, heuristically.

## 2. Proofs

### Lemma 1
At a tame p | b* (p ∤ 6), any order-6 character of (O/p^k)^× factors through (O/p)^×, because 1+pO is pro-p. So it equals χ_p(n)^j, and by (B1.5) that is R(n,p)^j χ_n(p)^j. The factor R depends only on the class of n mod 4.

At λ (the prime above 3): primary n lie in 1+3O ≅ Z_3^2 (torsion-free, since ω ∉ 1+3O). Its order-3 characters form (Z/3)^2, and they are trivial on (1+3O)^3 = 1+9O, so they are constant mod 9. The cubic supplements χ_n(ω) and χ_n(λ) span this group. There are no order-2 characters here.

At 2 (inert, norm 4): order-3 characters factor through F_4^×, giving χ_n(2)^{±2}. For order 2, O_2^×/O_2^{×2} has order 8. The characters of conductor dividing 4 form a group of order 4, constant on mod-36 classes. The symbol χ_n(2)^3 = (2/n)_2 = (2, Nn)_{Q_2} has conductor 8: take u = 1+4ω, whose norm is ≡ 5 mod 8. It therefore supplies the missing coset.

So ξ(n) = (const on σ)·χ_n(e_ξ), and likewise for any admissible λ.

### Lemma 2
Start from the proof of (B0.22) at C = d = 1, where k = ε e r^6 and (B1.17) gives the summand c_{n1}c̄_{n2}χ_{n1}(k)χ̄_{n2}(-k)R(n1,n2).

1. **Phases.** Inside σ, R(n1,n2) = R(σ,σ). By (B1.4) and bimultiplicativity, R(n,n) = Π R(p,p) = χ_n(-1). So R(σ,σ)χ̄_σ(-1) = 1, and the k-term is |Σ_n c_n χ_n(e)·…|^2 up to the smooth coupling.
2. **r-sum.** Σ_{r≠0} Ω̂'(e r^6 √Q/(n1 n2)) equals (const)·q_e^{-1/6}(q_{n1}q_{n2})^{1/6}Q^{-1/6}κ(Ω') + O(1). This follows from a lattice-point count and the substitution w = z^6. The r = 0 term, of size O(1), is excluded and is lower order because P_a^{1/6} → ∞.
3. **Column sums.** The extra coprimality conditions (n,r) = 1 change only the density, by a factor 1 + O(q_r-divisor terms). They are then ordinary smooth sums of a character mod 36·b*. That gives Y·Col + O(Y^{1/2+ε}).
4. **Other C and d.** The paper's dyadic weight is C_0^{-1/6}D_0^{-1} after summing over each dyad. That is q_C^{-7/6} and q_d^{-2} per individual modulus, so the full exceptional constant is the C = d = 1 term times an Euler product Π_{p∉S}(1+O(q_p^{-7/6})) = 1+O(η_S).
5. **Positivity.** w_e > 0 and κ(Ω') > 0, so the form is >= w_{e_0}κ|Col(e_0)|^2.

The non-exceptional terms are exactly the paper's O((Q/Y)(1+P_a^2/Y)Z^ε). At the optimum P_a^2/Y = Z^{-0.29}, which is negligible.

### Lemma 3
Use the v-separated form, row weight Ω_v(q_m/Q). Poisson in m modulo a·b* has length Q and prefactor Q/(q_a q_{b*}). The complete sum is τ·q_a·ξ̄(·)χ_a(H)·(twist⁻¹); by construction the twist cancels c_a's ψ. For H = r^6 (u = 1) we have χ_a(r^6) = 1_{(a,r)=1}, and

Σ_r Ω̂_v(r^6√Q/|a|) ≈ (q_a/Q)^{1/6}·κ(Ω_v).

So the a-weight is W_1(t)t^{-1+iv}·t^{1/6}: the same column power and the same Mellin point as in Lemma 2.

Cross-check with B0.28: X^{1/2-1/6} = X^{1/3}, and κ(Ω_v) corresponds to 𝓜(1/6) > 0, which the paper proves by the same Riesz-kernel positivity.

### Proposition O
For each piece (σ, v), Lemmas 2 and 3 give

||A_{σ,v}||_{Ω'} >= c(Q/Y)^{1/2}P_a^{1/12}|Col_{σ,v}(e_0)| = c(Q/Y)^{1/2}P_a^{1/12}|Sig_{σ,v}|/κ(Ω_v).

Integrating ∫|Ŵ_0(iv)|·…dv and summing over σ, the triangle inequality gives a bound >= c' P_a^{1/12}X^{1/2}·|Σ_{σ,v} signal|. This needs |κ(Ω_v)| <= κ(Ω') whenever Ω' >= |Ω| pointwise, which holds because the Riesz kernel is positive.

Each lever is then closed as follows.

* **S-class column weights or a change of ξ.** Lemma 1 makes the e_0 term equal to |signal|^2 class by class.
* **Radial W_0 and W_1.** These enter both quantities through the identical Mellin value at t^{-5/6+iv}.
* **Row weight.** κ(Ω') > 0 for every Ω' >= 0, by the Riesz kernel.
* **Twist by χ̄_a(u_0), u_0 ∉ S.** The signal moves to row u_0 and gains q_{u0}^{-1/6} from q_u^{-z}. The exceptional term needs k = u_0 e r^6, giving (P_a/q_{u0})^{1/6}. So bound/signal = q_{u0}^{1/6}(1+(P_a/q_{u0})^{1/6})^{1/2}, which is increasing in q_{u0}. The optimum is u_0 = 1, and the target also turns into a moving-conductor character.
* **Recombining ray classes, Ã_τ = Σ_σ R(τ,σ)A_σ.** This is again an S-class column weight, so Lemma 1 applies (R is a mod-4 bicharacter).

## 3. Uncertain steps

* **U1.** Lemma 2 is a main-term asymptotic. I did not write out the error terms in steps 3 and 4 beyond the paper's own bounds. The Euler-product step treats the paper's dyadic upper bounds as per-modulus bounds. This needs a large-cutoff S; enlarging S is allowed in the paper.
* **U2.** The claim in Lemma 3 that the Mellin points match was checked by scaling (powers of q_a, X and Q), not by a full residue computation.
* **U3.** The bound/signal comparison sums the pieces (σ, v) with absolute values. Cancellation between v-pieces in the signal only makes the obstruction stronger, so this direction is safe.
* **U4.** Item (iv) (y < x is worse) was not rechecked.

## 4. What was not shown, and the precise obstruction

**(a) Flaw in (iii).** Invariance under m → m r^6 is shared by the whole vector A, not only by A^{exc}. For (n,r) = 1, g(n,-m r^6) = χ̄_n(r^6)g(n,-m) = g(n,-m). Invariance therefore says nothing about alignment between T and A^{exc}, and CS sharpness for the actual inner product <A,T> is unproved. Only "CS cannot be beaten using norms alone" is established (Prop. O).

**(b) Heuristic: CS is not sharp for the true sum.** Suppose all zeros of L(s,·) satisfied Re s <= 0.85. The admissible range is κ >= 7/10, and κ = 2β−1 = 0.7 is the edge. Then the paper's own B0.51 machinery gives |I| << Z^{0.85-k} = Z^{0.169}, which is below X^{1/2} = Z^{0.187}. So there is no true-size barrier at P_a^{1/12}. The loss is purely a method loss.

**(c) Identification of A^{exc}.** The exceptional piece is the Patterson (sextic metaplectic Eisenstein) residue of the column Gauss-sum series: A^{exc}_m ≈ Y^{-1/3}ρ_ψ(m). The mean of |ρ|^2 is about Q^{-1/6}, so A^{exc} dominates exactly when Y^2 > Q (P_a > 1), and this reproduces (Q/Y)P_a^{1/6}. The Gram r^6-sum is the analogue of B0.28's ζ_K(6z) pole. Both are the same Riesz-1/6 phenomenon.

**(d) The remaining route needs a new estimate.** The projection route needs

|Σ_m Ω ξ(m) A^{exc}_{m,σ,v} T_σ(m)| << Q^{1/2}X^{1/2}P_a^{c/2}Z^ε with c < 1/6.

This is a correlation between sextic theta coefficients and the cubic theta sum. The coefficients for n = 6 are not explicitly known, so it is unavailable.

The equivalent dual route is CS over the columns a. It gives |I| <= Y^{-1/2}Q^{-1/2}E_a^{1/2}, where E_a = Σ_{a~Y}|Σ_m ξ(m)χ_a(m)T(m)Ω|^2.
* The sixth-power diagonal of E_a is harmless. Using (B0.26) dyadically, its exponent is at most M (numerical maximum 0.833 = M).
* E_a << Q(Q+Y)Z^ε would remove the exceptional term entirely (c = 0).
* Generic sextic large sieves do not reach this. A Heath-Brown-type (YQ)^{2/3} term gives probe exponent 0.2010, worse than 0.1938.
* By duality, the operator norm is necessarily >= Q P_a^{1/6}. So any such gain must use the arithmetic of T, for example an a-twisted version of the B2 reflected energy, which is not in the paper.
