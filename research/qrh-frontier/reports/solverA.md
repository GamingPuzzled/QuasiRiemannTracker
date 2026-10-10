# Solver A: lowering the slot cost 6κ in MP (B0.5)/(B3.40)

## Outcome

I did **not** prove MP with slot cost c·κ for any c < 6. This is a failure report. I do prove an
optimality statement for the paper's method: no local improvement to the B3.5–B3.10 induction can
lower the 6. That holds whichever of routes (a), (b) or (c) is tried. In the application range
(n1 = n2 = n ≤ 1/2) the bound 2n + 6κz ≤ M is already the best the framework allows, **even with
plain reflection at every child and with unused allowance carried from node to node**. Any gain
needs a new analytic input. Section 4 lists three candidate inputs.

---

## 1. Statements

### 1a. Claimed (negative / structural) proposition

Consider any argument built only from these parts:

* (P1) the first Poisson step (B3.52) followed by Cauchy–Schwarz into the positive Gauss-row norms
  (B3.54);
* (P2) amplification by row multipliers m that leave every Gauss row invariant,
  G(u, h m) = G(u, h) for every admissible column u (B3.58);
* (P3) the second transform (B3.60)–(B3.63);
* (P4) whole-plain reflection (B3.44)–(B3.45), which may be used at any node, including
  positive-slot nodes;
* (P5) pointwise slot removal through (B3.47), at cost κ per unit of slot length;
* (P6) allowances that are carried, so unused allowance at one node may be spent at its
  descendants (the paper discards it; allowing it only helps).

Such an argument proves Σ_k |S_k(Z^n)|^4 |Q_{k,I}|^2 ≪ Z^{M+ε}, with no loss, only if

    2·min(n, M−n) + 6κz ≤ M.

So for n ≤ M/2 the paper's hypothesis n1 + n2 + 6κz ≤ M is sharp for this framework. For
n > M/2 the only improvement is the trivial one: reflect the plain once at the start, replacing
n by M − n.

### 1b. What would have been needed

To reach c < 6, at least one of the following must fail. Each holds as an identity or is
unavoidable inside the framework (proofs in §2):

1. **Generic step.** A generic amplified child has A′ = A, M′ = M − 6ℓ and allowance
   Δ = ℓ, with equality in (B3.64) and (B3.70).
2. **Amplifier.** Every invariance amplifier satisfies (row growth) ≥ 6 × (density gain).
3. **Diagonal.** Every node needs A_eff − M ≤ (its current allowance), from (B3.61) and (B3.80).
4. **Slots.** Slots cannot be reflected. Their only removal mechanism costs κ per unit, (B3.47).

---

## 2. Proofs

### 2.1 Fact 1: the generic child saturates (B3.70) and (B3.64)

Take the configuration with no common support at either Poisson step:

* c = d = R = E = 0 and p = s_0 = 0 in §B3.6, so q̃ = q and K = K_0 = 2A − m (B3.53);
* main amplified norm, so w = w_o = 0, ℓ = ℓ_* = σ/3 and g = J_+ + 2σ (§B3.7);
* D_2 = E_2 = 1, so c_2 = d_2 = b_2 = p_2 = g_2 = t_2 = V = 0 (§B3.8).

This is the dominant (typical) part of the sum: coprime columns, no shared primes.

With these values:

* J = d − c + (K_0 − K) − 2w + w_o = 0, so g = 2σ = 6ℓ.
* a_0 = A − c − w = A.
* From (B3.63): m′ = 2a_0 − K − g = m − 6ℓ, q′ = q, and therefore M′ = M − 6ℓ.
* The child column is the same product l_1·l_2·(slots), so A′ = A − c − w − c_2 = A.
* (B3.70) becomes (A′ − M′) − (A − M) = g = 6ℓ, which is **equality** in the bound 6(w + ℓ).
* (B3.64) gives Δ = b_2 − p_2 + w + B_c + ℓ = ℓ, which is **equality** in Δ ≥ w + ℓ.

**Consequence for route (c).** Tighter accounting cannot help. Every slack term in (B3.70) is a
nonnegative common-support length: c, w, d + (K_0 − K) when J < 0, and c_2 − g_2 + t_2. The
extra terms b_2 − p_2 and B_c in Δ are of the same kind. All of them vanish on the generic part.

The defect growth is not caused by the slots. It comes from the row shrinking (M′ = M − 6ℓ) while
the column length stays fixed (A′ = A).

* For plain factors this growth does no harm, because reflection (B3.44)–(B3.46) shortens them
  whenever n > M′/2. This is why MP0 holds for all lengths.
* Slots have no functional equation. Their length therefore stays in A_eff, and only they turn
  the growth into a cost.

So the question "is the full 6(w + ℓ) charged to the slot part?" has the answer: it is charged to
A − M as a whole, and the slots are the only part of A that cannot discharge it.

### 2.2 Fact 2: the sixth-power amplifier is optimal among invariance amplifiers

Suppose X(hm) = X(h) for every admissible column polynomial X(h) = Σ_u B(u) G(u, h). The columns
u include every good prime π in the plain range. By (B3.58) and (B1.13), G(π, hm) is
conj(χ_π(m)) G(π, h). So we need χ_π(m) = 1 for every good prime π of the relevant norm.

The sextic symbols χ_π(·), π ranging over primes, separate K*/(K*)^6 away from S. This is Kummer
theory / Chebotarev applied to K(m^{1/6})/K. Hence m is a unit times a sixth power, up to S-units.

There are Y^{1/6+o(1)} such m with q_m ≤ Y. A multiplier set with density gain Z^ℓ therefore needs
multipliers of norm up to Z^{6ℓ}. The amplified rows hm then fill a ball enlarged by 6ℓ.

Since G(u, ·) is not invariant under any non-sixth-power m, a "cheaper" amplifier would have to
compare *different* polynomials X_m. Such polynomials cannot be merged into one ball norm without
losing exactly the factor |T| that was gained.

**Route (b): using the slot primes themselves.** The slot primes lie in the column, not the row,
and are already summed. Using them as amplifier primes changes nothing in this count.

**Cubes or squares as multipliers.** χ_u(p^3) is the quadratic symbol (p/u)_2, a character in u
modulo p. Multiplying h by p^3 therefore twists the column by a moving factor. That twist is
admissible in the τ-family (width +ℓ), but it differs for each p. Recombining the twisted copies
costs Z^ℓ, which cancels the gain. I checked this explicitly and the net gain is zero.

### 2.3 Fact 3: every node needs A_eff ≤ M + L

**Diagonal term.** At a node with current allowance L, the diagonal term j = 0 of (B3.60) has
excess (A − M) + g − ℓ − … over (B3.57), by (B3.61). The exceptional term (B3.80) has excess
(A − M)_+.

**Why it cannot be avoided.** After the Cauchy step (B3.54), N_C is a positive norm over a row
ball of length K ≈ 2A − M. Its column has length ≈ A − c. If A > M, then K exceeds the column
length, and Σ_h |X(h)|^2 ≥ (#h)·‖B‖^2 up to constants for any coefficients. That is the
large-sieve diagonal, which no later step recovers. The bound then gives Z^{A} rather than Z^{M}.

**Zero-frequency terms.** The zero-frequency and principal terms are O(Z^{M}) (§B3.6 and (B3.21)),
so L ≥ 0 at every node.

### 2.4 Fact 4: slots cannot be reflected; removal costs κ per unit

* Q_{k,i} is a prime sum. The only row-uniform bound available from the hypothesis (a zero-free
  region to (1+κ)/2) is (B3.47), |Q_I|^2 ≪ Z^{κz+ε}.
* This bound is sharp row by row under that hypothesis alone: a zero at β = (1+κ)/2 near the real
  line is not excluded for any individual ψ_k ν_i.
* Averaged removal (large values / zero density) would require controlling Σ_{k∈V} |S_k|^4 over
  sparse sets V. Through Hölder this needs Σ_k |S_k|^{4p} with p > 1. The family does not contain
  this moment: it has only two plain factors. Section 4 gives details.

### 2.5 The optimization (proof of 1a)

**Model.** Follow one branch of generic children. Let M run down from M_0 = 1, using reflection
whenever it helps. By Facts 1 and 2:

* the cumulative allowance available at width M is at most G(M) = (1 − M)/6;
* if r(M) is the slot length removed so far (nondecreasing), the remaining allowance is
  L = G − κr ≥ 0;
* by Fact 3, every node needs

      2·min(n, M − n) + (z − r) − M ≤ G(M) − κ r(M).                   (★)

Non-generic children have Δ > ℓ or smaller A′. The adversary can always choose the generic branch,
since it carries the main term. So (★) along the generic branch is necessary.

**Case M ∈ [2n, 1].** Here (★) reads

    (1 − κ) r ≥ 2n + z − M − (1 − M)/6.

* The right side increases as M decreases, because its derivative in M is −5/6.
* The budget gives κ r ≤ (1 − M)/6.

**Case M < 2n.** Both plains are reflected, and the right side becomes M − 2n + z − (1 − M)/6.
This decreases as M decreases. So the binding width is M = 2n.

**At M = 2n.** Take r_0 ≤ min{z, (1 − 2n)/(6κ)}. Then (★) requires

    z ≤ (1 − 2n)/6 + (1 − κ) r_0 ≤ (1 − 2n)/6 + (1 − κ)(1 − 2n)/(6κ) = (1 − 2n)/(6κ).

This is exactly 2n + 6κz ≤ 1.

**Sufficiency.** Removing all slots by width 2n meets the terminal condition κz ≤ G. This is the
paper's schedule.

**Numerical check.** I ran a small script (`lp.py` in the same scratchpad directory) that searches
over removal schedules in this model, with reflection allowed. It returns
z_max = (1 − 2n)/(6κ) to five decimals:

| n | κ | z_max |
|---|---|---|
| 0.30 | 0.70 | 0.09524 |
| 0.30 | 0.7498 | 0.08891 |
| 0.42 | 0.70 | 0.03810 |
| 0.42 | 0.7498 | 0.03557 |
| 0.48 | 0.70 | 0.00952 |
| 0.48 | 0.7498 | 0.00889 |

**General ρ.** With amplifier ratio ρ in place of 6, the same computation gives
z_max = (1 − 2n)/(ρκ). So c = ρ exactly, which ties the record bound directly to the degree of the
character family.

---

## 3. Route (a), examined in detail

### (a1) Putting the plains-and-slots column into canonical form (B3.6)/(B3.7)

**Why this would be very valuable.** In the application the problem is subcritical
(A = 2n + z < 1). The first-Poisson dual has row K = 2A − 1 < A, so it sits in a long-polynomial
regime. If the dual column were canonical, (B3.7) would give:

* first margin: F − M − z_0 = (1 − 2n − z) − z, i.e. 2n + 2z < 1;
* second margin: 4F − 3M − 6z_0 = 4 − 3·… > 0, easily.

That is slot cost 2, the same as MI: c = 2/κ ≈ 2.7.

**Why it fails.**

1. **The coefficients are not theta coefficients.** For squarefree u, the dual coefficient of a
   plain variable is γ_1(u) = μ(u)G(u)·conj(a_0(u)). This follows from (B0.2):
   γ_1 = μαG/γ_2 and γ_2^3 = μα. So the dual of a plain sum is an *inverse-theta* polynomial.

   By contrast, MI's inverse polynomial dualizes to a_0, the theta coefficient. The canonical
   estimate needs that coefficient exactly ("no further coefficient on n is permitted", after
   (B3.6)), because its terminal step completes and reflects the theta sum (B3.8) via RE. An
   inverse-theta column has no completion and no reflection. That is precisely why §B3.6–B3.8
   recurse back to a plain moment instead of terminating.

2. **The plain product is a twisted divisor sum, not a theta series.** The column l_1 l_2 carries
   a coefficient of divisor type. Since a_0(l_1 l_2) = a_0(l_1) a_0(l_2) χ_{l_2}(l_1)^4, one
   could treat l_2 as the label f. But Cauchy over the label costs Z^{V} = Z^{n_2}. In
   normalized terms, Σ_h |X|^2 ≤ Z^{V}·𝓔 with 𝓔 ≪ Z^{F}, a fatal loss.

### (a2) Separating the slots by Cauchy–Schwarz or Hölder, MI-style

Writing |S|^4 |Q|^2 as |S|^2 · |SQ|^2, or moving |Q|^2 to one Cauchy side as
(S S Q Q̄) · conj(S S), changes nothing useful:

* The first option needs max_k |S|^2 ≤ Z^{n} or the moment Σ|S|^4|Q|^4. Both are worse.
* The second makes one Cauchy side carry slot length 2z. That side is again a plain moment with
  slot length 2z. Balancing the two sides (§B3.6 ledger) is optimal when they are symmetric, so
  this is strictly worse.

### (a3) Expanding |Q_k|^2 into twisted fourth moments

Expanding |Q_k|^2 gives Σ_{p≠p′} Σ_k χ_{p p′^5}(k) |S_k|^4. To use this we would need a
square-root saving T^{1/2}, with T = Z^{2z}, in each twisted fourth moment. That is equivalent to
the original problem: Poisson leads back to the same Gauss-row norm.

---

## 4. Precise obstruction, and what a new proof would need

**Obstruction.** In the generic configuration of §2.1:

* (B3.70) holds with equality, A′ − M′ = A − M + 6ℓ;
* (B3.64) holds with equality, Δ = ℓ;
* the ratio 6 is the minimum possible by the sixth-power count (§2.2);
* slots cannot be reflected and cost κ per unit (B3.47);
* every node must satisfy A_eff ≤ M + L ((B3.61), (B3.80), forced by the positivity step
  (B3.54)).

The schedule optimization (§2.5) then gives exactly 2n + 6κz ≤ M. The binding width is M = 2n,
where the plains have stopped being reflectable "for free".

**Inputs that would beat it.** Any one of these would suffice, but none is in the paper:

1. **A bilinear (non-Cauchy) bound** for the first-Poisson form (B3.52) when A > M, using that
   the slot columns are prime sums with fixed ray characters. It would have to beat the
   large-sieve diagonal Z^{A}.
2. **A reflection-type identity for the slot factor.** Equivalently, a way to give the
   prime-sum part of the column a dual of length < z.
3. **A sparse-row fourth-moment bound** of the form

       Σ_{k∈V} |S_k(Z^n)|^4 ≪ |V|^{θ} Z^{M(1−θ)+ε}   for some θ > 0,

   or an eighth moment of S. Either would allow averaged rather than pointwise slot removal
   (zero-density style). Below the theoretical limit 1/κ, the gain in c depends on θ.

---

## 5. Steps I am not sure of

* **Fact 2 (Kummer separation).** The step "χ_π(m) = 1 for all good π in the plain range implies
  m is a unit times a sixth power" uses Chebotarev for K(m^{1/6}) with m of polynomial norm.
  This holds for primes of norm ≥ Z^{δ} once m is fixed. I did not check how it interacts with
  the masks and S-units. Those only change bounded factors and cannot change the exponent 1/6.
* **The proposition in §1a is a statement about the framework.** It assumes the argument
  follows a generic branch, and its necessary condition (★) uses the heuristic that the
  large-sieve diagonal is attained after Cauchy. That is a lower bound for the *method*, not for
  the true moment. The true moment is plausibly ≪ Z^{M} up to 2n + z ≤ M (c = 1/κ) under GRH-type
  heuristics. With arbitrary slot coefficients, generic large-sieve heuristics suggest the same
  limit.
* **Non-generic children.** I did not prove that a cleverer argument could not deliberately send
  mass into non-generic children (c, c_2, w > 0) to gain extra Δ. The generic part has full
  measure and is bounded separately, so I believe this is impossible.
* **The approximate-Q bound (B3.47).** I took it as the only pointwise input and did not
  re-derive it.
