# Solver I: B4, simultaneous zeros and mean values of inverse products

**Verdict.**
- **Nothing new is proved unconditionally, and θ does not move.**
- **B4 is now priced exactly.** The decisive finding is that the second character is irrelevant. Every gain a joint
  bound for M_u(ψ̄)·M_u(ηψ̄) could give at the binding bin is already available from a **single-character,
  two-length** product moment Σ_u |M_u(φ;U^{r1}) M_u(φ;U^{r2})|². That version also avoids an off-diagonal loss.
- **The product moment has a sharp break-even slope.** Write e = 1 + λ(r1+r2−1) beyond the conductor.
  - It pays **iff λ < 5/6**, provided it is combined with a Poisson-range "MI2" (an MI for degree-2 Möbius
    polynomials).
  - Without MI2 the break-even is 0.6743 (UM2 only) or 0.5908 (nothing in the Poisson range).
- **The 5/6 is structural.** It is the slope of the paper's sixth-power lift (B0.9).
- **Route (a), even if completed, ends in an exact tie.** Extending UM/MI to degree-2 Möbius coefficients and lifting
  gives λ = 5/6 on the nose, so the gain is zero.
- **Routes (b) and (c) give provably nothing.**
- **B4 is therefore another face of (★).** It is worth roughly half of (★) at equal λ, behind the same barrier.

Scripts are in `research/qrh-frontier/scripts/simul/`. No existing repository file was edited.

---

## 0. Premise: binding rows carry both zeros (checked)

In (B0.53) the s-line is moved to a+16e using **only** the denominator buffer. The paper says: "The denominator's
zero-free buffer permits this last move". The main-slot amplitude bound g_i ≤ δ/2 is B0.16, which applies to the
numerator presentation ψ̄_u only.

Refine each bin by two levels:
- δ_D, the zero level of L(ηψ̄_u);
- δ_N, the zero level of L(ψ̄_u).

A common buffered height for both comes from the same pigeonhole as (B0.13), using about 2/e heights. Then:
- **δ_D < δ_N.** Use the numerator's zero as witness. F equals F_paper(δ_N, q/δ_N) − (δ_N−δ_D)/2.
- **δ_D > δ_N.** Use the denominator's witness with p = q/δ_D ≤ δ_N/(2δ_D) < 1/2. Then
  F = F_paper(δ_D, p) < b, by ∂F/∂p = +0.0324.

So the bound is attained only at δ_N = δ_D = δ0 with p = 1/2, which confirms solverH. Off the diagonal the paper's
margin is only about 0.0345·|δ_D−δ_N|; this matters in §1.6.

## 1. Pricing

### 1.1 The count system
Lemma 1W gives witnesses (r_i, m_i) with r_i ≤ t_i, r_i+m_i ≥ t_i and m_i ≤ 1/2. Its proof works for any class-fixed
t0 ∈ (0, 3/2], not only [1, 3/2] (Claim 1). Apply it twice: either once per character (B4 as posed), or **twice to
the same zero with t1 ≠ t2** (single-character version). This gives two inverse witnesses |M(U^{r_i})|² ≥ U^{δ r_i}.

The single lines are the paper's lines, with the supply cap z ≤ L = ℓ/d made explicit (it matters for short r):
- A_I(r) = 1 − δr − 2q·min((1−r)/2, L), for r ≤ 1;
- long(r) = (1+5r)/6 − δr, for r > 1;
- A_P(m) = 1 − 2δm − 2q·min((1−2m)/(6κ), L).

The hypothetical **joint line** is e(r1,r2) − δ(r1+r2), with x = r1+r2. Then

    R_new = inf_{t1,t2} sup_{adversary (r1,r2)} min(all lines).

**Leakage pitfall (found and removed).** A naive family e = 1 + λ(x−1)_+ lets the optimizer take t1 → 0. That turns
the "joint" line into an improved *single* long bound, which is just (★) in disguise. All results below use

    e(r1,r2) = max( E(r1), E(r2), 1 + λ(r1+r2−1)_+ ),   E(r) = max(1, (1+5r)/6).

With this family the hypothesis implies nothing new for one factor alone.

There are three Poisson-range variants for x ≤ 1:
- **P0**: no joint bound there.
- **P1**: UM2, meaning Σ|M1M2|² ≪ U.
- **P2**: MI2, meaning Σ|M(U^{r1})M(U^{r2})Q_I|² ≪ U whenever r1+r2+2z ≤ 1 (and 2x+8z ≤ 3). Each Möbius factor
  costs 1 per unit length, as in MI.

### 1.2 Exact solution (closed forms, verified against an exact LP)
`joint_lp.py` solves the max-min exactly. The lines are piecewise affine, so the adversary problem is a maximum of
small LPs, with an outer Nelder–Mead over (t1, t2). `closed_forms.py` gives the equilibria.

At a bin, put A = 1−q, s = δ−q, B = 1−q/(3κ), w = δ−q/(3κ). Define the level-v crossings:
- r_v = (A−v)/s;
- m_v = (B−v)/(2w);
- t_v = (v−1/6)/(5/6−δ);
- x_b(v) = 1 + (v−1+δ)/(λ−δ).

The equilibria are:
- **P0** (best play is t1 = t2 ↑ t_v): 2r_v − 1 = x_b(v) − 1. The adversary's longest joint length is 2r_v.
- **P2**: the optimal play is asymmetric, with t2 ↑ t_v and t1 ∈ (m_v + t_s − t_v, x_b − r_v), where t_s = r_v + m_v.
  The adversary is squeezed between two lines:
  - the MI2 line at small x, which is below v iff x > r_v;
  - the λ-line at large x.

  This gives the equation 2(r_v+m_v) − t_v = x_b(v).
- **P1**: (1−v)/δ + 2m_v − t_v + r_v = x_b(v).

**LP versus closed form at the binding bin.** The table shows R_new. The LP and the closed form agree to all printed
digits (`out_lp_tmin0.txt`, `out_checkP1.txt`).

| λ | P0 (t ≥ 1) | P2 (closed) | P2 (LP) |
|---|---|---|---|
| 5/6 | R | R | R (ΔR = 0) |
| .80 | R | .589362 | – |
| .75 | R | .587562 | .587562 |
| .65 | R | .581819 | .581819 |
| .60 | R | .576816 | .576816 |
| .58 | .587604 | .574046 | .574046 |
| .55 | .578454 | .568523 | .568523 |

At λ = .50 the LP also gives exactly v_P2 = .552100. An earlier figure of .548397 came from a stale, buggy run and
is superseded by the regenerated `out_lp_tmin0.txt`.

### 1.3 Break-even slopes at the binding bin (`breakevens_exact.py`, 40 digits)
At the binding bin, R0 = .590332353526, t* = 1.163953576509, r_R = .745699942724 and m_R = .418253633785, with
r_R + m_R = t*.

| variant | λ_be |
|---|---|
| P0 | δ0 + (R0−1+δ0)/(2r_R−1) = **0.5907881353** |
| P1 | **0.6743482017** |
| P2 | δ0 + (R0−1+δ0)/(t*−1) = **5/6 exactly** |

**Why exactly 5/6.** From (B0.40), R−(1−δ) = δCV/(2(V+δC)) and t*−1 = δC/(2(V+δC)), so the ratio is V = 5/6 − δ in
every bin.

**Structural reading.** The P2 configuration builds a **virtual long witness** of length x = t* out of a short
factor (≈ m_R ≈ .42) and a medium factor (≈ r_R ≈ .75). A product bound then beats the paper iff it beats the long
bound E(t*) at the same total length. The full-minimax break-even agrees: there is no gain at λ = .592 for P0
(`out_breakeven.txt`) and none at 5/6 for P2.

### 1.4 θ as a function of λ (full minimax re-solve, harness model `slack2.solve`, grid error ≈ 3e-8)
Two single-character tables follow. The third column of the first table is the gain for (★) at equal λ, taken from
FINDINGS for comparison.

**Single-character P2 (MI2 plus slope λ)**, from `price_P2.py` / `out_P2.txt`:

| λ | θ | gain | (★) at equal λ |
|---|---|---|---|
| 5/6 | 0.8748866157 | 0 | 0 |
| .83 | 0.8748742394 | 1.24e-5 | – |
| .82 | 0.8748360194 | 5.06e-5 | 1.0e-4 |
| .80 | 0.8747543020 | 1.32e-4 | 2.6e-4 |
| .75 | 0.8745129021 | 3.74e-4 | 7.1e-4 |
| .70 | 0.8742005251 | 6.86e-4 | – |
| .60 | 0.8732952935 | 1.59e-3 | 2.4e-3 |

Near 5/6 the slope is dθ/dλ ≈ 3.7e-3, about half of (★)'s 7.5e-3. At the P2 equilibrium the factor lengths are
r1 ≈ 0.43–0.55 and r2 ≈ 0.75–0.80, so x ≈ 1.18–1.35.

**Single-character P0**, from `price_theta_minimax.py` / `out_theta_minimax.txt`:

| λ | gain |
|---|---|
| ≥ .592 | 0 |
| .58 | 1.03e-4 |
| .55 | 4.28e-4 |
| .50 | 1.07e-3 |

**A caveat on P3 (ΔR = −1e-3 ⇒ Δθ ≈ −1.4e-4).** This holds only for a *uniform* ΔR. A joint line reduces R only in
bins with δ > δ_c(λ), since λ_be(δ) increases with δ (`price_theta.py`). So the δ-saddle slides left, and the
realized gain is curvature-limited. For P0 at λ = .58, ΔR(δ0) = −2.7e-3 but Δθ = 1.0e-4, not 3.7e-4.

### 1.5 The second character is irrelevant
On diagonal classes (δ_N = δ_D, the same slots) the two-character and single-character count systems are identical
line by line. Off the diagonal the two-character joint line is weaker, at −δ_N r1 − δ_D r2, and the paper's slack
there is only about 0.0345·|δ_D−δ_N|.

The class-space minimax over (δ_N, δ_D, q) (`minimax3.py`, `out_minimax3.txt`, P0 family) puts the argmax
off-diagonal for λ ≤ .585. Its gains, for example 1.9e-4 at λ = .55 against 4.3e-4 single-character, are smaller.
They are at the level of that solver's ~5e-5 grid error for λ ≥ .58.

**Recommendation.** If anything in B4 is pursued, state it for one character with two lengths. Then the target η and
"simultaneity" play no role.

## 2. Attack

### (a) Extending UM/MI (B3.31, B3.39) to μφ1 * μφ2
- **Even full success yields nothing.**
  - Suppose UM2 holds in the Poisson range: Σ_{q_v≤H} |M_v(D1)M_v(D2)|² ≪ H^{1+ε} for H ≥ (D1D2)^{1+ϑ}.
  - Apply (B0.7) to each factor with the same a, take a two-scale Sobolev supremum, and average over q_a ≤ P with
    H = UP⁶.
  - This gives exponent max(1, (1+5x)/6), which is λ = 5/6 exactly (Claim 4).
  - Together with MI2 this is the exact P2 tie (Claim 2): **zero gain**.
  - To beat it, the product's lift slope must drop below 5/6. That is the (★) barrier, now for products.
- **Even UM2 is not a corollary of the canonical estimate.**
  - In the B3.31 proof, Poisson in the row turns μ(c)·(Gauss sum of χ̄_c) into the theta coefficient
    a_0 = ᾱγ₂, via (B0.2) μ(n)γ_{−1}(n) = χ_n(−1)G(n)^{−1}ᾱ(n)γ₂(n).
  - For c = μφ1 * μφ2, after (B3.3) separation the squarefree coefficient is μ(n)g(n), with g multiplicative and
    g(p) = ν1(p)q_p^{−it1} + ν2(p)q_p^{−it2}. The dual coefficient becomes a_0(z)g(z), but (B3.6) admits "no further
    coefficient on n".
  - Using a_0(ab) = a_0(a)a_0(b)χ_b(a)⁴ rewrites g as a label sum: f = a enters through χ_n(f)⁴. However, the a-sum
    is **coherent**, inside |·|², while the canonical estimate controls only Z^{−V}Σ_f w(f)|·|².
  - Cauchy–Schwarz over a costs #a ≈ D1. This returns the trivial pointwise-type bound.
  - Squareful collisions (p | a, p | b) also leave the squarefree family.
  - What would be needed is a **label-coherent (bilinear) RE**.
- **No inverse analogue of MP0.**
  - MP0's all-length statement rests on whole-factor reflection (B3.44), which brings each plain length down to ≤ M/2
    and then uses the uncentered Poisson core.
  - 1/L has no reflection, so only the Poisson core has an inverse analogue. By the above, that core is worth at most
    a tie.
- **Cost-2 version of MI2.**
  - MI2 needs cost 1 per unit of Möbius length. Taking the rough part of the short factor as MI prime slots costs 2
    per unit (arbitrary slot coefficients).
  - That line is 1 − q − (δ−q)r2 − (δ−2q)r1, which equals A_I(r2) at p = 1/2. It is useless (Claim 5(ii)).
  - The P2 mechanism needs cost 1 even though its short factor at the small-x end is tiny (r1 ≈ 0.02–0.1).

### (b) Cauchy–Schwarz, Hölder and pointwise bounds
- **Pointwise times moment gives nothing new.** Since π_j ≥ w_j on the class,
  Σ|M1M2|² ≤ U^{π(M2)}Σ|M1|² gives a line e1 + π2 − w1 − w2 ≥ e1 − w1.
- **Hölder gives a convex combination of existing lines** (P4 averaging lemma).
- **Cauchy–Schwarz to fourth moments** needs Σ|M|⁴ with M² of length x > 1. That is the same open problem; for
  unequal lengths it is also lossy (the average of e4(r1) and e4(r2) exceeds e(r1+r2)).
- **Conclusion: no gain** (Claim 5).

### (c) η is a fixed character
- If η restricted to ideals prime to S equals n ↦ χ_n(e) for an S-unit e (order dividing 6; see Tick 4), then
  ηψ̄_u = ψ_{u⁵e} on (n, eu) = 1.
- The joint moment then becomes a shifted-row correlation Σ_v |M_v(D1)M_{ev}(D2)|². It is again a fixed-coefficient
  degree-2 object of length x > 1, so no new mechanism appears, and Cauchy–Schwarz reduces it to fourth moments.
- For general finite-order η, ηψ̄_u lies outside the sextic family.
- Most decisively, §1.5 shows that correlating the two characters cannot add anything at the binding bin: one
  character already realizes every diagonal gain.

### Precise obstruction
B4 reduces to the single-character statement

    Σ^{(6)}_{q_u≍U} |M_u(U^{r1}) M_u(U^{r2})|² ≪ U^{1+λ(r1+r2−1)+ε},  r1≈.42–.55, r2≈.75–.80, 1 < r1+r2 ≤ 1.4,

with **λ < 5/6**. It also needs MI2 (cost-1 degree-2 Möbius MI) in the Poisson range. The method's best possible
output, a Poisson-range UM2/MI2 plus the sixth-power lift, is exactly λ = 5/6. Beating it is a mean-value theorem for
a degree-2 inverse polynomial beyond the conductor, which is the same barrier as (★).

**Advice for the KB.** Record B4 as failed (F12) with this pricing, and keep the (★) family as the single target.

## 3. Scripts (`research/qrh-frontier/scripts/simul/`)
- `base.py`: copy of smallgain/base.py (exact record point).
- `joint_lp.py`: exact max-min LP. With an argument t_min it prints the λ table.
  Outputs: `out_lp_tmin0.txt`, `out_lp_tmin1.txt`.
- `closed_forms.py`: the v_P0, v_P1, v_P2 equilibria and regime checks.
- `breakevens_exact.py`: break-evens to 40 digits.
- `check_P1.py`: LP check of P1. Output: `out_checkP1.txt`.
- `price_P2.py`: full θ minimax for P0/P1/P2, single character. Output: `out_P2.txt`.
- `price_theta.py`: helper functions, and λ_be(δ) across bins. (Its `__main__` smooth re-solve is unreliable because
  it jumps branches; use the minimax files instead.)
- `price_theta_minimax.py`: P0, single character. Output: `out_theta_minimax.txt`.
- `breakeven_scan.py`: global P0 break-even. Output: `out_breakeven.txt`.
- `minimax3.py`: two-character class space (δ_N, δ_D, q). Output: `out_minimax3.txt`.

---

## CLAIMS

**Claim 1 (Lemma 1W for every truncation; two truncations).**
- *Hypotheses:* the setting of Lemma 1W (B0.17)–(B0.18).
- *Statement:* the conclusion holds for every class-fixed t0 ∈ (0, 3/2]. It may be applied twice to the same zero and
  presentation with two class-fixed truncations t1, t2, giving witness pairs (r1, m1) and (r2, m2) simultaneously.
- *Proof sketch:* t0 enters only through (i) the vanishing of μ*1 on 1 < q ≤ U^{t0}, (ii) the bound
  C_φ(ρ+z) ≪ U^{O(t0)}, beaten by U^{−10} for t0 ≤ 3/2, and (iii) the dyadic split. The step m ≤ 1/2 uses (B0.15)
  only. The remark after Lemma 1W allows class-fixed truncations. For two truncations, take the product of the two
  finite dyad label sets.
- *θ consequence:* none by itself; used by Claims 2 and 3.

**Claim 2 (exact break-evens of a two-witness product line).**
- *Hypotheses:* fix a bin (δ, p, κ) with r_R ∈ (1/2, 1], (1−r_R)/2 ≤ L, r_R ≥ m_R and V = 5/6−δ > 0. All of these hold
  at the binding bin (δ0, 1/2, κ0).
- *Count system:* the paper's lines (§1.1, with supply caps) for witnesses from Claim 1, plus the joint line
  e(r1,r2) − δ(r1+r2) with e = max(E(r1), E(r2), 1+λ(x−1)) for x > 1. For x ≤ 1, use one of P0, P1 or P2.
- *Statement:*
  - (i) Under P2, R_new < R iff λ < 5/6.
    - The "no gain" direction for λ ≥ 5/6 is proved by an explicit adversary strategy:
      - if some t_i ≥ t*, the adversary takes a long witness r_i ≥ t*, so J ≥ long(x) ≥ R;
      - otherwise, if the maximal x is at least t*, use it;
      - else the minimal x is at most r_R, and the MI2 line is ≥ R.
    - The gain direction for λ ∈ (δ, 5/6) comes from the explicit choice t2 ↑ t_v, t1 ∈ (m_v + t_s − t_v, x_b − r_v),
      which gives R_new ≤ v_P2 (closed form, §1.2).
  - (ii) Under P0, the break-even is δ + (R−1+δ)/(2r_R−1), which is 0.59079 at binding.
  - (iii) Under P1, the break-even is 0.67435 at binding.
  - The equality R_new = v and the "no gain" directions for (ii) and (iii) are verified by the exact LP to ≤ 1e-9.
- *Key identity:* (R−(1−δ))/(t*−1) = 5/6 − δ, from (B0.40).
- *Proof location:* §1.2–1.3; `closed_forms.py`, `joint_lp.py`, `breakevens_exact.py`.
- *θ consequence:* a product hypothesis gains iff its slope beyond the conductor is below 5/6 (with MI2), or below
  0.674 / 0.591 without it.

**Claim 3 (θ prices, numerical, harness model).**
- *Hypotheses:* MI2 together with H2_λ holds in all bins, single character.
- *Statement:* θ ≤ the values in the §1.4 table, for example λ = .80 ⇒ θ ≈ 0.8747543 and λ = .75 ⇒ θ ≈ 0.8745129.
  Near 5/6, dθ/dλ ≈ 3.7e-3. P0 values are in the second table.
- *Status:* slack2 grid minimax (baseline grid error 3.1e-8); not an interval certificate.

**Claim 4 (lifting a Poisson-range product moment gives exactly 5/6).**
- *Hypothesis UM2:* for every fixed ϑ, ε > 0, Σ_{0<q_v≤H} |M_v(φ;D1)M_v(φ;D2)|² ≪ H(HD1D2)^ε whenever
  H ≥ (D1D2)^{1+ϑ}, uniformly for common profiles. Here φ_v = νψ_v^{±1}, and the two factors may carry different
  fixed ν.
- *Statement:* Σ^{(6)}_{q_u≍U} |M_u(U^{r1})M_u(U^{r2})|² ≪ U^{max(1,(1+5(r1+r2))/6)+ε}.
- *Proof sketch:* apply (B0.7) to each factor with the same a. The zero extension φ_{ua⁶} = φ_u·1_{(·,a)=1} holds for
  any fixed ν. The coefficient mass is d(a)². Take a two-scale Sobolev supremum (as in B0.8 / B3.44). Average over
  ≫ P ideals a with H = UP⁶ = max(2U, (D1D2)^{1+ϑ}), and use the injectivity of (u,a) ↦ ua⁶. This gives
  U^{1/6}H^{5/6}.
- *θ consequence:* combined with Claim 2(i), UM2 and MI2 together, if proven, give **zero** gain (an exact tie).

**Claim 5 (no-go for Hölder, pointwise and slot routes).**
- *Statement (i):* every count line obtained from products or Hölder combinations of UM, MI, MP, MP0 and (B0.9),
  together with the pointwise bounds (B0.15)–(B0.16), is at least a convex combination of the paper's lines at the
  same witness data. Proof sketch: §2(b); P4 plus π ≥ w on the class.
- *Statement (ii):* inserting the rough part of a Möbius factor into MI as prime slots (cost 2 per unit length) gives
  the line A_I(r2) − (δ−2q)r1, which equals A_I(r2) at p = 1/2.
- *θ consequence:* 0.

**Claim 6 (two characters dominated by one).**
- *Hypotheses:* the refined bins (δ_N, δ_D) of §0.
- *Statement:*
  - Binding requires δ_N = δ_D = δ0 and p = 1/2.
  - On diagonal classes, the two-character joint system equals the single-character two-length system with the same e.
  - Off the diagonal, the joint line is weaker (−δ_N r1 − δ_D r2), while the paper's slack is only about
    0.0345|δ_D−δ_N|.
  - Numerically, the two-character prices are smaller than the single-character prices (`minimax3.py`).
- *Proof location:* §0 and §1.5.
- *θ consequence:* B4 as posed can never beat its single-character analogue.

**Claim 7 (obstruction; explanatory, not an impossibility theorem).** The Poisson dual of μφ1*μφ2 in the B3.31
computation has coefficient a_0(z)g(z), with g(p) = ν1(p)q_p^{−it1} + ν2(p)q_p^{−it2}. That lies outside (B3.6).
Expanding g as a label sum puts a coherent label inside |·|², and Cauchy–Schwarz over it costs ≈ D1. Therefore UM2
and MI2 need a label-coherent (bilinear) version of RE. Location: §2(a).
