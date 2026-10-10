# QRH slack analysis — findings log

Target: beat cycle25 quartic, θ = 683505193/781250000 ≈ 0.87488664704 (b0 = 11/12 − ℓ0/4).
Source: proofs/cycle25-quartic-20261010/manuscript/paper.md (B0.3–B0.7, B2.6).

## Model (slack2.py)
- Principal probe (B0.27 + B0.36): b0 = x/2 + (y−x)/12 + (1+3ℓ−2M)_+/8 + 5/6 − x/3 − ℓ/6, M = x+y.
  With M = 1−ℓ this is 11/12 − ℓ/4. Bound improves iff ℓ grows.
- Constraint: nonprincipal-row exponent (B0.55) F(d;δ,p,κ) < β for all bins δ∈(0,κ], p∈[0,1/2],
  d∈[d_min,h], κ = 2β−1, h = 1+ℓ−x. Binding at d=h, p=1/2, δ0≈0.4693. z0 cancels at d=h.
- Optimizer reproduces b0 = 0.87488662 (paper 0.87488665; diff = grid error). Geometry x is
  already minimax-optimal (paper's B0.48 kills the h-coefficient at the binding point).

## Sensitivities (db0 per unit change of the constant)
| Constant (location) | per-unit | note |
|---|---|---|
| Hölder 1/4 in numerator 1/4+3R/4 (B0.41) | +0.075 | biggest lever; needs higher moment than 4th of plain sums |
| M+ℓ≤1 relaxation, *free* | −0.063 | but see below: not free |
| V0 = 5/6 in long inverse envelope (B0.40/B0.9) | +0.0074 | from sixth-power averaging exponent (1+5r)/6 |
| amplitude cap p≤1/2 (B0.37) | +0.0075 | |
| probe (y−x)/12 (B0.22 exceptional P_a^{1/6}) | ~−0.0005 per unit of divisor | small |
| plain capacity 1/(6κ), damping coeff | ≈0 | negligible |
| z0, supply 4/21, other endpoints | 0 | non-binding |

## Dead end 1: relaxing M+2z ≤ N_* (B2.44)
Hypothesis only forces the row branch max(H, v+ℓ_b)=H in (B2.38). Past it, the column branch
gives energy ≈ T_d = 2M+2z−N_* = M + s (s = M+ℓ−1), i.e. probe cost s/2.
Break-even energy cost is < 0.25·s (cost 0 → gain 1.25e-3 at s=.02; cost .25s → gain 0).
With the honest cost s, optimum stays at s=0. Would need a column-branch estimate saving ≥ ~0.8s
over (B2.38) — not plausible as a tweak.

## Next to examine
- Higher-moment numerator (B0.41): sixth moment of reflected plain sums needs 3(1−m) ≤ 1-type length;
  at binding m=y/h≈0.58 the dual length ≈0.42 → 3·0.42>1. Check mixed/partial Hölder or using
  marks to shorten.
- Long envelope V0: any improvement to (B0.9) sixth-power-free inverse moment.

## Tick 1 results (2026-10-10)
- **Saddle structure.** Binding point (δ0, p=1/2, d=h0) is stationary in both δ (max) and h (min).
  So improving estimates *away* from δ0 gives zero first-order gain. Verified: using
  min(Hölder, pointwise (B0.15)) for the numerator — pointwise wins for δ>0.49 — gain 7e-14 (none).
  Hölder/pointwise interpolation is linear in the exponent ⇒ endpoints only.
- **Conversion factor:** lowering F at the binding point by Δ lowers θ by ≈0.23·Δ.
  Last record step was ~3e-5 in θ ⇒ need ΔF ≳ 1.3e-4 at (δ0≈.4693, p=.5, κ0≈.7498).
- **Row count R ingredients** (Rgen.py reproduces paper R exactly):
  MI capacity (1−r)/2 → a_I(1−r): −0.008 per unit a_I; MP capacity: −0.0004; V0: +0.0075.
  All weak. R is not where the slack is.
- **Remaining first-order lever:** numerator exponent at binding, d(1/4+3R/4) − q·dz_A/2
  (≈0.5467 in Z-units). Needs an estimate for Σ_{u in bin}|S_u(reflected length ≈ U^{.42})|
  beating fourth-moment Hölder on a sparse row set of size U^{.59}.

## Tick 2 results
### Probe Gram exceptional term (B0.22) — largest lever found so far
Gram bound Σ|A_m|² ≪ (Q/Y)(1 + P_a^{c} + P_a²/Y), paper c = 1/6, P_a = Y/X = Z^{y−x} ≈ Z^{.087}.
Probe exponent = x/2 + c(y−x)/2 + ... (paper: (y−x)/12). Re-optimizing geometry (gram.py):

| c | b0 | gain vs record |
|---|---|---|
| 1/6 | 0.874886616 | 0 |
| 1/8 | 0.872448789 | 2.4e-3 |
| 1/12 | 0.869588586 | 5.3e-3 |
| 0 | 0.865152069 | 9.7e-3 |

(At c≤1/8 the constraint y−2x+ℓ<0 binds, x→1/3.) Even c = 1/6 − 0.01 would beat the last
record step (~3e-5) by ~20×.

**Assessment: probably a genuine main term, not a lossy step.** The exceptional frequencies
k = e·r⁶ (B0.22 proof) make the column character χ_{n1}(k) = χ_{n1}(e) on (n1,r)=1, so the column
sums give full mass with no r-dependent phase (L_C table B1.17 depends on k only via p|k). Summing
over r ≤ (RP_a/C0)^{1/6} is then a positive count. Dominant dyad is C0=D0=1 with e S-supported; since
ξ is locally quadratic, e with valuation 2 at each S-prime makes the twisted column character
principal, so the twist cannot cancel it either. The paper itself notes the bound fails for
arbitrary column coefficients.
**Ways around (need new ideas, not checked):** (i) choose column weights W_1/twist so that every
exceptional e gives a mean-zero column sum (needs a different ξ compatible with B0.4's phase
requirements ξ^6 principal on units); (ii) weighted Cauchy–Schwarz in (B0.27), since the excess Gram
energy lives on structured m and T(m) may be small there; (iii) avoid P_a>1 — the optimizer already
balances y−x>0 against −y/2 in F, so no free gain.

## Tick 3 results (literature)
- Dunn–Radziwiłł (columbia.edu/~goldfeld/Dunn.pdf): Heath-Brown's cubic large sieve is sharp under GRH,
  via dispersion of cubic Gauss sums. Supports "exceptional sixth-power term is a genuine main term":
  no generic large-sieve improvement will remove it.
- Sextic large sieve over Z[ω] has shape (M+N+(MN)^{2/3})(MN)^ε (cited in arXiv:2201.01885, Thm 1.3 ref).
  Nothing sharper found for the sextic Eisenstein case.
- **Only structural escape left:** the dominant exceptional dyad (C0=D0=1) needs the twisted column
  character n ↦ χ_n(e)·ξ̄(n)χ_n(b*)·R(n,n2) to be principal for some S-supported e. The paper notes
  "fixed finite phases need not have order dividing six". If the probe twist could be chosen with
  local order outside the group {n ↦ χ_n(e)} (e.g. order-4 or 1+3O unit component) the C0=1 exceptional
  term would vanish; remaining exceptional dyads have C0>1 with weight C0^{-1/6}. Risk: ξ's quadratic/
  cubic choice is used to make ξ^6 principal on units and to keep the principal signal (B0.35) nonzero.
  OPEN QUESTION for a human expert: does B0.3/B0.4 survive a twist outside the sextic group? If yes,
  re-derive B0.22 exceptional exponent; payoff up to ~1e-2 in θ (gram.py).

## Tick 4 — twist escape is blocked (dead end 2)
- ξ̄(u) weights the dual rows in (B0.28). Writing frequency = u·a⁶ factors ζ_K^S(6z) only if ξ(a)⁶=1;
  that ζ pole at z=1/6 is the principal signal (B0.35, c_S ∝ M(1/6)). So ξ must have order | 6.
- Locally at p∈S (Q≡1 mod 6), characters of order | 6 are exactly powers of the local sextic character,
  and n ↦ χ_n(p) (via reciprocity) generates that subgroup. Hence for any admissible ξ some S-supported
  e (valuations 0..5) makes n ↦ χ_n(e b*)ξ̄(n) principal ⇒ C0=1 exceptional Gram term is forced.
- Conclusion: the P_a^{1/6} term is structurally tied to the main term. Removing it needs a different
  probe design, not a parameter change.

## Summary after 4 ticks
Within the cycle25 framework, every first-order lever is either (a) tight at the saddle point,
(b) a genuine main term (Gram exceptional, sharp per Dunn–Radziwiłł-type evidence), or
(c) worth ≤0.01 per unit of an estimate whose improvement would be a new theorem (numerator
fourth-moment Hölder at a U^{.59}-sparse row set; MI capacity; long envelope V0).
No parameter-level improvement to θ exists. Candidate research directions, ranked by payoff:
1. New probe avoiding Cauchy–Schwarz in (B0.27) / the sixth-power Gram main term (up to ~1e-2).
2. Large-values estimate for reflected plain sums on sparse sets beating 4th-moment Hölder
   (0.23 θ-gain per unit F at δ0≈.469, p=1/2).
3. Higher-moment MP (needs sextic large sieve beyond length U) — blocked by (MN)^{2/3} term.

## Research session (after "do the research")
### Tried: joint large-values numerator bound using mixed MP (joint.py, branches.py)
(B3.40) proves MP for two *independent* plain sums S(Z^{n1})S(Z^{n2}) (zero-slot version for all
lengths, unconditional). So rows with a plain witness |S(U^m)|²≫U^{δm} and numerator level U^σ obey
the extra count 1−2σ−δm−2q·z_mix. Result at δ0: gain ≈1e-5 (grid noise); +2e-3 at δ=.55, +1e-2 at
δ=.70 — irrelevant because of the saddle. Reason: the bin max is attained jointly by the *long*
rows (r≈t∈(1,3/2], m≈0) whose only witness is the inverse one; no mixed moment exists for
M_u(U^r)·S_num with r>1 (MI needs r+2·len≤1; a B0.9-type product bound helps only if 2σ>5m'/6≈.35).
Also tested y<x geometry (P_a<1): worse (ylessx.py).

### Key structural observation
- At the binding point the Hölder numerator bound and the pointwise bound (B0.15) agree to 1e-4,
  and bin rows plausibly saturate the pointwise bound (a zero near a forces |S|~N^{δ/2}).
  ⇒ the numerator step is essentially sharp *given R*.
- R is a zero-density bound for the sextic family: #rows with a zero of real part ≥ a=(1+δ)/2 is ≤U^R.
  At binding: σ=a=.7347, R=.5903 ⇒ N(σ) ≪ F^{2.225(1−σ)} (density hypothesis: exponent 2).
  For comparison Guth–Maynard for ζ near σ=.73 gives ≈2.25: the paper's estimate is already strong.
- The whole excess R−(1−δ) = δCV/(2(V+δC)) comes from the *long envelope* slope V=5/6−δ, i.e. the
  sixth-power-free inverse mean value (B0.9):
      Σ_{q_u≍U, 6-free} |M_u(U^r)|² ≪ U^{(1+5r)/6+ε},   1<r≤3/2.
  If this held with exponent 1+s(r−1), s<5/6, V shrinks; s≤δ0≈.47 gives R=1−δ.

### Payoff of improving R (density.py): R = 1−δ + λ·excess
| λ | θ | gain |
|---|---|---|
| 1 (paper) | 0.874886616 | 0 |
| 0.9 | 0.874066724 | 8.2e-4 |
| 0.75 | 0.872832011 | 2.1e-3 |
| 0.5 | 0.870770421 | 4.1e-3 |
| 0 (density hyp.) | 0.869047620 | 5.8e-3 |
Even λ=0.99 (1% of the excess) beats the last record step (~3e-5) by ~3×.

### Precise open problem (the research target)
Improve the 5/6 in the sixth-power-free mean value of sextic Möbius polynomials beyond the conductor
range: for D=U^r, 1<r≤3/2,
    Σ_{q_u≍U, u sixth-power-free} |D^{-1/2} Σ_n μ(n)ν(n)χ_n(u) W(q_n/D)|² ≪ U^{1+s(r−1)+ε}, s<5/6,
or any zero-density estimate for the family {νψ_u} at σ≈.735 with exponent <2.225.
Notes: general-coefficient large sieves cannot do it (Heath-Brown/sextic (MN)^{2/3} term gives
exponent 2(1+r)/3 >(1+5r)/6; even the sharp U+D gives slope 1). Halász–Montgomery fails: the
off-diagonal n↦χ_n(u/v) has conductor ≈U² so Pólya–Vinogradov gives U > large value U^{δr}.
Must use Möbius/zero structure (Guth–Maynard-style large values for this family?) — genuinely new math.

## Attack session 2 ("keep attacking")
### Product-witness mean value (produm.py) — no gain
Every bin row has |M_u(U^r)S_u(U^m)|² ≫ U^{δ(r+m)}, r+m ≥ t0. Even granting a hypothetical UM for
products M·S (count ≤ e(r+m)−δ(r+m), e(x)=max(1,(1+5x)/6)), R: 0.59033 → 0.59026 (noise).
### R is a triple point
Fine grid: worst row (r≈.746, m≈.418, x=r+m≈1.164) satisfies A_I(r) = A_P(m) = E(x)−δx, independent
of t0 ∈ [1,1.15]. So R is set by three inequalities: MI line, MP line, sixth-power long bound E.
Conversion: ΔR = −1e-3 at δ0 ⇒ Δθ ≈ −1.4e-4.
### Literature calibration
Guth–Maynard (Annals 2026): N(σ,T) ≤ T^{30(1−σ)/13}; region formula 15/(3+5σ) ≈ 2.247 at σ=.735.
Chen arXiv:2507.08296 (GM for Dirichlet characters): exponent 7/3. Paper's sextic-family exponent at
σ=.735 is already 2.225 (prime-slot amplification), so importing GM/Chen does not beat it directly.
### The missing lemma (★)
Canonical estimate (B3.7) + UM give Σ_{q_k≤H}|M_k(D)|² ≪ max(H,D); the sixth-power trick turns this
into E(x)=(1+5x)/6. Needed:
  (★) Σ_{q_k≤H} |D^{-1/2}Σ μ(n)ν(n)χ_n(k)W(q_n/D)|² ≪ H^{1−λ+ε} D^{λ},  H ≤ D ≤ H^{1.2},  λ < 5/6.
Then E slope = min(λ,5/6). Payoff (lambda.py):
| λ | θ | gain |
|---|---|---|
| 5/6 | 0.874886616 | 0 |
| 0.82 | 0.874786254 | 1.0e-4 |
| 0.80 | 0.874627543 | 2.6e-4 |
| 0.75 | 0.874181591 | 7.1e-4 |
| 0.60 | 0.872523458 | 2.4e-3 |
| 0.50 | 0.871446061 | 3.4e-3 |
Obstructions/notes:
- The D-term is sharp only for adversarial coefficients; under GRH truth is ≈H. Sixth-power rows
  k=a⁶ contribute ≈H^{1/6}·D^{κ} (κ=2β*−1) — compatible with λ ≥ ~0.33, so no obstruction at λ≈.8.
- UM (B3.39) itself, applied to rows a⁶, implies Σ_{(n,a)=1}μνW(n/D) ≪ D^{11/12} on average ⇒ a
  zero-free strip Re s>11/12: UM is the quasi-RH engine (matches OpenAI's reported 11/12). (★) is a
  genuine strengthening of that engine beyond the Poisson range D<H, where the dual length D²/H > D.
- Heath-Brown identity route (M(D) → M(D1)M(D2)S(N), D_i ≤ D^{1/2}) needs mean values for M·M·S;
  only S·S (MP0, all lengths) and single M (UM) exist. Missing: an MI-type bound with a plain factor.
- MI's 2z (r+2z≤1) comes from RE row-branch energy M+z_a at T_d≈0, which is attained: no easy slack.

## Attack on (★) — three routes, all fail (with reasons)
1. **Poisson / duality alone.** Σ_{q_k≤H}|M_k(D)|² → (Poisson in k, conductor ≈D²) → same-shape sum with
   rows R=D²/H and Gauss-twisted Möbius coefficients; applying it twice is the identity. UM on the dual
   (R ≥ D^{1+θ} iff D ≥ H^{1+θ}) returns exactly D. Barrier max(H,D) is intrinsic to this method;
   beating it needs arithmetic (zero) input.
2. **Bootstrap via zero counts (bootstrap.py).** Mellin shift: rows with no zero right of σ have
   |M_k(D)|²≪D^{2σ−1}; so Σ ≲ U^{E'(x)}, E'(x)=max_δ[R(δ)+xδ] with the framework's own R(δ).
   Result: E'>E for x≤1.35 (x*=1.164: E'=1.160 (p=.5) / 1.172 (all p) vs E=1.137); max at δ≈.6–.7.
   Rows with zeros nearer β* dominate; framework is self-consistent, no self-improvement.
   (Also needs control of 1/L near zeros — B0.14 handles it via zero-free boxes.)
3. **Heath-Brown identity.** μ_{~D} = Σ_j ± M_{≤D^{1/J}}^{*j} * 1^{*(j−1)}. Tools: MP0 (Σ|S(n1)S(n2)|²≪U,
   ANY lengths, unconditional) + pointwise |M(D_i)|²≪D_i^κ, |S(N)|²≪N^{κ·min(n,1−n)_+} from β*.
   Terms with Möbius length d_M cost U^{1+κ d_M}; beating 5/6 at x*=1.164 needs d_M ≤ .18, but d_M up
   to x occurs. Bilinear grouping needs a general sextic large sieve: known (M+N+(MN)^{2/3}) gives
   ≥U^{1.44}; even the conjectural sharp (M+N) gives Type II ≪U^{x} > U^{(1+5x)/6}.
   ⇒ Decomposition discards the Möbius structure UM exploits (UM beats any general-coefficient sieve).

### Verdict
(★) is not reachable by Poisson, bootstrap, or combinatorial decomposition with current tools.
It would need a new idea that keeps Möbius structure past the Poisson range D>H — effectively a new
zero-density theorem for the sextic family at σ≈.73–.84. Genuine open research.

## New approach: MP slot cost 6κ (numerator capacity) — priced, then traced to a structural cause
Numerator Hölder regains q·z/2 per unit slot capacity; z_A=(1−2m')/(6κ)≈.035 while supply L=ℓ/d≈.21.
Replacing MP's 6κ by cκ everywhere (numerator z_A and A_P line), mpcost.py:
| c | θ | gain |
|---|---|---|
| 6 (paper) | 0.874886616 | 0 |
| 5.5 | 0.874783146 | 1.0e-4 |
| 5 | 0.874658236 | 2.3e-4 |
| 4 | 0.874309180 | 5.8e-4 |
| 3 | 0.873697462 | 1.2e-3 |
| 2 | 0.872226885 | 2.7e-3 |
Origin: (B3.70)–(B3.72): child affine defect grows ≤6(w+ℓ); removing a slot of length d lowers it by
6κd at pointwise cost κd (B3.47); closure needs cost ≤ child margin ⇒ c=6κ (paper: "cost-to-capacity
ratio 1/6"). The 6 in (B3.70) is the amplifier of B3.7: rows h→hp⁶, row growth 6ℓ for density gain ℓ.
Sixth powers are forced: G(u,hm)=χ̄_u(m)G(u,h), invariance for all u needs m a sixth power; a
non-invariant amplifier (h→hp) makes the u-coefficient row-dependent, outside the MP family.
The κ is the global zero-free cost; children are dual families, so the bin's sharper δ cannot be used.
Possible escape routes (not worked out): (a) handle slotted parts MI-style (MI's slots cost 2, no
amplifier, because its dual lands in the canonical long-polynomial regime); (b) an amplifier with better
row-growth/gain ratio for the sextic Gauss family. Both need re-proving the MP induction.

## Overall status after all attacks
Every lever reduces to a sextic-structural constant: Gram P_a^{1/6} (sixth-power frequencies),
long slope 5/6 (sixth-power averaging), MP slot cost 6κ (sixth-power amplifier). The framework is
locally optimal; progress needs a new mean-value theorem (★) or a new amplifier/slot mechanism.

## Parallel solver round (subagents; graders run only on positive claims)
- Solver B ((★), report solverB.md): FAILED, no claim. Adds: all-rows vs 6-free forms of (★) are
  equivalent with same λ (5/6 = lifted max(H,D) at H=D); sixth-power rows not the obstruction (pointwise
  handles them for λ≥1/3); exact mollifier identity M_k(D)=Σ F_k(B)M_k(D/B), X<B≤XH^{1+ε}, but no
  Möbius×plain mixed moment exists and even a hypothetical sharp one optimizes to λ=1; MI marks in the
  lift gain only if q>5/6 (q≤δ/2<3/8); Halász–Montgomery gives U^{r(1−δ)}, worse. Independently
  confirms: UM on rows a⁶ ⇒ ζ≠0 for Re s>11/12 — UM/RE carry the quasi-RH strength.
- Solver A (MP slot cost, solverA.md, lp.py): FAILED, no claim; argues c=6 is optimal *within the method*.
  Generic child (no common support) saturates (B3.70) and (B3.64) with equality (A'=A, M'=M−6ℓ, Δ=ℓ),
  so tighter accounting gains nothing. Invariance amplifier needs m = unit·sixth power ⇒ ratio 6 is a
  floor. Node constraint A_eff ≤ M+L from positivity (B3.54)/diagonal (B3.61),(B3.80). LP over
  reflection schedules reproduces z_max=(1−2n)/(6κ) exactly. MI-style route fails: dual of a plain sum
  is μ·Ḡ·ā₀ (inverse-theta), outside canonical estimate/RE. Needs: non-Cauchy bilinear bound for
  (B3.52) at A>M, a reflection identity for slot factors, or sparse-row/8th-moment bound for S.
  Self-flagged unverified: Kummer step with S-units/masks; method-vs-truth lower bound.
- Solver C (probe exceptional term, solverC.md): no improvement. Prop. O: per ray class/v-piece the
  Gram exceptional mass and the principal signal are driven by the SAME column functional at the same
  Mellin point (Lemma 1: |Col(e0)|=|Col'(1)|), so any Cauchy–Schwarz-based bound has bound/signal ≳
  P_a^{1/12}; column weights, other ξ, smooth profiles, non-S twists χ̄_a(u0) all fail (u0=1 optimal).
  Confirms (i),(ii) (closed the gaps at 2 and λ). **Refutes my (iii)**: m→mr⁶ invariance holds for all
  of A, so it says nothing about alignment with T; heuristically (zeros ≤.85) the true |I|≪Z^{.169} <
  X^{1/2}=Z^{.187}, so the P_a^{1/12} loss is a METHOD loss, not the truth.
  **New concrete target (best found):** an a-twisted reflected energy
      E_a = Σ_a |Σ_m ξ(m)χ_a(m)T(m)|² ≪ Q(Q+Y)Z^ε
  (bilinear, avoids CS) would give c=0 ⇒ θ gain up to ~9.7e-3 (gram.py). Sixth-power diagonal of E_a is
  harmless; generic large sieve gives exponent .2010 (insufficient) ⇒ needs a twisted version of B2 RE.
  Uncertain (per solver): exceptional asymptotic at main-term level only; Mellin-point match by scaling.

## Independent numeric audit of the record's endpoint (audit.py, audit2.py) — ALL PASS
High-precision (60 digits) + exact rationals where possible; grids, not interval arithmetic.
- T1.1 root interval & p'>0 on [1/6,1/5]; ℓ0 = 0.1671200785068101428882107.
- T1.5: θ−b0 = 3.587e-14 ∈ (0, 4e-14) (also via exact root interval); θ < 10499/12000 < ProofCouncil pin.
- B0.45 values match paper to 1e-14; ℓ0<x0<y0<2x0, y0<1/2<h0<2y0; κ0 = 2b0−1 exactly.
- F4(δ0,1/2,κ0)=b0 exactly; both B0.48 equations hold (residual <1e-58); R(δ0,1/2)=5/9+δ0/(18κ0).
- Grid max of F4 over δ∈(0,κ0], p∈[0,1/2] is b0 − 3.3e-11 (attained at binding point).
- B0.49 identity verified to 3e-58 on 200×51 grid; R-substitution exact; J∈[81.5,269.5]⊂(69,270];
  A1, A2 minorants hold; K*=179.168; 768T=117.638.
- T1.3: max ∂C/∂κ = .0415 < 1/18; max ∂F4/∂κ = .0085 < 1/10.
- Margins: other endpoint .8105 < b0−.06; floor .8599 < b0−.01; B0.56 = −.0859 < −1/20;
  ℓ0/h0 = .2105 > 1/5; y0 < 4ℓ0; 0 < h0−4ℓ0 < .13, h0 < .8.
Conclusion: the arithmetic layer of the record is correct and tight (no free rounding slack beyond 3.6e-14).
- Grader of Solver C (graderC.md): Lemma 1 correct w/ fixable gaps (b* primes must be primary & in S;
  φ(n)/q_n density factor missing but cancels). Lemma 2 as LOWER bound UNPROVEN: C>1, d>1 terms only
  controlled by averaged upper bounds with Z^ε and signed L_C; step S56 turns averaged dyadic bound into
  per-modulus. Lemma 3: m-side Mellin point verified, column side by scaling only. Prop. O "every CS split":
  UNPROVEN — covers only radial row weights, per-(σ,v) splits, S-class×radial column weights; not
  non-radial/arithmetic row weights, general column weights, CS without v-separation, or dual CS over
  columns. Honest content: B0.22 sharp for these coefficients, conditional on Lemma 2.
  Refutation of my (iii): CORRECT. Z^{.169} heuristic: INVALID (B0.51 needs β*>b0; at β=.85 margin
  becomes +.02) ⇒ only Z^{.1938} > X^{1/2}=Z^{.187}; u≠1 rows unhandled.
  ⇒ Net: neither "CS loss is forced" nor "true probe < CS bound" is established. E_a target stays open,
  with weaker heuristic support.
- Solver D (twisted energy E_a, solverD.md, solverD/expo.py): NOT PROVED; best available c = 1/3 (worse
  than paper's 1/6). Established (some steps sketched, U1–U3): a-route CS is equality on the u=1 piece
  (no built-in loss; target consistent with β≤.8677); twist χ̄_a(m)=χ_{a⁵}(m) shifts the theta index,
  not absorbable as a B2 label (V=y column branch 1.126>M) and m-Poisson is circular; reflected T on
  squarefree rows = fixed char·ᾱ(R)²γ₁(R)·short quadratic sum; Poisson in a cancels that phase and
  kills the P_a^{1/6} mechanism (diagonal ≤Y·Q^{7/6}; sixth-power frequencies have no main term due to ᾱ²).
  Precise obstruction: a sharp sextic large sieve for ψ(m)ᾱ(m)²R(m) with dual length K=Q²/Y>Q, i.e.
  non-correlation of T with sextic Patterson theta coefficients (not explicitly known for n=6).
  Probe exponents: HB-type sieve c=1/3 → .2010; classical c=1 → .2299; paper c=1/6 → .1938; target c=0 → .1866.

## Round summary
All four targets (MP slot cost, (★), Gram exceptional term, twisted energy) fail with precise
obstructions. Two reduce to the same deep input: sextic Patterson-theta / large-sieve non-correlation.
Record arithmetic independently verified. No improvement to θ.

## Numerical experiment: sextic large sieve (SLS_T) — evidence it is TRUE (exp/sieve*.py)
Setup: Z[ω], split-prime proxy; moduli m = sf products of split primes, N(m)∈[Q,2Q]; arguments k primary,
N(k)≤K, (N(k),6)=1; ratio = Σ_k|Σ_m b_mχ_m(k)|² / ((#k+#m)‖b‖²) (sharp sieve ⇒ ratio ≲1 flat).
Symbols validated by supplementary laws & multiplicativity; Gauss sums |g̃|=1.
- POSITIVE CONTROL works: cubic sieve with cubic Gauss-sum vector shows Patterson excess 1.13–1.56 vs
  random ~0.5–0.94; excess grows by ≈2^{1/3} when Q doubles (Q=1500→3000) = predicted (MN)^{2/3} loss.
  (Sextic Gauss vector: small excess ~8%.)
- FAITHFUL THETA PROXY b_m = ᾱ(m)²Σ_{N(μ)≤D} g̃₃(μ)α(μ)N(μ)^{-1/2}(μ/m)₂ (structure of Solver D's
  𝓡, ρ_μ from B2.8 cubic theta coefficients), D≈Q^{0.8}: ratios = random within 1–2% at Q=1500 and 3000,
  sextic and cubic sieves, both conj conventions. No loss.
- MECHANISM: angular twist kills the loss. ᾱ²·Gauss, α²·Gauss, ᾱ·Gauss all drop to random level
  (vs 1.2–1.56 untwisted), both scales. Consistent with: HB's (MN)^{2/3} term comes from the Patterson
  pole of the Kubota series; nontrivial infinity type removes the pole (intertwining L(ᾱ^{12})≠1).
- Analytic route: Poisson in k turns (SLS_T) into a dual sieve over h ≤ Y=Q²/K (Q^{1/2}<Y<Q) with
  vector b_m·g̃₆(m); HB there loses (QY)^{2/3}/Q = Z^{.029} = P_a^{1/3} (Solver D's c=1/3). Removing
  the pole-driven term would give the sharp bound ⇒ c=0 ⇒ θ ≈ 0.865152 (gain 9.7e-3), conditional.
Caveats: split primes only (inert primes, S-units, branches, ϑ twist, smooth V♯ omitted); small Q;
sharp cutoffs; evidence, not proof.
**Conjecture (ATS, angular-twisted sextic sieve):** for b_m = η_∞(m)·c_m with η_∞ a nontrivial angular
Hecke character and c_m of Gauss/theta type, Σ_{N(h)≤H}|Σ_{N(m)~Q} b_mχ_m(h)|² ≪ (H+Q)‖b‖²Z^ε,
in the range Q^{1/2} ≤ H ≤ Q.
- Literature (literatureE.md, litE/): no published ATS. Dunn–Radziwiłł (Annals 2024, arXiv:2109.07463)
  Prop 5.2: Σ g̃(c)(r/c)₃(c/|c|)^ℓ N(c)^{-s} ENTIRE for ℓ≠0, simple pole at 5/6 for ℓ=0 — confirms the
  mechanism for n=3. Thm 1.4: (AB)^{2/3} optimal under GRH, extremal β_b = g̃(b)W. Props 9.1–9.2/(1.13):
  sieve minus explicit correction ≪ (AB)^{2/3−δ}; correction vanishes for angular-twisted vectors but
  Prop 9.2 hypotheses exclude them (to examine). de Faveri–Dunn–Hoffstein arXiv:2607.07911 conj (1.8):
  Ξ_n = A+B+A^{1−1/n}B^{2/n}+A^{2/n}B^{1−1/n}; untwisted n=6 exceeds Q iff H≥Q^{1/2} ⇒ ATS must use twist.
  Blomer–Goldmakher–Louvel IMRN 2014 Thm 1.3 = source of sextic (M+N+(MN)^{2/3}). DDHL arXiv:2306.11875
  Prop 4.3: quartic twisted Kubota series entire for ℓ≠0. n=6, ℓ≠0 pole structure unverified.
- Solver E (ATS proof, solverE.md, solverE/psr.py, expo.py): NOT PROVED for our vector; unconditional c=1/6.
  DR Prop 9.2: the only essential hypothesis is Axiom (3.2) with η≤1/4 — pointwise √ cancellation for every
  dual frequency k ≤ B²/A (|β|≤1, roughness, ℓ≠0 clause inessential/unused, flagged U-E1). Our 𝓡 sums
  ~D=Q^{0.8} quadratic chars with ℓ¹-mass D^{1/2}: GRH gives η≈.40 ⇒ E_a exponent 2.33 ≫ 1.666.
  Theorem E1 (sketch, GRH-conditional): sharp ATS for pure twisted Gauss vectors ψᾱ^ℓμ²g̃₆W, ℓ≠0 — the twist
  does remove the Patterson obstruction; does NOT cover 𝓡.
  Conditional: Hypothesis PSR (pointwise √Q bound for ψᾱ²𝓡 against all sextic chars up to K) ⇒ c=0 ⇒
  θ≈0.865152. Numerically PSR holds (max over k ~ log #k, flat in D=60..400, Q=1500), but PSR is stronger
  than GRH (incoherence across ~D quadratic twists); mean-square version circular under Poisson.
  HB recursion (one Poisson in m) gives c=1/6 = de Faveri–Dunn–Hoffstein conjectured Ξ₆ at this geometry ⇒
  1/6 conjecturally optimal for GENERIC vectors; beating it needs non-correlation of T with sextic theta
  coefficients τ₆(h'm⁵), h'≤P_a (cubic analogue works via Patterson's explicit τ₃; τ₆ unknown).
  Twisted sextic Kubota series not the bottleneck (lose Z^{.29} pointwise even with Lindelöf; entireness n=6 = inference U-E4).

## FINAL STATUS
Unconditional θ unchanged. Best conditional: PSR (numerically supported, stronger than GRH) ⇒ θ≈0.865152.
The precise missing input: non-correlation of the cubic-theta sum T with sextic theta coefficients τ₆.

## Further thinking (no agents) — angles examined
1. Angular-twist the probe's columns to kill the Gram exceptional term: kills the principal signal too
   (u=1 row then has L(w,ᾱ^j) with no pole at w=1) — consistent with Solver C Lemma 1. Dead.
2. Linear combinations of probes / per-bin probe choice: one probe must control all bins; no mixing. Dead.
3. Prime-sum witnesses (explicit-formula detector) counted by MI with r=0: capacity 2z≤1 (and 8z≤3)
   ⇒ count ≤ U^{1−δ·3/8}≈U^{.82} ≫ U^R. Dead.
4. Cube-shift trick u→ua³ (rows ≤ U·P³, slope 2/3 instead of 5/6): shifted rows carry different polynomials
   (twist by (a/n)₂); no invariance ⇒ cannot transfer |M_u|². Dead.
5. Functional equation for 1/L beyond conductor: 1/L(1−s) sits deep in the strip, no expansion. Dead.
6. Structure of (SLS_T), worked out: Σ_k|B(k)|², B(k)=Σ_μρ_μF(kμ³), F(j)=Σ_m ψᾱ²(m)χ_m(j)w.
   - m=m' part: K·Σ_m|Σ_μρ_μ(μ/m)₂|² ≪ K(Q+D)Σ|ρ|² by Heath-Brown's SHARP quadratic large sieve. Fine.
   - m≠m' part: Poisson in k (dual length Y=Q²/K<Q) ⇒ vectors Σ_μρ_μ ψᾱ²(m)g̃₆(m,μ³): per-μ these are pure
     twisted Gauss vectors (Solver E Thm E1 applies, GRH-cond.), but summing over μ loses ℓ¹-mass D^{1/2};
     orthogonality across μ would come from Σ_m(μμ'/m)₂ (quadratic sieve) but extracting it requires
     another Poisson ⇒ circular. Both sides (m: ᾱ², μ: α(μ)) carry angular twists ⇒ no theta main terms,
     matching numerics. The gap is a genuinely trilinear (m, μ, h) estimate.
7. Conditional results on GRH are vacuous here (GRH for the family gives θ=1/2 outright); PSR is the only
   meaningful conditional input.
8. Quadratic-family redesign (HB quadratic sieve is sharp; square trick slope 1/2): but Gram exceptional
   becomes P^{1/2} and the principal pole moves to ζ(2z) at 1/2 — likely why OpenAI chose sextic.
   Cannot evaluate without rebuilding the theory.

## Theta-residue identification experiment (exp/theta_id.py, theta_id2.py)
Idea: dangerous probe piece is Σ_m ρ₆(m)γ₁(m)ᾱ²(m)ψ(m)(μ/m)₂ (μ small). If ρ₆γ₁ were Hecke-character-like at
squarefree m, each μ-term = smooth sum of an angular Hecke char over a range > conductor ⇒ rapid decay ⇒ c=0.
Measure ρ₆ via A_m = Σ_{a~Y} g̃(χ_a)·conj χ_a(m) (Kubota sum at index m; Patterson residue ∝ ρ(m)).
- VALIDATION (cubic): A_m ∝ conj g̃₃(m), corr .75–.86 vs baseline .06–.08; energy 2–2.7× random. Recovers
  Patterson's τ₃ ∝ conj g̃₃.
- SEXTIC: residue present (energy +30% at Y=20000,Q=1000, 'lam' convention = paper's e(x/(λa))), but NOT a
  single Gauss sum or pairwise product; only identifiable component ≈ constant in m (corr .195 vs .10);
  top-3 candidates explain ≤3.8% (random 0.3%). ⇒ ρ₆ at split-squarefree m is not a simple Gauss expression.
