# Solver G: direct numerical test of the T ⟂ τ₆ non-correlation (HB-recursion route, solverE §3.2(a))

**Verdict.**
* **Numerics: no T-specific correlation is detected** (one caveat at the largest Q, see §3.2). At every accessible size, the probe vector T(m) = ᾱ(m)²γ₁(m)𝓡(m) is at the random-phase level against the untwisted sextic Kubota sums A_m(h′). This also holds for its 36 ray-class twists and for every one of its μ-components. The level is flat in Q (#m = 131 → 2100) and in Y (12.5k → 200k).
* **The test detects real main terms.** Validated on the cubic analogue, where the answer is known:
  * the aligned vector reaches full alignment, r ≈ #m·f, growing linearly in Q;
  * the ᾱ²-twisted version drops to the null.
* **No proof follows.** Within the multiple-Dirichlet-series (infinity-type) framework, the ℓ ≠ 0 mechanism gives a qualitative "no main term" for each fixed μ. Quantitatively, the route needs near square-root cancellation on average over μ ≤ D = Q^{0.8}, which no convexity-level argument reaches. The obstruction is again the ℓ¹ versus ℓ² mass of 𝓡, the same wall as PSR, now seen from the pole side.

Scripts: `exp/noncorr/` (`lib.py`, `validate_lib.py`, `noncorr.py`, `grid.py`, `pooled.py`, `tube_hull.py`, `summarize.py`; outputs in `exp/noncorr/out/`). Run from that directory.

---

## 1. What is measured

Split-prime proxy over Z[ω], reusing the conventions of `exp/eis.py` and `theta_id.py`. `lib.py` is a vectorized rewrite (discrete-log tables per rational prime). `validate_lib.py` checks it against the slow reference code: 0 symbol mismatches, Gauss sums agree to 4·10⁻¹⁴, |g̃| = 1, and multiplicativity holds.

**Objects.**
* a: squarefree products of split primary primes, N(a) ∈ [Y, 2Y].
* m: the same, with N(m) ∈ [Q, 2Q].
* h′ ∈ {1} ∪ {small primary elements, N ≤ 400}, split and inert.

  A_m(h′) = Σ_{a~Y} g̃₆(a) χ_a(h′) χ̄_a(m)     (untwisted Kubota sum; 'lam' additive character e(Tr x/(λa)), the paper's convention)
  C(h′)   = Σ_m v(m) A_m(h′)

The probe vector is T(m) = ᾱ(m)²·γ₁(m)·𝓡(m), where:
* γ₁ = g̃₆(m);
* 𝓡(m) = Σ_{N(μ)≤D} ρ_μ (μ/m)₂, with ρ_μ = g̃₃(μ)α(μ)N(μ)^{-1/2}, the same proxy as `sieve3.py`;
* D = Q^{0.8}.

The cubic validation uses g̃₃ and χ² throughout (the "chk" convention, in which Patterson's τ₃ ∝ conj g̃₃ is recovered, as in theta_id).

**Statistics.** All are normalized so that the null value is about 1.

* **r(h′) = |C|² / Σ_m |v_m|²|A_m|².**
  * For random-phase v, E r = 1 exactly, with an Exp(1) distribution.
  * A vector fully aligned with the pole gives r ≈ #m·f.
  * f is the pole energy fraction of A, estimated from two independent random halves of the a's: f = 4Re⟨A₀,A₁⟩/‖A‖².
* **#m·ρ̂², the pole-projection estimator.**
  * Two-half version: #m·Re(C₀C̄₁)/(‖v‖²Re⟨A₀,A₁⟩).
  * Pooled version: over 5 disjoint a-windows [Y_i, 2Y_i], which share the same residue direction ϱ₆(dm) but have independent noise (`pooled.py`).
  * It estimates #m·|⟨v, conj pole⟩|²/(‖v‖²‖pole‖²). Null ≈ 1; full alignment = #m.
  * It isolates the pole part and is insensitive to the non-pole noise.
* **κ_a = |C(h′)|² / (Σ_a|g̃(a)χ_a(h′)|² · E_a), with E_a = Σ_a |Σ_m v_m χ̄_a(m)|².**
  * This is the Cauchy–Schwarz ratio over a, i.e. the "Y part" bound of solverE §3.1. It is ≤ 1 always.
  * The generic route loses exactly when κ_a is O(1).
  * Random-like behaviour is κ_a ≈ 1/#m or smaller.
  * The m-side Cauchy–Schwarz ratio |C|²/(‖v‖²‖A‖²) is ≈ r/#m and behaves the same way.

**Comparison vectors.**
* Null: random phases with |v| = |T| (200–400 draws).
* Positive controls:
  * untwisted aligned vectors: g̃·ν (cubic); the constant / ray-class family ν (sextic);
  * the data-driven pole proxy conj(A from the other half).
* Infinity-type tests: the same vectors times ᾱ².
* Components of T:
  * T·ν for the 36 ray-class characters ν = χ_m(λ)^iχ_m(−1)^jχ_m(ω)^k (standing in for ν₁ and ψ);
  * the per-μ pieces ᾱ²γ₁(μ/m)₂;
  * the untwisted γ₁𝓡ν.

## 2. Cubic validation (answer known: τ₃ ∝ conj g̃₃)

`grid.py --n 3 --conv chk` (Y = 5k–40k, Q = 500–4000) and `pooled.py --n 3` (5 windows Y = 2.5k–40k):

| Q | #m | pole z (pooled) | aligned g̃₃·ν: max #m·ρ̂² (pooled) | T: #m·ρ̂² | T·ν (36): mean±se | T-components (#μ): mean±se | ᾱ²·g̃₃·ν mean | RAND mean±se |
|---|---|---|---|---|---|---|---|---|
| 500 | 131 | 34 | **131** | 2.75 | 1.36±0.17 | 1.48±0.19 | 0.59 | 1.02±0.05 |
| 1000 | 272 | 48 | **274** | 0.11 | 0.73±0.12 | 1.04±0.18 | 0.52 | 1.07±0.05 |
| 2000 | 525 | 67 | **529** | 0.21 | 0.78±0.12 | 1.01±0.10 | 0.41 | 1.02±0.05 |
| 4000 | 1069 | 93 | **1086** | 1.10 | 1.19±0.15 | 0.94±0.07 | 0.27 | 1.00±0.05 |

Per-cell grid (16 cells), at h′ = 1:
* The aligned control has r = 115 → 1000 ≈ #m·f, linear in Q, with κ_a ≈ 0.2–0.67 (O(1): it saturates Cauchy–Schwarz over a).
* T has r averaging 0.95 over 16 cells (null 1 ± 0.25), max 3.6, and κ_a·#m ≈ 0.05–1.5, the same as RAND (0.4–1.4).
* ᾱ²·(pole proxy) has r ≈ 0.0–1.4, against 63–764 for the untwisted pole proxy. This is the infinity-type cancellation of Σ ᾱ²|τ₃|².

So the setup detects a genuine pole main term at full strength and correctly reports its absence for T, which is the expected answer in the cubic case (solverE §3.2(a)).

## 3. Sextic results

*(The per-cell h′ ≠ 1 statements in §3.1–3.2 refer to the first runs at index m h′⁵ and are superseded by §3.2b.)*

### 3.1 The pole is weak at accessible sizes

f(h′=1) = 0.1–0.45, with per-cell z = 2–5.6 (`grid.py`, Y = 25k–200k). Pooled over 5 windows, the pole z is given in the table below. This matches theta_id (+30% energy).

The pole is concentrated at h′ = 1. For prime h′ it is ≈ 0 (f(h′ ≠ 1) mean ≈ 0), consistent with Kazhdan–Patterson-type vanishing of τ₆ at h′⁵m for most local exponents. The cubic case behaves similarly: f is about 0.4 at h′ ≠ 1, against 0.85 at h′ = 1.

**Sextic positive control.** Untwisted sextic Gauss vectors do **not** align with the pole (g̃₆·ν family at null, as theta_id found). The constant vector does:
* r(h′=1) = 10–92;
* pooled #m·ρ̂² ≈ 0.3–0.4·#m, i.e. Σ_m ρ₆(m) has a main term.

Its ᾱ²-twist drops to r ≈ 0.0–0.2. Untwisted ray-class twists ν: max r grows with Q (≈15 → 46 on average over Y).

### 3.2 T and all its pieces are at the null

Grid (16 cells, h′ = 1):

| Q | #m | mean over Y: r(T) | T·ν (36) ⟨r⟩ | T-comps ⟨r⟩ | RAND ⟨r⟩ | +ctrl max r | #m·f | κ_a(T)·#m | κ_a(RAND)·#m | κ_a(+ctrl) |
|---|---|---|---|---|---|---|---|---|---|---|
| 500 | 131 | 0.65 | 1.10 | 1.32 | ≈1 | 14.7 | 9–59 | 0.009 | 0.014 | 1.6e-3 |
| 1000 | 272 | 2.14 | 0.79 | 0.95 | ≈1 | 27.9 | 52–75 | 0.052 | 0.027 | 2.0e-3 |
| 2000 | 525 | 0.35 | 0.90 | 0.98 | ≈1 | 38.6 | 64–152 | 0.010 | 0.043 | 1.8e-3 |
| 4000 | 1069 | 1.27 | 1.23 | 1.07 | ≈1 | 45.9 | 96–232 | 0.141 | 0.105 | 1.5e-3 |

r(T) over all 16 cells: mean 1.10 (null 1 ± 0.25).

**One outlier.** At (Y=50k, Q=1000), r = 6.65 and #m·ρ̂² = 23.
* Bonferroni p ≈ 0.02 over 16 cells.
* It does **not** recur in the three independent a-windows at the same Q: r = 0.00, 1.72, 0.18.
* A real pole correlation would recur, because the residue direction is the same for every window.
* The pooled estimator below, which uses all windows, is the decisive check.

**Pooled cross-window estimator**, h′ = d = 1 (`pooled.py`, 5 disjoint windows Y = 12.5k–200k, `out/pooled_sextic*.txt`). Entries are #m·ρ̂²: null ≈ 1, full alignment = #m. Families are shown as mean±se / max.

| Q | #m | pole z | T | const (untw) | ν untw (36) | ᾱ²ν | T·ν (36) | T-comps (#μ) | g̃₆·ν untw | ᾱ²g̃₆·ν | RAND (400) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 500 | 131 | 13 | 0.3 | 41 | 2.2±1.2/41 | 1.11 | 0.73±0.19 | 1.98±0.38 (39) | 1.09 | 1.11 | 0.75±0.07 (sd 1.4) |
| 1000 | 272 | 16 | −1.7 | 88 | 3.5±2.5/88 | 0.84 | 0.44±0.11 | 1.06±0.16 (68) | 1.03 | 1.74 | 0.80±0.08 |
| 2000 | 525 | 17 | 0.1 | 107 | 6.3±4.2/114 | 0.25 | 0.47±0.18 | 1.08±0.19 (115) | 0.73 | 1.54 | 1.08±0.10 |
| 4000 | 1069 | 21 | −1.6 | 33 | 7.8±5.4/197 | 0.84 | 1.75±0.46 | 1.31±0.18 (203) | 1.14 | 1.91 | 0.82±0.10 |
| 8000 | 2128 | 33 | 3.8 | −30 | 10.8±11/400 | 0.51 | **2.47±0.37** | 0.81±0.09 (350) | **3.49±0.73** | 1.91 | 1.04±0.10 (sd 2.0) |

The large-window-only rerun (Y = 50k–200k) gives: Q = 4000: T −2.2, T·ν 1.64±0.41; Q = 8000: T 9.6 (RAND sd 2.4), T·ν 3.07±0.53, g̃₆·ν 4.08±0.86, g̃₆𝓡ν 2.32, ᾱ²ν 0.56, T-comps 0.80.

* **Positive controls (untwisted const / ν) grow with Q** (max 41 → 400, i.e. ≈ 0.2–0.3·#m). That is a main term.
* **ᾱ² kills them** (ᾱ²ν at 0.25–1.1).
* **T and its μ-components are at the null through Q = 4000.**
* **Q = 8000 anomaly.** Every γ₁-carrying family is mildly elevated, at ≈ 2–4× null, i.e. ρ² ≈ 1.5–4.5·10⁻³:
  * T·ν, untwisted g̃₆·ν, g̃₆𝓡ν, ᾱ²g̃₆·ν, and T itself in the large-window run (≈ 3.6σ nominal).
  * The excess is the same with or without ᾱ² and absent for the T-components (ᾱ²g̃₆(μ/m)₂). It is therefore not the angular-sensitive Patterson mechanism. It looks like a weak γ₁-correlated component common to the windows when Q is large relative to the window range, possibly a proxy artefact (split-prime-only a's, sharp cutoffs).
  * Nominal significance is inflated: the two runs share data, and family members are not independent.
  * I cannot resolve it at these sizes. **It is flagged as the one place where the Q^{-1} decay is not clean.**

### 3.2b Index correction and the binding frequencies d ≠ 1

**The first runs used the wrong index for h′ ≠ 1** (found independently by the literature reader and by solverK).
* g̃(a)χ_a(h′)χ̄_a(m) = g̃(a, m·h′⁵), so those runs probed τ₆(m h′⁵). There the pole vanishes identically for prime h′, consistent with the measured f ≈ 0.
* The h′ = 1 results above are unaffected.
* Fix: `--index md`, now the default in `grid.py` / `pooled.py`, uses χ̄_a(d), i.e. A_m(d) = Σ_a g̃₆(a)χ̄_a(d)χ̄_a(m), whose pole is ∝ ϱ₆(dm), the index in (★★).
* Cubic check: with index md the pole is strong at every d (z ≈ 30–70, against ≈ 0–5 before), and T stays at null (#m·ρ̂² 0.7–1.2).

**Sextic, index md.** The pole at d ≠ 1 is now present: grid f ≈ 0.1–0.3; pooled z = 6–28.

T at d ≠ 1, pooled mean of #m·ρ̂² over the 7 nontrivial d: 0.97, 0.14, 1.04, 1.68, 2.05 for Q = 500…8000 (null ≈ 1).

Per N(d) (`grid.py --index md`, Y = 50k–200k, Q = 1000–4000; 9 cells; conjugate pairs d, d̄ give identical statistics, so ≈ 4 independent d per cell):

| N(d) | ⟨r(T)⟩ | ⟨κ_a(T)⟩·#m | ⟨κ_a(T·ν)⟩·#m | ⟨κ_a(RAND)⟩·#m | ⟨f⟩ |
|---|---|---|---|---|---|
| 1 | 1.45 | 0.052 | 0.035 | 0.040 | 0.22 |
| 7 | 1.09 | 0.026 | 0.030 | 0.034 | 0.16 |
| 13 | 0.84 | 0.025 | 0.027 | 0.030 | 0.16 |
| 19 | 1.71 | 0.085 | 0.038 | 0.033 | 0.17 |
| 25 | 1.10 | 0.034 | 0.032 | 0.032 | 0.18 |

* Over all 63 d ≠ 1 instances, mean r(T) = 1.2.
* The largest single values (r = 5.05 at N(d) = 19; 4.88 at N(d) = 25) are within the expected maximum of ~45 independent Exp(1) draws.
* No trend in N(d) is visible.
* **Caveat on d-range:** N(d) ≤ 25 only. In the experiment P_a = Y²/Q ~ 10⁶–10⁸, so this covers only N(d) ≤ P_a^{6s} with s ≲ 0.03, which is the low end of solverK's (Π_s) range. Larger d (N(d) up to 400 is cheap) was not run for lack of time.

### 3.2c Explicit size of the pole-part correlation relative to Cauchy–Schwarz

The pole-part Cauchy–Schwarz ratio is
  ρ² = |⟨T, conj ϱ₆(d·)⟩|² / (‖T‖²‖ϱ₆(d·)‖²)
(equal to 1 for the aligned worst case, which is the c = 1/6 situation). `summarize.py` gives 2σ upper bounds from the pooled estimator, d = 1:

| Q | 500 | 1000 | 2000 | 4000 | 8000 |
|---|---|---|---|---|---|
| ρ²(T) ≤ | 2.3e-2 | 1.2e-2 | 7.7e-3 | 3.6e-3 | 3.6e-3 (6.9e-3 large-window run) |
| ⟨ρ²(T·ν)⟩ ≤ | 8.5e-3 | 2.4e-3 | 1.6e-3 | 2.5e-3 | 1.5e-3 |
| 1/#m | 7.6e-3 | 3.7e-3 | 1.9e-3 | 9.4e-4 | 4.7e-4 |

* Log-log slope in Q: T −0.70, T·ν family −0.49. Restricted to Q ≤ 4000: ≈ −0.9 (T). The no-main-term value is −1; a main term is 0. The cubic analogue gives −1.17 / −1.05.
* **Rough exponent statement.** Through Q = 4000 the pole part of |C(d)|² is ≈ Q^{-1} times its Cauchy–Schwarz bound, i.e. full square-root cancellation. Including Q = 8000, the 2σ-safe statement is ρ² ≲ Q^{-1/2}, with a possible floor ρ² ~ 10⁻³ (γ₁-related, not ᾱ²-sensitive, see the anomaly in §3.2).
* The needed saving is only Z^{-0.029} on |C|² (Z^{-0.0144} on the amplitude), i.e. ρ² ≤ P_a^{-1/3}.
  * In paper geometry Q^{-1/2} is Z^{-0.42}.
  * In the experiment the analogous threshold P_a^{-1/3} is ≈ 1e-2 (Y = 10⁵, Q = 8000).
  * Even the flagged Q = 8000 values (≤ 7e-3) are below it.
  * Extrapolating a possible 10⁻³ floor to the paper's sizes is not justified either way.

### 3.3 Coherence across μ (the ℓ¹ versus ℓ² issue)

Let c_μ = Σ_m ᾱ²γ₁(μ/m)₂A_m(1). Compare
* |Σ_μ ρ_μ c_μ|² (what T sees),
* Σ_μ |ρ_μ|²|c_μ|² (incoherent sum),
* (Σ_μ |ρ_μ||c_μ|)² (the worst case, i.e. Cauchy–Schwarz / triangle in μ).

Observed:
* incoherence ratio |Σρc|²/Σ|ρc|²: 0.0–3.0 in all cells, cubic and sextic;
* worst-case ratio: 25–130, growing roughly in proportion to #μ (39 → 203).

So numerically the μ-pieces are mutually incoherent and 𝓡 behaves like its ℓ² mass O(1), not its ℓ¹ mass D^{1/2}. This is the same phenomenon PSR asserts on the dual side.

### 3.4 Scaling answer

A main term would make r and #m·ρ̂² grow ∝ #m (cubic aligned control: ×8 from Q = 500 to 4000), and κ_a stay O(1). For T we see:
* r flat at ≈ 1;
* κ_a·#m = O(1), i.e. κ_a decays like 1/Q exactly as for random vectors;
* pooled ρ̂² consistent with 0.

In exponent terms this is square-root cancellation: |C(h′)|² ≈ Σ_m|T|²|A_m|². That is a factor ≈ #m below the Cauchy–Schwarz "Y part", far more than the P_a^{1/6} = Z^{.0144} needed for c = 1/6 → 0 (§4.1).

## 4. Does a provable statement follow?

### 4.1 What is actually needed (pole part, h′ ≈ 1)

Write X = |Σ_m T(m)ϱ₆(dm)|, so that C_pole ≈ Y^{2/3}X.

Removing the P_a^{1/6} loss of the one-step HB recursion needs |C(h′)|² ≪ Y·P_a^{-1/6}·E_a with E_a = Q². That is,

  X ≪ Q·Y^{-1/6}·P_a^{-1/12} = Z^{0.7490}.

The **true trivial bound** is X ≤ ‖T‖₂‖ρ₆‖₂ ≈ Q^{1/2}·Q^{5/12} = Q^{11/12} = Z^{0.7635}. This uses:
* Σ_m|ρ₆|² ≍ Q^{5/6}, which is exactly the Gram exceptional term of B0.22;
* ‖𝓡‖₂² = O(log).

So only a saving of Z^{-0.0145} over Cauchy–Schwarz is needed. The cubic-type situation (exact orthogonality) is not required: any power saving below Q^{-0.0175} on X suffices for (a).

The non-pole part (b) is measured in the same r and is also at the random level numerically. Proving it is the Type-II problem of solverE §3.2(b), which is not addressed here.

### 4.2 The infinity-type (Dunn–Radziwiłł Prop 5.2, ℓ ≠ 0) argument

Consider the triple series

  Z(s,w,u) = Σ_{a,m,μ} g̃₆(a)χ̄_a(d)χ̄_a(m) · ψᾱ²(m)g̃₆(m)(μ/m)₂ · g̃₃(μ)α(μ) · N(a)^{-w}N(m)^{-s}N(μ)^{-u}.

Then Res_{w=2/3} Z = Σ_m T_u(m)ϱ₆(dm)N(m)^{-s}.

Using (μ/m)₂ = χ_m(μ)³, the m-node is an **ℓ = 2 twisted sextic Kubota series at index aμ³**. The DR mechanism (U-E4) says it is **entire**. The μ-node is an ℓ = 1 twisted cubic Kubota series, also entire. So the only polar divisor is the a-side Patterson line w = 2/3. The twist removes exactly the s-polar line that, in the untwisted case, intersects it and produces the main term. That is the precise sense in which ᾱ² "forces" cancellation of the pole part.

**Quantified by Bochner's tube theorem.** `tube_hull.py` runs an exact LP over the holomorphy regions:
* R-m: m outermost;
* R-a: a outermost;
* R-μ: K continued via ℓ ≠ 0 entireness.

It uses convexity growth exponents, all heuristic and flagged:
* Kubota series: N(n)^{(1−σ)/2};
* Θ_m(u): conductor N(m)^2.

The LP gives the following.

* **Per fixed μ** (u large, μ-growth ignored): the residue series is holomorphic for Re s > **5/6** on w = 2/3. So
  Σ_{m~Q} ψᾱ²γ₁(m)(μ/m)₂ϱ₆(dm) ≪ Q^{5/6+ε}(N(μ)N(h′))^{O(1)},
  a power saving over the per-μ trivial Q^{11/12}.
  * **This is the provable-in-principle statement**: no main term per μ.
  * It is conditional on (i) U-E4 (n = 6, ℓ ≠ 0 entireness with polynomial bounds), (ii) Patterson continuation of Ψ(w;n) with uniform convexity in n, (iii) the standard several-variable Phragmén–Lindelöf / Bochner step.
  * It is the sextic analogue of what Patterson's explicit τ₃ gives for free in the cubic case.
* **For T itself it fails, by a lot.**
  * With the μ-sum inside the hull: X ≪ Z^{1.30} (conductor m²) or Z^{1.06} (conductor m).
  * Even with every growth exponent halved (a subconvexity-type input): Z^{0.97}.
  * Requirement: Z^{0.749}.
  * Reason 1: the framework bounds ρ₆ pointwise via convexity, at N^{1/6} at w = 2/3, whereas its true average size is N^{-1/12}.
  * Reason 2: the μ-sum is charged its ℓ¹ mass D^{1/2} = Z^{0.333}.
* **Per-μ bounds plus the triangle inequality in μ cannot work at all.** Suppose each μ-piece is bounded by Q^σN(μ)^θ. The requirement becomes σ ≤ 0.4997 − 0.80θ, while square-root cancellation per μ is σ = 5/12 = 0.4167.
  * So one needs essentially square-root cancellation per μ (within Q^{0.083}), uniformly in μ ≤ Q^{0.8} with conductor exponent θ ≲ 0.1. That is Lindelöf/GRH-strength.
  * With Cauchy–Schwarz over μ instead, the same threshold must hold on average over μ ≤ D. The available mean-value tool (Heath-Brown's quadratic large sieve) gives (D+Q)Σ|ρ₆|², which is exactly the trivial bound because D < Q.

**Net.** The infinity-type argument is the right mechanism. It explains the numerics and plausibly yields a rigorous qualitative per-μ statement, conditional on U-E4. It does **not** yield the quantitative bound the route needs.

The missing input is an incoherence-across-μ statement for the family μ ↦ Σ_m ψᾱ²γ₁(m)(μ/m)₂ρ₆(m) (μ ≤ Q^{0.8}). Equivalently, it is a beat-the-quadratic-large-sieve bound for the vector b_m = ψᾱ²γ₁ρ₆(m) in the range D < Q. This is the pole-side mirror of PSR. The numerics (§3.3) say it is true, with the incoherence ratio O(1).

## 5. Caveats (numerics are evidence, not proof)

* **Split primes only.** a, m and μ are squarefree products of split primes. Inert primes, S-units, non-squarefree a (where the Kubota series has its full Euler structure), branches, the ϑ twist and the smooth V♯ are omitted. Sharp cutoffs throughout. The Kubota sum is therefore not the true Patterson series; the pole is present but diluted (f ≤ 0.45).
* **Small sizes.** #m ≤ 1069 (2100 in pooled), Y ≤ 2·10⁵, D ≤ Q^{0.8} ≤ 1330. The sextic pole has z ≈ 3–6 per cell, so sensitivity to a *partial* correlation is limited. The pooled estimator bounds ρ² ≲ few/#m (see table). P_a in the experiment is Y²/Q ≈ 10⁵–10⁷, far from the paper's regime, and the exponent claims are extrapolations.
* **Proxy for T.** ν₁ and ψ are represented by small ray-class characters, and ρ_μ by the sieve3 proxy (normalization N^{-1/2}, flagged in solverE). The real T(R) of solverD Claim 3 also carries marks and branch sums.
* **Positive control in the sextic case.** The sextic positive control is the constant / ray-class family plus the data-driven pole proxy, since τ₆ is not explicit and no explicit sextic analogue of g̃₃ exists. The test therefore shows that T is orthogonal to the empirically estimated pole direction, not to an identified τ₆.
* **Theory (§4.2).**
  * The convexity exponents for the twisted Kubota series and for Θ_m are heuristic.
  * U-E4 (ℓ ≠ 0 entireness for n = 6) is an inference, as in solverE.
  * The several-variable convexity step is assumed.
  * The LP is bookkeeping, not a proof.

## 6. Bottom line

* At the correct index ϱ₆(dm), for d = 1 and N(d) ≤ 25:
  * T shows no detectable correlation with the sextic Patterson pole. It is at the random level in Y, and in Q up to 4000, in a setup that detects the cubic main term at full strength.
  * The pole-part Cauchy–Schwarz ratio is ρ² ≈ Q^{-1} through Q = 4000; it is ≲ Q^{-1/2} (2σ) including Q = 8000.
  * At Q = 8000 a weak excess (ρ² ~ 10⁻³) appears in every γ₁-carrying family, twisted or not. It is unresolved and possibly a proxy artefact.
* Twisting by ᾱ² kills every untwisted main term tested, in both the cubic and the sextic case.
* The infinity-type mechanism proves (conditionally) only per-μ absence of main terms. The quantitative c = 1/6 → 0 step needs near square-root cancellation coherent over μ ≤ Q^{0.8}, the same μ-incoherence wall as PSR.
* Unconditional θ is unchanged.
