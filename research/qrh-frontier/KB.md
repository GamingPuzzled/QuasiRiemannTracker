# QRH knowledge base (harness memory)

Structured memory for the advisor → solvers → verifier loop, modelled on the UCLA First Proof batch-2 harness
(github.com/1stproof/batch-2, batch-2-submissions/ucla). The advisor reads this file before each round, gives each
solver only the context it needs, and records results here. `FINDINGS.md` remains the detailed log.

Target: lower θ below 683505193/781250000 ≈ 0.87488664704 (cycle25 quartic; θ = 11/12 − ℓ/4).
Paper: proofs/cycle25-quartic-20261010/manuscript/paper.md. Exponent model: scripts/slack2.py.

## Verifier protocol
Every positive claim (lemma, estimate, θ gain) goes to an independent verifier agent that has not seen the solver's
reasoning beyond the write-up. The verifier returns SCORE 1–10, MAJOR_GAPS, MINOR_GAPS.
- ≥ 9: accepted → `proven_results`.
- 7–8: one targeted patch round on the minor gaps, then re-verify (max 2 patches per gap).
- ≤ 6: → `failed_attempts`, with the gap as the reason. Future solver prompts receive it.

## proven_results (verified; usable as lemmas)
- P1 Record arithmetic is correct; θ − b0 = 3.6e-14 is the only rounding slack (audit.py, solverH).
- P2 Active constraints at the optimum: M+ℓ=1, p ≤ 1/2, the δ/h saddle, d=h. All other certificate constants have positive margin (solverH).
- P3 Conversion: ΔF at the binding point (δ0≈.4693, p=.5, κ0≈.7498) gives Δθ ≈ 0.23·ΔF. ΔR = −1e-3 gives Δθ ≈ −1.4e-4.
- P4 Hölder mixing of moment lines never beats the minimum line. The mixed MP line is the average of the pure lines (solverH).
- P5 The amplitude cap p ≤ 1/2 is the sharp pointwise bound (B0.16), saturated by single-zero rows (solverF).
- P6 Twist escape: ξ must have order dividing 6, which forces the C0=1 exceptional Gram term (tick 4).

## failed_attempts (do not repeat without a new idea)
- F1 Relaxing M+2z ≤ N_* (column-branch cost s, break-even needs ≤ .25s).
- F2 Parameter retuning or geometry (x is minimax-optimal; y<x is worse).
- F3 Improving estimates away from δ0 (saddle point; zero first-order gain).
- F4 (★) via Poisson/duality (involution), via a zero-count bootstrap (E′ > E), or via the Heath-Brown identity (loses Möbius structure).
- F5 Product-witness mean value M·S (R changes only at the noise level).
- F6 MP slot cost: c = 6 is optimal within the method (generic child saturates B3.70).
- F7 Gram exceptional term via column weights, other ξ, or non-S twists (Solver C Prop. O, partially proven).
- F8 Twisted energy E_a: the generic sieve gives c = 1/3; DR §9 needs pointwise PSR (stronger than GRH) (solverD/E).
- F9 Identifying ρ₆ at squarefree m as a Gauss-sum expression (theta_id).
- F10 Amplitude cap: lowering p_max, slot-based zero detectors, height averaging (solverF).
- F11 Certificate-constant sharpening, triple-point mixing, bin splitting (solverH).
- F12 ε-Hölder for the numerator: higher moments via pointwise input, MP0 factorization, Halász–Montgomery or
  Huxley-type large values (solverJ; adversary bin saturates the toolkit).

## bottlenecks (precise missing inputs, with payoff)
- B1 Twisted energy E_a ≪ Q(Q+Y)Z^ε for T = ψᾱ²γ₁𝓡 (equivalently τ₆ non-correlation or PSR): c 1/6 → 0, θ ≈ 0.865152.
- B2 (★) Möbius sextic mean value beyond the Poisson range with λ < 5/6: λ = .8 gives θ ≈ 0.874628.
- B3 MP slot cost below 6κ (new amplifier with a better row-growth to gain ratio): c = 5 gives θ ≈ 0.874658.
- B4 NEW (solverH): a mean value for 1/(L(ψ̄_u)L(ηψ̄_u)), counting rows with simultaneous low zeros. This adds a new count line at the triple point. Not yet priced.
- B5 Any change of ≥ 1.3e-6 in the Hölder 1/4 or the Gram 1/12 coefficient gains 1e-7.
- B6 NEW (literatureG): (★★), our specific vector run through de Faveri's proof (arXiv:2610.04045, Lemma 6.1):
  Σ_{d≤P_a}(P_a/N d)^{1/6}|Σ_{a~Y}F_a g̃₆(a)χ̄_a(d)|² ≪ Z^ε·Y·Σ|F_a|².
  Needed: beat Cauchy–Schwarz by Z^{.0144} in the pole part (correlation of c_m with ρ₆(dm)) and by Z^{.0072} in the
  off-pole part (needs Lindelöf-on-average for sextic Kubota series). ANY partial saving gives c < 1/6, which beats
  the record. Normalizations are heuristic and unverified.

- B6 update (solverK, PENDING VERIFICATION, see verifierK): a perfect pole gives c = 1/7 (θ ≈ 0.873641), not 0.
  A saving s must hold uniformly for N(d) ≤ P_a^{6s}, which gives c = max(1/6 − s, 1/7). The missing theorem is the
  averaged non-correlation (Π_s) in solverK §5. Angular-twist handle (i) fails (non-abelian torus).

- B5 update (solverJ, unverified negative claims): there is no ε-improvement of the numerator Hölder step using the paper's
  toolkit. The construction scripts/epsholder/adversary.py builds an exponent-level virtual bin that satisfies every row
  inequality to within 2.2e-16 (UM, B0.9, MI, MP/MP0/mixed, pointwise, the cap, Lemma 1W) and saturates the fourth moment.
  So any gain needs an estimate outside the toolkit. The precise target is a sixth moment on bin rows beating the trivial
  bound by U^{η·n} with η > 0.00237 (e.g. η = 0.0025 gives θ = 0.8748847499834). General-coefficient sieves give exponent
  ≥ 3n = 1.263 > trivial 1.1976, and the short-length range (3n ≤ 1) is never available because n > 8/21.

## corrections
- C1 (literatureG): the pole-part theta index is τ₆(h′m), not τ₆(h′m⁵) as solverE §3.2 and FINDINGS state. With m⁵ the
  term would vanish trivially (τ₆(Mp⁵)=0); with h′m it is needed at unknown squarefree indices.
- C2 (literatureG): de Faveri arXiv:2610.04045 Thm 1.1 PROVES the n-th order sieve A+B+A^{1−1/n}B^{2/n}+A^{2/n}B^{1−1/n}
  (n ≥ 3), superseding BGL (M+N+(MN)^{2/3}). This gives c = 1/6 unconditionally (= record, no gain) and makes 1/6 the
  proven generic barrier. Its (★) and zero-density consequences are weaker than the record's.

## frontier (open tasks / in flight)
- Non-correlation numerics DONE (solverG, evidence only, split-prime proxy). At index m·d the pole-part correlation
  of T with ϱ₆ is at the random level through Q = 4000 (ρ² ≈ Q^{-1}; needed ρ² ≤ Z^{-0.029}). Only N(d) ≤ 25 was
  tested. At Q = 8000 every γ₁-carrying vector shows ρ² ~ 1e-3 (≈3σ), unresolved, so the safe statement is
  ρ² ≲ Q^{-1/2}. The cubic validation passes. No proof: coherence across μ ≤ Q^{0.8} hits the same wall as PSR.
  Open check: re-run Q = 8000 (and larger N(d)) to decide whether the excess is real or a proxy artefact.
- arXiv reading: done (literatureG). Dunn 2403.13151 and Castillo–de Faveri–Dunn 2604.01316 do not transfer.
  Bröker–Hoffstein Conj. 5.7 (τ₆(π)² = cubic Gauss sum × patternless unit) supports PSR heuristically.
  Unread lead: Diaconu 2004 (general-order Kubota series; likely "entire for ℓ ≠ 0" at n = 6).
- Round 2 solvers in flight: B4 simultaneous zeros (solverI), B5 ε-Hölder (solverJ).
- Mathbox verification session (branch research/qrh-mathbox) in flight.

## advisor_notes
- Round 1 (2026-10-11): solverF (amplitude cap), solverH (small gains) and literatureF all returned negative; they are recorded above.
- Round 2: literatureG returned B6 and corrections C1–C2. Correction C1 was relayed to the running noncorr agent.
  Plan for round 3: a solver on B6 partial savings (pole part first), informed by the noncorr numerics.
