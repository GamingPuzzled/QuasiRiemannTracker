# Literature search F: new inputs that could lower θ ≈ 0.87488664704 (search date 2026-10-11)

Scope: items NOT already in FINDINGS.md / literatureE.md (DR, DFDH 2607.07911, BGL, DDHL 2306.11875,
Guth–Maynard, Chen 2507.08296, BB06, DIPP 2607.27131, Gao–Zhao 2201.01885 are excluded).

**Access caveat (important).** In this session the egress proxy allowed only github.com and
raw.githubusercontent.com. arxiv.org, export.arxiv.org, semanticscholar, ems.press, chinta.ccny.cuny.edu,
*.github.io and journal sites all failed (403 from the proxy, or DNS failure in WebFetch). So:
- **[O]**: I opened the source and read the quoted text myself (GitHub-hosted only).
- **[S]**: I saw it only as a web-search result summary. I did NOT open it, and the statement is
  unverified. These are listed in §3 as leads to open later, not as evidence.

---

## 1. Priority 1: new QRH records or claims (all [O])

**Bottom line: as of 2026-10-11 there is no claim below 0.87488664704 anywhere I could find.**

| Source | Claimed θ | Status |
|---|---|---|
| Tracker catalogue, `GTimG/QuasiRiemannTracker` main `catalogue/results.json` (as_of 2026-10-10) [O] https://raw.githubusercontent.com/GTimG/QuasiRiemannTracker/main/catalogue/results.json | best entry: cycle25-quartic-20261010 = 0.87488664704 | the other entries are OpenAI 7/8, Argonaut 0.87499975, Liu 0.87495707, ProofCouncil/nielstron 0.874957019…, Levy 10499/12000 ≈ 0.8749167 |
| Tracker PR list [O] https://github.com/GTimG/QuasiRiemannTracker/pulls?q=is%3Apr | 11 PRs, all dated Oct 9–10; the latest record is #9 (quartic, c ≈ 0.12511335) | — |
| PR #7, Gaozhongpai, "exact slack retuning" [O] https://github.com/GTimG/QuasiRiemannTracker/pull/7 | 0.874957019420098946… | closed; "recovers numerical slack…introduces no new arithmetic estimate"; superseded |
| PR #8, Hailey Collet, "single compensating prime", 29/33 [O] https://github.com/GTimG/QuasiRiemannTracker/pull/8 | 29/33 ≈ 0.8788 | non-record alternative proof; open (see §2.1) |
| Forks (akashlevy, GamingPuzzled, Gaozhongpai, HaileyStorm, nielstron) [O] https://github.com/GTimG/QuasiRiemannTracker/forks | no new bound in their branch lists | — |
| openai/math `history.md` [O] https://raw.githubusercontent.com/openai/math/main/history.md | the Oct 7 2026 entry has withdrawals and fixes, none touching QRH/family 003 | — |
| LW56e2/quasi-riemann-improved README [O] https://raw.githubusercontent.com/LW56e2/quasi-riemann-improved/main/README.md | 4374785349/5·10⁹ ≈ 0.8749570698 | see note below |
| stevenecreech/math README [O] https://raw.githubusercontent.com/stevenecreech/math/main/README.md | proposes 437479/500000 = 0.874958 | "full dependency build and final axiom checks are pending" |

Note on LW56e2: its README proves an "exact two-bin obstruction". It shows that retuning the same two
parameters with the specified balanced counting estimate cannot lower the boundary by 10⁻¹². This
matches our FINDINGS conclusion that the framework is at a parameter-level saddle.

Kintali, "A Short Proof of the Quasi-Riemann Hypothesis" (shivakintali.github.io/papers/QRH.pdf) [S] is a
weaker claim (47/48 per the search summary), so it is irrelevant here.

## 2. Mathematical items I opened

### 2.1 Collet, "A zero-free half-plane from a single compensating prime" (29/33), PR #8 manuscript [O]
https://raw.githubusercontent.com/HaileyStorm/QuasiRiemannTracker/single-prime-29-33-20261010/proofs/single-prime-29-33-20261010/manuscript/zeta_zero_free_4_33_self_contained.md

- **Thm D.1 (mixed plain moment), verbatim shape.** ψ_k(n) = τ(n)χ_n(k)1_{(n,R)=1} with moving twist
  τ(n) = ν(n)∏_{p|A}χ_n(p)^{e_p}, where q_{rad A} ≤ Z^q, and S_k(N;W) = N^{-1/2}Σψ_k(n)W(q_n/N). Then
  Σ_{0<q_k≤Z^m, ψ_k^*≠1} |S_k(Z^{n1};W1) S_k(Z^{n2};W2)|² ≪ Z^{m+q+ε},
  "no lower conductor bound or restriction on n1+n2". Corollary (D.4): Σ_{q_k≤H}|S_k(N)|⁴ ≪ H^{1+ε}.
- **(3.2)/(3.4) and the sixth-power averaging.** The paper uses Σ_{q_u≤H}|M_u(D)|² ≪ H(HD)^ε for H ≥ D^{1+ϑ}, plus
  the reconstruction M_u = Σ_{d|rad a} μ(d)φ_u(d)q_d^{-1/2} M_{ua⁶}(D/q_d), which gives U^{1/6}H^{5/6}. This
  is the same 5/6 slope as our (B0.9). There is no new input toward (★).
- **Lemma 4.1 (Gram).** The bound is (1 + P_a^{1/6} + P_a²/Y)Z^ε, with exceptional frequencies (k) = e·a⁶. This is
  identical to the record's B0.22 exceptional term.
- **Compensating prime.** It inserts 1_{p | c n³} on the whole theta index. "The factor Z^{ℓ/6} … is the gain
  from compensation" (line 297). The record's b0 = 11/12 − ℓ/4 already contains this ℓ.
- **Which target it addresses.** D.1 bears on the MP slot cost (6κ) and the MP0 family. However, the record's
  own (B3.40) already proves the MP family with "full moving supports" and moving-radical lengths
  (cycle25 paper.md lines 2358 and 2656). **D.1 is therefore the same zero-slot MP0 already in the record, not a
  new input.** It does not address (★), E_a/ATS, PSR or τ₆.
- **Assessment: no gain.** It is an independent re-derivation of components already in the record, and it gives a
  weaker θ.

### 2.2 Tracker forks, LW56e2 obstruction [O]
These are covered in §1. They add evidence that parameter retuning is exhausted, which is consistent with FINDINGS.
They add no new estimate.

## 3. Unopened leads [S]: titles and URLs from search summaries only (statements unverified)

These are the items most worth opening when arXiv access is available. The relevance column is my inference
from the search summary alone.

| Lead | What the search summary said | Target | Likely use |
|---|---|---|---|
| Bröker–Hoffstein, *Fourier coefficients of sextic theta series*, Math. Comp. 85 (2016); arXiv:1312.0568 https://arxiv.org/abs/1312.0568 | partially proven, partially conjectured values of τ(π)², τ(π²), τ(π⁴) for primes π on the 6-fold cover of GL2; n ≥ 4 coefficients "very mysterious" | τ₆ non-correlation | Only prime-power data. Squarefree τ₆(m) (which we need) is apparently not given, matching our theta_id experiment (ρ₆ is not a simple Gauss expression). Low. |
| Chinta–Friedberg–Hoffstein, *Double Dirichlet series and theta functions* https://chinta.ccny.cuny.edu/publ/patterson.pdf | conjecture for sextic theta coefficients via a non-Eulerian double Dirichlet series identity; proved over rational function fields | τ₆ | The conjectured relation (sextic ↔ cubic/quadratic) could give a heuristic model for τ₆ to test PSR against. It is conjectural over number fields. Medium for heuristics, nil for proof. |
| Friedberg–Ginzburg, *Descent and theta functions for metaplectic groups*, JEMS https://ems.press/content/serial-article-files/32310 | descent integrals conjecturally realize quartic/sextic theta; some descent coefficients arithmetic for infinitely many primes | τ₆ | No explicit τ₆ at squarefree indices. Low. |
| Stable wave front set of theta representations, arXiv:2411.02073 | says GL2 theta coefficients for n ≥ 4 remain open | τ₆ | Confirms that τ₆ is open. |
| Dunn, *Metaplectic cusp forms and the large sieve*, ANT 19:9 (2025); arXiv:2403.13151 | new large sieve for a bilinear form with kernel ρ_f (cubic metaplectic cusp-form coefficients), proved by the Browning–Vishe circle method, Voronoi and Gauss–Ramanujan sums; linear estimate beats level 2/3 | E_a/ATS, τ₆ non-correlation | The closest template for "bilinear forms with metaplectic-coefficient kernels". Cubic and cuspidal only (no pole), so it is not directly applicable to sextic theta. **Highest-priority lead to open.** |
| Castillo–de Faveri–Dunn, *Non-vanishing for quartic Hecke L-functions…*, arXiv:2604.01316 (Apr 2026) | has a section on "Twisted sums of quartic Gauss sums" and a mollifier | E_a (n=4 analogue) | It may contain angular-twisted quartic Gauss-sum bounds. n=4 uses the sharp quadratic sieve (DDHL trick), which does not transfer to n=6 (literatureE §5). Low–medium. |
| Hamdar, *Hecke L-functions away from the central line*, arXiv:2509.06152 | first moment of cubic Hecke L over Q(√−3) across the strip; secondary term Q^{5/6} at s=1/2; metaplectic input "only known for cubic, not ℓ ≥ 4" | (★) heuristics | It confirms the Q^{5/6} secondary term (the Kubota pole) in the cubic family. It is a moment of L, not of 1/L. Low. |
| Gao–Zhao, *Mean squares of quadratic twists of the Möbius function*, JNT 247 (2023); arXiv:2212.09241; and *On the moments of averages of quadratic twists of the Möbius function*, arXiv:2405.18778 | asymptotic for Σ_{d≤X}^♭ |Σ_{n≤Y} μ(n)(8d/n)|²; the follow-up gives Y ≪ X^{1/3−ε} unconditionally, and Y ≪ X^{2/(3k)−ε} for the k-th moment | (★) (quadratic analogue) | This is the quadratic analogue of (★), but in the short range Y < X, whereas (★) needs D > H. Not applicable to the long range. Low. |
| Cossaboom, *Zero density estimates for GL2 automorphic L-functions*, arXiv:2609.27715 (Sep 2026) | Guth–Maynard large values with ℓ^q (q ≥ 8/3) coefficient control instead of ℓ^∞; Dirichlet polynomials of length up to T | (★) / R density | The ℓ^q relaxation is the type of GM extension needed for families with non-ℓ^∞ coefficients. But GM-type exponents (≈2.25 at σ≈.735) do not beat the record's 2.225 (FINDINGS "Literature calibration"). Low. |
| Dubon, *Zero-density concentration for Dirichlet polynomials*, arXiv:2609.17875 | deterministic concentration of the zeros of finite Dirichlet sums | — | Irrelevant. |
| arXiv:2604.23900, *Non-vanishing of cubic twists of GL_n(Q) L-functions* | uses μ on Z[ω] with cubic symbols | (★)-adjacent | Abstract not retrieved. Unknown. |

**Not found anywhere (searched, no hits):**
- an explicit formula for squarefree sextic theta coefficients;
- the residue of the n-th order Kubota series at 1/2+1/n for n ≥ 4 (still described as open, consistent with DIPP);
- a sextic or n-th order large sieve with angular-twisted coefficients (ATS);
- any improvement of BGL's (MN)^{2/3} for n ≥ 4;
- zero-density estimates for sextic Hecke families near σ ≈ 0.73–0.84;
- long-range mean values of μ twisted by n-th power residue symbols.

## 4. Assessment against FINDINGS targets

| Target | New input found? |
|---|---|
| (★) mean value λ < 5/6 | No. The only adjacent work is quadratic and short-range (Gao–Zhao [S]), or GM-type and too weak (Cossaboom [S]). |
| E_a / ATS | No. Dunn 2403.13151 [S] is the best template (a bilinear large sieve with a metaplectic kernel), but it is cubic and cuspidal. |
| PSR | No. |
| τ₆ non-correlation | No. τ₆ at squarefree indices remains unknown (Bröker–Hoffstein and CFH give only prime-power data or conjectures [S]). |
| MP slot cost 6κ | No. Collet's D.1 [O] is the same MP0 family already in the record's B3.40. |

Net: there is no published or preprint input (to 2026-10-11) that lowers θ. The record is unchallenged. The
recommended next reads, once arXiv is reachable, are Dunn 2403.13151 (bilinear metaplectic large sieve), then
Castillo–de Faveri–Dunn 2604.01316 (twisted quartic Gauss sums), then CFH (a heuristic τ₆ model to test PSR).
