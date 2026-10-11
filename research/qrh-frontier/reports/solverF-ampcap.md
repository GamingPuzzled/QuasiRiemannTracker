# Solver F — the amplitude cap p ≤ 1/2 in (B0.37)

**Verdict: NO IMPROVEMENT. The cap is forced.** It is the explicit-formula size of a prime sum when a zero
of real part a = (1+δ)/2 is present. The bin's own defining zero attains it, so the binding class p = 1/2 is
the generic "single zero saturates everything" configuration. It is not a lossy estimate. θ is unchanged.

Scripts: `scripts/ampcap/slack2_ampcap.py` (an unmodified copy of slack2.py) and `scripts/ampcap/pmax_scan.py`.

## 1. Where p ≤ 1/2 comes from

- **Bins (§B0.2, B0.13).** A nonfloor bin row u has a zero of some member of 𝒳(u) = {νψ_u^{±1}} with
  real part σ ∈ [a, a+e), and no zero with real part ≥ a+2e up to height 3^{i+1}T. Here δ = 2a−1.
- **What p measures (§B0.5, B0.37).** Each physical main slot of (B0.32)–(B0.33) contains the centrally
  normalised prime sum P_i^{-1/2} Σ_𝔭 ψ̄_u(𝔭) w_i(q_𝔭/P_i)(q_𝔭/P_i)^{it}, with P_i = Z^{ℓ_i}. Its size
  is P_i^{g_i}. The quantity q = L^{-1} Σ_i (ℓ_i/d) g_i is the length-weighted mean amplitude exponent, so the
  slot product has size Z^{ℓq}. Then p = q/δ.
- **Why the cap is 1/2.** The cap comes from (B0.16): |P^{-1/2} Σ ψ̄_u(𝔭) w (q_𝔭/P)^{it}| ≪ P^{δ/2+O(e)+ε}.
  This bound follows from −L'/L, the zero-free region right of a+2e, and Borel–Carathéodory (B0.14).
  Hence g_i ≤ δ/2, which gives q ≤ δ/2 and p ≤ 1/2.
- **How p enters F (B0.46/B0.47).** The parameter p enters in three places:
  - +ℓq: the slot product is paid pointwise in (B0.41).
  - −q·min(d, 2y−d)/(12κ): the Hölder/MP mark gain qz_A/2.
  - R_κ(δ,p): marks reduce the row count through the MI and MP lines (B0.38) and C_κ(p).

## 2. Size of the lever (pmax_scan.py)

The derivative dF/dp at the binding point (δ0 = .469345, κ0 = .749773, h0 = .793948) splits as follows:

| piece | value |
|---|---|
| amplitude cost ℓ0·δ0 | +0.07844 |
| numerator mark gain −δ0(h0−4ℓ0)/(12κ0) | −0.00655 |
| row count (3h0/4)∂R/∂p | −0.03947 |
| **total** | **+0.03242** (matches the finite difference) |

The total is positive, so the worst case sits at the cap. Re-optimising the geometry (x, ℓ) with p_max moved
gives these results:

| p_max | b0 | gain |
|---|---|---|
| 0.51 | 0.8749595500 | −7.3e-5 |
| 0.50 | 0.8748866157 | 0 (record, up to grid error) |
| 0.499 | 0.8748791367 | 7.5e-6 |
| 0.495 | 0.8748488706 | 3.8e-5 |
| 0.49 | 0.8748101951 | 7.6e-5 |
| 0.45 | 0.8744689784 | 4.2e-4 |
| 0.40 | 0.8739600099 | 9.3e-4 |

So dθ/dp_max ≈ 0.0075. Reaching 1e-6 would need Δp ≈ 1.3e-4, i.e. every slot exponent g_i at most
δ/2 − 6e-5 at the binding point. Reaching 3e-5 would need Δp ≈ 4e-3 (Δg ≈ 1.9e-3). Any such gain
must be a power saving that holds uniformly over the binding class.

## 3. Why no alternative estimate helps (each route checked)

**(a) The pointwise bound (B0.16) is sharp.** Take a zero ρ = σ + iγ with σ ∈ [a, a+e) of the slot member,
and let |γ − t| = O(1). The explicit formula then gives a term ŵ_i(ρ − 1/2 − it)·P_i^{ρ−1/2}, of size
≍ P_i^{δ/2}, at every slot length simultaneously. Any power saving below δ/2 is therefore impossible for
such rows.

What remains heuristic: I did not prove that no cancellation among zeros occurs. This does not matter
for the verdict. Upper bounds cannot be improved for rows where one zero dominates, and nothing in the
framework excludes such rows.

**(b) No consistency constraint excludes p = 1/2 rows from the bin count.** The same zero that puts u
in bin δ produces the Lemma 1W witnesses and saturates the slots. So the class p = 1/2 is the natural
class, not an extreme one. Excluding it would require counting rows with a zero at a better than U^{R}. That
is the zero-density lever (R, V0, and (★) in FINDINGS), not the cap.

There is a minor caveat. If the bin zero belongs to ηψ̄_u but not to the slot member ψ̄_u, then p can be
smaller. The worst case, with all members having zeros at a, cannot be ruled out rowwise.

**(c) Treating the extreme-amplitude rows with the pointwise numerator gives nothing.** For rows with p near
1/2 one could replace Hölder (B0.41) by the pointwise numerator (B0.15). FINDINGS tick 1 (pointwise.py)
already took the minimum of the two over all (δ, p). The gain was 7e-14, because the two bounds agree to
1e-4 at the binding point and the saddle kills any off-point gain.

**(d) Using the saturated slots as stronger zero detectors fails.** At p = 1/2 a mark of length z has
|·|² ≥ U^{δz}, the same strength per unit length as the zero witnesses. Three ways to use this all fail:

- **More supply does not help.** Mark capacity is limited by MI, with r + 2z ≤ 1 and z ≤ 4/21 at the
  crossing, not by supply. The supply L = ℓ0/h0 = .2105 already exceeds 4/21 (§B0.5: "if L > 4/21 the
  supply suffices"). Extra amplitude cannot be converted into count.
- **Trading inverse length for marks loses.** Each unit of r traded for marks costs δ and gains q = δ/2,
  which is a net loss. This is the positive slope δ − q of A_I.
- **Marks as the sole witness are weaker.**
  - MP with m = 0 gives count U^{1−δL} ≈ U^{0.90}.
  - The k-th power of the slot product of length ≤ U under a general sextic large sieve
    (M+N+(MN)^{2/3}) gives U^{4/3−δ} ≈ U^{0.86}.
  - Even the conjectural de Faveri–Dunn–Hoffstein Ξ_6 at A = B = U gives U^{7/6−δ} ≈ U^{0.70}.
  - All of these exceed R = 0.590.

  Only a sharp (A+B) sextic sieve for prime products of length U would give 1−δ = 0.53. That sieve is
  false for generic coefficients according to the Ξ_6 conjecture, and unproved for prime products.

**(e) Averaging over heights gives no power saving.** Saturation needs Im z within O(1) of a zero ordinate.
However, the Gaussian e^{(s+z−1)²} and the Mellin decay confine the effective heights to (log)-ranges, and
T = U^{τ} with τ arbitrarily small. Restricting heights therefore saves only U^{o(1)}.

**(f) Moving Re z for the slots gives no first-order gain.** The variable z is shared by ζ(6z), q_u^{-z}, X^{1/2−z}
and the slots. The slot amplitude is linear in Re z, and z0 cancels in F at d = h. So moving Re z gives
no first-order gain at the binding point.

**(g) Changing the slot design does not help.** Prime sets with sector or Hecke-angle restrictions, other
weights w_i, or thinner prime sets all obey (B0.16) with the same δ/2. This holds after enlarging Θ by the
angular characters, and their principal signal changes only by constants. No exponent changes.

## 4. Where the p-dependence could actually move

The ℓq cost (+0.078) is only partly offset (−0.046) by the gains from marks. The gains are limited by MP slot
capacity, with z_A = (1−2m)/(6κ), which is the 6κ amplifier (Solver A, mpcost.py), and by MI capacity, with
z ≤ (1−r)/2 (Rgen.py: −0.008 per unit). Improving the p-trade-off is therefore equivalent to those
already-priced levers. It is not a separate lever.

## Rigor status
- Proved / verified:
  - The source chain (B0.13 → B0.14 → B0.16 → B0.37).
  - The F-decomposition: the finite difference matches the analytic pieces exactly.
  - The p_max scan (same optimiser as slack2; grid error ≈ 3e-8).
- Heuristic, but not needed for the verdict:
  - Sharpness of (B0.16), meaning no systematic cancellation among zeros at aligned heights.
  - The "generic saturation" picture in (b).
- Not established: no improvement is claimed.
