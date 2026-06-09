# Phase 40 — Resolution of the Phase 37 Fourier anomaly

*Written after Phases 40 and 40b (the decisive controls).*

## The question

Phase 37 found a persistent degree-1 (`U²`) Fourier anomaly in the
p-adic digit sequence of the Teichmüller lift error
`s(k) = z(δ(k)) = z([k]τ(G) − τ(kG))`, elevated 5–11σ over a random
`S¹` baseline, universal across CM/non-CM curves. The Phase 37 report
called it a *new structural property* "not explained by any universal
identity." That claim was never actually controlled against the two
most likely mundane explanations. Phase 40/40b supply those controls.

## Control 1 (Phase 40) — is it a prefix-sum artifact?

`s` is the cumulative sum (mod N) of its own first differences
`d(k) = s(k+1) − s(k)`. Prefix-summing is a low-pass filter
(`ŝ(ξ) = d̂(ξ)/(1−e^{2πiξ/n})`, a `1/ξ` blow-up at low frequency), so a
generic concern is that the `U²` elevation is just an artifact of
summing *any* increments with this multiset.

Test: shuffle `d`, rebuild `s'` by prefix-summing in the new order,
re-measure `U²`. Result across p = 31…521:

- `z(incShuffle vs random) ≈ 0` everywhere (−0.26 … +1.04).
  Shuffled-increment prefix sums sit **at** the random baseline.
- `z(real vs incShuffle)` is large and **grows with n** (≈ +2 at p=67
  to ≈ +10–12 at p=521).

**Verdict:** NOT a prefix-sum artifact. The digit-`j` nonlinearity
(extract digit `j`, exponentiate) destroys the random walk's
low-frequency mass, so a shuffled walk looks random. The signal lives
in the *specific ordering* of the cocycle increments, and its
significance grows with data — the signature of a real effect.

## Control 2 (Phase 40b) — is it just the Phase 21b antisymmetry?

The one universal identity that constrains the ordering is Phase 21b:

    δ(k) + δ(n−k) = [n]·τ(G).

Verified here to be **exact in the z-coordinate**: on p=97,
`s(k) + s(n−k) ≡ C (mod N)` takes a *single* value across all 78
pairs (carries included, not just leading order). A rigid reflection
about `k = n/2` couples Fourier mode `α` with `−α`, which is exactly
what elevates `U²`.

Test: restrict to the first half `k = 1…⌊(n−1)/2⌋`. Within one half the
reflection maps each `k` *outside* the window, so it imposes no internal
constraint. Re-run the Phase 40 increment-shuffle control on the half.

Result across p = 67…521:

- `z(real vs incShuffle)` collapses to noise (−2.44 … +1.57), **no
  growth with n**. At p=257 the full sequence (n=256) gave z ≈ +8; the
  half (n=128) gives z ≈ +0.1–0.6. Same ballpark n, opposite result —
  not a statistical-power loss.

**Verdict:** the entire Phase 37 anomaly is the Fourier-domain shadow
of the Phase 21b antisymmetry. Remove the `k ↔ n−k` reflection (by
halving, or by shuffling increments across the full range) and the
`U²` elevation vanishes completely.

## Consequence

The Phase 37 "degree-1 Fourier anomaly" is **real but not new**. It is
not an independent structural property of the Teichmüller cocycle; it
is the spectral signature of a single, exact, *universal* identity
(`δ(k) + δ(n−k) = [n]τ(G)`) that was already known, already documented
(Phase 21b), already understood to be a **public** consequence of the
curve's geometry, and already shown to give a mere 2-to-1 search
reduction with **zero** cryptanalytic leverage.

This corrects the Phase 37 report's three central claims:
- "Not explained by any universal identity" → **false**; it is exactly
  the Phase 21b identity.
- The varying decay exponent `β ∈ [0.28, 0.47]` "depending on curve
  family" → consistent with finite-size noise around a single
  reflection-induced spectrum; no curve-family-specific structure
  survives the half-sequence control.
- "Most pronounced on non-CM curves, ruling out CM as the source" →
  irrelevant; the source is the antisymmetry, which is universal
  (holds for all `j`, confirmed at 40-bit non-CM in Phase 28).

## Status of the project

The last open lead is closed. Every structural probe across Phases
1–40 reduces to one of:
1. a homomorphism (canonical lift) → transports DLP without weakening,
2. an Aut(E)-equivariance → public, `(g−1)/|Aut|` rank deficit,
3. the antisymmetry `δ(k)+δ(n−k)=[n]τ(G)` → public, 2-to-1 only,
4. or certified pseudorandomness (DFT, Berlekamp–Massey, Mahler,
   valuation, multiplicative/additive fits).

No efficiently-computable-from-`Q` functional on the lift depends
linearly on `k`. The p-adic lifting attack surface on ordinary
elliptic curves is inert. This is the clean, complete negative result
the paper should report — now with **no unexplained anomaly** left
dangling.
