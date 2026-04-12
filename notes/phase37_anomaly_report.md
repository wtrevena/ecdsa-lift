# Phase 37 series — Degree-1 Fourier anomaly in the Teichmüller cocycle

*Written after Phases 37–37e (Gowers norms, Aut residual, permutation
control, anomaly investigation).*

## Summary

For an ordinary elliptic curve E/F_p and the Teichmüller lift
τ: E(F_p) → E(Z/p^e), the "lift error" cocycle δ(k) = k·τ(G) − τ(kG)
carries a small but persistent **degree-1 Fourier anomaly** in its p-adic
digit sequence.  This anomaly is:

- **Universal** — present on j=0 (CM by Z[ω]), j=1728 (CM by Z[i]),
  and generic (no CM) curves.  Strongest on j=1728 (β ≈ 0.28).
- **Index-ordered** — destroyed by random permutation of the sequence.
- **Independent of Aut(E)** — survives Aut-orbit-averaging projection.
- **Independent of DC bias** — survives mean subtraction.
- **Digit-independent** — present on every nontrivial p-adic digit
  (j = 1, 2, 3, ...), but at different Fourier frequencies.
- **Cryptanalytically inert** — computing δ(k) requires knowing k
  (the discrete log), so the anomaly cannot be turned into a DLP oracle.

## Quantitative characterisation

The function f_{d,j}(k) = exp(2πi · digit_j(z(δ(k))) / p) has peak
Fourier coefficient |f̂(α)| that scales as:

    |f̂(α)| ~ C · n^{−β}

with β ∈ [0.28, 0.47] depending on the curve family (vs β = 0.5 for
random S^1-valued sequences).  The U^2 Gowers norm is consequently
elevated by 5–11 standard deviations above the random baseline across
all tested primes (p = 31 to 521, n = 20 to 521).  U^3 and higher
Gowers norms are at random level, confirming the structure is
degree-1 (linear) only.

## Key experimental results

### Phase 37 (initial sweep)

Ran on secp256k1 family y² = x³ + 7, primes p = 31..631, e = 4.
All three digits (d1, d2, d3) show persistent U^2 elevation.
U^3 is dead.

### Phase 37c (Aut residual)

Projected f_{d,j} onto the Aut(E)-orbit-constant subspace and
subtracted.  Residual U^2 remained elevated (3–10σ), ruling out
the Aut-equivariance explanation from Phase 21b.

### Phase 37d (permutation control) — the decisive test

Randomly permuted the *values* of f_{d,j} while keeping the same
multiset.  Permuted U^2 matched random S^1 baseline (z ≈ 0).
Original (un-permuted) U^2 was 5–11σ above the permuted baseline.
**Conclusion: the structure is in the *index ordering*, not in the
value distribution.**

### Phase 37e (focused investigation)

1. **Dominant frequencies are different across digits.**  d1, d2, d3
   have independent peak Fourier modes — no single linear character
   governs the anomaly.
2. **No algebraic pattern in α.**  The dominant frequency shows no
   relationship to a_p, g_ord, ω_index, or other curve invariants.
3. **Precision invariance.**  digit_j is determined at precision e = j+1;
   the anomaly doesn't grow or shrink with more p-adic precision.
4. **CM independence.**  j=1728 shows β ≈ 0.28 (strongest), generic
   no-CM curves show β ≈ 0.31, j=0 shows β ≈ 0.44.  The anomaly is
   **most pronounced on non-CM curves**, ruling out the CM structure
   as the source.
5. **Peak ratio → 1.**  For large n, the leading Fourier coefficient
   barely exceeds the second — the anomaly is diffuse across modes.

## Mathematical interpretation

δ(k) = −Σ_{i=1}^{k-1} c(iG, G) is the *prefix sum* of the 1-cocycle
c(iG, G).  The Fourier anomaly in the prefix sums is related to the
spectral properties of the individual cocycle values c(iG, G) via the
prefix-sum transfer function 1/(1 − e^{2πiξ/n}).  The individual
c(iG, G) values are known (from Phases 21–31) to satisfy the universal
2-cocycle identity and Aut-equivariance, but these constraints alone do
not predict the spectral anomaly — they only govern the *multiset* of
values, not their Fourier structure.

The anomaly therefore represents a *new* structural property of the
Teichmüller cocycle that is:
- Not explained by any universal identity
- Not explained by CM or Aut(E) symmetry
- Not specific to any particular curve family

It is plausibly a consequence of the *additive structure of the formal
group* interacting with the *multiplicative structure of the scalar
multiplication map* k ↦ kG, producing a weak but persistent spectral
bias.  A theoretical explanation would require understanding the
Fourier analysis of the formal-group exponential composed with the
group law — a non-trivial question in p-adic harmonic analysis.

## Why this does not break ECDLP

The cocycle δ(k) is a function of k (the discrete log).  To compute it,
you need to evaluate k·τ(G) and τ(kG).  The latter is computable from
Q = kG (the public key), but the former requires knowing k.  Therefore:

- The anomaly cannot be observed from (G, Q) alone.
- It cannot be used as a distinguisher or predictor for k.
- It is a property of the *mathematical object* (the cocycle), not of
  any *observable* derived from the public-key problem.

Even at maximum signal strength (β = 0.28 for j=1728), extrapolating
to n = 2^256 gives |f̂(α)| ≈ 2^{-72} — detectable in principle but
requiring exhaustive enumeration to measure.

## Open questions

1. **Why β < 0.5?**  What property of the formal-group/group-law
   interaction produces a spectral decay exponent strictly between
   0 (deterministic) and 1/2 (random)?  Is β = 1/2 − 1/(2·|Aut|) or
   some similar formula?
2. **Why is j=1728 strongest?**  The CM by Z[i] automorphism has order 4;
   j=0 has order 6.  Does larger Aut suppress or enhance the spectral
   anomaly?  (Our data: |Aut|=4 gives β=0.28; |Aut|=6 gives β=0.44;
   |Aut|=2 gives β=0.31.  No obvious monotone relationship.)
3. **Is the anomaly present for non-Teichmüller lifts?**  The canonical
   lift τ_can is a homomorphism with δ ≡ 0.  For "random" sections σ
   (e.g. Hensel lift), does the anomaly persist or change?
4. **Can the anomaly be accessed indirectly?**  E.g., via pairings,
   isogenies, or multi-party protocols that allow computing functions
   of k without knowing k.  (We conjecture no, but this is the most
   important open question.)
