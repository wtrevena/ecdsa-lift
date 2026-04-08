# Graveyard of dead routes

If a route fails, **write down here** why it failed, with enough
specificity that we don't accidentally try it again in a different
disguise three weeks from now.

## D1 — Naive isogeny-graph CM accumulation
**Claim attempted.** Walk the isogeny graph from secp256k1, accumulate
one extra linear relation on the secret scalar per isogeny, reduce the
DLP to a short-vector problem in a CM lattice.

**Why it's dead.** Deuring's theorem. For an ordinary elliptic curve
with CM by the maximal order `O_K`, every `ℓ`-isogeny can only map to
a curve whose endomorphism ring is either the same maximal order or
the non-maximal suborder `Z + ℓO_K`. You never *gain* endomorphisms
by walking — at best you stay put, at worst you lose some. The
"linear relations accumulate" story was wrong: there was only ever
one CM relation, and it was already there at step 0.

**Escape routes that survive.** (a) Cross into the supersingular
isogeny graph. This is a fundamentally different object and
genuinely adds endomorphisms, but connecting ordinary to supersingular
requires leaving the domain of the naive ℓ-isogeny walks.
(b) Pairing lifts to higher-genus Jacobians — `secp256k1` embeds into
`Jac(C)` for `C` a cover, which might carry extra structure.

## D2 — Teichmüller-lift `δ(k)` as a recoverable encoding
**Claim attempted.** The Teichmüller lift error
`δ(k) = [k]G̃ - τ([k]G)` at precision `p^4` is injective in `k`,
so maybe it's a *computable* function of `k` (polynomial, linear,
multiplicative, CM-equivariant) and therefore invertible.

**Why it's dead.** It's injective but pseudorandom. Every algebraic
structure we tested failed. Under the Teichmüller lift specifically,
inverting `k ↦ δ(k)` requires a lookup table and is equivalent to
brute force.

**Escape routes that survive.** It's the *lift* that's pseudorandom,
not the `p`-adic approach. A *homomorphic* section (Phase 4) would
replace the Teichmüller scrambling with a linear structure in one
move.

## D4 — Homomorphic section does NOT directly give the attack
**Claim attempted.** Construct a homomorphic section
`s : E(F_p) → E(Z/p^e Z)` via the cohomological retraction
`s(g) := τ(g) - μ([n]·τ(g))` where `τ` is the Teichmüller lift and
`μ` is the inverse of multiplication-by-`n` on the kernel of
reduction `A = Ê(pZ/p^e)`. Then use the homomorphism to recover the
discrete log by, e.g., applying `[n]` and taking the formal logarithm.

**Why the construction works.** The group-cohomology argument
`H^2(E(F_p), Ê(pZ_p)) = 0` when `gcd(n, p) = 1` guarantees the
short exact sequence splits, and the explicit retraction above
provides the splitting. Implemented in `src/lifts.py::
splitting_section_proj`, verified as a homomorphism on the
cases where the projective arithmetic is robust (p = 13, 31, 67
with prime cyclic n, precision `p^3`, 100% of pair-checks pass).

**Why it nonetheless doesn't give an attack.** If `s` is a
homomorphism, then
```
[n] · s(P) = s([n] · P) = s(O) = O    for every P,
```
so `[n] · s(P)` is the identity in `E(Z_p)` for every `P` and the
formal-logarithm trick yields zero regardless of the secret. The
whole point of Smart–SSSA on anomalous curves is that `n = p`, so
`[n] · τ(P)` lands in the formal group (order `p^{e-1}`) in a way
that depends linearly on `P` via the non-homomorphic correction of
`τ`. A *homomorphic* `s` has no such correction, which is exactly
why it's useless for this attack vector.

**What it still buys us.** The explicit `s` is, essentially, a
*computable canonical lift* from `E(F_p)` to `E(Z_p)` (the
analogue of the Serre–Tate canonical lift for `j = 0` curves).
This is a genuinely useful object — for Hida theory, `p`-adic
modular forms, height pairings — but none of these are directly
cryptanalytic. The attack would have to come from a different
use of the lift, such as:
  (a) Exploiting algebraic relations between `s(P)` for different
      `P` via the CM endomorphism (the `λ`-eigenspace of the GLV
      endomorphism acts on `Im(s)` and might linearize something).
  (b) Pairing `s(P)` into a higher-dimensional object (abelian
      variety, Jacobian of a cover) where a different DLP holds.
  (c) Using `s` to populate a lattice in `Q_p^k` on which LLL /
      BKZ might find short vectors related to the secret.

None of these are free. They all require new ideas.

## D3 — Smart–SSSA for non-anomalous curves
**Claim attempted.** Multiply a lifted point by `n = |E(F_p)|` to
push it into the formal group `Ê(pZ_p)`, apply the formal log,
recover `k` from the ratio.

**Why it's dead for non-anomalous curves.** When `n = p` (anomalous),
the ambiguity from the choice of lift vanishes mod `p` and the method
is clean. When `gcd(n, p) = 1`, the ambiguity is a unit in `Z_p` times
an arbitrary formal-group element, and different lifts give different
answers. The method is only well-defined up to that ambiguity, and
the ambiguity swamps the signal.

**Escape route that survives.** If we had a *canonical*, homomorphic
choice of lift, the ambiguity disappears by construction — because
any two homomorphic sections differ by a homomorphism `E(F_p) → Ê`,
and that homomorphism is forced to be zero when `gcd(n, p) = 1`. So
Phase 4 doesn't just solve its own problem, it also rescues Smart–SSSA.
