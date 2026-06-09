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
(But see D4: the homomorphic section *does* exist, and it does
*not* yield the attack.)

## D5 — CM-endomorphism lift error `ψ(P) = [λ]·τ(P) − τ(λ·P)`
**Claim attempted.** The GLV endomorphism `φ(x, y) = (ζ x, y)` is a
Z-linear map on `E(F_p)` with eigenvalue `λ`. Lift `φ` naively by
applying `[λ]` on the lifted side and compare to the Teichmüller
lift of `[λ]·P`. The error `ψ(P)` lands in the formal group. Ask
whether `k ↦ ψ([k]G)` is additive (a homomorphism `Z/n → Ê`) or
CM-equivariant (`[λ]·ψ(P) = ψ([λ]·P)`).

**Why it's dead.** Phase 5 measured both. On every j = 0 test
curve (`p ∈ {31, 43, 67, 73, 79, 97, 103}`, precision `p^3`):
- `ψ` is injective on `k` (good).
- Additivity fails hard: ≤ 24/576 pairs are consistent, and the
  first failure is already at `(j, k) = (1, 2)`.
- CM-equivariance fails completely: 0/39.
So `ψ` is just another pseudorandom injection into `Ê`, no better
than the Phase 1 `δ`.

## D6 — Hidden Number Problem on the Teichmüller cocycle
**Claim attempted.** Tabulate the formal-group parameter
`z_k = z(τ(G) + τ([k]G) − τ(G + [k]G))` of the Teichmüller
cocycle, fit `z_k ≈ k · z_1 + small` modulo `p^e`, and run an
HNP-style LLL reduction on `(k_i, z_{k_i})` pairs to recover the
secret scalar.

**Why it's dead.** Phase 6 shows `z_k` is not linear in `k`: the
"residues" `z_k − k · z_1 mod p^e` are uniformly distributed in
`[0, p^e)` with average ratio-to-`p^e` near 0.5 across every tested
curve, and LLL on the natural HNP lattice finds only short vectors
that are just modular reductions of the input rows — no leakage.

## D7 — MOV embedding degree for ordinary j = 0 curves
**Claim attempted.** For the test curves, compute the embedding
degree `k = min{k : n | p^k − 1}` and hope it's small enough to
make the DLP in `F_{p^k}` tractable.

**Why it's dead.** Phase 7 measured `k` for every test curve. It's
generically large (order of `n`, not small). For `secp256k1` it is
astronomical — this is precisely why MOV was ruled out as a threat
when the curve was chosen. The j = 0 CM structure does NOT reduce
the embedding degree. There's also no genus-2 cover with a smaller
embedding degree (Galbraith et al.). Known.

## D8 — Polynomial fit of `δ(k)` mod `p^e`
**Claim attempted.** Even if finite differences of `δ(k)` don't
vanish (Phase 2), maybe a Lagrange-interpolated polynomial of
modest degree fits `δ(k)` globally mod `p^e`. If yes, one can
evaluate the polynomial at an unknown `k` and read off the secret.

**Why it's dead.** Phase 8 tried degrees 1 through 7 on every test
curve via modular Vandermonde. The fit count is always exactly
`degree + 1` — only the points used to solve the system. No
polynomial of degree ≤ 7 describes `δ(k)` on any curve. Given
Phase 2 showed even order-4 finite differences don't vanish, this
closes the polynomial-fit family.

## D9 — LLL on the lift coordinate lattice
**Claim attempted.** Put `(k_i, X_i, Y_i, Z_i)` rows with projective
`τ(k_i G)` coordinates into a lattice, add modular rows scaled by
`p^e`, and LLL-reduce to find a short vector relating `k_i` to the
lift coordinates.

**Why it's dead.** The lattice is rank-deficient (many rows in 4
columns), so LLL trivially finds a zero vector as its shortest
element. Any non-trivial short vector it finds is a `p^e`-multiple
of an input — i.e. a reduction artifact, not a relation involving
`k`. Phase 9 confirmed this on every curve. No linear structure
between `k` and the projective lift coordinates is visible to LLL
at our lattice rank.

## D10 — Ordinary ↔ supersingular isogeny bridge
**Claim attempted.** Find an isogeny from the target ordinary curve
to a supersingular curve, exploit the latter's maximal endomorphism
ring to solve DLP, and transfer back.

**Why it's dead (textbook).** Trace of Frobenius is an isogeny
invariant mod `p`. Ordinary curves have `t ≢ 0 (mod p)`,
supersingular have `t ≡ 0 (mod p)`. No F̄_p-rational isogeny
connects them; they sit in different connected components of the
`ℓ`-isogeny graph. No correspondence via a higher-genus curve
crosses the partition either, since the Jacobian of any rational
cover of an ordinary curve is itself ordinary. Phase 10 verified
`t ≢ 0 (mod p)` for every test curve.

## Synthesis across Phases 1–10

The p-adic lifting attack surface on j = 0 ordinary curves has been
comprehensively mapped and found inert under every structural lens
we could aim at it:

| Phase | Target                                      | Verdict |
|-------|---------------------------------------------|---------|
| 1     | Injectivity of `δ(k)`                       | Yes, but pseudorandom |
| 2     | Polynomial / additive structure of `δ`      | None detectable |
| 3     | Smart–SSSA on non-anomalous curves          | Blocked (D3) |
| 3.5   | Cohomological existence of a homomorphic `s`| Yes, and constructive |
| 4     | Homomorphic section → ECDLP attack          | Construction works, attack doesn't (D4) |
| 5     | CM-endomorphism lift error `ψ`              | Pseudorandom (D5) |
| 6     | HNP on cocycle `z_k`                        | No bias (D6) |
| 7     | MOV on j = 0 curves                         | Not applicable (D7) |
| 8     | Low-degree polynomial fit of `δ`            | No fit (D8) |
| 9     | LLL on lift coordinates                     | Rank-deficient / no leak (D9) |
| 10    | Ordinary ↔ supersingular bridge             | Impossible (D10) |

**What the session actually produced.** Not an ECDSA break. A
computable canonical lift (Phase 4) — essentially a j = 0 analogue
of the Serre–Tate construction, with explicit projective-coordinate
arithmetic over `Z/p^e`. That's a legitimate number-theoretic
object of independent interest (p-adic modular forms, height
pairings, overconvergent cohomology), just not a cryptanalytic
weapon.

**What would have to be true for p-adic lifting to break ECDSA.**
Some functional on the canonical lift would have to depend
*linearly* on the scalar `k` modulo `p^{e-1}` while being
*efficiently computable* from the reduction mod `p` alone. Every
natural such functional we tried is either trivially zero (by
homomorphism, D4), or pseudorandom (δ, ψ, cocycle z). The
obstruction isn't a single theorem; it's the accumulated failure
of every explicit construction to distinguish `k · z_1` from
uniform noise in `Z/p^{e-1}`. An ECDSA break from this direction,
if it exists, lies behind a structural insight that none of the
Phase 1–10 probes were sensitive to.

## Phases 11–22: deeper pseudorandomness tests and reformulations

After Phase 10 the "pseudorandom" label was audited: previous
phases only ruled out linear/low-degree-polynomial fits. Phases
11–22 added genuine randomness tests and algebraic reformulations.

| Phase | Target                                              | Verdict |
|-------|-----------------------------------------------------|---------|
| 11    | Genus-2 Jacobian cover (Cantor/Mumford)             | Blocked: compose bug, open |
| 12    | GLV CRT decomposition `k = k₁ + k₂ λ`               | No joint linear fit |
| 13/13b| Weil pairing via Miller w/ auxiliary divisor        | Blocked: affine pivots non-unit, open |
| 14    | Multi-base δ matrix rank in Z/p^e                   | Trivially 0 (every δ ∈ Ê) |
| 14b   | Normalized rank of `δ/p` at each p-adic layer       | **Layer 0 rank-deficient, layers ≥1 full** |
| 15    | Triple cocycle `γ = δ²c`                            | ~zero (as expected) |
| 15b   | p-adic depth of γ                                   | Confirms δ²c = 0 |
| 16    | Frobenius error `F(τ(P)) − τ(P)`                    | Identically zero (sanity) |
| 17    | DFT spectrum of `δ(k)`                              | Flat (peak/mean < log L) |
| 18    | Berlekamp–Massey linear complexity of `δ` digits    | Complexity = L/2 (max) |
| 19    | p-adic valuation distribution of δ                  | Uniform on Ê |
| 20    | Multiplicative structure `δ(jk)` vs `δ(j) δ(k)`     | None |
| 21    | Symmetry `δ(n-k) = ±δ(k)`                           | Wrong form |
| 21b   | **`δ(k) + δ(n-k) = [n] τ(G)`**                      | **100% verified** |
| 22    | δ invariance under isomorphic Weierstrass models    | Not invariant |

### Phase 21b — the first exact algebraic relation in δ

Derivation: `δ(n-k) = [n-k]τ(G) − τ((n-k)G) = [n]τ(G) − [k]τ(G)
+ τ(kG) = [n]τ(G) − δ(k)`. Verified 100% on every curve.
Consequence: `δ` is antisymmetric around `k = n/2` up to the
constant `[n]τ(G) ∈ Ê`. This gives a 2-to-1 reduction of the
scalar search space but does **not** invert δ — pseudorandomness
within each half is preserved under Phases 17, 18, 20.

### Phase 14b — leading-digit rank deficit (fully resolved in 14c/d/f)

Matrix `M[H][k] := δ_H(k) / p mod p`. Layer 0 (the leading
p-adic digit) has rank 5–10 instead of full rank 12 across all
test curves. Layers 1 and 2 are full rank.

**Resolution (Phases 14c → 14d → 14f):**

- **Phase 14c** — scaled rows/cols up to 32. Rank *saturates* on
  every curve far below min(rows, cols). Not a small-matrix
  artifact; genuine low-rank structure.
- **Phase 14d** — hypothesized Aut(E)-orbit explanation. For j = 0
  with **prime** g_ord, verified on every curve:
  - Every orbit of ⟨λ, −1⟩ ⊂ (Z/g)* has size 6.
  - Number of orbits = (g−1)/6 = observed rank exactly.
  - Within each orbit, every row is an F_p-scalar multiple of a
    single representative (kernel basis vectors have exactly two
    nonzero entries).
- **Phase 14f** — extended to j ≠ 0. Universal rule confirmed:

      rank = min(size, (g_ord − 1) / |Aut(E)|)

  holds for **every** curve with prime g_ord across j = 0
  (|Aut| = 6) and generic j (|Aut| = 2). j = 1728 not independently
  verified due to difficulty finding prime-order |E(F_p)| at tested
  primes, but theory predicts rank = (g−1)/4 with |Aut| = 4.

**Conclusion.** The leading-digit rank deficit is **exactly** the
Aut(E) equivariance of the Teichmüller lift. For j = 0 curves the
action is μ₆ and the deficit is dramatic; for generic curves it is
μ₂ = {±1} and the deficit is just a factor 2. The structure is a
property of the public curve, visible to any attacker, and gives
**zero cryptographic advantage**. The Phase 14b "last open lead"
is closed as negative.

### What was certified as pseudorandom (properly, this time)

- `δ(k)` has flat DFT spectrum (Phase 17)
- `δ(k)` digit sequences have maximal linear complexity (Phase 18)
- No multiplicative/XOR/additive collision (Phase 20)
- Uniform valuation distribution on Ê (Phase 19)

### Still open (blocked on implementation, not theory)

- **Phase 11 (genus-2)**: FIXED in `phase11b_genus2_fixed.py` —
  we Hensel-lift each affine point on `C(F_p)` to `C(Z/p^e)` and
  build the Mumford representative from the lifted coordinates,
  guaranteeing `v² ≡ f (mod u)` exactly and avoiding the non-exact
  reduction failure in Cantor's algorithm. Runs cleanly through
  24 Cantor compositions on 5/6 test curves. **Finding**: no
  linear-in-k structure in any Mumford coordinate (u₀, u₁, u₂,
  v₀, v₁). The genus-2 cover adds no new lever.
- **Phase 13 (Weil pairing on the lift)**: attempt
  `phase13c_miller_fraction.py` tracks Miller numerator and
  denominator as separate ring elements to skip divisions. It
  still fails because the affine scalar ladder for `P` of prime
  order `g_ord` necessarily reaches `T ≡ -P (mod p)`, and `T + P`
  is a kernel-of-reduction point that affine addition cannot
  represent. A full fix would rewrite Miller on ProjCurve with
  projective line-evaluation formulas — deferred. Note however
  that for `G` of prime order, the Weil pairing is trivial:
  `e_n(G, kG) = e_n(G, G)^k = 1`. The interesting object would
  be the Tate pairing of `δ(k)` against something, but Phase 26
  already ruled out nontrivial linear structure in
  `log_Ê(δ(k))`, which is the Tate local pairing.

### Phase 26–29 — "key new experiments" (all negative)

- **Phase 26 — Formal log of δ**. Taking `log_Ê(z(δ(k)))` with
  `log_Ê` the formal logarithm (power series from
  `pari.ellformallog`) should linearize the formal group. On
  `y² = x³ + b`, `log_Ê` has support only at `z^{6i+1}` (sparse,
  a direct trace of the μ₆ symmetry). Sampling `k = 1..40` and
  fitting: no linear or low-degree polynomial structure.
- **Phase 27 — Mazur–Tate height**. Bilinearity
  `h_p(kP, kP) = k² h_p(P, P)` verified numerically via
  `pari.ellpadicheight` on `y² = x³ + 3` at primes
  `p = 7, 13, 19`. But the global Mazur–Tate height is defined
  on `E(Q)`, not on `E(F_p)`; the LOCAL contribution at p for
  points in `Ê` is exactly `−log_p(z)`, i.e. Phase 26. Mazur–Tate
  adds nothing beyond the formal log angle.
- **Phase 28 — Non-CM curves at (near-)cryptographic sizes**.
  On random curves with `j ∉ {0, 1728}` and prime group order
  at prime sizes `p` of 20, 24, 28, 32, 40 bits:
  - Leading-digit rank `= min(size, (n-1)/2) = (n-1)/|Aut|` on
    every case. The Aut-orbit rank rule extends to cryptographic
    scale and non-CM curves.
  - Phase 21b z-additive antisymmetry
    `z(δ(k)) + z(δ(n-k)) ≡ z([n]τ(G)) (mod p^e)` holds 8/8 on
    every case. This relation, originally observed only on j=0
    CM curves, is **universal**.
- **Phase 29 — Mahler expansion of δ**. Mahler coefficients
  `a_n = (Δⁿf)(0)` of `k ↦ z(δ(k))`, computed via finite
  differences over `Z/p^e` with `e = 5`. On every tested j=0
  curve: `a_0 = a_1 = 0` (trivially, since δ(0) = δ(1) = 0),
  and **`v_p(a_n) = 1` uniformly for all n ≥ 2**. Nonzero
  coefficients are evenly distributed across residue classes
  mod 6. Conclusion: δ is continuous but "maximally
  p-adically rough" — Mahler coefficients do not decay, so
  δ is not locally analytic, and the μ₆ symmetry is invisible
  at the Mahler level. A nontrivial low-rank Mahler support
  would have been the last possible structural handle; it is
  absent.

### Phase 14b/c/d/f closes the last structural lead

Every surviving lead from Phases 1–22 was a small-matrix artifact
or an Aut(E)-equivariance. The canonical-lift attack surface on
j = 0 (and by extension any j) curves is inert.

## Final verdict — the p-adic lifting attack on ECDSA is dead

What every phase certifies, in aggregate:

1. **A computable canonical lift exists** (Phase 4 — essentially a
   Serre–Tate construction for j = 0, with explicit projective
   arithmetic on `Z/p^e`).
2. **The Teichmüller error δ(k)** is the natural
   scalar-dependence measure on the lift and is provably:
   - pseudorandom at the DFT level (Phase 17),
   - maximal in Berlekamp–Massey linear complexity (Phase 18),
   - uniform in p-adic valuation on Ê (Phase 19),
   - free of polynomial / additive / multiplicative fit (Phases
     2, 8, 20),
   - coordinate-dependent (Phase 22).
3. **The only exact algebraic relations** in δ are:
   - `δ(k) + δ(n−k) = [n]·τ(G)` (Phase 21b) — antisymmetry.
   - `δ_{ζH}(k)/p ≡ χ(ζ) · δ_H(k)/p mod p` (Phase 14d/f) —
     Aut(E)-equivariance of the leading p-adic digit.
4. **Both relations are "free" from the curve's public structure**
   and give no advantage against the ECDLP.
5. **No natural functional on the lift depends linearly on k
   modulo p^{e−1} while being efficiently computable from the
   reduction mod p alone** — every candidate we constructed either
   vanishes by homomorphism (Phase 3.5, 4) or is pseudorandom.

The obstruction remains what it was after Phase 10, now with much
higher confidence: *the canonical lift is a legitimate
number-theoretic object, but not a cryptanalytic weapon against
ECDSA on j = 0 or any tested curve*.

With Phases 26–29 added, every "one more thing to try" has been
tried:

| Lever                          | Result                            |
|--------------------------------|-----------------------------------|
| Formal logarithm linearization | No low-degree structure (Ph 26)   |
| Mazur–Tate p-adic height       | Reduces to formal log (Ph 27)     |
| Non-CM curves at 40 bits       | Same rank rule / antisymmetry     |
| Mahler expansion               | Uniform `v_p = 1`, no decay (Ph 29) |
| Genus-2 Jacobian cover         | No linear Mumford structure (Ph 11b) |
| Weil pairing                   | Trivial on prime-order G          |

Every structural probe — whether coming from formal groups,
p-adic heights, p-adic Fourier analysis (Mahler), cohomological
automorphism equivariance, higher-genus Jacobians, or pairings —
confirms the same verdict: the Teichmüller (or any efficiently
computable) lift-error `δ(k)` is algebraically inert.

The only exact relations are:

1. `δ(k) + δ(n−k) = [n]·τ(G)` — universal (all j, any curve,
   confirmed at 40-bit non-CM).
2. `δ_{ζH}(k)/p ≡ χ(ζ)·δ_H(k)/p mod p` — Aut(E)-equivariance
   of the leading digit, giving rank deficit `(g−1)/|Aut|`.

Both are **public** consequences of the curve's geometry and
provide zero cryptanalytic leverage.

## Phase 32b — explicit Serre–Tate canonical lift (j = 0)

For ordinary `E/F_p` with `j = 0` and `p ≡ 1 mod 3` we constructed
the Serre–Tate canonical lift `τ_can : E(F_p) → E(Z_p)` algorithmically
via Newton iteration on the Frobenius equation `π·Q − Q = 0`, where
`π = a + b·ω` is the explicit element of `Z[ω]` of norm `p` and
trace `a_p`.  The Newton scaling is the tangent-space inverse
`(π_t − 1)^{-1}` with `π_t = a + b·ω_lift ∈ Z_p`.

Pure-Python projective-coordinate version (Phase 32) failed on 4 of 7
test primes due to Z-coordinate overflow when scalar-multiplying
deep kernel-of-reduction points.  The Sage version (Phase 32b),
running on `EllipticCurve(Zp(p, prec))`, **converges on all 7 test
primes** (31, 43, 67, 73, 79, 97, 103) at `prec = 10`.

**Validation.**  For each prime, we computed
`δ_can(k) := [k]·τ_can(G) − τ_can(kG)` for `k = 1..12` and found
**12/12 zeros** (vs `1/12` for the naive Hensel lift).  This is
the strongest possible verification that `τ_can` is a homomorphism
on `⟨G⟩`.

**Cryptanalytic verdict.**  Because `τ_can` is a homomorphism, the
DLP in `⟨τ_can(G)⟩ ⊂ E(Z_p)` is *exactly* as hard as the DLP in
`⟨G⟩ ⊂ E(F_p)`.  The lift moves the problem without weakening it,
exactly as Lubin–Serre–Tate theory predicts.  Phase 32b is therefore
an **engineering deliverable** — the first explicit, working
canonical-lift algorithm in this codebase, validated by Sage — but
**not** a cryptanalytic lever.

What it *does* close: any speculation that there might be a
"slightly wrong" canonical lift that's still computable and yet
secretly exploitable.  Newton on the Frobenius equation produces
*the* canonical lift (uniquely characterised in `ker(π−1)`), and it
is unconditionally a homomorphism.  Non-canonical sections (like
the naive Hensel lift) are exactly the ones whose `δ` is non-zero,
and Phases 21–29 already showed that *those* `δ` are
algebraically inert.
