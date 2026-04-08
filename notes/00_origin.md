# Origin

This repository came out of a long exploratory conversation about
whether a frontier LLM pointed at the mathematical foundations of
cryptographic primitives (as opposed to implementation bugs in the
code that uses them) could find something real.

The thesis the conversation settled on:

> Cryptographic hardness is always stated inside a specific
> mathematical framework. Historically, "impossible" problems become
> tractable when lifted into a higher-dimensional / richer mathematical
> space where additional structure becomes visible. Transformer models
> operate natively in very high-dimensional representation spaces and
> might have non-obvious inductive biases for exactly this kind of
> embedding search.

The specific target: **secp256k1**, because `j = 0` gives it CM by
`Z[ω]` and therefore more algebraic structure to work with than a
generic curve.

## What was already established before the repo existed

### Phase 1 — Information survives the lift
At `p`-adic precision `≥ p^4`, the Teichmüller lift error
`δ(k) := [k]G̃ - \widetilde{[k]G}` (computed in `Z/p^4 Z`)
is **injective** as a function of the secret scalar `k` across every
small prime tested. The discrete-log information is provably present
in the `p`-adic data.

### Phase 2 — But scrambled
Under the Teichmüller lift specifically, `δ(k)` is:
- not polynomial in `k` (finite differences don't vanish),
- not additive (`δ(a+b) ≠ δ(a) + δ(b)`),
- not multiplicative,
- not CM-equivariant (`δ(λk)/δ(k)` isn't constant).

So with *this* lift, inverting `k ↦ δ(k)` is equivalent to brute
forcing the DLP.

One structural oddity worth keeping: when `p = n` (anomalous curves),
`δ(k) = δ(n - k)` exactly — a `Z/2`-symmetry induced by the
interaction of the Teichmüller `p`-th-power construction with the
group order. This symmetry *breaks* as soon as `n ≠ p`, which is
suggestive.

### Phase 3 — The Smart–SSSA wall
The `[n]`-multiplication trick (`ψ(P) := log_Ê([n]·P̃) / p`) which
solves the DLP in linear time for anomalous curves (`n = p`) does
**not** generalize: for `gcd(n, p) = 1` the map depends on the choice
of lift, and the ambiguity is a unit in `Z_p` times an arbitrary
formal-group element. The "obvious" attack dies exactly there.

### Phase 3.5 — The cohomological twist
Here is where the conversation left off, and where this repo picks
up. The short exact sequence of abelian groups
```
0 → Ê(pZ_p) → E(Z_p) → E(F_p) → 0
```
has associated obstruction class in `H^2(E(F_p), Ê(pZ_p))`. When
`gcd(n, p) = 1` — i.e. for non-anomalous curves, including
secp256k1 — this `H^2` vanishes because `|E(F_p)| = n` is invertible
on `Ê(pZ_p)`. So the sequence **splits as an extension of abelian
groups**: a homomorphic section
```
s : E(F_p) → E(Z_p),   s(P + Q) = s(P) + s(Q)
```
provably exists.

If we had such an `s` in our hands, the game is over:
```
ψ(P) := log_Ê( [n] · s(P) )   ∈   Z_p
```
is a group homomorphism `E(F_p) → (Z_p, +)` (because `[n] · s(P) ∈ Ê`
by the sequence, and `log_Ê` linearizes the formal group), and then
```
k = ψ(Q) / ψ(G)   (mod p-adic precision, then CRT'd up)
```
recovers the secret scalar in a constant number of group operations.

So Phase 4 is the whole ballgame: **the section exists. Can we
construct it explicitly?** The Teichmüller lift is a set-theoretic
section but not a homomorphism — that's exactly the obstruction Phase
2 ran into. The question is whether there's a *computable* choice of
section, or whether the existence proof is purely non-constructive.

## Open questions this repo will chew on

1. **Phase 4 (primary).** Constructive splitting. Candidates:
   - Average the Teichmüller section over the group law:
     `s(P) := (1/n) · Σ_{Q ∈ E(F_p)} ( τ(P + Q) - τ(Q) )`
     where `τ` is the Teichmüller lift. This is the usual trick for
     constructing a cocycle trivialization when `|G|` is invertible.
     It requires `n` to be invertible in `Z_p`, which is exactly our
     `gcd(n, p) = 1` assumption. Plausibly works.
   - Push to the formal group first via `[n]`, then use the
     logarithm's linearity to define `s` implicitly.
   - Canonical Serre–Tate lift: for `j = 0` the canonical lift is
     explicit, and it might *already* be a homomorphic section.

2. **Phase 5.** If constructive Phase 4 fails (averaging has a
   hidden obstruction, etc.), move to SageMath and the canonical lift.

3. **Phase 6.** Independent angle — isogeny walks, crossing to
   supersingular territory, picking up extra endomorphisms en route.

4. **Phase 7.** The meta-concern: lattice PQC might be *more*
   vulnerable than ECDSA to the same style of attack.

## How to collaborate with a model on this

- One experiment per script, each small enough to read in one
  sitting.
- Every experiment writes JSON to `results/`. Models can diff runs.
- `notes/obstructions.md` is the graveyard. Any route that dies, it
  dies with a specific reason so we don't re-enter it.
