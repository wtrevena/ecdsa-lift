# The inertness theorem (statement and proof)

*An explicit instance, in canonical-lift / 2-cocycle language, of the
well-known principle (Silverman's "Four Faces of Lifting"; Gadiyar–Padma
for `F_p^*`) that a p-adic lift transports the ECDLP without weakening
it. Stated and proved here in the precise form our experiments verify.*

## Setup

Let `p ≥ 5` be prime and `E/F_p` an elliptic curve with `#E(F_p) = N`.
Fix a point `G ∈ E(F_p)` of order `n` and assume **`gcd(n, p) = 1`**
(the ECDSA case: `n` is a large prime `≠ p`). Let `Q = kG` be the public
key; the ECDLP is to recover `k ∈ Z/nZ` from `(E, G, Q)`.

For `e ≥ 1` write `E(Z/p^e)` for the points of a fixed Weierstrass lift
of `E` to `Z/p^e Z`. Reduction mod `p` gives the exact sequence of
abelian groups
```
        0 → Ê(pZ/p^e) → E(Z/p^e) --red--> E(F_p) → 0,            (∗)
```
where `Ê(pZ/p^e)` is the kernel of reduction, a `p`-group of order
`p^{e-1}`, with the formal-group parameter `z = -X/Y` giving an
isomorphism `Ê(pZ/p^e) ≅ (pZ/p^e Z, +_F)` (`+_F` the formal group law),
and the formal logarithm `log: Ê → (pZ/p^e, +)` an isomorphism
satisfying `log([m]u) = m·log(u)`.

A **section** is any map `τ: E(F_p) → E(Z/p^e)` with `red ∘ τ = id`
(not necessarily a homomorphism). Its **lift error** (a 1-cochain) is
```
        δ_τ(k) := [k]·τ(G) − τ(kG)  ∈  Ê(pZ/p^e).
```
(`δ_τ(k)` lies in `Ê` because `red([k]τ(G) − τ(kG)) = kG − kG = O`.)

## Lemma 1 (canonical lift; existence, uniqueness, homomorphy)

Because `gcd(n,p)=1`, the subgroup `⟨G⟩` of order `n` lifts uniquely to
a subgroup `⟨G̃⟩ ⊂ E(Z/p^e)` of order `n` with `red(G̃)=G`, and
`τ_can: jG ↦ [j]G̃` is a group homomorphism (the canonical lift on
`⟨G⟩`), i.e. `δ_{τ_can} ≡ 0`. Explicitly `G̃ = [m]·τ(G)` for any section
`τ`, where `m ≡ 1 (mod n)` and `m ≡ 0 (mod p^{e-1})` (CRT, valid since
`gcd(n, p^{e-1}) = 1`).

*Proof.* `E(Z/p^e)` is a finite abelian group of order `N·p^{e-1}`; its
`n`-torsion injects under reduction into `E(F_p)[n]` (an element of `Ê`
killed by `n` is killed by `gcd(n, p^{e-1}) = 1`, hence trivial), and
surjects onto `E(F_p)[n] ⊇ ⟨G⟩` because `n`-torsion is étale/coprime to
`p` (Hensel). So `red: E(Z/p^e)[n] → E(F_p)[n]` is an isomorphism and
`G̃ := (red|_{[n]})^{-1}(G)` is the unique order-`n` lift; `⟨G̃⟩` is its
span and `τ_can` is a homomorphism by construction. For the formula:
write `τ(G) = G̃ + c` with `c ∈ Ê` (possible since `red(τ(G))=G=red(G̃)`).
Then `[m]τ(G) = [m]G̃ + [m]c = [m mod n]G̃ + [m mod p^{e-1}]c
= [1]G̃ + [0]c = G̃`, using `ord(G̃)=n` and `exp(Ê)=p^{e-1}`. ∎

## Theorem (inertness of the lift error)

Let `τ` be any section and write `c := τ(G) − G̃ ∈ Ê`. Then:

**(a) Decomposition.** For every integer `k`,
```
        δ_τ(k) = [k]·c − c_k,      where c_k := τ(kG) − [k]G̃ ∈ Ê.
```
The first term `[k]c` is the only part depending on the scalar through
multiplication; `c_k` depends only on the *point* `kG` (public).

**(b) The secret part is governed by `k mod p^{e-1}`.** `[k]c` depends on
`k` only through `k mod p^{e-1}`; equivalently `log([k]c) = k·log(c)`
in `pZ/p^e`, of additive order `p^{e-1}/gcd(...)`, a power of `p`.

**(c) The public residue does not give `[k]c` (computational form).**
Since `gcd(n, p^{e-1})=1`, the CRT isomorphism
`Z/np^{e-1} ≅ Z/n × Z/p^{e-1}` makes `k mod n` (fixed by the public `Q`)
and `k mod p^{e-1}` (which governs `[k]c`) independent coordinates of a
`k` ranging over `Z/np^{e-1}`. The public data fixes only `k mod n`.
**Caveat (do not overstate):** once a canonical representative
`k ∈ [1,n-1]` is fixed — as every experiment does — `k mod p^{e-1}` and
hence `[k]c` are *determined*; the correct claim is **computational**, not
information-theoretic: no efficient algorithm computes `[k]c` from
`(E,G,Q)` short of solving the ECDLP (this is the contrapositive of (d)).

**(d) Well-definedness obstruction.** `δ_τ` is a well-defined function of
the *instance* `(E,G,Q)` (i.e. of `k mod n`) if and only if `c = 0`,
i.e. iff `τ` is canonical on `⟨G⟩` — in which case `δ_τ ≡ 0`. For
`c ≠ 0`, replacing the integer representative `k` by `k+n` shifts `δ_τ`
by the fixed nonzero kernel element
```
        [n]·τ(G) = [n]·c  ≠ O,
```
which is exactly the Phase 21b antisymmetry constant
(`δ_τ(k) + δ_τ(n−k) = [n]τ(G)`). **Note:** this antisymmetry requires
the section to be odd, `τ(−P) = −τ(P)`, so that
`τ(kG) + τ(−kG) = O`; it holds for the x-Teichmüller / Hensel sections
used here (the two Hensel `y`-roots are negatives) but is a property of
the section, not of every conceivable lift.

**Corollary (no lift oracle).** No section's lift error provides a
function of the public data that depends on the secret: for `gcd(n,p)=1`
every `δ_τ` is either identically zero (canonical) or not a function of
the instance at all (non-canonical, ambiguity `[n]c ≠ O`). The p-adic
canonical lift moves the ECDLP into `E(Z/p^e)` without weakening it.

*Proof.* (a) `δ_τ(k) = [k]τ(G) − τ(kG) = [k](G̃+c) − ([k]G̃ + c_k)
= [k]G̃ + [k]c − [k]G̃ − c_k = [k]c − c_k.`
(b) `c ∈ Ê` has order dividing `exp(Ê)=p^{e-1}`, so `[k]c=[k mod p^{e-1}]c`;
the log statement is the isomorphism `log` with `log([k]c)=k log(c)`,
and the order of `k ↦ k·log(c)` in `(pZ/p^e,+)` is
`p^e / gcd(log(c), p^e)`, a power of `p` dividing `p^{e-1}` since
`v_p(log(c)) ≥ 1`.
(c) Immediate from `gcd(n, p^{e-1})=1` and CRT.
(d) `δ_τ` depends on the integer `k` only via `[k]c` (the term `c_k`
depends on the point `kG`, a function of `k mod n`). `k ↦ [k]c` factors
through `k mod n` iff `[n]c = O`. Since `c ∈ Ê` and `gcd(n,p^{e-1})=1`,
`[n]c=O ⇔ c=O`. When `c≠O`, `δ_τ(k+n)−δ_τ(k)=[n]c=[n]τ(G)` (as
`[n]G̃=O`), the antisymmetry constant. ∎

## The boundary is sharp (anomalous case)

If instead `n = p` (anomalous, `#E(F_p)=p`), the hypothesis
`gcd(n,p)=1` fails and every step degenerates in the attacker's favour:

- `k mod n = k mod p` now *coincides* with the residue the formal group
  sees, so (c) collapses — there is no independent hidden coordinate;
- the obstruction `[n]c = [p]c` **vanishes**, because `c ∈ Ê(pZ/p^2)`
  has order dividing `p`, so `[p]c = O`; thus `δ_τ` becomes well-defined
  on the instance;
- the formal logarithm of `[p]·(lift)` is then a homomorphism
  `E(F_p) → F_p` recovering `k` — Smart's attack.

So the single inequality `gcd(n,p)=1` is exactly the line between
"lift is inert" (ordinary ECDSA curves) and "lift breaks the DLP"
(anomalous curves). This is verified numerically in Phase 41 (the
obstruction `[n]τ(G)` is a nonzero `p`-power-order kernel element on
every ordinary test curve) and Phase 42 (th