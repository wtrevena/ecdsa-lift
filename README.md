# ecdsa-lift

Experimental cryptanalysis of secp256k1 via higher-dimensional lifts.
Can structure invisible on `E(F_p)` become visible when the curve is
lifted into `E(Z_p)`, `E(Q_p)`, formal groups, isogeny neighborhoods,
or higher-genus Jacobians?

## Thesis

secp256k1 has `j = 0`, so it has complex multiplication by `Z[ω]`.
That algebraic structure is currently used *defensively* (GLV). The
question this repo exists to answer is whether the same structure can
be turned offensive by lifting the discrete-log problem into a space
where the group law linearizes.

## Roadmap

- [x] **Phase 1** — Teichmüller lift error `δ(k)` is injective at
      precision `≥ p^4`. DLP information provably survives the lift.
- [x] **Phase 2** — `δ(k)` under the Teichmüller lift is pseudorandom:
      not polynomial, linear, multiplicative, or CM-equivariant.
      Inversion equivalent to brute force.
- [x] **Phase 3** — `[n]`-multiplication trick (Smart-SSSA) only works
      for anomalous curves (`n = p`). For `gcd(n, p) = 1` the choice
      of lift introduces a unit ambiguity that kills the reduction.
- [x] **Phase 4** — Group-cohomology splitting. Constructed `s`
      explicitly via `s(g) = τ(g) - μ([n]·τ(g))` where `μ` is the
      inverse of mult-by-`n` on `Ê`. Verified as a homomorphism on
      several small j=0 curves in projective coordinates. **But**
      `[n]·s(P) = s(O) = O` trivially, so the Smart-SSSA-style
      formal-log attack does *not* follow. See
      `notes/obstructions.md::D4` for the full story.
- [x] **Phase 5 (CM-linearization).** `ψ(P) = [λ]·τ(P) − τ(λP)` is
      injective but neither additive nor CM-equivariant — just
      another pseudorandom injection into `Ê`. See D5.
- [x] **Phase 6 (HNP on cocycle).** No linearity, residues
      uniformly distributed, LLL finds no non-trivial short vector.
      See D6.
- [x] **Phase 7 (MOV embedding degree).** Generic large embedding
      degrees for ordinary j=0 curves; secp256k1's `k` is
      astronomical. MOV inapplicable. See D7.
- [x] **Phase 8 (polynomial fit).** Modular Vandermonde fit for
      `δ(k)` degree ≤ 7 produces only trivial matches (exactly the
      points used to solve the system). See D8.
- [x] **Phase 9 (LLL on lift coordinates).** Rank-deficient
      lattice; LLL finds the zero vector. No linear leakage from
      projective lift coordinates. See D9.
- [x] **Phase 10 (ordinary↔supersingular bridge).** Impossible by
      trace-of-Frobenius invariance under isogeny. See D10.

**Result: ECDSA not broken.** See `notes/obstructions.md` for the
full synthesis — every phase was a negative result, and the
accumulated obstructions describe a comprehensive (though not
provably complete) wall around the p-adic lifting attack surface
on j = 0 curves.

## Layout

- `notes/` — background, derivations, prior-art summary, the full
  chat transcript that motivated the repo (`00_origin.md`).
- `src/` — reusable primitives: curve arithmetic over `Z/p^k Z`,
  Teichmüller lifts, formal-group logarithm, Hensel lifting.
- `experiments/` — numbered, reproducible scripts. Each one prints a
  clear pass/fail and writes a JSON result under `results/`.
- `results/` — raw outputs of each experiment run.

## Ground rules

1. **Every experiment is reproducible.** Fixed seeds, primes written
   into filenames, JSON output.
2. **Negative results are valuable.** If a route dies, document *why*
   it died in `notes/obstructions.md` so we don't re-enter it.
3. **Small primes first, always.** Nothing scales to secp256k1 until
   it works on `p = 17, 19, 23, 29, 31, 37, 41, 43, 47, 61`.
4. **SageMath for anything touching canonical lifts, isogeny graphs,
   or `p`-adic fields beyond raw `Z/p^k`.** Python stdlib only for
   low-level exploration.

## Phases 40–42 and the paper

After the Phase 1–39 negative sweep, three phases closed the loop and
became the basis of a short paper (`paper/ecdsa_lift_paper.tex`, builds
to a 7-page PDF):

- **Phase 40 / 40b** — *mechanism-derived controls* that overturn the
  Phase 37 "Gowers anomaly." The increment-shuffle null (40) refutes the
  prefix-sum-artifact hypothesis; the half-sequence null (40b) localizes
  the entire signal to the exact public identity
  `δ(k)+δ(n−k)=[n]τ(G)` (Phase 21b). The anomaly is real but not new and
  carries no advantage. See `notes/phase40_resolution.md`.
- **Phase 41** — numerical verification of the **inertness theorem**
  (`notes/theorem.md`): for `gcd(n,p)=1`, the secret-bearing part of any
  section lift error is governed by `k mod p^{e−1}`, independent of the
  public `k mod n`; the well-definedness obstruction is the nonzero
  kernel element `[n]τ(G)` of maximal order `p^{e−1}`. All core claims
  pass on every test prime (pure-Python `z`-coordinate arithmetic).
- **Phase 42** — the boundary is **sharp**: on an anomalous curve
  (`#E(F_p)=p`, `gcd(n,p)=p`) the obstruction vanishes and Smart's attack
  recovers `k` with 100% success (12/12). The single inequality
  `gcd(n,p)=1` is the exact line between inert and broken.

The paper's thesis is methodological: in structure-hunting cryptanalysis
a control is valid only against a *named* alternative hypothesis, and
nulls must be derived from mechanisms. The inertness theorem (an explicit
instance of Silverman's "Four Faces of Lifting" and Gadiyar–Padma, not
new mathematics — see `notes/literature_check.md`) is what certifies,
after the fact, that the suspicion was structurally forced.
