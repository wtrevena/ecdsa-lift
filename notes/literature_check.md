# Literature check — how novel is the inertness theorem?

*Done before writing the paper, to calibrate claims honestly.*

## What is genuinely prior art (the theorem is NOT new mathematics)

The conceptual core of our "inertness theorem" is established:

1. **Silverman, "The Four Faces of Lifting for the ECDLP"** (ECC 2007
   talk) and **"Lifting and Elliptic Curve Discrete Logarithms"**
   (Springer LNCS 5747, 2009). This is *the* reference that
   systematically frames the four lifting strategies (local/global
   field × torsion/non-torsion points) and explains why each fails.
   Our setting — local field, p-adic canonical/Teichmüller lift — is
   one of his four faces. The principle that the lift *transports* the
   DLP without weakening it, and that the only local case that works is
   the anomalous (Smart) one, is his.

2. **Gadiyar & Padma, "The discrete logarithm problem over prime
   fields: the safe prime case. The Smart attack, non-canonical lifts
   and logarithmic derivatives"** (arXiv:1702.07107, 2017), building on
   Gadiyar–Maini–Padma, Indocrypt 2004. This is the *multiplicative
   group* analogue and states our core idea almost verbatim:
   - the Teichmüller (canonical) lift has vanishing logarithm and gives
     **no** information,
   - **"solving the discrete logarithm problem is equivalent to the
     construction of a non-canonical lift,"**
   - the carry/Iwasawa-log obstruction, and the Smart/Satoh–Araki/
     Semaev anomalous case as the degenerate point where it dissolves.

3. **Smart (1999), Satoh–Araki (1998), Semaev (1998)** — the anomalous
   `n = p` attack itself, which our Phase 42 reproduces.

## What (at most) is ours

The theorem, as we state it for elliptic curves, is a **precise,
explicit instance** of the above principle in the canonical-lift /
2-cocycle language, with three specifics I did **not** find stated
verbatim in the sources above:

- the well-definedness obstruction identified concretely as the kernel
  element `[n]·τ(G)` (= our Phase 21b antisymmetry constant);
- its additive order being a pure power of `p`, so the secret-bearing
  part `[k]c` is governed by `k mod p^{e-1}`, independent of the public
  `k mod n` (the clean `gcd(n,p)=1` ⇒ CRT-independence statement);
- the sharp-boundary framing: the obstruction is nonzero **iff**
  `gcd(n,p)=1`, vanishing exactly on anomalous curves.

**Honest assessment of novelty: low.** This is folklore-made-precise.
It is squarely within — and arguably a corollary of — Silverman's
framework, and conceptually identical to Gadiyar–Padma transposed from
`F_p^*` to `E`. It should be presented as a clean explicit instance and
unifying lens, **not** as a new theorem, with Silverman and
Gadiyar–Padma cited prominently.

## Consequence for the paper

This confirms the earlier strategic call: **Path 1 (the theorem) cannot
carry the paper as a novel result.** Its role is expository/unifying —
to explain *why* all 40 empirical phases had to come out inert, and to
retroactively justify the methodological suspicion. The distinctive
contribution is **Path 2**: the methodology of systematic negative
exploration plus mechanism-derived adversarial controls, with the
Phase 37 → 40 false-positive-and-kill as the worked example. The
theorem is the "why it had to be a false positive" that closes that
loop.

## Other references worth citing

- J. H. Silverman, *The Arithmetic of Elliptic Curves*, 2nd ed.,
  GTM 106 — formal groups (Ch. IV), kernel of reduction / filtration
  (Ch. VII), the exact sequence `0 → Ê → E(Z_p) → E(F_p) → 0`.
- Lubin–Serre–Tate (1964) — canonical lift existence/uniqueness.
- Satoh (2000) — canonical-lift point counting (the computational
  instance of the same object; our Phase 32b is a small version).
- Silverman, "On Partial Lifting and the ECDLP" (ASIACRYPT 2004) and
  "Lifting Elliptic Curves and Solving the ECDLP" — partial-lift
  obstructions.
