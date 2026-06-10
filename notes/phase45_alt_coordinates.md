# Phase 45 — Alternative coordinates / uniformizers for the lift error

*Falsification experiment. Prior on success ~5%. Result: clean null.
This closes a standing referee question ("did you check whether a
non-Weierstrass presentation exposes the structure?").*

## The question

The inertness theorem (`notes/theorem.md`) says the abstract object
`delta(k) = [k]tau(G) - tau(kG) in Ehat(pZ/p^e)` is cohomologically
inert. Prior phases certified the *Weierstrass* parameter
`z(delta(k)) = -X/Y` pseudorandom: no polynomial / Mahler /
multiplicative fit, flat DFT, near-maximal Berlekamp-Massey complexity.

But "efficiently computable low-degree relation" is *presentation
sensitive*. The abstract object being inert does not, by itself, forbid
a different coordinate system from making the *same* sequence
algebraically simple (the SEA / Edwards / AGM precedent: a change of
model can turn an opaque computation into a transparent one). Phase 45
re-expresses `delta(k)` in alternative presentations and re-runs the
structure probes on the re-expressed sequence.

## What was tried

Code: `experiments/phase45_alt_coordinates.py`
Data: `results/phase45_alt_coordinates.json`
Curves (j = 0, fixed list, all non-anomalous, gcd(n,p)=1, prime order):
`(p,b) in {(31,3),(43,7),(67,2),(73,13),(79,3),(97,5),(103,5)}`.
Precision `e = 8` (chosen so the `z^7` formal-log term is active:
`z in pZ` => `v_p(z^7) = 7 < 8`; below `e = 8` the formal log collapses
to the identity and the uniformizer change would be vacuous).

No PARI on this VM (`cypari2` is absent — `phase26` imported it and would
not run here). The formal-group power series are computed **natively**
with exact rational arithmetic and cross-checked two independent ways
against the documented j = 0 coefficients:
`log(z) = z + (3b/7) z^7 + (15 b^2/13) z^13 + ...`
(note: the `src/formal_group.py` docstring's "2b^2/13" for the z^13 term
is a low-precision approximation and is *wrong*; the correct coefficient
is `15 b^2 / 13`, confirmed by both the `omega = dx/2y` integration and a
direct `w`-recursion expansion. This does not affect the verdict.)

### Presentation A — uniformizer reparametrizations `t = f(z)`

For each of six presentations we recomputed `t(delta(k))` and ran four
probes: (1) modular Vandermonde polynomial fit, degrees 1..8 (phase08);
(2) first-difference test (linear / quadratic in k); (3) additive DFT
flatness, peak/mean of nonzero bins (phase17); (4) Berlekamp-Massey
linear complexity over F_p of the top p-adic digit (phase18).

| id | uniformizer | rationale |
|----|-------------|-----------|
| A0 | `t = z` | control (the certified-pseudorandom baseline) |
| A1 | `t = log_Ehat(z)` | the canonical linearizing parameter (`+_F` -> `+`) |
| A2 | `t = z/(1+3z)` | Mobius reparametrization |
| A3 | `t = z + 2 z^2` | quadratic reparametrization |
| A4 | `t = z + 2 z^2 + 2 z^3` | cubic-ish reparametrization |
| A5 | `t = z(1 + w)`, `w = -Z/Y` | honest *non-z* parameter built from both kernel coordinates |

A1 was verified non-degenerate: `log(z) != z` for every sampled k, with
`v_p(log(z) - z) = 7` exactly — the `z^7` term is genuinely contributing
a nonzero p-adic digit, so the probe is testing a real alternative
uniformizer rather than silently re-testing `z`.

### Presentation B — Edwards / Montgomery / theta

**INFEASIBLE for every candidate, and this is itself a finding.** A short
Weierstrass curve maps to Montgomery form iff it has a rational point of
order 2, and to (untwisted) Edwards form iff it has a rational point of
order 4. Every candidate has *prime* group order (43, 31, 73, 67, 97,
79, 97), hence **no nontrivial rational 2- or 4-torsion at all**. So
neither model exists over these fields. This matches the task's a priori
hint that j = 0 curves generally lack the order-4 point Edwards needs;
here it is even stronger (no order-2 either).

To guard against a *silent null from a dead pipeline*, a **positive
control** is included: a planted linear sequence `t(k) = A k + B mod p^e`
is fed through the exact same `probe_sequence`. The pipeline correctly
returns **STRUCTURE FOUND** (perfect degree-1 polynomial fit, constant
first difference). So a null in Presentation A is a real null, not a
broken detector.

## Results

**Presentation A: clean null in all 6 presentations x 7 curves.**

- Polynomial fit: best is always the trivial "degree-8 through 9 points"
  baseline (`9/total` exact hits, `poly_perfect_fit = False`). A genuine
  low-degree relation would give `hits = total`. No presentation does.
- First differences: never constant (no linear), never second-constant
  (no quadratic), for any presentation.
- DFT peak/mean: ~1.7–3.6 across the board (flat-spectrum regime;
  threshold for "peak" is 4.0). No reparametrization concentrates the
  spectrum. The formal-log presentation A1 is no flatter or peakier than
  the identity — it just shuffles the same noise.
- Berlekamp-Massey complexity: `L ~= total/2` everywhere (e.g. 21/42,
  15/30, 30/60) — the maximal-complexity signature of a pseudorandom
  sequence. Identical (to within +/-1) across all presentations.

The presentation does not move the needle on any probe. Crucially the
*formal logarithm itself* (A1) — the one map that genuinely linearizes
the formal group law — does not linearize `delta(k)` as a function of
**k**. That is the expected consequence of the theorem: `log` linearizes
`+_F` (the group structure on the kernel), but `delta`'s opacity is in
its *dependence on the scalar k*, which `log` leaves untouched.

**Presentation B: infeasible (documented), positive control passes.**

## Verdicts

| presentation | verdict |
|--------------|---------|
| A0 identity `z` | still pseudorandom (baseline, as expected) |
| A1 formal log | still pseudorandom |
| A2 Mobius `z/(1+3z)` | still pseudorandom |
| A3 quadratic `z+2z^2` | still pseudorandom |
| A4 cubic `z+2z^2+2z^3` | still pseudorandom |
| A5 xy-norm `z(1+w)` | still pseudorandom |
| B Edwards | INFEASIBLE (no rational 4-torsion; prime order) |
| B Montgomery | INFEASIBLE (no rational 2-torsion; prime order) |
| B theta | INFEASIBLE (subsumed by B: needs 2-torsion structure) |
| positive control (planted linear) | STRUCTURE FOUND (pipeline live) |

**Overall: ALL presentations still pseudorandom — clean null.**

## Why this is the *expected* outcome (theory check)

A reparametrization `t = f(z)` with `f` an invertible analytic change
near 0 is, on the kernel, just a coordinate relabelling of the same
group `Ehat`. The probes ask whether `t(delta(k))` is low-degree *in k*.
But `f` acts on the *value* `delta(k)`, not on its argument `k`; a
fixed invertible `f` cannot convert a sequence with maximal BM
complexity into a low-complexity one unless `f` itself encodes the
hidden recurrence — and a *fixed, k-independent, low-degree* `f` has no
room to do that (that is exactly what the BM-invariance under A1..A5
shows numerically). The only `f` that could help is one tuned to k,
which is begging the question. So the theorem's inertness is robust to
this whole family of moves, as observed.

## Honest prior-updated assessment

Going in, the prior on success was ~5%. Nothing here beat the baseline
on any probe, the formal log (the strongest single candidate) was inert,
and there is a clean structural reason why a fixed reparametrization
*cannot* expose k-structure that BM says is not there. Posterior on
"some elementary coordinate change exposes a low-degree relation in k"
drops to roughly **1%**.

**Recommendation: close this lane.** Specifically:

- Do NOT spend more time on fixed invertible reparametrizations of z,
  other formal-group uniformizers, or twist/automorphism relabellings.
  These are now covered (A1..A5 here; phase22 for twists; phase26 for the
  formal log at high precision). Add to `notes/obstructions.md`.
- The ONLY residual, and it is a long shot, is a model that is *not* an
  invertible reparametrization of the same kernel — i.e. a genuinely
  different curve model (Edwards/Montgomery/Jacobi-quartic) where the
  Teichmuller section itself is redefined, so `delta` is a *different
  cochain*, not a relabelling of this one. That is rigorously infeasible
  for the prime-order j = 0 curves studied here (Presentation B). It
  could in principle be probed on a *sanity curve* engineered to carry
  4-torsion, but (i) that is no longer secp256k1-like, (ii) phase22
  already showed different Weierstrass models give pseudorandom-related
  deltas, and (iii) the theorem applies model-independently. Low value;
  pursue only if a referee explicitly demands the Edwards sanity curve.

## Reproduce

```
python3 experiments/phase45_alt_coordinates.py
# writes results/phase45_alt_coordinates.json
```

Deterministic (fixed candidate list, fixed `e = 8`, no randomness in the
math; `SEED = 45` is recorded for provenance only).
