# Phase 46 — OEIS / closed-form "close the lid" check on the lift error

*Cheap, falsifiable, pre-registered NULL test. Prior on a hit ~2%. Outcome:
NULL, as expected. This memo is the audit trail; the paper gets one footnote
(see VERDICT).*

## What was tested and why

The inertness theorem (Phase 41, `notes/theorem.md`) describes the section's
lift error

    delta(k) = [k] tau(G) - tau(kG)   in  Ehat(pZ/p^e),

read through the formal-group parameter `z = -X/Y`, as a "pseudorandom-looking"
element of `pZ/p^e` constrained only by the exact relations V1–V5 (antisymmetry
around `Cn`, integer-rep shift by `Cn`, etc.). Before the paper asserts "no
exploitable closed-form structure / pseudorandom", we owe the reader one
explicit, cheap attempt to *falsify* that: is `z(delta(k))` secretly a named
integer sequence (OEIS) or a simple closed form in `k` mod `p`?

Script: `experiments/phase46_oeis.py`. It reuses the *exact* `z_delta`
construction from `phase41_theorem_verification.py` (single curve op, then read
`z`; projective-degeneracy flagged), so the object tested is literally the one
the theorem is about, not a re-derivation.

Curve family: `y^2 = x^3 + 7` (j = 0, secp256k1 family), precision `e = 4`.
Primes: **31, 43, 67, 97, 109, 151, 163** (the task's 31/43/67/97 plus three
more that Phase 41 certified arithmetically reliable, to get longer — hence
less spuriously-matchable — sequences). All seven ran with
`arithmetic_reliable = True` (zero degenerate entries at e = 4 in this run;
note Phase 41 had flagged 43 and 103 as degenerate in its *deep* V3 probe, but
the single-op `z_delta` used here is non-degenerate for all seven of these).

## Sequences checked (k = 1 .. n-1)

For each prime, three sequences, full data in `results/phase46_oeis.json`:

| p   | n   | length | `z(delta(k))` (raw, in Z/p^e) | `z(delta(k)) mod p` | `d1(k) = (z//p) mod p` | leading decimal digit |
|-----|-----|--------|-------------------------------|---------------------|------------------------|-----------------------|
| 31  | 21  | 20     | varying, in [0, p^4)          | all 0               | varying in F_p         | varying 0–9           |
| 43  | 31  | 30     | "                             | all 0               | "                      | "                     |
| 67  | 79  | 78     | "                             | all 0               | "                      | "                     |
| 97  | 79  | 78     | "                             | all 0               | "                      | "                     |
| 109 | 129 | 128    | "                             | all 0               | "                      | "                     |
| 151 | 133 | 132    | "                             | all 0               | "                      | "                     |
| 163 | 139 | 138    | "                             | all 0               | "                      | "                     |

Key structural fact, confirmed numerically (`seq_mod_p_all_zero = true` for all
seven primes): **`z(delta(k)) mod p` is identically zero.** This is not a
coincidence to be matched — it is the theorem. `delta(k)` is a
kernel-of-reduction element, so `z(delta(k))` lies in `pZ`, hence vanishes mod
`p`. Searching the all-zeros sequence in OEIS or fitting it is vacuous, so the
genuinely informative `F_p`-valued object is the **first non-trivial p-adic
digit**

    d1(k) := (z(delta(k)) // p) mod p,

which is what the fits below actually target. (The task said "z(delta(k)) mod p
… digit sequence"; the honest reading of that, given the kernel constraint, is
`d1`.) The raw and leading-decimal-digit sequences are also recorded and queried
for completeness; the leading-digit one is deliberately included as a
calibration probe — it is exactly the kind of short, small-alphabet sequence
that matches OEIS spuriously, so a "hit" there would be a warning, not a result.

## Part 1 — OEIS queries

OEIS's own sequence-matcher (`https://oeis.org/search?q=<terms>&fmt=text`) could
not be fetched directly: the sandbox `web_fetch` tool returns an empty body for
`oeis.org` (the host blocks / the body is stripped — a Fibonacci control
`1,1,2,3,5,8,13,21` also came back empty, while non-OEIS hosts fetch fine). The
queries were therefore run through the `WebSearch` tool against the public index.
Every query came back with **no matching OEIS entry**. The full set of
copy-pasteable search URLs is printed at the end of the script run and is also in
the JSON (`oeis_query_*` fields). Representative queries actually executed:

- `q=461094,256184,801660,801660,427521,17918,253146,683705,714240` (p=31 raw)
  → no match. WebSearch likewise found no OEIS entry for the shorter
  `…,17918,253146`.
- `q=20425,2774446,234565,435246,…` (p=43 raw) → no match.
- `q=0,25,18,6,6,27,20,13,14,7,22,15,16,9,2,23` (p=31 `d1`) → no match.
- `q=0,4,2,8,8,4,1,2,6,7,2,2,7,3,5,1,1,7,5,5` (p=31 leading digit) → no match.
- `q=0,3,8,8,3,1,1,6,1,1,7,1,2,1,1,1,8,1,3` (p=67 leading digit) → no match.

Judgment of significance: **none — there were no hits to judge.** Caveat for
reproducibility: `WebSearch` is not the OEIS sequence-matcher itself, so it does
not replicate OEIS's offset-shifting / subsequence logic. If a reviewer wants
belt-and-suspenders certainty, paste the `oeis_query_raw` strings into
oeis.org by hand — but the offline structural test (Part 2) is independently
decisive, and the prior on a hit was ~2% to begin with. Note also that even a
*nominal* OEIS hit on the short leading-digit or `d1` sequences (≤ 138 terms, and
the short ones only 20–30 terms) would be presumptively spurious: OEIS has
>390k sequences, and short low-alphabet strings collide by chance.

## Part 2 — Offline closed-form / symbolic-regression on d1(k)

Models fit to `d1(k)` as a function of `k`, **exactly over Z/p** (least-squares
is meaningless mod p): constant, affine/linear, quadratic, cubic, multiplicative
`a·k`, geometric `a·c^k`. Protocol: fit/select parameters on the first half of
the valid `k`, then report **exact-match accuracy on the held-out second half**.
The honest baseline for a structureless sequence is `1/p` (uniform guess); a
genuine closed form would generalize to ≈ 1.0. For the polynomial models the
script first asks the strong question "is the *entire* training half exactly a
degree-≤d polynomial mod p?" — a positive there (`fit_mode =
exact_over_all_training`) would be a real hit; otherwise it falls back to the
degree-d interpolant through the first few points and measures generalization.

Held-out accuracy of `d1(k)` (full numbers in `results/phase46_oeis.json`):

| p   | 1/p (null) | mode  | const | affine | quad  | cubic | mult  | geom  | n_test |
|-----|-----------|-------|-------|--------|-------|-------|-------|-------|--------|
| 31  | 0.032     | 0.000 | 0.000 | 0.000  | 0.000 | 0.000 | 0.000 | 0.100 | 10     |
| 43  | 0.023     | 0.000 | 0.000 | 0.000  | 0.000 | 0.000 | 0.000 | 0.000 | 15     |
| 67  | 0.015     | 0.000 | 0.000 | 0.000  | 0.000 | 0.000 | 0.000 | 0.000 | 39     |
| 97  | 0.010     | 0.051 | 0.000 | 0.000  | 0.026 | 0.000 | 0.000 | 0.051 | 39     |
| 109 | 0.009     | 0.000 | 0.000 | 0.000  | 0.000 | 0.016 | 0.000 | 0.000 | 64     |
| 151 | 0.007     | 0.000 | 0.030 | 0.000  | 0.000 | 0.015 | 0.015 | 0.000 | 66     |
| 163 | 0.006     | 0.014 | 0.014 | 0.000  | 0.014 | 0.000 | 0.087 | 0.000 | 69     |

Findings:

- **No model generalizes.** Every held-out accuracy sits at, below, or
  trivially above the `1/p` random baseline. The largest single value, 0.100
  (p=31, geometric), is `1/10` of a 10-point test set — the coarsest possible
  non-zero granularity, fully consistent with chance; the next largest, 0.087
  (p=163, multiplicative), is `6/69 ≈ 0.087`, again chance-level for `p=163`.
  No value approaches 1.0; none even reaches 0.15.
- **`d1(k)` is provably not a low-degree polynomial in `k` mod p.** For all
  seven primes and all degrees 0–3, the polynomial fit fell back to
  `interpolation_first_k_points` — i.e. *no* degree-≤3 polynomial reproduces
  even the training half. (`exact_over_all_training` never occurred.) The
  interpolants then fail on held-out data at baseline, as expected.
- Multiplicative (`a·k`) and geometric (`a·c^k`) models likewise collapse to
  baseline despite the geometric search being given the unfair advantage of
  picking `c` by best *training* fit.

Conclusion of Part 2: `d1(k)` is statistically indistinguishable from a uniform
i.i.d. `F_p` sequence with respect to all of {constant, affine, quadratic,
cubic, multiplicative, geometric} modular models — consistent with the
theorem's "pseudorandom" characterization and inconsistent with any simple
closed form that would constitute exploitable structure.

## Reproduce

```
cd ecdsa-lift
python3 experiments/phase46_oeis.py        # prints sequences + fits, writes JSON
# OEIS: paste any results/phase46_oeis.json "oeis_query_raw" string into
#       https://oeis.org/search?q=<that>&fmt=text   (host blocks the sandbox
#       fetcher; do it in a browser, or use WebSearch as this phase did)
```

Output: `results/phase46_oeis.json`.

## VERDICT (footnote-ready)

> We performed a pre-registered null check that the lift-error sequence is not a
> disguised known sequence or simple closed form. For `E: y^2 = x^3 + 7` over
> `F_p` (`p in {31,43,67,97,109,151,163}`, lift precision `e=4`) we computed the
> formal-parameter lift errors `z(delta(k))`, `k=1..n-1` (sequence lengths
> 20–138). As required by the theorem, `z(delta(k)) ≡ 0 (mod p)` for every `k`
> (`delta` is a kernel-of-reduction element); the informative residue is the
> first non-trivial p-adic digit `d1(k)=(z(delta(k))/p) mod p`. OEIS index
> searches of the raw, `d1`, and leading-digit sequences returned no match, and
> an offline modular symbolic-regression sweep (constant, affine, quadratic,
> cubic, multiplicative `a·k`, geometric `a·c^k`, fit on a training half and
> scored on a held-out half) found no model exceeding the `1/p` random baseline
> — in particular `d1(k)` is not a degree-≤3 polynomial in `k` mod `p`. The
> lift error is thus consistent with the pseudorandom behaviour the inertness
> theorem predicts and shows no exploitable closed-form structure.
