# Phase 44 — Multi-instance / cross-curve statistical aggregation

*The one attack lane the inertness theorem (notes/theorem.md) does not
directly cover. Result: DEAD. A clean, expected null.*

## The question

The inertness theorem proves the lift error `δ(k) = [k]τ(G) − τ(kG)` is
*per-instance* inert for `gcd(n,p)=1`: the only public structure is

- **(H1)** an `Aut(E)`-equivariance of the leading base-`p` digit of `z(δ)`;
- **(H2)** the antisymmetry pairing `z(δ(k)) + z(δ(n−k)) = Cn = z([n]τ(G))`
  (Phase 21b / 40b; the *entire* Phase-37 Fourier anomaly was its shadow).

The theorem is about a *single* instance. It says nothing, on its own,
about whether the **joint distribution of `δ`-features across a large
ensemble of curves** carries a secret-correlated signal that only emerges
in aggregate — the differential/linear-cryptanalysis analogy, where single
traces look random but pooled statistics leak. Phase 44 tests exactly that.

## Design (pre-specified, fixed seed 20260609, e = 4)

- **Ensembles.** (A) `j=0` CM family `y²=x³+7`, `p≡1 mod 3`, `p∈31..200`,
  ordinary only (exclude `N=p` anomalous and `N=p+1` supersingular) — 16
  curves. (B) generic non-CM random `(a,b)`, `j≠0,1728` (trivial `Aut`),
  drawn deterministically — 8 curves. Ensemble (B) is the clean cross-curve
  test because **(H1) is vacuous there**.
- **Features of `z(δ(k))`** (kernel element ⇒ `z≡0 mod p`): base-`p` digits
  1,2,3; `z mod p²`; `p`-adic valuation.
- **Labels.** full secret `k`; coarse residues `k mod 2,3,5`; normalized
  secret bucket `knorm = ⌊8k/n⌋` (curve-independent, used for pooling).
- **Public-identity quotient (the critical control).**
  - **(H0)** drop the trivial `k=1` (`δ(1)=0` identically — see below).
  - **(H2)** restrict to the first half `k ≤ ⌊(n−1)/2⌋` (the Phase-40b
    control; within one half the reflection imposes no internal constraint).
  - **(H1)** report a `beyond_H1` feature set that **drops the leading
    digit** `f_digit1`; the headline uses first-half + beyond_H1, and is
    corroborated on the non-CM ensemble where (H1) is absent.
- **MI estimation, two independent estimators.** (i) histogram / plug-in MI
  (exact discrete MI, the workhorse, per-curve and pooled, null = 300 label
  permutations); (ii) sklearn `mutual_info_classif` **continuous kNN**
  (Kraskov) as a genuinely-different pooled cross-check (null = 120 perms).
  z-score = `(MI_real − mean(MI_null)) / std(MI_null)`.
- **Falsifier.** if pooled MI beyond the public identities sits at the
  shuffle null floor, the lane is dead.

## Headline numbers (first-half, beyond_H1)

| ensemble | pooled max |z| hist | pooled max |z| kNN | per-curve max |z| |
|----------|--------------------:|-------------------:|------------------:|
| j0 (16)  | 1.74                | 1.26               | 2.26              |
| noncm (8)| 2.16                | 0.63               | 2.09              |

Global maxima: **pooled |z| = 2.16, per-curve |z| = 2.26.** Over ~32 pooled
+ ~120 per-curve tests, the expected null max |z| is ≈ 2.6–3.0. Everything
sits at or below the null floor; the largest value (noncm `kmod3`, z=−2.16)
is a *negative* fluctuation (MI below the null mean), i.e. anti-signal.

**VERDICT: DEAD** — pooled MI beyond the public identities is at the
shuffle null floor. Falsifier triggered.

## The one false alarm, and why it was false (kept for the record)

The *first* run (before the (H0) quotient) reported a WEAK-SIGNAL: the
`j0 / knorm / histogram` pooled cell at **z ≈ 3.4**, stable across re-seeds
and surviving a within-curve blocked shuffle. Skeptical follow-up killed it:

1. The independent **kNN estimator on the identical data gave z ≈ 0.9** — a
   flat null. Real MI shows in both estimators; plug-in-only elevation
   signals plug-in bias, not signal.
2. Decomposing by feature, the entire bump was **`f_val` (the valuation),
   z ≈ 13** — while every other feature was at the floor.
3. The valuation distribution is near-constant: every `k` has `v_p(z(δ))=1`
   except **exactly one `k` per curve with `v=e=4`**, and that point is
   **always `k=1`, on every single curve** — because `δ(1) = [1]τ(G) −
   τ(G) = 0` *by definition*. Pooling this fixed, secret-free `k=1` point
   across 15 curves is what the plug-in MI was reading.
4. **Excluding `k=1` collapses all three estimators to the floor**:
   hist z 3.54→0.37, blocked-null hist z 3.82→0.62, kNN z 0.33.

So `δ(1)=0` is just the most trivial public identity of all. It is now
quotiented out (H0). This episode is the phase's main methodological
lesson: a plug-in MI bump that the kNN estimator does not see is almost
always a degenerate/near-constant feature, here a single public fixed point.

## Assessment (prior-updated)

The repo's prior was a strong null (30+ years of literature; Phases 1–43
ruling out the per-instance angle). Phase 44 extends that null to the
**cross-curve aggregation** lane, the last elementary statistical avenue the
theorem did not formally close, on two ensembles (CM and non-CM), two
estimators, with the known public identities (H0 trivial, H1 leading-digit,
H2 antisymmetry) properly quotiented out. No secret-correlated structure
survives pooling. Prior essentially unchanged: the lane is dead.

## Cheapest follow-up if one wanted to push further (none indicated)

Nothing here justifies follow-up. If a future referee insists: (a) scale the
ensemble to `p≲1000` and ~100 curves to tighten the null-max bound (cheap,
hours); (b) replace the per-feature MI sum with a *joint* multivariate
estimator (sklearn kNN already does the continuous joint and saw nothing);
(c) test higher-order pooled statistics (pairwise feature products) — but
the theorem's CRT direct-sum argument (Phase 38) predicts these are
independent across curves, so the expected outcome is again null.

## Reproduce

```
pip install scikit-learn --break-system-packages   # adds sklearn + scipy
python3 experiments/phase44_aggregation.py          # ~30 s, writes the JSON
# -> results/phase44_aggregation.json
```
Fixed seed 20260609; deterministic ensembles; the JSON carries every
per-curve and pooled MI/null/z plus the verdict.
