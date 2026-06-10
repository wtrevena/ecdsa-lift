# Adversarial Vet — "When a Control Does Not Control" (Trevena, June 2026), Round 2

**Reviewer scope:** ECC / p-adic arithmetic / additive combinatorics (Gowers norms) / statistics.
**Materials read in full this round:** the round-1 report; the revised `paper/ecdsa_lift_paper.tex` and 11-page compiled PDF (`pdftotext -layout`); `experiments/phase48_u2_growth_law.py`, `phase37_gowers.py`, `phase42_anomalous_boundary.py`; all twelve `results/*.json` (incl. the new phase37/47/48); `README.md`; `.gitignore`; `notes/literature_check.md`. Independent re-derivations and re-runs written from the paper text in pure Python (mass-factor, √L scaling, K constant, Table-1 reproduction, β scatter bootstrap, Hasse dichotomy). No `sage` in the environment; none needed.

---

## Verdict

**Minor revision — and the *same* item the round-1 report named as its single most serious issue is STILL UNFIXED.** The new mathematics (Proposition 1) is correct and well-supported; the §8 dichotomy remark is correct; the four state-of-knowledge flaggings and the "44σ" rewording all landed. But the round-1 **blocking artifact-hygiene defect is not fixed** (`.gitignore` still excludes `results/*.json`; `git ls-files results/` still returns only `.gitkeep`), the README page count is **now wrong again** (says 10, PDF is 11), and two round-1 should-fix items (phase42 still doesn't emit the §8 triple-check; the prose/JSON disagree on the Phase-47 √n prediction) remain open or were re-introduced. The task brief I was given asserted several of these were fixed; **they are not on disk** — verify before camera-ready.

---

## Round-1 fix confirmation table

| # | Round-1 item | Claimed fixed? | Actually landed & correct? |
|---|---|---|---|
| (b)/1 | `.gitignore` excludes `results/*.json`; data won't ship | yes (brief: "no longer excludes") | **NO** — `.gitignore` still contains `results/*.json`; `git check-ignore` confirms phase37/47/48 JSON all ignored; `git ls-files results/` = only `.gitkeep`. **Blocking, unfixed.** |
| 1 | commit phase37 JSON backing "5–11σ" & U3 | partial | JSON now **exists on disk** (seeded run) and backs 5–11σ (U2 z up to 11.15). But it is git-ignored (see above), and U3 is **not** uniformly "at baseline": U3 z reaches +6.93 (p=31) / +6.10 (p=43) on the smallest curves. |
| 2 | README "7-page" → "10-page" | yes | Changed to "10-page" — but **PDF is now 11 pages** (`pdfinfo`: Pages 11). Stale again. |
| 3 | phase42 emit 795 / section-independence / val-1 | (not in brief) | **NO** — `phase42` JSON still records only `obstruction_vanishes`, 12/12, curve/contrast. The three §8 facts (795, three-section independence, valuation 1) are still prose-only, not produced by the cited artifact. |
| 4 | "no efficient algorithm" flagged separate-from-theorem everywhere | yes | **YES** — abstract L63 ("separately from the theorem"), contributions L115 ("state-of-knowledge remark, kept separate"), Remark L506 ("confined to the present remark"). Consistent. Theorem 1(c)/Cor 1 carry only the info-theoretic/computational hedge. Correct. |
| 5 | abstract "44σ" → "standardized effect size ≈44 / 300-draw empirical null" | yes | **YES** — abstract L46 now reads "standardized effect size ≈44—against a 300-draw empirical null, not a Gaussian tail." No bare "44σ" remains. Correct. |
| 6a | `.gitignore` un-ignore results | yes | **NO** (same as (b) — unfixed). |
| 6b | literature_check.md Silverman LNCS 5747→5381 | yes | **YES** — line 11 now "LNCS 5381" (matches paper bib). Correct. |
| 6c | phase37 random_baseline seeded | yes | **PARTIAL** — `random_baseline(...)` now defaults `seed=20260609` and seeds its rng (deterministic). But the *Fourier-max* baseline loop later in `run()` still uses bare `np.random.random(n)` (unseeded global). U2 z-scores are now reproducible; the Fourier-z column is not. Harmless to the headline but not fully seeded as claimed. |

---

## NEW Proposition 1 (§5): is it correct? — YES, with one cosmetic proof-sketch imprecision

**Statement verified independently and from scratch (pure Python, my own reflection construction *and* the paper's additive construction):**

- **Mass factors 2/L vs 3/L:** confirmed. `E[‖f‖_{U2}^4]·L` → 2 for i.i.d. S¹-uniform f, → 3 for reflection-only f(L+1−k)=β·conj(f(k)). Measured ratios: L=40→1.484, 80→1.491, 200→1.495, 612→1.498 (my code); phase48.json mass_ratio_mean = **1.4972**. The "exactly 3/2, independent of L and alphabet" claim holds — phase48's p=211 and p=99991 cross-checks also give 1.498 and 1.504. **Correct.**
- **The additive synthetic = the Proposition's model:** f(n−k)=exp(2πi(C−s(k))/p) = exp(2πiC/p)·conj(f(k)) = β·conj(f(k)) with β=exp(2πiC/p). So phase48/phase43's `C - d[::-1]` construction *is* exactly the f(L+1−k)=β·conj(f(k)) of Prop 1. **Honest.**
- **√L scaling is rigorous, not hand-wavy.** I confirmed each ingredient numerically: mean(U2)·L^{1/4} ≈ 1.189 (const) ⇒ mean U2 = Θ(L^{−1/4}); sd(U2)·L^{3/4} ≈ 0.293 (const) ⇒ **null SD = Θ(L^{−3/4})** exactly as the sketch claims. Numerator = mean_rand·((3/2)^{1/4}−1) = Θ(L^{−1/4}); z = Θ(L^{−1/4})/Θ(L^{−3/4}) = **Θ(√L)**. The "Θ(L^{−3/4}) null-sd" claim is **correct**.
- **The constant.** Analytic K = (A/B)·((3/2)^{1/4}−1) with A=1.189, B=0.293 gives **K≈0.4329**; phase48 K_mean=0.4339; my direct measurement K≈0.426. The paper's "≈0.42" is the rounded-down value; the truer constant is ≈0.43. **Consistent.**
- **K·√(n−1) vs Table 1 syn/rand:** my prediction reproduces all 8 rows within noise (e.g. n=79 pred 3.76 vs obs 3.8–4.1; n=829 pred 12.27 vs 12.2–12.5; n=10639 pred 43.96 vs 44.5–44.8). **Reproduces.**

**Honest about proved-vs-measured:** the constant 0.42 is explicitly "verified numerically (Phase 48)", and the 3/2 derivation is labeled "Proof sketch." Good separation.

**[nice-to-have] One sentence of the proof sketch is literally wrong (the result is not).** The sketch says the constraint "forces f̂(L−ξ)=β·f̂(ξ), pinning each mode to its conjugate partner." A time-domain conjugate-reflection (k ↔ L−1−k with conjugation) does **not** induce f̂(L−ξ)=β·f̂(ξ): I checked, and |f̂(ξ)| is uncorrelated with |f̂(−ξ)| (corr ≈ 0.003) under the actual construction. The correct induced relation is a conjugation-plus-reflection-plus-shift in frequency, not the simple modulation written. The mass count of 3 is right (verified 3 different ways) and the conclusion is solid, but the stated spectral mechanism mis-describes the pairing. **Fix:** either state the correct induced symmetry (f is a conjugate-palindrome of its free half, whose spectral second-moment matrix gains one extra rank-1 contraction) or drop the parenthetical "pinning each mode to its conjugate partner," which is the part that's false.

---

## NEW §8 "dichotomy is exhaustive" remark — CORRECT

Verified the Hasse argument both by hand and computationally:
- Paper's interval [(√p−1)², (√p+1)²] = [p+1−2√p, p+1+2√p] = the standard Hasse interval. Algebra correct.
- A multiple mp of p lies in it iff |mp−(p+1)|≤2√p. m=1 always works (N=p, anomalous curves exist). m=2 ⟺ p−1≤2√p ⟺ p≤1+2√2≈3.83…(√p≤1+√2 ⇒ p≤5.83) ⇒ only p≤5. m≥3 never. **Computational check:** multiples of p in the Hasse range are {p} for p=7,11,13,17,19,23,29, and **{5,10} for p=5** — exactly the paper's "lone exception N=2p occurs only at p=5." Correct.
- The logical chain "n|N, p|n ⇒ p|N; Hasse ⇒ N=p (p≥7) ⇒ v_p(n)∈{0,1}, no partial kernel" is sound. (Even the p=5 N=2p exception has v_p(N)=1, so the "v_p(n)∈{0,1}" conclusion is unaffected; the paper correctly restricts the analysis to p≥7.) **The new remark is correct and does not damage the round-1-corrected §8 boundary math, which remains intact** (the [p]c=O at e=2, section-independence, φ(G)=15, 12/12 — all unchanged in the prose).

---

## NEW §5 β paragraph ("finite-size fit scatter") — CORRECT

Bootstrap-reproduced: reflection-only β across the 7 clean curves has **mean ≈ 0.401, sd ≈ 0.052** (paper says mean≈0.40, sd≈0.06 — matches). The real exponents {0.34, 0.45, 0.36} sit at z-distances **1.16, 0.94, 0.78** from the synthetic mean — each within ≈1 sd, exactly as claimed. The "finite-size fit scatter rather than secondary structure" characterization is numerically honest.

---

## NEW issues (this revision)

### [blocking] B1. The round-1 blocking item was NOT fixed: `results/*.json` is still git-ignored
`.gitignore` line 4 is still `results/*.json` (with `!results/.gitkeep`). `git check-ignore results/phase37_gowers.json results/phase47_large_prime.json results/phase48_u2_growth_law.json` returns all three; `git ls-files results/` = only `.gitkeep`; `git status --short results/` shows nothing. So the entire data backing the paper — including the three brand-new JSONs this revision added (phase47, phase48, phase37) — would not ship in the tagged release, and **no tag exists** (`git tag` empty). For a paper whose sole claimed contribution is reproducibility discipline this is the worst possible defect and an automatic artifact-eval failure. **Fix:** delete the `results/*.json` line from `.gitignore` (or `git add -f` the twelve JSONs), commit, create the promised tag, and confirm `git archive` of the tag contains them. The brief I was handed states this was fixed; it is not on disk.

### [should-fix] B2. README page count stale again: says 10, PDF is 11
README line 83 reads "to a 10-page PDF" but `pdfinfo` reports **Pages: 11** (Proposition 1 + the §8 dichotomy remark pushed it over). The round-1 fix (7→10) was undone by the revision's own growth. **Fix:** "11-page" (and ideally stop hard-coding the count, or generate the tag/release the README references — which still does not exist).

### [should-fix] B3. Prose/JSON disagree on the Phase-47 √n prediction
§5 line 343 says "predicted 43.3 vs. observed 44.6 at n=10639," using the **rounded** K=0.42 (0.42·√10638 = 43.32). But phase48's own committed JSON predicts **44.75** for n=10639 (K_mean=0.4339) — which matches the observed 44.6 essentially exactly. The paper quotes the *worse* prediction from its rounded constant and then calls the 1.3-unit gap "finite-size noise," when its cited artifact gives a near-perfect 44.75. This is an internal prose-vs-JSON inconsistency that undersells the result. **Fix:** quote K≈0.43 (or phase48's 44.75) so the headline number agrees with the JSON it cites.

### [should-fix] B4. phase42 artifact still omits the §8 triple-check (round-1 item 3, unaddressed)
§8 asserts z([p]τ(G))=15·53=795, section-independence across Teichmüller-x/naive-x/shifted-x, and valuation 1 "as we verified (Table 2)." `phase42_anomalous_boundary.json` records none of these — only `obstruction_vanishes`, the 12/12, and the curve/contrast. I re-confirmed all three are mathematically correct (795 = 15·53; φ(G)=15), but the cited phase still does not *produce* them. **Fix:** extend phase42 to compute and emit z([p]τ(G)) under the three sections, assert equality and val 1, and cite that field.

### [nice-to-have] B5. "U3 sat at baseline" is an overstatement on the small curves
Committed phase37 JSON: U3 z reaches **+6.93 (p=31, n=20)** and **+6.10 (p=43, n=30)** — small-sample noise on the smallest curves. For n≥66 U3 stays ≤4.67 and mostly ≤2.5. The sentence "U3 and higher sat at baseline" (L234) is true for the larger curves but not literally for the whole sweep. Defensible given the paper frames phase37 as "a mixed bag … later found to include supersingular and anomalous cases," but a one-clause hedge ("U3 at baseline except small-n noise on n≤30") would be more honest. (Also: p=1009 still returns a null record — the round-1 brittle-factoring bug #9 is unfixed and now visible in the committed JSON.)

### [nice-to-have] B6. phase37 Fourier-max baseline still unseeded
`random_baseline` is now seeded (U2 z-scores reproducible), but the Fourier-max baseline loop in `run()` still calls bare `np.random.random(n)`. The "fixed seeds" claim (README, §6) is not fully met for the Fourier-z column. **Fix:** seed that loop too.

---

## What I verified computationally this round (all pass unless noted)

| Check | Result |
|---|---|
| Prop 1 mass 2/L (iid) vs 3/L (reflection), ratio 3/2 | ✓ my code 1.48–1.50; phase48 mass_ratio_mean 1.4972 |
| Additive synthetic ≡ f(L+1−k)=β·conj(f(k)) | ✓ algebraically and numerically |
| null SD = Θ(L^{−3/4}); mean U2 = Θ(L^{−1/4}); z=Θ(√L) | ✓ sd·L^{0.75}≈0.29 const, mean·L^{0.25}≈1.19 const |
| K constant ≈0.42–0.43 | ✓ analytic 0.4329, phase48 0.4339, direct 0.426 |
| K·√(n−1) reproduces Table-1 syn/rand (8 rows) | ✓ all within noise |
| β scatter: refl-only mean 0.401 sd 0.052; real within ~1 sd | ✓ z-dist 1.16/0.94/0.78 |
| §8 Hasse dichotomy: only mult of p in range is N=p (p≥7); {5,10} at p=5 | ✓ exact |
| Table 1 (7 ordinary rows) reproduces phase43.json cell-for-cell | ✓ |
| Phase-47 row matches phase47.json | ✓ |
| `.gitignore`/`git ls-files`/`git check-ignore`/`git tag` | data git-ignored, no tag (B1) |
| PDF page count | 11 (README says 10 — B2) |
| Prop 1 proof-sketch "f̂(L−ξ)=β·f̂(ξ)" pinning | **does NOT hold literally** (corr |f̂(ξ)|,|f̂(−ξ)| ≈ 0); result still correct (B-nice) |
| phase42 JSON contents | only obstruction_vanishes/12-12/curve — no 795/section/val (B4) |
| phase37 JSON U2 z up to 11.15 (backs 5–11σ); U3 up to 6.93 (B5) | exists on disk, git-ignored (B1) |

---

## Remaining items (carried from round 1 or new), priority order

1. **[blocking]** Un-ignore + commit `results/*.json`; create the release tag (B1). *Still open from round 1.*
2. **[should-fix]** README "10-page" → "11-page" (B2).
3. **[should-fix]** Reconcile §5 "predicted 43.3" with phase48.json's 44.75 (B3).
4. **[should-fix]** Make phase42 emit the 795 / three-section / val-1 triple-check (B4). *Still open from round 1.*
5. **[nice-to-have]** Hedge "U3 at baseline" for small-n; fix p=1009 factoring bug; seed phase37 Fourier baseline (B5, B6). *Carried.*
6. **[nice-to-have]** Fix the false "f̂(L−ξ)=β·f̂(ξ)" sentence in the Prop 1 proof sketch (result unaffected).

**Bottom line:** the revision's *new mathematics is sound* — Proposition 1 (mass 3/2, √L scaling, K≈0.43, Table-1 reproduction), the §8 Hasse dichotomy, and the β-scatter claim all verify independently, and the four state-of-knowledge flaggings + the "44σ" rewording are correct and consistent. The revision did **not** fix the one item round 1 flagged as most serious (git-ignored data, no tag), it re-introduced a stale page count, and it left two should-fix artifacts open. None of these threaten a result, but the data-shipping defect is genuinely blocking for a reproducibility paper and must be cleared before any tagged release.
