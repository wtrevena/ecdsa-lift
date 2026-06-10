# Adversarial Vet — "When a Control Does Not Control" (Trevena, June 2026), Round 1

**Reviewer scope:** ECC / p-adic arithmetic / additive combinatorics (Gowers norms) / statistics.
**Materials read in full:** `paper/ecdsa_lift_paper.tex` (and the 10-page compiled PDF), `notes/theorem.md`, `notes/literature_check.md`, `notes/phase37_anomaly_report.md`, `notes/phase40_resolution.md`, and the scripts `phase37/40/41/42/43/47` plus all nine committed `results/*.json`. Independent re-derivations and recomputations were written from the paper text alone in pure Python (`/tmp/verify_*.py`).

---

## (a) Verdict

**Minor revision.** (Venue-dependent: *accept after minor revision* at a methodology/negative-results venue — IACR Communications in Cryptology SoK track, CFAIL, JCEN; *reject on novelty/fit* at a results-oriented venue such as J. Cryptology or DCC, exactly as the paper itself concedes.)

The mathematics is correct, the central conceptual error flagged by the previous review (§8/Table 2 "the obstruction vanishes") is now **fully and correctly repaired**, and the paper's headline numbers reproduce to two decimals from the committed scripts and JSON. The remaining issues are presentational and artifact-hygiene, none of them threatening a result. There is **no blocking mathematical problem.** The single item I would insist on before camera-ready is an artifact-packaging defect (the data is git-ignored), not a scientific one.

---

## (b) The single most serious remaining issue

**The repository's `.gitignore` excludes every `results/*.json` file, so the data the paper repeatedly points to would not ship in the tagged release.**

`.gitignore` contains `results/*.json` with only `!results/.gitkeep` un-ignored; `git ls-files results/` returns *only* `.gitkeep`. Yet:
- the author footnote promises "**Code, data**, and a per-phase index: the `ecdsa-lift` repository (release tagged with this paper)";
- the Appendix caption (line ~630) states "Full scripts, fixed seeds, and **JSON outputs are in the `ecdsa-lift` repository**";
- the entire thesis of the paper is reproducibility discipline.

For a paper whose contribution *is* methodological rigor, shipping a release in which the result data is silently excluded is the worst possible look and will be caught by any artifact-evaluation committee. **Fix:** remove the `results/*.json` ignore (or force-add the nine JSONs and check them in), create the promised tag, and verify `git archive` of the tag contains the JSON. This is a 2-minute fix but it is load-bearing for the paper's only claimed contribution. Tagging it the *most serious* item is a comment on how little else is wrong, not on its intellectual depth.

---

## (c) Ranked issues

### [should-fix] 1. Original-sweep numbers in §4 have no committed artifact
§4 (line ~204) asserts "$5$–$11\sigma$ on the original exploratory sweep (Phase 37 …)" and "$\Uk^3$ and higher sat at baseline" (line ~208). **There is no `results/phase37_gowers.json`** in `results/` (only phase40 onward exist), and even when present it would be git-ignored (issue (b)). I *ran* `phase37_gowers.py`: it reproduces the claim qualitatively (U2 z ≈ 8–11 on p=577/631; U3 z ≈ −1.4 to +2.2, i.e. at baseline), so the claim is true — but it is presently unverifiable from the shipped artifact, and the prior review's point M3 ("include the U3 numbers") is therefore still open. **Fix:** commit a `phase37` JSON with the per-curve U2/U3 z-scores backing both the "5–11σ" and "U3 at baseline" sentences. (Note also: `random_baseline()` in `phase37` is *not* seeded inside `run()`, so phase37 z-scores are nondeterministic run-to-run — fine for a qualitative claim, but seed it before committing the JSON.)

### [should-fix] 2. README is stale and contradicts the submission
`README.md` line 83 still says the paper "builds to a **7-page PDF**"; the compiled PDF is **10 pages**. The prior review already flagged "7 pages" (it was then 9). It is now wrong by three pages and describes a release/tag that does not yet exist. **Fix:** update the page count and create the tag (ties into issue (b)).

### [should-fix] 3. §8 numerical claims are stated in the prose but not recorded in `phase42`'s JSON
§8 (lines ~480–484) asserts three specific facts about the anomalous curve that I verified by hand but that the committed `results/phase42_anomalous_boundary.json` **does not contain**: (i) `z([p]τ(G)) = 15·53` (= 795 mod 53²); (ii) section-independence across the Teichmüller-x, naive-x, and shifted-x sections; (iii) valuation 1. The JSON records only `obstruction_vanishes: true`, the 12/12 recovery, and the curve. The script `phase42` likewise computes `obstruction_vanishes()` ([p]c = O) but **never computes or emits the section-independence triple-check or the value 795** that the paper reports as verified. So three load-bearing sentences in the corrected §8 are asserted-as-verified but are not actually produced by the cited phase. I independently confirmed all three are correct (see (d)), but the artifact should produce them. **Fix:** extend `phase42` to compute and emit `z([p]τ(G))` under the three sections and assert equality, then cite that field.

### [should-fix] 4. The abstract carries a state-of-knowledge claim adjacent to the theorem
The revision correctly moved "no efficient algorithm … is known" *out of* Theorem 1 clause (c) and into the well-written Remark (lines 447–454, explicitly "confined to the present remark"). Good. But the **abstract** (line 61, "for which no efficient algorithm from $(E,G,Q)$ is known") and the **Contributions** bullet (line 112, "which no known efficient algorithm computes from it") still state the computational/state-of-knowledge claim in the same breath as the theorem's content, with no flag that this part is *not* part of the proven theorem. A careful referee re-reading the abstract after the Remark will note the asymmetry. **Fix:** one clause in the abstract/Contributions marking the "no known algorithm" half as a state-of-knowledge statement (mirroring the Remark), keeping the theorem's *proven* content (CRT-independence over the natural lift) cleanly separated.

### [nice-to-have] 5. Code inconsistency: synthetic-reflection draws differ between Phase 43 and Phase 47
In `phase43`, `synthetic_reflection_u2` draws `r = rng.integers(0, N)` — **uniform over all of $\Z/p^e$**, so the synthetic "kernel" elements can have a nonzero units digit (digit 0 ≠ 0), unlike every real $z(\delta(k))$ which is $\equiv 0 \bmod p$. In `phase47` the synthetic draws are correctly kernel-restricted (`RNG.integers(0, pe//p, h) * p`). I tested whether this matters: for digits $j\ge 1$ the two constructions give U2 indistinguishable to within $\pm0.2\sigma$ (see (d)), because the digit-$j$ ($j\ge1$) marginal is uniform either way. **So this is harmless to every reported number**, but it is a latent inconsistency in a paper preaching protocol discipline, and a referee who reads both scripts will ask. **Fix:** make the two synthetics identical (kernel-restricted) and note that digit-$j\ge1$ statistics are invariant to the choice.

### [nice-to-have] 6. "44σ" / σ language vs. 300-draw empirical nulls
The abstract says "$\sim$$44\sigma$" and Table 1's last row reports +42.7..44.8. §6 (lines 365–368) already adds the correct caveat ("$z$-values are standardized effect sizes against $300$-draw empirical nulls, not calibrated Gaussian tail probabilities; we use them only comparatively"). This addresses prior point S2. The one residual tension: a "44σ" deviation is *literally impossible* to certify from a 300-sample empirical null (the null's own SD is estimated from 300 draws, and the empirical max is finite), so quoting a bare "44σ" in the **abstract** — before the §6 caveat — slightly oversells. The number is real *as an effect size* (real/inc U2 sits ~44 SD above the increment-shuffle mean), but "σ" invites the tail-probability reading the paper disclaims. **Fix:** in the abstract write "an effect size of $\sim$$44$ (standardized against the matched null)" or similar, not "$44\sigma$."

### [nice-to-have] 7. Lemma 1 surjectivity step is terse
The proof (lines 387–394) asserts $E(\Z/p^e)[n]$ "surjects onto it (étale $n$-torsion lifts by Hensel)." This is correct (the multiplication-by-$n$ map is étale over $\Z_p$ since $\gcd(n,p)=1$, so $n$-torsion lifts uniquely), but the one-word "Hensel" is doing a lot of work for the load-bearing lemma. **Fix:** one extra clause — "since $[n]$ is étale at $p$ ($\gcd(n,p)=1$), the smooth/Hensel lifting of $E$ lifts each $n$-torsion point uniquely" — costs a line and removes a referee question. (Not blocking; the claim is standard and I verified the conclusion numerically.)

### [nice-to-have] 8. `notes/literature_check.md` mis-cites the very paper it is calibrating against
`literature_check.md` line 10 cites Silverman "Lifting and Elliptic Curve Discrete Logarithms" as "**LNCS 5747, 2009**." The correct volume is **LNCS 5381** (SAC 2008; DOI 10.1007/978-3-642-04159-4), which is exactly what the **paper's** bibliography gets right. So the *paper* is correct and the *internal note* is wrong — harmless to the submission, but worth fixing so the note doesn't mislead a future editor. (I verified LNCS 5381 via Springer.)

### [nice-to-have] 9. Minor robustness bug in `phase37_gowers.py`
`phase37` throws `ValueError("Could not factor 1036488922561 …")` on p=1009 (a brittle ad-hoc prime-extraction in `z_of_delta`). It does not touch any reported curve, but it will trip anyone re-running the sweep. **Fix:** pass `p` explicitly rather than re-deriving it by factoring $N=p^e$.

---

## (d) What I verified computationally, and the result

All checks pass. Pure-Python reimplementations written from the paper text, different code path from the repo where feasible.

| Claim | Location | Independent result |
|---|---|---|
| Table 1 curve orders ($y^2=x^3+7$, $p\equiv1\bmod3$, $N$ prime, $\ne p,p{+}1$) | Table 1 | All seven exact: $N=$ 79,139,199,313,397,613,829, all prime, gcd$(N,p)=1$. ✓ |
| Phase 47 curve = smallest $p\ge10^4$, $p\equiv1\bmod3$, $N$ prime | §5/§9/Tab.1 | $p=10477$, $N=10639$ prime. ✓ |
| §3 example $p=97$ gives $n=79$ (distinct curve from $p=67$) | §3 | $N(97)=79$. ✓ |
| Antisymmetry $z(\delta(k))+z(\delta(n-k))$ single-valued $=z([n]\tau(G))$ | (3), §3 | $p=97$: single value over all 39 pairs $=z([n]\tau(G))$, val 1. $p=67$: same, $C_n=12151388$ (matches `phase41` JSON). ✓ |
| **§8 anomalous boundary**: $z([53]\tau(G))=15\cdot53=795\ne0$, val 1 | §8 | Exactly 795. ✓ |
| **§8 section-independence** of $[p]\tau(G)$ at $e=2$ | §8/Tab.2 | Teichmüller-x, naive-x, x+17p all give 795. ✓ |
| **§8 ambiguity dies**: $[p]c=O$ for kernel $c$ at $e=2$ | §8/Tab.2 | $z(c)\in p\Z$, $[p]c=O$. ✓ |
| §8 homomorphism $\varphi(P)=\frac1p\log([p]\tau(P))$ linear in secret | §8 | $\varphi(G)=15$, $\varphi(kG)=15k\bmod53$, additive $\varphi(P{+}Q)=\varphi(P){+}\varphi(Q)$. ✓ (the $\varphi(G)=15$ matches the paper exactly) |
| **Does §8 need $e=2$?** | §8 ("take $e=2$") | **Yes for the literal statement.** $[p]c=O$ holds for *every* kernel $c$ only at $e\le2$; at $e=4$ a kernel $c$ of order 2 survives ($[p]c\ne O$). So "take $e=2$" is mathematically *necessary* for "$[p]c=O$ for every $c$," not cosmetic. (Subtlety the paper need not but could note: the *attack invariant* $\varphi$ — the first digit of $\log[p]\tau$ — stays section-independent and $=15$ at $e=3,4$ too, because the surviving $[p]c$ lives in $p^2\Z$; so Smart's attack itself is robust above $e=2$, but the clean "$[p]c=O$" phrasing genuinely requires $e=2$.) The paper's restriction is correct. ✓ |
| Lemma 1: $\widetilde G=[m]\tau(G)$, $m\equiv1\bmod n$, $m\equiv0\bmod p^{e-1}$ | Lem.1 | Order $n$, reduces to $G$, $\mathrm{red}([2]\widetilde G)=2G$, map $jG\mapsto[j]\widetilde G$ a homomorphism so $\delta_{\mathrm{can}}\equiv0$. ✓ |
| Theorem 1(a),(b),(d) algebra | Thm.1 | $\delta_\tau(k)=[k]c-c_k$; $[k]c$ via $k\bmod p^{e-1}$; shift $=[n]c=[n]\tau(G)\ne O$; iff chain $[n]c=O\Leftrightarrow c=O$. All hold. ✓ |
| Table 1 all five $z$-columns, all 7 curves, 3 digits | Table 1 | Reproduce `phase43` JSON **cell-for-cell** (e.g. $p=823$ real/inc $+10.5..12.5$, syn/rand $+12.2..12.5$). ✓ |
| §4 range "$z$(real vs rand) ranges 3.3 to 13.6" | §4 | JSON min/max = 3.31 / 13.62. ✓ |
| §5 ranges (half/inc $-1.3..1.5$; syn/rand $+3.8..12.5$; real/syn $-1.4..1.5$) | §5 | JSON: $-1.34..1.49$; $3.78..12.46$; $-1.35..1.51$. ✓ |
| Fourier decay $\beta_{\mathrm{real}}=\{0.34,0.45,0.36\}$, $\beta_{\mathrm{syn}}=\{0.41,0.41,0.40\}$ | §5 | Refit from JSON peak-Fourier: $\{0.34,0.45,0.36\}$ and $\{0.41,0.41,0.40\}$ — exact to 2 dp. ✓ |
| Phase 47 Table-1 last row vs JSON | Tab.1 | real/inc $+42.7..44.1$, inc/rand $-0.0$, half/inc $-1.1..1.8$, syn/rand $+44.5..44.8$, real/syn $-2.0..0.7$ — all match JSON. ✓ |
| Phase 41 order list $\{21,31,79,67,79,111,129,133,86,139\}$, 5 composite, 2 arith-degenerate primes | §9 | Exact; composite $=\{21,111,129,133,86\}$; deep-multiply degenerate on $p=43,103$ only (= "two primes"). ✓ |
| Smart 12/12 on $y^2=x^3+4x+7/\F_{53}$, $G=(0,22)$ | §8/Tab.2 | `phase42` JSON: 12/12, curve/G match. ✓ |
| H1 prefix-sum transfer function $1/(1-e^{2\pi i\xi/n})$, low-pass pole at $\xi=0$ | §5 H1 | Standard summation filter; correct. ✓ |
| Phase 43 vs Phase 47 synthetic-draw discrepancy is benign | (issue 5) | U2 of full-range vs kernel-restricted synthetics differ by $<0.2\sigma$ for digits 1,2,3. ✓ |
| Citations: Silverman SAC'08 = LNCS 5381; Gadiyar–Padma CMJ 68 (2018) 1115–1124, arXiv:1702.07107 | refs | Both correct in the **paper** (note in `literature_check.md` is wrong, issue 8). ✓ |

**Net:** every checkable mathematical and numerical claim in the paper is correct and reproduces. The only things that *don't* reproduce from the shipped artifact are claims whose backing JSON is absent/git-ignored (issues 1, 3) — the underlying computations, when re-run, confirm the prose.

---

## (e) What the paper gets impressively right

1. **The §8/Table 2 correction is exactly right, and subtle.** The previous reviewer's "major point" — that "$[n]\tau(G)$ vanishes" at the boundary is false — is now repaired everywhere (abstract, §8, Table 2 caption "What vanishes … is the *ambiguity*, not the invariant"). And the corrected statement is genuinely correct at the level of detail few authors would get right: $[p]\tau(G)\ne O$ (=795), section-independent, secret-linear, while $[p]c=O$ for kernel $c$ — all four verified independently. The paper even correctly restricts to $e=2$, which is *necessary* (not optional) for the literal "$[p]c=O$ for every kernel element."

2. **Reproducibility is exceptional.** Cell-for-cell agreement of Table 1, two-decimal agreement on the fitted $\beta$ exponents, exact match of the Phase-41 order list and the "two degenerate primes," and the Phase-47 row all reproduce from the scripts. This is far above the norm.

3. **The positive control is the right move and is honestly bounded.** Synthesizing reflection-only sequences and showing they reproduce the U2 elevation *and its growth with $n$* — then explicitly **declining** to claim the digit-to-digit $\beta$ variation is matched (and naming carry-nonlinearity as the likely residual) — is precisely the discipline the paper preaches, applied to itself. The $\beta_{\mathrm{real}}$ vs $\beta_{\mathrm{syn}}$ mismatch is disclosed rather than buried.

4. **The "control that cannot separate its hypotheses" framing is genuinely instructive** and correctly argued: the permutation test collapses H1/H2/H3 simultaneously, so it cannot attribute the signal; the H1 (increment-shuffle), H2 (half-sequence), and H2-positive (synthetic-reflection) nulls each isolate one mechanism. The logic "controls are designed to *remove* signal, not to find it, so no multiplicity correction is needed for a discovery we are not claiming" (§6) is sound.

5. **Honest novelty calibration.** The theorem is repeatedly and correctly framed as an explicit instance of Silverman and Gadiyar–Padma, not new mathematics; `literature_check.md` reaches the same conclusion before drafting. The citations are accurate.

6. **The boundary-sharpness demonstration is a clean, falsifiable hinge:** the *same code path* is inert for $\gcd(n,p)=1$ and a total break at $n=p$, with the hypothesis $\gcd(n,p)=1$ as the literal dividing line — and the 12/12 recovery backs it.
