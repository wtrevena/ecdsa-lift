# Referee-Grade Assessment: "When a Control Does Not Control" (ecdsa_lift_paper.pdf)

Prepared 2026-06-09. Framed as the review you should expect to receive, plus how to get ahead of it.

## 1. Verdict in brief

The mathematics is correct (I re-derived Lemma 1 and Theorem 1 line by line), the methodological lesson is real, and the writing is unusually honest about novelty. The two things most likely to sink it are (a) **venue fit** — the paper self-declares its theorem as folklore, so a mainstream crypto journal will likely reject on novelty regardless of quality — and (b) **one genuine internal inconsistency in Table 1** plus several statistical-reporting gaps that a careful referee will treat as symptomatic in a paper whose entire thesis is statistical discipline. All are fixable. Expected outcome at a results-oriented journal: reject. At a methodology/SoK-friendly venue with the fixes below: major revision → accept.

## 2. What the paper claims

A negative exploration (40 phases) of p-adic canonical-lift attacks on ECDLP for ordinary curves; a case study of a 5–11σ Gowers U² anomaly falsely certified by a permutation control; mechanism-derived nulls localizing it to the antisymmetry identity δ(k)+δ(n−k)=[n]τ(G); an inertness theorem (acknowledged folklore, after Silverman and Gadiyar–Padma) explaining why all such anomalies must be artifact-or-identity when gcd(n,p)=1; and the sharp boundary at anomalous curves (Smart's attack).

## 3. Correctness of the mathematics

**Verified sound** (re-derived by hand):

- Exact sequence (1), kernel of order p^{e−1}, additivity of z at e ≤ 4. ✓
- Lemma 1: injectivity (gcd(n, p^{e−1}) = 1), surjectivity (étale n-torsion Hensel-lifts), and the CRT idempotent G̃ = [m]τ(G). ✓
- Theorem 1(a),(b),(d) and Corollary 1: the algebra is exactly right, including δ(k+n)−δ(k) = [n]c = [n]τ(G) and the iff chain [n]c=O ⇔ c=O. ✓
- §7 anomalous degeneration: standard, correct (Smart/Satoh–Araki/Semaev). ✓

**Three precision issues a strong referee will catch:**

1. **Theorem 1(c) overclaims "information-theoretic independence."** The ECDLP instance fixes k ∈ Z/nZ; once a canonical representative k ∈ [1, n−1] is chosen (which is what every experiment in the paper does), k mod p^{e−1} — and hence [k]c — *is* mathematically determined by (E,G,Q). The clause "no value of [k]c is determined by (E,G,Q)" is therefore false as stated. The correct statement is computational, not information-theoretic: an oracle for [k]c on canonical representatives yields k mod p^{ord(c)} and hence (with k mod n) recovers k — i.e., computing [k]c is exactly DLP-hard. Alternatively, state independence for k uniform on Z/(n·p^{e−1}). The proof of (d) is unaffected; only the gloss in (c), the abstract, and Corollary 1 need restating.
2. **The antisymmetry (3) needs an explicit oddness hypothesis.** δ(k)+δ(n−k) = [n]τ(G) − (τ(kG)+τ(−kG)), so (3) holds iff τ(−P) = −τ(P). That holds for the x-Teichmüller/Hensel sections used (the two Hensel y-roots are negatives), but it is a property of the section, not universal. One sentence fixes this; without it, "exact, universal, public identity" is too strong.
3. **Formal log coefficient.** For y² = x³ + b: ω(z) = 1 + 3b z⁶ + O(z¹²), so log(z) = z + (3b/7)z⁷ + …, not (b/7)z⁷. Harmless to the argument (only the valuation matters at e ≤ 4), but it's checkable and wrong.

## 4. Empirical and statistical assessment

1. **Table 1 is internally inconsistent — this is the most important fix.** Rows p = 67, 97 are consistent with j=0 curves of subgroup order n = 79 (I verified by enumeration that both fields admit j=0 curves with N = 79; L = n−1 = 78 ✓). But 167, 257, 521 are all ≡ 2 (mod 3), so every y² = x³ + b over those fields is supersingular with N = p+1 (verified for 167, 257 by point counts) — there are no ordinary j=0 curves there. And if L = n−1, those rows force n = 167, 257, 521 = p, i.e., *anomalous* curves, contradicting the paper's own gcd(n,p) = 1 setting. Each row needs: the exact curve, n, and the definition of L. As printed, a referee can argue three of five rows are outside the paper's regime.
2. **No statistical methodology.** "5–11σ" and the z-ranges are given with no permutation counts, no description of how σ was estimated, no definition of the S¹-random baseline, and no multiplicity correction across 40 phases × digits × encodings. In a paper whose thesis is "controls must be derived from mechanisms," a referee will be merciless about this.
3. **Post-hoc exclusions.** Two smallest primes excluded ("dominated by baseline noise"), two more dropped from the shift identity ("deep multiply degenerates"). Both may be legitimate, but the criteria read as post hoc. State them as pre-specified rules with the flagging probe described — the irony of unprincipled exclusions in this particular paper will not be lost on reviewers.
4. **Missing positive control for H2 — the single best addition.** The half-sequence control is absence-of-signal evidence. Synthesize random sequences with the reflection constraint imposed and show they quantitatively reproduce the U² elevation, its growth with n, and the β ∈ [0.28, 0.47] range. That converts "the signal vanishes when reflection pairs are removed" into "reflection alone predicts exactly what we saw," and it directly practices the paper's own preached discipline.
5. **No artifact.** "The ecdsa-lift repository" has no URL or DOI, and 38 of 40 phases exist only as a summary paragraph. Venues that would accept this paper (see §5) require data/code availability.

## 5. Novelty and venue

The paper preempts the novelty objection honestly, but honesty doesn't change the bar: Journal of Cryptology, Designs Codes & Cryptography, etc. will likely reject a self-declared folklore theorem plus a methods lesson. Realistic targets:

- **IACR Communications in Cryptology** — explicitly welcomes SoK and practitioner papers, requires correctness, relevance, and available code/data. Plausible home *if* the artifact is published and §4's issues are fixed; pitch as methodology/SoK.
- **CFAIL** (Conference for Failed Approaches and Insightful Losses in Cryptology) — thematically a near-perfect fit; check current proceedings status.
- Journal of Cryptographic Engineering, or an ePrint/arXiv note if priority of the methodological framing is the goal.

The "AI-assisted exploration pipelines are false-positive engines" angle (§9) is the most timely and distinctive claim in the paper; consider promoting it from a discussion paragraph to a framing device — it is the strongest answer to "why does this need to be published?"

## 6. Minor issues

- Ref [10] is malformed (Gowers and Green–Tao lumped, no venue/year); the U²↔Fourier correspondence is classical and needs only a textbook citation (e.g., Tao–Vu, *Additive Combinatorics*).
- Ref [3] should cite the published version: Czechoslovak Mathematical Journal 68 (2018), 1115–1124, not just the arXiv preprint. Ref [1] is a talk; give a retrievable source or fold into [2].
- The automorphism-equivariance "rank law" in §3 is asserted, never defined or used again — explain or cut.
- "All 78 pairs" for p = 97: n = 79 gives 39 unordered reflection pairs; clarify the count.
- digit_j(·) is never precisely defined (base-p digits of z(δ(k)) indexed from which end?).
- Phase numbering (37, 40, 40b, 41, 42) appears without a map; a half-page phase table or appendix would ground §3's "forty phases."

## 7. The referee report you should expect

*Summary:* technically sound, well written, honest; methodological point is valuable; theorem is acknowledged folklore. *Major concerns:* Table 1 curve specification (three rows apparently outside the stated regime); no statistical methodology for the headline σ claims; no artifact despite the contribution resting on 40 unreported experiments; Theorem 1(c) misstated as information-theoretic. *Recommendation:* major revision at a methodology-friendly venue; reject at a results venue with advice to resubmit elsewhere.

## 8. Recommended actions, in priority order

1. Fix Table 1: per-row curve equation, n, definition of L; resolve the n = p contradiction for rows 3–5.
2. Publish the repository with a DOI and a phase index; CI-style script per table.
3. Add the synthetic positive control for H2 (reflection-constrained random sequences reproducing U² elevation and β).
4. Restate Theorem 1(c)/Corollary 1/abstract computationally ("DLP-hard to evaluate," not "information-theoretically independent").
5. Add a statistics subsection: permutation counts, σ estimation, baselines, multiplicity handling, pre-specified exclusion rules.
6. State the oddness hypothesis τ(−P) = −τ(P) for identity (3); fix the log coefficient to 3b/7; reference hygiene per §6.
7. Re-target the submission (CiC as SoK/methodology, or CFAIL), and consider leading with the AI-assisted-exploration framing.
