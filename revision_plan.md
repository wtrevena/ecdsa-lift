# Unified Revision Plan

**Paper:** *When a Control Does Not Control: A Mechanism-Derived Refutation of a False Positive in p-adic Cryptanalysis of the ECDLP* (v2, June 2026)
**Sources merged:** Referee Report 1 (`referee_report_ecdsa_lift.md`: M1–M5, minor points, next-steps §7, strengthening appendix A–F) and Reviewer 2 (`reviewer_feedback.txt`, 14 sections).
**Consensus verdict:** Both reviews say the same thing — technically sound, unusually honest, correctness verified where checked; the problems are framing, scope, calibration, and reproducibility, not mathematics. R1: minor-to-major revision. R2: major revision. Plan below targets the union of both.

---

## 0. Reviewer 2 calibration — points already satisfied in the manuscript

R2 reviewed the paper *as reconstructed from Referee Report 1*, not the manuscript itself (their own opening line). The following R2 asks are already done in the paper; for each, the fix is a cheap visibility edit, not new work:

| R2 ask | Already in paper | Remaining action |
|---|---|---|
| Separate computational claim from theorem (§4) | Remark after Thm. 1 does exactly this | Adopt R2's suggested wording at each *invocation* (see item 2) |
| Don't read 44σ as a Gaussian tail (§6B) | §6: "standardized effect sizes against 300-draw empirical nulls, not calibrated Gaussian tail probabilities" | Add raw exceedance (0/300) per R1-B8; move the disclaimer earlier |
| Phase 47 not cherry-picked (§6) | Reported separately, marked †, labeled post-hoc | Reframe as successful out-of-sample prediction of the √n law (R1-F6) — stronger than bracketing |
| "Not contradicting Smart" (§4) | §8 *is* the Smart boundary; Table 2 | Add the one explicit sentence R2 suggests |
| Multiplicity disavowal (§6C) | §6: "no multiplicity correction… do not claim discovery" | Quantify it (item 8) rather than just disavow |
| ECDSA-scale caveat exists (§9) | §7/§9: "relevance… through the theorem's generality" | Move to abstract/intro (item 6) |

Do not concede these as missing in any response letter — cite the manuscript locations.

---

## 1. Decision to make first: venue and identity

Everything else depends on this (R1-M1/E2, R2 §3/§12). Both reviews agree the paper fails at a venue expecting new cryptanalysis and succeeds at one accepting negative results/methodology (IACR CiC, or a methodology/AI-for-science venue). Decide the canonical venue **before** editing, then:

- Retitle toward the methodology contribution (R2: "mechanism-resolved controls for automated cryptanalytic false positives" register, not "refutation of p-adic cryptanalysis").
- Rewrite abstract + intro to claim exactly: (i) worked localization of a false positive to one public identity; (ii) positive-control attribution; (iii) closed-form U² artifact signature (3/2 mass, √n law). Nothing more (R1-M1, R2 §13.1).
- Decide whether to state the AI-assisted nature of the exploration explicitly (R1-E1, R2 §12). If yes, add a short subsection: what the pipeline did, what it concluded, where the human audit intervened. Both reviews flag this as the timeliest angle.

## 2. Required: theorem language (R1-M2, R2 §4, §13.2)

- At every use of Theorem 1, say whether it leans on the proved well-definedness/CRT clause or the state-of-knowledge Remark; soften "forecloses."
- Adopt R2's formulation nearly verbatim where the paper currently says "Theorem 1 forecloses its being cryptanalytic": *"…has no well-defined secret-dependent content available from public data, except in the anomalous boundary case. Under the standard assumption that ECDLP remains hard on ordinary non-anomalous curves, observed structure should be presumed identity-or-artifact unless it survives section-invariance and mechanism-derived controls."*
- Add the explicit "recovers Smart as the boundary case, does not contradict it" sentence in §8.

## 3. Required: Proposition 1 tightening (R1-M3/A1/A2/A12, R2 §5, §13.3–4)

1. Replace the empirical "coefficient 3, Phase 48" step with the two-line analytic derivation: reflection makes each mode an improper complex Gaussian with |ρ| = σ², so E|Z|⁴ = 2σ⁴ + |ρ|² = 3σ⁴.
2. Derive the 0.43 constant in closed form, stating exactly what is standardized (mass vs norm; analytic vs empirical null variance). R1 reproduced scaling but the constant only to ~20% — purely a standardization ambiguity; fix by specification.
3. State explicitly that the β·conj(f) model is the carry-free idealization; the match to real digit sequences is empirical evidence (Table 1 cols. 7–8), not a corollary. Optional upgrade: the noisy-reflection lemma (R1-A10) giving mass (2+(1−q)²)/L with per-digit carry probability q — would convert the residual β scatter from "consistent with noise" to "predicted."
4. Handle edge modes: ξ=0, ξ=L/2, odd-L midpoint, μ_p vs S¹ alphabet (exact: 2/L² − 1/L³).

## 4. Required: scope section (R1 §7, R2 §10, §13.9)

Add new §11, "What this refutes, what it does not refute, and how future searches should be designed," per the eight numbered recommendations in Referee Report 1 §7 (placement, 1–1.5 pages, cross-references, one-row schema-killed table, merge "Honest scope" into it, resolve Phases 44–46 first, non-claims phrasing, name the protocol). R2 §10 independently demands the same section with the same killed/not-killed split — treat R1 §7 as the implementation spec. Add "extension-field cases" to the not-killed list (R2's addition; matches R1-A11 — the Hasse dichotomy argument is prime-field-specific and must say so).

## 5. Required: reproducibility, elevated to major (R1-M4/B7/C6, R2 §8, §13.7)

R2's framing is right and should be internalized: for a paper whose thesis is "controls can lie," the reproducibility bar is higher than normal. Ship all of:
pinned commit/release (ideally Zenodo DOI), exact Python/Sage versions, one-command reproduction for Table 1 and for Prop 1, CI running a deterministic subset, frozen pre-registration file with seed and timestamp, and a hard separation of submitted vs post-submission phases (44–46 either in or out — see item 4).

## 6. Required: ECDSA framing (R1-M5, R2 §9)

Move the scale caveat into the abstract/introduction, using R2's suggested sentence: experiments are diagnostic, not an empirical security evaluation of deployed curves; the extension to ECDSA-type curves is through the structural theorem, not experiment.

## 7. Strongly recommended experiments (both reviews converge)

Priority order:

1. **Asymmetric-section kill test** (R1-C1, R2 §7A). Break τ(−P)=−τ(P) deliberately; prediction: elevation collapses on the same curve/secret/digits. Both reviews call this the single strongest addition — it supplies the converse to the positive control and answers "of course reflection-built synthetics reproduce the reflection statistic."
2. **Equivalence testing + power analysis** (R1-B1/B2, R2 §6A, §13.5). TOST or minimum-detectable-effect bound for real-vs-syn ("we exclude any residual H3 component larger than X% of the elevation"); explicit power statement for the half-sequence control (√L scaling ⇒ a real signal would still appear at ≈ elevation/√2).
3. **Look-elsewhere quantification** (R1-B3, R2 §6C, §13.6). Expected max-z under a global null across the actual probe × encoding × digit × curve battery (analytic or by rerunning the battery on i.i.d. surrogates). Cite Gelman–Loken; optionally specification-curve (Simonsohn et al.).
4. **Non-CM curves** (R1-C2, R2 §7B). Add j=1728 and 2–3 random-j curves to the pre-registered table; Phase 14's Aut(E)-equivariance shows CM demonstrably touches digit statistics, so the all-j=0 Table 1 looks cherry-picked even though it isn't.
5. **Large-L synthetic convergence** (R1-C4, R2 §7D). Prop 1 constant at L = 10⁵–10⁶; near-zero cost.
6. **Smart verification: all k** (R1-C7, R2 §7C). At p=53, run all k ∈ [1, p−1]; label it implementation validation.
7. **Fixed-shift discriminator demo** (R1-C3). One synthetic run showing the non-reversing constraint leaves mass at 2/L.

## 8. Strongly recommended presentation (R1-D1/D2, R2 §11, §13.10)

Four artifacts, identical lists in both reviews:

1. Hypothesis × control matrix (H1/H2/H3 × four controls; cells kills/passes/predicts).
2. Null-distribution plot: rand, incShuf, syn nulls with the real statistic marked — "44σ from random yet inside the synthetic null" in one panel.
3. Growth-law plot: z(real/rand) and z(syn/rand) vs n with the predicted √n curve.
4. Attack-schema-killed table (schema | why dead | boundary | not covered | theorem/control).

Plus (R1 only, cheap): boxed practitioner checklist (D3), notation table (D8), name the methodology (D4).

## 9. Optional / third tier

Generalize Thm. 1 to abelian-group extensions (A5); cohomological remark (A6); log-coordinate version removing e≤4 (A7); G_m contrast (A8); negation-map subsumption sentence (A9); U³ corollary (A4); pseudo-variance classification theorem (A3); literature positioning — negative controls (Lipsitch et al.), exchangeability (Phipson–Smyth, Winkler et al.), cluster-failure (Eklund et al.), NIST SP 800-22 calibration norms (B4/B5/E3, R2 cites the same sources); digit-uniformity chi-square (B6); β-fit CIs and joint digit regression (B9); Phase-37 contamination quantification (C8); encoding-robustness row (C5); structured abstract (D6); cross-field precedents table (E4); "what standard practice would have said" run (F1).

## 10. Suggested execution order

1. Venue + title + AI-disclosure decision (item 1) — gates all prose edits.
2. Resolve Phases 44–46 in/out (gates items 4, 5).
3. Math edits: items 2, 3 (self-contained, no new experiments).
4. New experiments: item 7 (kill test first; it may also yield a figure).
5. Scope section + figures: items 4, 8 (write after experiments so the matrix/table reflect final content).
6. Reproducibility packaging: item 5 (last, against the frozen revision).
7. Abstract/intro rewrite: item 1's edits finalized once everything above is fixed.

**Bottom line:** no reviewer disputed any result. The revision is: narrow the claim, harden Prop 1 into a fully analytic statement, draw the kill boundary explicitly, add the asymmetric-section experiment plus equivalence/power/look-elsewhere calibration, ship airtight reproducibility, and add four figures. Done well, both reviews indicate acceptance at a negative-results/methodology venue.
