# Assessment: "When a Control Does Not Control: A Mechanism-Derived Refutation of a False Positive in p-adic Cryptanalysis of the ECDLP" (Trevena, June 2026)

## Summary of the submission

The paper reports, from a long negative exploration of p-adic lifting attacks on the ECDLP, a worked false positive: a Gowers U2 probe of the Teichmüller lift error δ(k) showed a 3–14σ elevation that a standard permutation control certified as real. The paper shows the permutation control could not separate the three candidate mechanisms, builds one mechanism-derived null per candidate, and localizes the entire signal to the exact public antisymmetry δ(k)+δ(n−k)=[n]τ(G) via a quantitative positive control. An inertness theorem (an acknowledged instance of Silverman and Gadiyar–Padma) explains why any such anomaly on ordinary curves must be artifact-or-identity, with the boundary at anomalous curves (Smart's attack). Claimed novelty is methodological only.

## Independent verification performed

All checks below were done from scratch (pure Python/numpy, written from the paper's text alone, different generator points, different seeds). Scripts: `verify.py`, `verify2.py`.

**Mathematics (checked by hand):** Lemma 1 and Theorem 1(a),(b),(d) are correct as stated and proved. Identity (3) follows in two lines from τ(−P)=−τ(P). The formal-log claim (first correction at degree 7 for j=0, hence z-additivity of kernel arithmetic at e≤4) is correct. The Remark on clause (c) correctly distinguishes the information-theoretic from the computational claim.

**Computations (all reproduce):**

- All seven Table 1 curve orders are exactly correct (n = 79, 139, 199, 313, 397, 613, 829; all prime; all p ≡ 1 mod 3; gcd(n,p)=1). The §3 example (p=97, n=79) and the anomalous curve (#E(F53)=53 for y²=x³+4x+7) are correct.
- The §9 Phase-41 order list matches exactly, in increasing-p order, for p = 31, 43, 67, 79, 97, 103, 109, 151, 157, 163 (n=86 is a non-maximal-order subgroup of N=172 at p=157 — see minor point M4). The two anomalous members of the family in that range (p=61, p=127, where N=p) were correctly excluded.
- Antisymmetry (3) holds exactly: z(δ(k))+z(δ(n−k)) takes a single value, equal to z([n]τ(G)), across all reflection pairs (p=67 checked exhaustively; all 7 curves checked in the Table 1 reproduction).
- Obstruction [n]τ(G): nonzero, valuation 1, additive order p^(e−1) — as claimed (Theorem 1(b),(d), Table 2 ordinary row).
- **Full Table 1 reproduces within Monte Carlo noise** in all five columns, on all seven curves, all three digits — including the growth of the elevation with n (real/inc from ~+3.4 at n=79 to ~+12 at n=829, syn/rand tracking it, real/syn and half/inc and inc/rand all at noise).
- **The Fourier decay exponents reproduce to two decimals**: β_real = {0.34, 0.45, 0.36} for digits 1,2,3 vs paper's {0.34, 0.45, 0.36}; β_syn = {0.41, 0.41, 0.40} vs paper's {0.41, 0.41, 0.40}.
- Smart's attack: 12/12 recovery reproduced on the anomalous curve.

This degree of reproducibility from the text alone is exceptional and should be weighed heavily in the paper's favor.

## Major point (must fix): §8 / Table 2 misstate what vanishes on the anomalous curve

The paper claims (abstract, §8, Table 2) that on the anomalous curve "the obstruction [n]τ(G) … vanishes." This is false as a value statement, and it contradicts the paper's own (correct) 12/12 recovery:

- Exact computation: z([53]τ(G)) = 15·53 ≢ 0 (mod 53²); valuation exactly 1 at any working precision. [p]τ(G) ≠ O at e=2.
- Smart's attack *requires* this to be nonzero: φ(G) = z([p]τ(G))/p must be invertible mod p to compute k = φ(Q)/φ(G). If [n]τ(G) were O, the attack would fail (the canonical-lift failure mode of Satoh–Araki).

What actually vanishes is the **section ambiguity**: for any kernel element c, [p]c = O at e=2 (verified), so [p]τ(G) is independent of the choice of section — confirmed exactly: three different sections (Teichmüller-x, naive-x, x+17p) all give z = 795 mod 53². The correct sharp dichotomy, cleaner than the one stated, is about the action of [n] on the kernel of reduction:

- gcd(n,p)=1: [n] is an automorphism of the kernel, so the section ambiguity survives at full strength; δ is not a function of the instance and the only well-defined object is the public constant (3).
- n=p: [n] annihilates the kernel, so the ambiguity dies; [p]τ(P) becomes a well-defined, generically **nonzero**, secret-linear invariant — Smart's homomorphism.

Likewise "the lift error becomes well-defined" (§8) is incorrect as stated: δ(k+p)−δ(k) = [p]τ(G) ≠ O. What becomes well-defined is [p]τ(P), not δ.

Required changes: rewrite the §8 narrative and the Table 2 anomalous-row entry/caption ("vanishes" → "section-independent and nonzero"; or report [n]c for kernel c, which does vanish), and adjust the abstract sentence "when it fails (n = p, anomalous curves) the obstruction vanishes." The sharp-boundary thesis survives intact — indeed strengthened — but the object that changes at the boundary must be stated correctly. Also note the conflation's origin: in the ordinary case [n]τ(G) = [n]c by Lemma 1, but for n=p no order-n lift G̃ need exist (and on this curve none does, since [p]τ(G) ≠ O), so the identification silently fails exactly at the boundary the section is about.

**Confirmed by the source repository.** The `ecdsa-lift` repo's own `experiments/phase42_anomalous_boundary.py` is correct precisely because the *code* keeps the two objects separate: `obstruction_vanishes()` computes [n]c = [p]c for the kernel element c (docstring: "the obstruction [n]c = [p]c with c in the kernel of reduction") and returns True, while `psi_log()` recovers k from `[p]·lift` = [p]τ(G), which is nonzero. The `notes/theorem.md` "boundary is sharp" section likewise says only "[n]c = [p]c vanishes." So the prose error is confined to the paper's §8/Table 2/abstract, where [n]c is rewritten as [n]τ(G); the experiment underneath is right. This also means the error was introduced in the write-up and was not caught by the repo's prior self-review (`feedback/ecdsa_lift_paper_referee_assessment.md`), which checked the ordinary-case identity δ(k+n)−δ(k)=[n]c=[n]τ(G) — true on ordinary curves — but never probed the anomalous-case labeling.

## Secondary points

S1. **Theorem 1(c) mixes a theorem with a state-of-knowledge claim.** "…and no efficient algorithm computing it from (E,G,Q) is known" is not a mathematical statement and should move from the theorem to the (already well-written) Remark. Clause (c) should retain only the CRT/independence content.

S2. **σ language vs 300-draw nulls.** z-scores up to ±13.6 are standardized deviations against an empirical null of 300 draws; Gaussian tail language ("5–11σ") overstates what the null can certify. The paper already disclaims discovery claims (§6); a sentence stating that z-values are effect sizes, not calibrated tail probabilities, would close the gap.

S3. **Residual β mismatch is real but honestly disclosed.** I confirm both sets of exponents; reflection-only synthetics reproduce the mean decay (~0.4) but not the digit-to-digit variation. The paper's disclaimer is accurate; if possible, one sentence on candidate origin (e.g., digit/carry nonlinearity) would help, but leaving it open is acceptable.

S4. **Artifact availability.** The `ecdsa-lift` repository exists and is well organized (numbered, reproducible per-phase scripts; fixed seeds; JSON outputs; src primitives), and I cross-checked it: its `results/phase41_theorem_verification.json` reports V1_Cn = 12151388 for p=67, identical to my independent reimplementation. For a paper whose entire contribution is methodological discipline, this must be made public (with DOI) and available to referees — I could not confirm a public URL. Two housekeeping items: the README still says the paper "builds to a 7-page PDF" (submission is 9 pages) and describes a tagged release that should be created; and `feedback/ecdsa_lift_paper_referee_assessment.md` is a prior self-review of the 7-page version — the submission already addressed most of its points (Table 1 now uses ordinary p≡1 mod 3 curves with n prime, the §5 positive control and §6 statistical protocol were added, the log coefficient was fixed to 3b/7, the oddness hypothesis is stated, and Theorem 1(c) was restated computationally with the Remark). The one substantive item that prior review did **not** raise, and that the revision did not fix, is the §8/Table 2 obstruction labeling (Major point above).

S5. **Scale.** All experiments are at p ≤ 823, e=4. §9 already concedes one large-prime-order ordinary curve "would strengthen the demonstration"; a referee will likely ask for it. The theorem's regime-agnosticism is correctly argued, so this is presentation, not soundness.

## Minor points

M1. Appendix A consists only of a heading; Table 3 floats after the references. Fix layout.
M2. Ref [1] is slides (gray literature); fine given [2], but consider citing [2] as primary throughout.
M3. "U3 and higher sat at baseline" (§4) is asserted without numbers; include them in the artifact.
M4. §9: state explicitly that n is the order of G and may be a proper divisor of #E(Fp) (n=86 at p=157, where #E=172); Table 1's "n = #E(Fp) prime" applies only to the clean set.
M5. Two curves with n=79 (p=67 in Table 1, p=97 in §3) may momentarily confuse; a parenthetical would help.

## Overall evaluation

Soundness: high — every checkable claim checked out except the §8/Table 2 mischaracterization, which is a localized (and fixable) conceptual error that does not threaten the main results; the experiment it describes (12/12 recovery) is correct.
Reproducibility: exceptional (full independent reproduction from the text, including two-decimal agreement on fitted exponents).
Novelty: methodological only, as honestly claimed; the worked example is genuinely instructive and timely for automated/AI-assisted exploration pipelines.
Fit: depends on venue — this is a methodology/negative-results paper; at a venue that accepts such (e.g., a communications-style crypto journal or a reproducibility-oriented track), it merits acceptance after revision. At a venue requiring new cryptanalytic results or new mathematics, the fit question dominates, and the paper itself concedes the theorem is not new.

**Recommendation: minor-to-moderate revision.** Required: the §8/Table 2 correction (Major point), S1, S4. Suggested: S2, S3, S5, M1–M5.
