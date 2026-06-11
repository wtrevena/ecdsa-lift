# Referee Report

**Manuscript:** *When a Control Does Not Control: A Mechanism-Derived Refutation of a False Positive in p-adic Cryptanalysis of the ECDLP*
**Author:** W. T. Trevena
**Date of review:** 10 June 2026

---

## Note to the editor (please read first)

The corresponding-author address on the manuscript (`william.todd.trevena@gmail.com`) is a close match to my own name and contact details. If I am in fact connected to the author, this is a conflict of interest and I should be recused from a formal decision-making role. I have written the report below as a substantive technical assessment regardless, so that it is useful either as a referee report or as critical feedback to the author, but the editor should confirm reviewer independence before relying on it for an accept/reject decision.

---

## 1. Summary of the submission

The paper is a methodological case study built around a single negative cryptanalytic result. In an automated, multi-phase exploration of p-adic / Teichmüller lifting attacks on the ECDLP for ordinary curves (the ECDSA regime), a higher-order Fourier probe (a Gowers U² statistic on the digits of the Teichmüller lift error δ(k)) returned a large apparent anomaly (5–11σ, growing to an effect size ≈44 on a large-prime curve). A standard value-permutation control endorsed it as genuine index-ordered structure.

The paper's thesis is that this endorsement was invalid because the permutation control cannot separate the competing explanations for the signal. The author enumerates three candidate mechanisms (H1 prefix-sum/low-pass artifact, H2 antisymmetry shadow, H3 genuine new structure), builds one mechanism-derived null per hypothesis, and — decisively — a *positive* control: synthetic sequences carrying only the exact public antisymmetry δ(k)+δ(n−k)=[n]τ(G) reproduce the statistic and its growth with n. A closed-form result (Proposition 1) shows the antisymmetry inflates the U² mass by exactly 3/2 and forces a √n growth of the standardized gap, i.e. the "growth = accumulating signal" reading is inverted. An inertness theorem (Theorem 1) then argues that for gcd(n,p)=1 the only secret-dependent part of any section lift error, [k]c, is governed by k mod p^{e−1}, which the public data does not pin down over the natural lift; the boundary case n=p recovers Smart's attack, and §8 argues the dichotomy is exhaustive (no intermediate "graded leak").

The author is explicit that the cryptanalytic verdict is the expected negative one and that Theorem 1 is not new mathematics — it is an instance of Silverman's and Gadiyar–Padma's principle, recast in lift-error language. Novelty is claimed only for the methodology and for the explicit identification of the obstruction with the antisymmetry constant.

## 2. Overall assessment

This is a careful, honest, and well-written paper. The central methodological point — that a control is a null model testing exactly one named alternative, and that a permutation test which collapses several mechanisms at once cannot attribute the signal to any of them — is correct, clearly argued, and well illustrated. The positive control (synthesizing the surviving mechanism and showing it quantitatively reproduces the statistic) is the strongest part of the paper and is genuinely good practice that deserves wider adoption, especially for automated/AI-assisted "fan-out many probes" pipelines.

My reservations are not about correctness, which I checked and found sound (Section 3), but about the **weight and novelty** of the contribution relative to a journal bar (Section 5). The technical core is, by the author's own framing, a repackaging of known results, and the statistical methodology, while correctly applied and nicely localized, is close to standard practice in applied statistics (permutation nulls test a specific exchangeability hypothesis; mechanism-specific nulls and positive controls are established tools). The paper's value is in the *synthesis and the worked example*, which is real but should be judged against the target venue.

## 3. Correctness — independent verification

I independently reproduced or checked the key quantitative and theoretical claims. All held up.

**Proposition 1 (the 3/2 mass ratio and √n law) — reproduced numerically.** Simulating S¹-uniform i.i.d. sequences versus reflection-constrained sequences (free half i.i.d., second half fixed by f(L+1−k)=β·conj f(k)) and computing the U² mass Σ|f̂(ξ)|⁴, I obtained, averaged over many draws:

| L | E[‖f‖⁴]·L, i.i.d. | E[‖f‖⁴]·L, constrained | ratio |
|---|---|---|---|
| 78 | 1.99 | 2.90 | 1.46 |
| 312 | 2.00 | 2.98 | 1.49 |
| 828 | 2.00 | 2.99 | 1.50 |
| 10638 | 2.00 | 3.00 | 1.50 |

The 2/L → 3/L jump and the 3/2 ratio are confirmed and converge as L grows, and the standardized synthetic-vs-random gap scales as Θ(√L). The precise constant (the paper's 0.43√n) depends on whether one standardizes the U² *norm* or the *mass* and on the null's variance estimate; I reproduced the scaling and the constant to within ~20%, which is consistent with that standardization choice. The qualitative claim that matters — *growth-with-n is the signature of the artifact, not of accumulating signal* — is solidly supported.

**The antisymmetry identity (3) — confirmed analytically.** δ(k)+δ(n−k)=[n]τ(G) is an immediate algebraic consequence of the section property τ(−P)=−τ(P): δ(k)+δ(n−k) = [n]τ(G) − (τ(kG)+τ((n−k)G)), and since (n−k)G=−kG, the bracket is τ(kG)+τ(−kG)=O. The identity is therefore elementary and not in doubt. That it holds *in the z-coordinate as a single value* relies on the kernel being additive in z at precision e≤4; I verified symbolically that the formal group law for y²=x³+b has its lowest cross-term at total degree 7 (the −3b z₁⁶z₂ … terms), so for kernel elements of valuation ≥1 the correction has valuation ≥7 and vanishes mod p⁴. The precision argument is correct.

**Theorem 1 and Lemma 1 — logic checked.** Parts (a),(b),(d) are direct algebra and are correct. Lemma 1 (unique order-n lift, δ_can≡0) is the standard fact that reduction is an isomorphism on n-torsion when gcd(n,p^{e−1})=1, argued correctly. Clause (c) is the delicate one and is discussed below as a major point.

**The boundary dichotomy (§8) — confirmed numerically.** I checked that for E/F_p the only multiple of p in the Hasse interval [(√p−1)², (√p+1)²] is N=p itself for all p≥7, with the N=2p exception occurring only at p=5, exactly as the paper states. Hence v_p(n)∈{0,1} and there is no intermediate regime; the "graded inertness theorem has empty interior" claim is correct.

What I did **not** do: re-execute the full p-adic lift pipeline that produces the δ(k) sequences and Table 1's seven-curve numbers. Given that identity (3) is elementary and Proposition 1 explains the entire elevation in closed form, I did not consider this necessary for the report, but if the editor wants end-to-end reproduction the author's tagged repository should be run independently (see major point M4).

## 4. Strengths

- The positive control is the right tool and is executed convincingly: a synthetic model of the surviving mechanism reproduces the statistic *and its growth law*, which is much stronger evidence than the negative controls alone.
- Proposition 1 turns an empirical observation into a closed form and correctly reframes the √n growth as an artifact signature. This is the paper's most quotable technical contribution and gives auditors a concrete discriminator (mass 3/L for index-reversal coupling vs 2/L for a shift).
- Intellectual honesty throughout: the author repeatedly distinguishes what is proved from what is a state-of-knowledge remark (notably the Remark after Theorem 1), declines to claim a reduction beyond what holds, and foregrounds that the theorem is not new mathematics.
- The statistical protocol (§6) is pre-registered, seeded, and reports effect sizes against empirical nulls with explicit disavowal of Gaussian tail interpretation and of multiplicity-based discovery claims. This is exemplary for a negative-result paper.

## 5. Major points

**M1 — Novelty and venue fit.** The paper should state more sharply, up front, what is new. By the author's own account the cryptanalysis is negative and expected, Theorem 1 is an instance of [2,1,3], and the methodological principles (a permutation null tests one exchangeability hypothesis; build mechanism-specific nulls; use a positive control) are individually standard. The genuine contributions appear to be (i) the *worked localization* of a specific, seductive anomaly to a single public identity, and (ii) the closed-form 3/2 / √n artifact signature. I recommend the introduction be rewritten to claim exactly these and no more, and that the editor weigh whether (i)+(ii) clear the bar for a research article at this venue versus a methodological note/SoK-style contribution. I think the material is publishable, but the framing currently oversells via the abstract's "5–11σ" headline before deflating it.

**M2 — The status of Theorem 1(c) as a "theorem."** The security-relevant content of the inertness result is computational (no efficient algorithm produces [k]c from (E,G,Q)), and the author correctly confines this to a Remark rather than the theorem. But this means the *theorem proper* establishes essentially a well-definedness/CRT-independence statement over an artificially extended index set (k uniform on Z/np^{e−1}Z), while the part that would actually certify security is a conjectural state-of-knowledge assertion no stronger than "this lift does not help." This is fine and honestly presented, but the paper occasionally leans on "Theorem 1 forecloses its being cryptanalytic" (e.g. §5) in a way that slightly overstates what the theorem alone gives. I'd ask the author to (a) make explicit, at each invocation, whether they are using the proved well-definedness clause or the conjectural computational remark, and (b) soften "forecloses" to reflect that the foreclosure is modulo the standard ECDLP-hardness state of knowledge.

**M3 — Proposition 1's model vs the digit-encoded sequence.** Prop 1 analyzes the idealized constraint f(L+1−k)=β·conj f(k). The actual statistic is computed on the phase encoding f_j(k)=exp(2πi·digit_j(z(δ(k)))/p), and digit extraction does not commute with the additive reflection s(n−k)=C−s(k) (base-p carries). The paper acknowledges the finer digit-dependent Fourier decay is only partially reproduced and attributes the residual to finite-size scatter — which is reasonable and supported by the resampling argument — but Prop 1 should state explicitly that its "β·conj f" model is the carry-free idealization and that the match to the real digit sequence is empirical (Table 1, cols. 7–8), not a corollary of the proposition. As written, a reader could infer Prop 1 proves the match for the real sequences, which it does not.

**M4 — Reproducibility.** The paper rests heavily on numerics in Table 1 and on the 47-phase exploration indexed in Appendix A. The repository is referenced but the manuscript does not pin a commit/release hash, does not state software versions (the "pure-Python lift" and the Sage validation in Phase 30–36), and Phases 44–46 are described as "post-submission, see repository," which is awkward for a static reviewed artifact. I recommend: a pinned release tag/commit in the paper; a one-command reproduction script for Table 1 and for Proposition 1's constant; and removal or clear bracketing of any phase whose results are not fixed at submission time.

**M5 — Scope/over-generalization to ECDSA-size curves.** The numerical verification of Theorem 1 is on small curves (p=31–163) plus one curve at p=10477, and the author is careful to say the relevance to real ECDSA curves is "through the theorem's generality, not through these particular orders." That caveat is appropriate, but because the abstract and title invoke "ECDSA," the paper should add one explicit sentence in the introduction clarifying that nothing is claimed or tested at cryptographic sizes and that the theorem (not the experiments) is what carries the regime-agnostic conclusion.

## 6. Minor points

- Abstract: "an apparent 5–11σ effect is cheap" and "standardized effect size ≈44 … not a Gaussian tail" — the rhetorical setup is effective but the abstract is dense and front-loads sensational numbers before the deflation. Consider trimming.
- §4: "z(real vs rand) ranges 3.3 to 13.6" vs the abstract's "∼3–14σ" — make the two consistent (round identically).
- §5, Proposition 1 proof sketch: state explicitly that the ξ=0 mode contributes at O(1/L²) and is asymptotically negligible, so that the "L−1 nonzero modes" accounting is rigorous.
- §6: "We make no multiplicity correction and therefore do not claim discovery" is good; consider moving this sentence earlier (it preempts the obvious objection to the σ language).
- §8: the parenthetical that δ itself stays ill-defined while [p]τ(P) becomes well-defined is important and easy to miss; consider promoting it from a parenthetical.
- References: [1] is an invited-talk slide deck "available via the author's page" — give a stable URL or archive link. [7] Lubin–Serre–Tate is cited as lecture notes; give a canonical reference. Check that [3] (Gadiyar–Padma, arXiv:1702.07107) page range and the INDOCRYPT 2004 companion are both cited consistently.
- Table 1: define "real/inc," "inc/rand," etc. in the caption on first use without requiring the reader to reconstruct them from §5; the column semantics (which null is the denominator of each z) should be unambiguous.
- Typography: ensure σ/z usage is uniform (the paper sometimes writes "σ" and sometimes "z" for the same standardized quantity).
- "Honest scope" paragraph (§10) is excellent and could be referenced from the introduction so the novelty boundary is set early.

## 7. Recommended next steps (concrete)

The author has received feedback elsewhere suggesting a dedicated scope-and-research-map section, and I endorse a refined version of that suggestion as the single highest-value revision beyond M1–M5. The paper currently disposes of the false positive and proves inertness, but it leaves the *boundary of the refutation* implicit, scattered across §7's Remark, §10's "Honest scope" paragraph, and the phase index. That gap invites exactly the two misreadings most likely to sink the paper in review: an overclaim reading ("the author asserts all p-adic ECDLP research is dead" — which Theorem 1 does not and cannot establish, see M2) and an underclaim reading ("the author merely found and fixed their own false positive — why is that publishable?"). A section that draws the boundary precisely defuses both at once, because the correct message sits exactly between them: *this specific attack schema is dead for structural reasons, and the methodology that killed it is the reusable contribution.* Crucially, the section should be framed as a **scope boundary** — a statement of what Theorem 1 and the controls jointly entail — and not as a loose "future work" list, which would dilute the paper's disciplined tone.

I recommend the new section be titled along the lines of **"What this refutes, what it does not refute, and how future searches should be designed,"** and be built from three tightly-scoped paragraphs:

**(a) The killed schema, stated exactly.** One paragraph defining the attack class the paper terminates: extraction of the secret k from *statistics of section lift errors* δτ(k) — or any digit, phase, or other deterministic encoding thereof — on ordinary curves with gcd(n, p) = 1. This covers apparent pseudorandomness failures, Fourier/Gowers elevation, spectral-decay anomalies, digit nonuniformity, prefix-sum/walk artifacts, and machine-discovered structure in the ordered sequence (k, δτ(k)). The paragraph should state the resulting prior explicitly: by Theorem 1(c)–(d) and Corollary 1, any such signal computable from (E, G, Q) carries a *structural presumption of identity-or-artifact*, rebuttable only by surviving mechanism-derived controls of the kind in §5. This is the paragraph where the 47-phase exploration earns its keep: Appendix A is the empirical census of this schema, and the theorem is the reason the census had to come back empty.

**(b) What is not killed.** One paragraph listing, affirmatively, what lies outside the refutation: nonce-leakage attacks on ECDSA, side channels, invalid-curve and twist attacks, small-subgroup attacks, small-embedding-degree/MOV mechanisms (inert here only for the independent reason already noted in §10), quantum algorithms, deliberately weakened curve families, and — most importantly for the mathematical audience — *non-section-based global lift constructions*. The paper proves inertness for section lift errors; it claims no impossibility theorem for all p-adic or global lifts, and Silverman's "four faces" taxonomy [2, 1] should be invoked here to place the killed schema as one cell of that map, not the map. This paragraph is the overclaim antidote.

**(c) A protocol for future searches.** One paragraph (or a short numbered checklist, which I would prefer — it makes the methodology citable as such) prescribing what any new direction in this area must do before its anomalies deserve attention: (i) pre-specify the public invariant computed from (E, G, Q); (ii) show it is well-defined independent of section choice — the Theorem 1(d) test; (iii) show it is not forced to be canonically zero (Lemma 1) or a public identity (the fate of (3)); (iv) pair every proposed signal with mechanism-derived nulls covering walk/summation artifacts, algebraic symmetries (reflection foremost), and automorphism/coordinate effects, *plus* a positive control synthesizing the surviving mechanism, per §5; (v) accept as evidence only signals that survive all of the above, scale on large ordinary prime-order curves, and admit a mathematical explanation not reducible to the known inertness. This paragraph is the underclaim antidote: it converts the worked example into a reusable instrument.

Concretely, I recommend the author:

1. **Place the section between the current §10 and §11**, as a new §11 ("What this refutes, what it does not refute, and how future searches should be designed"), with the Conclusion renumbered to §12. Placing it after the Discussion lets it inherit the methodological lessons of §10(i)–(iv) by reference; placing it before the Conclusion lets the Conclusion's final sentence point at the protocol rather than restate it.

2. **Keep it to roughly one to one-and-a-half pages** for a venue of this mathematical register: three paragraphs of 6–10 lines each plus one small table (item 4). Longer and it becomes the loose future-work section the framing must avoid; shorter and the not-killed list degenerates into an unargued bullet dump.

3. **Anchor every claim in paragraph (a) to an existing result by explicit cross-reference**: the identity-or-artifact prior to Theorem 1(c)–(d) and Corollary 1; the symptom list to the phase encoding f_j of §4 and the probes of Appendix A; the artifact signatures to Proposition 1's 3/2-mass / √n discriminator (which deserves a sentence here as the auditor-facing tool it is). Per M2, each invocation must say whether it leans on the proved well-definedness clause or the state-of-knowledge Remark — paragraph (a)'s "strong prior" language is the honest formulation precisely because the computational claim is conjectural.

4. **Add a one-row "attack schema killed" summary table** in the style of Table 2 — columns roughly: *schema* (statistics of δτ(k) and encodings, ordinary gcd(n, p) = 1), *why dead* (Thm. 1: secret part [k]c not publicly determined; visible structure = identity (3) or artifact), *boundary* (n = p: Smart, §8), *not covered* (pointer to paragraph (b)). This gives skimming reviewers and future citers the paper's verdict in one glance and is cheap to typeset.

5. **Merge, do not duplicate, the existing "Honest scope" paragraph (§10).** Its content (theorem not new mathematics; section lift errors only; pairings separate) is a strict subset of the new section. Replace it with a one-sentence forward pointer ("the precise boundary of the refutation is drawn in §11"), and likewise add the single early cross-reference from the introduction that I asked for in §6 (minor points) of this report — that resolves M1's framing concern and the new section simultaneously.

6. **Resolve the status of Phases 44–46 before this section is written.** They are marked "post-submission, see repository" in Appendix A (see M4). Either bring their results into the fixed manuscript (in which case the aggregation/MI estimators and alternative-uniformizer probes belong in paragraph (a)'s symptom census) or bracket them out of the paper entirely and have the new section's protocol note that the repository hosts continuing post-publication probes under the same pre-registered discipline. A scope section that cites unfixed results would undercut its own purpose.

7. **State the not-killed list as the paper's non-claims, not as open problems.** The phrasing should be of the form "this paper establishes nothing about X" rather than "X remains interesting future work" — the former is a boundary, the latter is the loose framing to avoid, and the former also pre-empts any reviewer attempt to read an impossibility claim into Corollary 1.

8. **Consider naming the checklist of paragraph (c)** (even something plain like "the section-invariance protocol") and echoing that name once in the abstract's final sentence. If the methodology is the contribution — as M1 concludes and the author's own "Honest scope" agrees — it should be packaged so that subsequent papers can comply with it, cite it, or be measured against it. That, more than any individual numerical result here, is what would make this paper load-bearing for the field.

I want to be explicit about why this addition materially strengthens the submission rather than padding it. Under the overclaim reading, a reviewer skims the title and Corollary 1 and objects that no theorem of this shape can close off all lifting approaches; paragraph (b) and recommendation 7 make that objection unraisable because the paper will have said it first, with the taxonomy citation in hand. Under the underclaim reading, a reviewer grants correctness but sees only an internal post-mortem; paragraph (c) and recommendations 4 and 8 answer that the deliverable is a transferable falsification protocol with a closed-form artifact signature, of which the Gowers episode is the worked demonstration. The two readings cannot both be pre-empted by the current diffuse scoping; a single section that kills the schema sharply *and* hands the field the instrument that killed it does both, and aligns the paper's structure with what §2 of this report already identified as its real contribution.

## 8. Recommendation

**Minor-to-major revision** (technically sound; revisions are about framing, scope, and reproducibility, not correctness).

The mathematics and statistics are correct as far as I checked, and the paper is unusually honest about what it does and does not establish. My substantive asks are: sharpen the novelty claim and align it with the "Honest scope" paragraph (M1); separate the proved clause of Theorem 1 from the conjectural computational remark at every point of use (M2); clarify that Proposition 1 models a carry-free idealization (M3); pin reproducibility artifacts (M4); and add an explicit ECDSA-scale caveat (M5). With those addressed I would support publication, with the venue-fit judgment (research article vs methodological note) left to the editor.

I would be glad to see a revised version.

---

# Appendix: Exhaustive strengthening analysis

The following is a comprehensive, prioritized inventory of ways the paper could be improved or strengthened, going beyond the major/minor points above (checked for non-duplication with M1–M5 and the minor list). Items marked **[H/M/L]** importance, **(S/M/L)** effort.

## A. Mathematical strengthening

**A1. Replace the empirical step inside Proposition 1's proof with the actual two-line computation. [High] (Small)**
The proof sketch currently says the coefficient 3 was "verified by direct second-moment measurement, Phase 48" — an empirical step inside a proposition's proof is a genuine weakness a referee will pounce on. The computation is short: the reflection f̂(L−ξ) = β·conj f̂(ξ) makes each Fourier mode an *improper* (non-circularly-symmetric) complex Gaussian whose pseudo-variance ρ = E[f̂(ξ)²] attains |ρ| = E|f̂(ξ)|² = σ². For such a variable E|Z|⁴ = 2σ⁴ + |ρ|², hence 3σ⁴ = 3/L² exactly. Stating it this way makes the "3" a one-line consequence of impropriety, not a measurement.

**A2. Derive the 0.43 constant in closed form (and state the standardization). [High] (Small–Medium)**
The √n coefficient is currently empirical, and this report could only reproduce it to ~20% because the standardization (U² norm vs U² mass; analytic vs 300-draw sd) is ambiguous. All ingredients are Gaussian: Var of the i.i.d. mass is 20/L³ by Wick (E|Z|⁸ = 24σ⁸ for circular Z), the mean gap is 1/L, and a delta-method step converts to the norm. Publish the constant as a formula, then report the measured value as confirmation. This also fixes the awkwardness that a "closed form" proposition ends with a measured constant.

**A3. Unify the discriminator into a small classification theorem. [High] (Medium)**
The remark "index-reversal coupling → 3/L; fixed shift → 2/L" generalizes cleanly via A1: any linear constraint induces a pseudo-variance ρ(ξ) per mode, and the mass inflation factor is exactly 1 + ⟨|ρ|²⟩/(2σ⁴) ∈ [1, 3/2]. This gives auditors a *continuum* lookup: which symmetries inflate U², by how much, with reflection as the extremal case. One proposition replacing the parenthetical remark would be the paper's most reusable mathematical artifact.

**A4. Extend the analysis to U³ to explain (not just report) "U³ sat at baseline." [High] (Medium)**
The paper asserts empirically that U³ and higher showed nothing. That is a prediction Prop-1-style machinery can derive: compute the expected U³ mass under the reflection constraint and show the inflation is asymptotically negligible (or a smaller factor). Turning the unexplained observation into a corollary closes a loose end and preempts "why didn't the artifact show in U³?"

**A5. State Theorem 1 at its natural generality: any abelian-group extension. [High] (Small)**
The proof uses only: an extension 0 → K → E′ → E → 0 of abelian groups, K a p-group of exponent dividing p^{e−1}, and gcd(n, p) = 1. State the abstract version in one paragraph and derive the elliptic case as an instance. This makes "inertness" portable to Jacobians, algebraic tori, and any future lifting setting — and directly answers "is this just about one curve coordinate?"

**A6. Add the cohomological one-liner. [Medium] (Medium)**
Lemma 1 and Theorem 1(d) are the statements that H¹ and H² of Z/n with coefficients in the p-group Ê vanish (n invertible on Ê): sections = splittings, the unique order-n lift is the vanishing of the obstruction class, and [n]τ(G) is precisely the norm of the cochain c. A short remark casting δ as a 1-cochain whose failure to descend to Z/n is measured by [n]c would please number theorists and make "not new mathematics" precise rather than apologetic.

**A7. Work in log coordinates to remove the e ≤ 4 restriction. [Medium] (Small–Medium)**
The z-coordinate single-valuedness of identity (3) leans on log(z) ≡ z mod p^e for e ≤ 4. Restating all digit statistics on log(z(δ(k))) makes every identity exact for *all* e and all Weierstrass models, eliminating a hypothesis. At minimum, add a remark that the e ≤ 4 restriction is cosmetic and the log-coordinate version is unconditional.

**A8. Add the multiplicative-group contrast as a remark. [Medium] (Small)**
For F_p^× the Teichmüller lift *is* a homomorphism, so δ ≡ 0 and there is no anomaly to chase; the elliptic case is interesting precisely because the canonical section is not multiplicative pointwise on E(F_p) under non-canonical sections. One paragraph contrasting G_m (Gadiyar–Padma's setting) with E explains *why* this false positive could exist at all, and tightens the link to [3].

**A9. Quantify and subsume the "2-to-1 search reduction." [Medium] (Small)**
The reflection k ↔ n−k is exactly the Q ↔ −Q identification already exploited by the negation-map speedup in Pollard rho (a standard √2 gain). One sentence noting the antisymmetry yields *nothing beyond the textbook negation map* fully closes the "but it IS an advantage" door.

**A10. A "noisy reflection" lemma to bridge Prop 1 to the real digit sequences (beyond M3's caveat). [Medium-High] (Medium)**
M3 asks only for a caveat that carries break the idealization. Go further: if the phase-level reflection holds for a fraction 1−q of pairs (carries/borrows at digit j occurring with computable probability q), the mode pseudo-variance scales by (1−q) and the mass becomes (2 + (1−q)²)/L by A3's formula. Estimate q per digit from the data and check the predicted partial inflation — this could *quantitatively* explain the residual digit-to-digit β variation now attributed to scatter.

**A11. Scope the exhaustive dichotomy to prime fields explicitly, or extend it. [Medium] (Small)**
The Hasse-interval argument ("only multiple of p in range is N = p") is special to E/F_p. Over F_{p^m}, p | N with v_p(N) > 1 is possible and the anomalous-type analysis differs. Either add "over prime fields" to every statement of the dichotomy or include a remark on what changes for extension fields. Currently §8's "the dichotomy is exhaustive" reads more general than it is.

**A12. Edge cases in Prop 1. [Low] (Small)**
Handle the self-paired midpoint when L is odd (f(k)² = β forces a deterministic point), the ξ with L−ξ = ξ mode, and state in what sense "independent of the value alphabet" holds (the o(1/L) term is alphabet-dependent; for the actual μ_p-valued phases E|f̂|⁴ = 2/L² − 1/L³ exactly).

## B. Statistical strengthening

**B1. Formal equivalence testing for "real vs syn indistinguishable." [High] (Small–Medium)**
The decisive claim is currently absence-of-significance (z ∈ [−1.4, +1.5]). Apply TOST or report the minimum detectable difference: pre-specify an equivalence margin (e.g., 10% of the elevation) and show the data exclude any residual H3 component larger than that margin. This converts the paper's strongest empirical claim from "we saw nothing" to "we bound what could be there."

**B2. Power analysis of the half-sequence control. [High] (Small)**
A skeptic's best line: "the half-sequence collapse could be power loss — half the length, half the sensitivity." Since z scales as √L, an H3 signal of the observed size would have appeared at ≈ elevation/√2 on the half; report the detectable effect size at L/2 with 300 draws and note the observed half/inc values exclude it. Also state explicitly that the half-null is recomputed at the half-length.

**B3. Quantify the look-elsewhere effect of the 47-phase exploration. [High] (Medium)**
The "false-positive engine" claim is qualitative. Estimate the expected maximum standardized statistic under a global null across the actual number of probes × encodings × digits × curves run (a max-of-correlated-z calculation, or empirical: rerun the full probe battery on i.i.d. surrogates and report the max z observed). Showing "a 5σ headline is *expected* from this many looks" makes the methodological thesis quantitative, and connects to Gelman–Loken forking paths and specification-curve analysis (cite Simonsohn–Simmons–Nelson).

**B4. Position within the negative-/positive-control literature. [High] (Small)**
The methodology has named precedents the paper should cite: negative controls in epidemiology (Lipsitch–Tchetgen Tchetgen–Cohen 2010), placebo tests in econometrics, exchangeability requirements for permutation tests (Phipson–Smyth; Winkler et al.), null-model controversies in network science/ecology (degree-preserving rewiring choosing the wrong null), and GWAS/neuroimaging permutation pitfalls under dependence (Eklund et al.'s cluster-failure as the famous "control that did not control"). This both strengthens the scholarship and *answers* the "is the lesson new?" attack: the lesson is established elsewhere and demonstrably missing from cryptanalytic practice.

**B5. One sentence of permutation-test theory naming the failure. [Medium] (Small)**
The value-permutation null assumes exchangeability of values under H0; any index-coupled constraint (the reflection) violates exchangeability, so rejection is uninformative about *which* dependence exists. Saying this in the standard statistical vocabulary, with citation, shows the failure was an instance of known theory — which is precisely the paper's point that practice lags theory.

**B6. Justify (or test) the i.i.d.-uniform baseline at the digit level. [Medium] (Small–Medium)**
The rand null is uniform on the S¹ continuum; the real f_j takes values in μ_p, and digits of z(δ(k)) are only *pseudorandomly* uniform (Phases 1–2). Add a chi-square uniformity test of the digit multiset per curve, and a remark that continuum-vs-μ_p alphabet changes E|f̂|⁴ only at O(1/L³) (per A12), so the baseline choice is immaterial at the reported precision.

**B7. Externalize the pre-registration. [Medium] (Small)**
"Fixed before the clean re-run; seed 20260609" is self-attested. Commit the protocol file with a dated, hashed tag (or OSF/Zenodo timestamp) and cite that artifact. Distinct from M4's pinning: this is about *timestamping the design*, not the code.

**B8. Report exceedances and show the nulls. [Medium] (Small)**
z = 44 against a 300-draw null invites the objection that no 300-sample null can calibrate such a value. Report the raw exceedance (0/300) alongside z, and (with D1) show the null histograms with the real value marked — the picture "real sits *inside* the synthetic null but 44 sd's from the random null" is the whole paper in one panel.

**B9. Honest uncertainty on the β fits. [Medium] (Small)**
The decay exponents come from 7-point log-log fits; report per-fit standard errors/CIs, and replace the informal "within ≈ 1 standard deviation" with a joint regression including digit as a factor and a formal test for digit dependence (which should fail to reject). Seven points is thin; saying so explicitly is stronger than letting the reader notice.

## C. Experimental strengthening

**C1. The asymmetric-section kill test — a new falsifiable prediction. [High] (Medium)**
Identity (3) requires τ(−P) = −τ(P). Construct a deliberately asymmetric section (e.g., choose the Hensel y-root by a parity convention that breaks the negation symmetry); the reflection identity then fails, and the H2 mechanism predicts the U² elevation *collapses to baseline* on the same curve, same secret, same digits. Running this is the cleanest possible confirmation that the antisymmetry — and nothing else about δ — generates the anomaly. Arguably the strongest single experiment the paper could add.

**C2. Run the clean protocol on non-CM curves. [Medium-High] (Medium)**
Table 1 is exclusively j = 0 (CM by Z[ω], extra automorphisms — and Phase 14 found an Aut(E)-equivariance of the leading digit, so CM demonstrably *does* touch the digit statistics). The original sweep's "CM, j=1728, and non-CM alike" remark covers only the decay exponent. Add j = 1728 and 2–3 random-j curves to the pre-registered table to show 3/2-inflation and real-vs-syn equivalence are CM-independent.

**C3. Demonstrate the fixed-shift discriminator empirically. [Medium] (Small)**
The claim "a constraint sending mode ξ elsewhere leaves the mass at 2/L" is stated but never shown. One synthetic experiment (e.g., s(k + m) = C − s(k) for fixed m ≠ reflection) with a one-line result would let auditors trust the discriminator as tested, not asserted.

**C4. Push the Prop 1 constant to L = 10⁵–10⁶. [Medium] (Small)**
Phase 48 is curve-free (pure synthetic), so verifying the constant at L two orders larger costs almost nothing and forestalls "constant drift" worries. Separately, if the pure-Python lift permits, one real curve at n ≈ 10⁵–10⁶ would extend Table 1's confirmatory row.

**C5. Encoding robustness row. [Low-Medium] (Medium)**
Show the localization survives at least one alternative encoding (e.g., phase-encode z mod p² directly, or balanced digits): the elevation should transform exactly as the reflection predicts. Phases 44–46 hint at this but are bracketed as post-submission; promoting one stable result into the paper would close the "you only tested one encoding" objection.

**C6. CI + archival artifacts (beyond M4's pinning). [Medium] (Small)**
Add a CI workflow running a fast subset reproducing Table 1 and the Prop 1 constant on every push; archive the release on Zenodo for a DOI; target the venue's artifact-evaluation badge.

**C7. Strengthen the Smart-attack verification. [Low] (Small)**
"12/12 random secrets" undersells a deterministic attack: either state it is deterministic and the 12 trials test the implementation, or run all k ∈ [1, p−1] (trivial at p = 53) and say "p−1 of p−1."

**C8. Quantify the Phase-37 contamination. [Low-Medium] (Small)**
The original sweep "later found to include supersingular and anomalous cases" — report what the 5–11σ range becomes after removing the invalid curves, so the reader knows the false positive wasn't *itself* partly an artifact of bad curve selection.

## D. Presentation and structure

**D1. Add figures — the paper has none. [High] (Medium)**
Four concrete panels: (a) log-log z(syn/rand) and z(real/rand) vs n with the closed-form √n line through both — the paper's thesis in one plot; (b) overlaid Fourier spectra |f̂(ξ)|: real vs reflection-only synthetic vs i.i.d., showing the mode-coupling signature; (c) null-distribution histograms (rand, incShuf, syn) with the real value marked — the "z = 44 yet inside the syn null" picture (pairs with B8); (d) the hypothesis-control schematic (D2). For a methodology paper, (c) and (d) carry the argument better than any table of z-ranges.

**D2. A hypothesis × control matrix. [High] (Small)**
A 3×4 table — rows H1/H2/H3, columns {value-permutation, increment-shuffle, half-sequence, synthetic positive} — with cells "kills / passes / predicts elevation." It shows at a glance that the permutation column kills all three rows (the invalidity) while the other columns separate them (the fix). This single table is the paper's argument and would be the most-reproduced figure in talks.

**D3. A boxed practitioner checklist. [High] (Small)**
The Discussion's four lessons, recast as a numbered audit procedure: (1) enumerate mechanisms; (2) one mechanism-derived null each; (3) positive control on the survivor; (4) check the growth law against the artifact signature (A3's table); (5) seek a structural inertness reason. Boxed and titled, it makes the paper citable as a protocol, not just a story.

**D4. Name the methodology. [Medium] (Small)**
"Mechanism-derived nulls" appears but is not branded. A consistent, memorable name ("null attribution," "mechanism-resolved controls," "positive-control attribution") used in the title of §5, the checklist, and the conclusion materially helps adoption and citation.

**D5. Make the argument's logical structure explicit. [Medium] (Small)**
The paper's full argument is a conjunction: controls *localize* (empirical), Prop 1 *quantifies* (closed form), Theorem 1 *forecloses* (structural) — and no single leg suffices (the steelman of H3 in particular is killed only by the theorem, since "structure visible only in full sequences and mimicked by reflection" is logically possible). One short paragraph or a small diagram laying out which claim rests on which leg preempts several referee attacks at once.

**D6. Consider a structured abstract. [Low-Medium] (Small)**
The two audiences (crypto, methodology/statistics) would both be served by Context / The false positive / Method / Result / Transferable lesson structure, which also defuses the "sensational numbers first" problem structurally rather than by deletion.

**D7. Appendix A upgrades. [Low] (Small)**
Per-phase pointers to repo scripts and one-line descriptions of each phase's null; remove or version-freeze the "post-submission, see repository" rows.

**D8. Notation table. [Low] (Small)**
The chain δ_τ(k) → z(δ(k)) → s(k) → d(k) → f_j(k) is built up across three sections; a five-row table at the end of §2 would save every reader one backward pass.

## E. Positioning and impact

**E1. Be explicit about the AI-assisted nature of the exploration. [High] (Small)**
The paper says "automated, multi-phase exploration" and quotes an "internal report" that drew the wrong inference, but never says what generated it. If this was an LLM-agent pipeline, saying so plainly — with a short subsection on what the agent did, what it concluded, and where the human audit intervened — turns this into one of the first documented case studies of an AI-generated mathematical false positive *with a full post-mortem*. That is the timeliest thing in the paper and is currently buried. It also changes the citable audience: AI-for-math (cite Davies et al. Nature 2021, the Ramanujan Machine and its spurious-pattern critiques), automated-discovery norms, and agent-evaluation communities.

**E2. Dual-venue strategy. [Medium] (Small)**
The full paper fits IACR Communications in Cryptology (which accepts negative/systematization results) or a methodology-friendly crypto venue; a distilled 4-page version of §5–6 + checklist would fit Patterns / Harvard Data Science Review / an AI-for-science workshop. Decide which is the canonical version before revising, since the framing edits differ.

**E3. Borrow from symmetric-crypto distinguisher practice. [Medium] (Small)**
The block-cipher/TRNG community has mature norms for calibrating many statistical distinguishers (NIST SP 800-22's known issues with dependent tests, multiple-test corrections for randomness batteries). Citing it shows where calibrated-null practice already exists inside cryptography.

**E4. A cross-field "precedents" mini-table. [Medium] (Small)**
One table mapping each pitfall to a named precedent: running-sum low-pass artifact ↔ spurious regression/unit roots (Granger–Newbold); permutation under dependence ↔ GWAS/neuroimaging (Eklund); wrong null model ↔ network-science rewiring debates; forking paths ↔ Gelman–Loken. Demonstrates the lessons are convergent across fields and crypto is the outlier in lacking them.

## F. Defensive (anticipated attacks not covered above)

**F1. "A competent analyst would have caught this immediately." [Medium-High] (Medium)**
Counter concretely: (i) the permutation control *is* the textbook recommendation in most applied references; (ii) show that standard randomness batteries (NIST-style shuffling/spectral tests) also endorse the δ digit sequences — i.e., the false positive survives the tools an ordinarily diligent analyst would deploy. A half-page "what standard practice would have said" run is cheap and devastating to this objection.

**F2. "Exchangeability is textbook — the lesson isn't new." [Medium] (Small)**
Concede detection-of-invalidity is textbook (B5) and sharpen that the contribution is *attribution*: a constructive protocol (mechanism nulls + positive control + closed-form growth signature) that identifies *which* mechanism, with a quantitative match. One sentence distinguishing "this control is invalid" (known) from "here is exactly what the signal is, to the constant" (this paper).

**F3. "z = 44 from a 300-draw null is statistically meaningless." [Medium] (Small)**
Handled by B8 + the analytic null of A1/A2: once the rand-null mean and variance are closed-form, the standardized gap has an honest reference and the 300 draws become a check, not the definition.

**F4. "The half-sequence collapse is just power loss." [Medium] (Small)**
Handled by B2; flag it explicitly in §5 text since it is the most natural objection to the H2 negative control.

**F5. "Your positive control is unfalsifiable — any sequence with the reflection would match." [Low-Medium] (Small)**
Not quite, and the paper should say why: the positive control makes a *quantitative* prediction (elevation magnitude and √n growth, matched within ±2 z across two orders of magnitude in n), and the construction in C1 supplies the converse kill test. State the falsifiability conditions explicitly.

**F6. "Phase 47 is post-hoc; the headline 44σ row is cherry-picked." [Low] (Small)**
The paper already brackets it honestly; add that the Phase-47 value was *predicted in advance* by the 0.43√n law (predicted 44.8, observed 44.6) — that makes the post-hoc row a successful out-of-sample prediction, far stronger than "reported separately." This framing currently appears only inside Prop 1's discussion and deserves a sentence in §5 and the abstract.

## Top-10 priority shortlist (highest leverage per unit effort)

1. **A1** — remove the empirical step from Prop 1's proof (impropriety/pseudo-variance argument).
2. **D1(c)+D2** — null-histogram figure and hypothesis×control matrix.
3. **A2** — closed-form 0.43 constant with explicit standardization.
4. **C1** — asymmetric-section kill test (new falsifiable prediction).
5. **E1** — make the AI-assisted-exploration angle explicit.
6. **B1+B2** — TOST equivalence bound and half-test power analysis.
7. **A3** — pseudo-variance classification of constraint-induced U² inflation (the auditor's lookup table).
8. **B4+B5** — situate within negative-control/exchangeability literature.
9. **A5** — state Theorem 1 for general abelian-group extensions.
10. **D3+F6** — boxed practitioner checklist; reframe Phase 47 as an out-of-sample prediction.

