# Referee Report

**Manuscript:** *When a Control Does Not Control: A Mechanism-Derived Refutation of a False Positive in p-adic Cryptanalysis of the ECDLP* (v3)
**Author:** W. T. Trevena
**Reviewer recommendation:** **Major revision (leaning toward acceptance)** at a methodology- / negative-results-friendly cryptography venue, conditional on repairing the converse control and the proof of Proposition 2, and softening Corollary 1.

---

## 1. Summary of the submission

The paper is a methodological case study, not a new attack. During a 47-phase automated exploration of p-adic / formal-group lifting attacks on the ECDLP, a Gowers U²/Fourier probe on the Teichmüller lift error δ(k) produced a large, growing effect (≈3–14 standardized units on a clean set, ≈44 on a large-prime curve) that a standard value-permutation control endorsed as "real, index-ordered structure." The paper argues this was a *false positive against the wrong null*, and dissects it completely:

- It builds **mechanism-derived nulls** (increment-shuffle for a prefix-sum hypothesis H1; a half-sequence restriction and a reflection-only positive control for an antisymmetry hypothesis H2; the claimed "new structure" H3) and localizes the entire signal to one exact public identity, the antisymmetry δ(k)+δ(n−k)=[n]τ(G).
- It gives a **closed-form account** (Prop. 1–2): a pseudo-variance identity E|Z|⁴ = 2σ⁴+|ρ|² classifies U² mass inflation by index-coupling constraints, yielding the maximal 3/2 inflation for index reversal and a growth law z(syn/rand) = 2^{7/4}(3^{1/4}−2^{1/4})√L ≈ 0.427√L.
- It proves an **inertness theorem** (Thm. 1) at abelian-group-extension generality: for gcd(n,p)=1 the only secret-dependent part of a section lift error is [k]c, governed by k mod p^{e−1}, which is not determined by the public data; the visible δ is identity-or-artifact.
- It demonstrates the **boundary is sharp**: at n=p (anomalous curves) the ambiguity is annihilated and the surviving invariant recovers the secret — Smart's attack — verified 52/52.

The positioning (a reusable "mechanism-resolved falsification protocol" plus a U² artifact signature, with an honest scope map) is the right one, and the manuscript is substantially stronger than the description of the prior version suggests. The remaining issues are narrow but include one that is, as written, a genuine defect.

---

## 2. Verification performed by this referee

I independently re-derived the central mathematics and reproduced the key experiments (Python/NumPy/SymPy, matching the paper's stated stack). Results:

**Probability/Fourier core — confirmed.**
- Proposition 1, E|Z|⁴ = 2σ⁴+|ρ|² for a complex Gaussian: re-derived via Isserlis (E[X⁴]=3a², E[Y⁴]=3b², E[X²Y²]=ab+2c² ⇒ 2σ⁴+|ρ|²); Monte Carlo gives 8.998 vs. theory 9.000. ✔
- The inflation factor 1+⟨|ρ|²⟩/(2σ⁴)∈[1,3/2] and the maximal 3/2 for index reversal: simulated mass ratios across L∈{78,…,3000} land on 1.490–1.503. ✔
- The growth constant: simulated z(syn/rand)/√L = 0.42–0.44 across the same range, matching 2^{7/4}(3^{1/4}−2^{1/4}) = 0.4267. ✔
- U³: the reflection's U³ *excess* decays like 1/L (0.0009 → 0.0002 → 0.00004 at L=64,128,256) while the U² ratio holds at 3/2 — consistent with the paper's "U³ and higher sit at baseline." ✔

**Number theory — confirmed.**
- Curve orders for y²=x³+7: #E(F₆₇)=79, #E(F₁₆₃)=139, #E(F₂₁₁)=199, all prime (matches Table 1). ✔
- Lemma 1: with m≡1 (mod n), m≡0 (mod p^{e−1}), Gₑ=[m]τ(G) satisfies [n]Gₑ=O and red(Gₑ)=G. ✔
- Antisymmetry (3): with correct (additive) kernel arithmetic at e=2, z(δ(k))+z(δ(n−k)) is single-valued (constant 1340 across all reflection pairs on p=67). ✔
- §8 boundary: y²=x³+4x+7 over F₅₃ is anomalous (#E=53); Smart's attack with a proper odd section recovers **52/52** residues. ✔
- The Hasse dichotomy: brute-force over all curves shows exactly one family with N=2p at p=5 and **no** curve with p|N, N≠p for p=7,…,29 — confirming the "empty interior" claim. ✔

**Two findings from this verification feed directly into the comments below:**

1. *The converse control's section is not a section.* Implementing the "smaller integer representative" y-root rule, I find it fails red(τ′)=id on **exactly 50%** of affine points (39/78 on F₆₇, 99/198 on F₂₁₁, 26/52 on F₅₃): for half the points it lifts −P, not P. This is precisely why "half the lift errors leave the kernel of reduction," and it is the substance of Major Comment M1. (As a clean illustration of the same sign trap, when I first ran Smart's attack with a root-choosing rule that ignored the residue, recovery was 26/52 with every failure of the form k_rec = p−k — a pure sign flip; respecting the residue fixed it to 52/52.)

2. *The proof of Proposition 2 omits its one non-trivial step.* The asserted null standard deviation SD = 2^{−7/4}L^{−3/4} requires Var(U² mass) ≈ 4/L³. The naive "independent modes" estimate gives 20/L³ (since Var(|Z|⁴)=20σ⁸ for a proper complex Gaussian) — five times too large, which would give a constant ≈0.19√L, not 0.427√L. The factor of 5 comes from the Parseval constraint Σ_ξ|f̂(ξ)|²=1: modelling the periodogram as a uniform Dirichlet vector gives Var(mass) = 4(L+5)/[(L+1)(L+2)(L+3)] − 4/(L+1)² = 4/L³ + O(L⁻⁴), which my simulation confirms (var_rand = 8.1e-6 at L=78 vs 4/L³=8.4e-6, not 20/L³=4.2e-5). So the final constant is **correct**, but the proof as written neither derives this variance nor even keeps "mass" and "norm" straight (see M2).

---

## 3. Major comments

### M1. The asymmetric-section converse control is not, as constructed, a section — and this confounds the converse test. *(Must fix.)*

§5 replaces the x-Teichmüller section with τ′ that "selects the Hensel y-root by a negation-insensitive rule (smaller integer representative)," reports that the antisymmetry loses single-valuedness, and concludes the elevation collapses because the reflection was broken. The manuscript then says, revealingly, that "with τ′ exactly half the lift errors leave the kernel of reduction."

That sentence is the tell. For *any* genuine section τ′ (red∘τ′=id), one has δ_{τ′}(k)=[k]τ′(G)−τ′(kG), which reduces mod p to kG−kG=O — so δ_{τ′}(k) lies in the kernel of reduction for *every* k, regardless of whether τ′ is odd. (Oddness, τ′(−P)=−τ′(P), is what the *antisymmetry* (3) needs; kernel-membership of δ is automatic for any section.) If half the lift errors leave the kernel, τ′ is not a section: the "smaller representative" rule selects the root reducing to −y for the points where that root is the smaller integer, so τ′ lifts −P there. I confirmed this directly: the rule fails red(τ′)=id on 50% of points on every curve tested (§2 above).

Consequently the experiment, as written, does not demonstrate "break the reflection while preserving the section-lift-error attack surface." It demonstrates the weaker "a non-section lift rule that is half-wrong destroys the identity and the statistic" — which a hostile referee will reject as confounded (the collapse could be an artifact of points leaving the kernel and being filtered, not of breaking oddness).

**Fix (constructive).** Build a *genuine* asymmetric section by perturbing a valid odd section with a non-odd kernel offset:
τ′(P) = τ(P) ⊕ ι(u(P)), with u(P) ∈ Ê(pZ/p^e) chosen so that u(−P) ≠ −u(P) (e.g. u(P) a fixed nonzero kernel element for P in a chosen half of the orbit and 0 on the complement, or u(P)=φ(x(P)) for any non-even φ). Then red(τ′)=id holds by construction (offsets are in the kernel), δ_{τ′} stays in the kernel for all k, the usable length is not halved, and τ′(−P)=−τ′(P) genuinely fails. Re-running on this τ′ would make the converse control airtight. I would be glad to see this; if the elevation still collapses (as it should — the reflection is exactly the constraint generating the U² inflation, per Prop. 2), it is the single most convincing experiment in the paper. Until then the converse-control claim should be marked provisional.

### M2. The proof of Proposition 2 skips the only step that is not bookkeeping, and conflates "mass" with "norm." *(Must fix; the result itself is correct.)*

The proof writes "the i.i.d. mass has standard deviation Θ(L^{−3/4}); … its leading constant is SD = 2^{−7/4}L^{−3/4}." Two problems:

(a) **Units.** The U² *mass* M=Σ_ξ|f̂(ξ)|⁴ has mean 2/L and standard deviation Θ(L^{−3/2}); it is the U² *norm* M^{1/4} that has SD Θ(L^{−3/4}). The sentence attaches the norm's scaling to the word "mass." Please separate them explicitly.

(b) **The missing variance.** Obtaining SD(norm)=2^{−7/4}L^{−3/4} via the delta method requires Var(M)=4/L³. The natural first guess — independent proper-complex-Gaussian modes — gives Var(M)=L·Var(|Z|⁴)=20/L³ and hence a constant ≈0.19√L, contradicting the data. The correct value 4/L³ is forced by the Parseval constraint Σ|f̂|²=1 (a uniform-Dirichlet periodogram), which reduces the variance by exactly a factor of 5. This is *the* substantive analytic step behind the 0.427 coefficient, and the proof currently omits it. I verified both the analytic Dirichlet computation (Var(M)=4/L³+O(L⁻⁴)) and the simulation. Please add the one- or two-line Dirichlet/Parseval variance derivation; it converts the "fully analytic" claim from asserted to demonstrated, and it is genuinely the nicest part of the calculation.

Relatedly, "the reflection inflates the U² mass by **exactly** 3/2" coexists with o(1/L) remainders. State this as the *limiting* inflation factor (mass = 3/L+o(1/L) vs 2/L+o(1/L)), and handle the finite-L edge cases (zero mode; the self-paired ξ=L/2 mode for even L; the L-odd midpoint) explicitly — my simulations show their effect is O(1/L), so a sentence suffices, but a careful referee will look for it.

### M3. Corollary 1 states a computational-hardness claim as a corollary of a group-theoretic theorem. *(Should fix.)*

Theorem 1 (a)–(d) is a clean algebraic statement: δ_τ factors through k mod n iff c=0, and for c≠0 the secret-bearing part [k]c is governed by k mod p^{e−1}, which the public instance does not pin down over the natural lift. That is correct and well proved. But Corollary 1 ("no section lift error yields an *efficiently computable* function of (E,G,Q) that depends on the secret") is not a consequence of the theorem alone — it leans on ECDLP-hardness and on "no efficient algorithm is known," which the manuscript correctly flags as a state-of-knowledge remark elsewhere. Calling it a Corollary gives a computational claim the force of a theorem. Either (A) rename it ("Consequence, under the standard ECDLP-hardness assumption") or (B) split it into the formal equivalence (δ_τ factors through k mod n ⇔ c=0) plus a clearly labelled computational remark. Cryptography referees are acutely sensitive to exactly this algebraic-vs-computational boundary.

### M4. The abstract's "5–11σ is cheap" elides the more interesting mechanism. *(Should fix.)*

There are two distinct false-positive engines, and the paper conflates them. (i) *Look-elsewhere*: many probes make moderate |z| outliers likely; the manuscript's own calibration (expected max |z|≈3.2–3.6; an empirical battery hitting 5.0) supports this for ~5σ. (ii) *Wrong null*: a real public algebraic identity tested against a generic random baseline produces an effect that grows without bound (here ∝√L, reaching 44). The headline 11σ (and the 44) is overwhelmingly mechanism (ii), not (i) — multiplicity arithmetic cannot manufacture 44σ. The current phrasing risks suggesting the whole thing is multiple testing, which undersells the paper's actual and sharper point: the baseline was wrong because the statistic encodes a public identity. Suggest revising to distinguish the two explicitly (e.g., "apparently anomalous separations are cheap — they arise either from multiple probing or, as here, from real public structure tested against the wrong null").

### M5. The methodology needs statistics/causal-inference scholarship, not only crypto references. *(Should fix.)*

The transferable thesis — "a control is valid only against a named alternative; derive nulls from mechanisms; require a positive control" — is, as the paper honestly notes, textbook in spirit. Its contribution is the disciplined *import* into automated cryptanalysis plus the to-the-constant dissection. To protect against the "this is just textbook statistics" reflex, cite the relevant lineage: negative/positive controls and placebo tests (Lipsitch–Tchetgen Tchetgen–Cohen), exchangeability requirements for permutation tests, specification-curve / multiverse analysis (Simonsohn; Steegen et al.), garden-of-forking-paths and multiple-comparison problems (Gelman–Loken; Benjamini–Hochberg), null-model controversies in ecology/networks, and the cryptographic randomness-test-battery literature (NIST SP 800-22 and its known calibration pitfalls; recent battery-recommendation surveys). Then make the framing explicit: the principle is textbook; the contribution is showing, exactly, what a plausible-but-invalid control launders, and providing a structural reason it had to.

### M6. Reproducibility coordinates are incomplete for a paper whose thesis is reproducibility. *(Should fix.)*

A static paper cannot point to "the ecdsa-lift repository." Provide the full URL **and** the full commit hash **and** the release tag, plus an archival DOI (Zenodo/Software Heritage). Unify "commit 674750a" with the "release tagged with this paper." Give expected runtime and hardware, and a *minimal* reproduction command for the headline claims (Table 1, Prop. figures), not only `run_all.sh` (which, as written `run all.sh`, has a space — confirm the actual filename). Phases 44–46 are listed as "post-submission, see repository," which is awkward in a frozen artifact; either freeze them at a stated commit or move them to an explicit "not part of this paper" note.

### M7. Adjacent p-adic / formal-group ECDLP literature should be engaged. *(Should address.)*

Beyond Silverman and Gadiyar–Padma, there is directly relevant work on p-adic/formal-group approaches to the ECDLP and on cocycle/height-function framings (e.g., the line of work on "the ECDLP over the p-adic field and formal groups," and p-adic height functions of cryptanalytic interest). Since the paper's whole object — section lift errors as 1-cochains valued in the kernel of reduction, with a cohomology-vanishing inertness argument — sits squarely in this space, a short related-work paragraph situating the contribution against it (and confirming the inertness conclusion is consistent with, or sharper than, prior framings) would strengthen the novelty claim and pre-empt an obvious referee question. I verified all eight existing references are real and correctly characterized; this is about completeness, not correction.

---

## 4. Minor comments

1. **Name the protocol.** Give the falsification workflow a citable name (e.g. "Mechanism-Resolved Null Protocol" or "Section-Invariance Protocol") and use it in the abstract, §11, and conclusion.
2. **Add a figure.** A methodology paper should not rest on dense tables alone. The highest-value figure: (A) overlaid null distributions of the U² norm — random, increment-shuffle, reflection-synthetic, with the real value marked; (B) z(syn/rand) and z(real/rand) vs √L with the 0.427√L line; (C) a hypothesis × control matrix (H1/H2/H3 against permutation / increment-shuffle / half / positive / converse), making visually obvious that the permutation control is the only column that fails to separate.
3. **Demote the carry model to a model.** "Modelling the rest as independent gives…" reads as derivation; write "a simple independence model for borrow failures predicts…". State that q_j depends on C, the digit distribution, p, and the sample, so "≈0.2" is empirical, not universal. And: if the L≈10638 row already bears on the digit-1 ≈3.0 / digits-2,3 ≈2.7 prediction, report the measured per-digit masses rather than only forecasting them at L≳3000.
4. **U³ baseline: prove or soften.** The 1/L decay of the U³ excess is currently asserted as derived. Either add a short lemma (the extra difference operator couples reflection partners on an O(1/L) fraction of shift–position pairs) or write "heuristically, and confirmed numerically." My simulation supports the 1/L decay, so a brief argument should be attainable.
5. **TOST detail.** Specify the equivalence margin, whether the test is on mass or norm, whether the residual is defined in raw statistic / standardized z / fraction of elevation, whether digits are pooled, and whether 300 draws suffice for the stated interval.
6. **Half-sequence power.** Give the per-curve (or min/max) pair (z_expected-half, z_observed-half) so the 4–9σ retained-power claim is auditable from the paper, not only the repository.
7. **Abstract density.** The abstract carries the full effect-size catalogue, both controls, both propositions, the exact 0.427 constant, the theorem, the boundary, the caveat, and the lesson. Move the exact coefficient, "52/52," and the detailed theorem phrasing into the introduction; keep the abstract to the five load-bearing claims.
8. **Title.** "Refutation" can read as refuting an external program; since the target is an internally generated false positive, consider "Mechanism-Derived Nulls for a False Positive in p-adic ECDLP Cryptanalysis" if the venue is conservative. The current title is acceptable.
9. **Notation table.** A short table for δ, z(δ), digit_j, the phase encoding f_j, and U² mass vs norm would help non-specialist readers; the audience for a methodology paper is broader than the p-adic-ECDLP audience.
10. **AI-assisted disclosure.** §10 has an "automated and AI-assisted exploration" subsection but never says whether the 47-phase pipeline was LLM-assisted, agent-generated, or scripted. If AI was materially involved, say what it generated, what the internal report concluded, and where human audit intervened — that would make the paper a genuinely valuable case study in AI-assisted false positives. If not, soften the language so it does not read as opportunistic positioning.
11. **One-line search-reduction caveat.** Add a sentence making explicit that the reflection yields only the familiar k ↔ n−k negation symmetry (a 2-to-1 reduction), not any new advantage — so no reader mistakes "real public structure" for a cryptanalytic gain.
12. **"Indistinguishable."** Continue to use it only in the bounded TOST sense you define; flag once more at first use in §5 that it means "equivalent within the stated margin," not "no difference."

---

## 5. Opportunities / suggested directions (optional, for the authors' consideration)

These are not required for acceptance; they are where I think the work could reach beyond a single negative result toward something more significant.

**D1. Promote the pseudo-variance result into a standalone "improper-Gaussian audit."** Proposition 1 is not ECDLP-specific: any time an automated pipeline reports elevated 4th-moment / U² Fourier mass, the per-mode pseudo-variance ρ(ξ) diagnoses whether it is an index-coupling artifact, and the [1,3/2] classification reads the responsible *linear constraint* off the data. Packaged as a general diagnostic (estimate the pseudo-variance spectrum; identify which constraint — reflection, fixed shift, automorphism — generates the inflation; predict the growth exponent), this is more broadly citable than the ECDLP case study and could anchor its own paper. A natural extension: derive the inflation dictionary for higher Gowers norms U^k and for multidimensional Fourier (Gowers on (Z/L)^d), giving "constraint ↦ inflation factor" in general.

**D2. Hunt for the regime the theorem proves is empty *for elliptic curves over prime fields*.** The cleanest result in §8 is that the partial-degeneration regime (0 ≠ ker([n]|Ê) ⊊ Ê) has empty interior on E(F_p), by Hasse. But the manuscript itself notes this can fail for E/F_{p^m} and for higher-dimensional abelian varieties, where the kernel of reduction is a p-group of higher rank and [n] can have a *nontrivial but proper* kernel — a "graded leak." That is exactly where a genuine new lift attack could live, and the paper's machinery is already set up to look. A systematic search for near-anomalous structures in g≥2 Jacobians / abelian surfaces, or over small extension fields, guided by the theorem's one hypothesis, is the most promising route to a *positive* result and would convert this from "we explained why nothing was here" to "here is precisely where to look next."

**D3. Prove a completeness statement for public constraints on δ.** The exploration found exactly two surviving public relations: the antisymmetry and the automorphism-equivariance of the leading digit. Theorem 1 explains why anything public must be identity-or-artifact, but it does not enumerate the identities. Characterizing *all* public linear constraints on δ (and showing the U²/Gowers signature is fully determined by that finite list) would upgrade "we tried 47 things and they failed" to "here are all the things that can appear, and here is why the list is complete." That is a qualitatively stronger paper.

**D4. Use the growth law as an experiment-design tool.** The 0.427√L law lets the authors *solve for* the L needed to exclude an H3 component of a target size, replacing the post-hoc admission that small curves are underpowered with an a-priori power calculation. This would make the equivalence claims quantitatively tight and is essentially free given the closed form.

**D5. Build a small benchmark of cryptanalytic false positives.** Given the AI-assisted-exploration framing, formalize the protocol as a checklist and pair it with an open benchmark of known false positives (this U² one plus a handful of others), demonstrating that the protocol catches them while a naïve permutation test does not. This is the move that turns one case study into a reusable, citable methodology — likely the highest-leverage direction for impact.

**D6. Transfer the lens to other "lift-and-hope" primitives.** The section-invariance discipline and the pseudo-variance audit are not ECDLP-specific. A companion application to lattice/LWE structure, isogeny-based schemes, or discrete logs in class groups would multiply the methodological reach and is a natural follow-on.

---

## 6. Recommendation, itemized

**Must fix before acceptance**
- **M1** — repair the converse control: prove τ′ is a section or replace it with a genuine asymmetric section (kernel-offset construction), and re-run.
- **M2** — complete the Proposition 2 proof: add the Parseval/Dirichlet variance derivation (Var(mass)=4/L³) and separate "mass" from "norm"; state 3/2 as the limiting factor and dispatch the finite-L edge modes.
- **M3** — rename or split Corollary 1 so the computational claim is not presented as a theorem.

**Should fix**
- **M4** — distinguish look-elsewhere from wrong-null in the abstract.
- **M5** — add statistics/causal-inference references and make the "principle is textbook, contribution is the import + dissection" framing explicit.
- **M6** — full reproducibility coordinates (URL + commit + tag + DOI, runtime, minimal command); resolve Phases 44–46.
- **M7** — engage the adjacent p-adic/formal-group ECDLP literature.
- Minor items 1, 2, 3, 5 (figure; protocol name; carry-model demotion; TOST detail).

**Nice to have**
- Minor items 4, 6–12, and any of the directions in §5.

**Overall.** The mathematics I was able to check is correct, including the headline 3/2 inflation, the 0.427√L law, Smart's 52/52 recovery, and the Hasse dichotomy; the inertness theorem is a sound (and honestly credited) repackaging; and the pseudo-variance classification is a clean, genuinely reusable contribution. The paper has a coherent, defensible identity as a mechanism-resolved falsification protocol for automated cryptanalytic false positives. One construction (the converse control) is, as written, mathematically invalid, and one proof (Prop. 2) omits its key step; both are fixable without new ideas. With those repaired and the computational/algebraic boundary cleaned up, I would support acceptance.
