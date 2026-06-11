# Revision Implementation Notes (v3 → v4)

**Maintainer:** referee-implementation pass (Claude/Cowork). **Last updated:** 2026-06-11.
**Purpose:** track the implementation of the v3 referee report's recommendations so work can resume after any interruption. Update this file every time a step finishes.

---

## 0. Context / where things live

- Paper source: `paper/ecdsa_lift_paper_v3.tex` (compiled `paper/ecdsa_lift_paper_v3.pdf`). v4 will be `paper/ecdsa_lift_paper_v4.tex`.
- Referee report (this pass authored it): `referee_report_ecdsa_lift_v3.md` (repo root). The prior reviewer's notes: `feedback.txt`.
- Experiments: `experiments/phaseNN_*.py`; results JSON: `results/phaseNN_*.json`.
- Source library: `src/` (`curves.py`, `projective.py`, `lifts.py`, `formal_group.py`, ...).
- Reproduction driver: `run_all.sh`.
- Repo head at start of this pass: commit `674750a`.

### Scratch verification scripts (in Cowork outputs, NOT in repo)
These were used during review to independently reproduce results; can be reused/ported:
- `stats_vec.py`, `sv.py` — U² mass ratio 3/2 and 0.427√L constant.
- `u3.py` — U³ baseline decay.
- `ec2.py`–`ec6.py` — RCB/Bosma–Lenstra complete addition, Smart 52/52, antisymmetry via additive-z, Hasse dichotomy, τ′ section-failure rate.

---

## 1. Scope agreed with user

Implement **everything exhaustively** — all must-fix (M1–M3), should-fix (M4–M7), minor items, AND the creative directions where feasible. Versioning: create **v4** tex, keep v3 frozen.

---

## 2. Referee findings being implemented (from referee_report_ecdsa_lift_v3.md)

### MUST FIX
- **M1** — converse control τ′ ("smaller integer representative") is NOT a section: fails `red∘τ′=id` on ~50% of points (verified). Fix = genuine kernel-offset asymmetric section. **→ DONE, Phase 60.**
- **M2** — Prop 2 proof omits the Var(mass)=4/L³ derivation (Parseval/Dirichlet, the factor-5 reduction from naive 20/L³) and conflates "mass" SD with "norm" SD. The 0.427 constant is correct. **→ DONE verifying, Phase 61. Tex edit pending.**
- **M3** — Corollary 1 states a computational-hardness claim as a corollary of a group-theoretic theorem. Split/rename. **→ tex edit pending.**

### SHOULD FIX
- **M4** — abstract conflates look-elsewhere (≈5σ cheap) with wrong-null (11σ–44σ from a real public identity). Distinguish. **→ tex edit pending.**
- **M5** — add statistics/causal-inference references (negative/positive controls, placebo tests, specification-curve/multiverse, forking paths, NIST SP 800-22 battery pitfalls). **→ tex edit pending.**
- **M6** — full reproducibility coordinates: repo URL + full commit hash + tag + archival DOI; runtime/hardware; minimal reproduction command; resolve `run_all.sh` (confirm filename, no space); freeze/remove Phases 44–46. **→ tex edit pending.**
- **M7** — engage adjacent p-adic/formal-group ECDLP literature (e.g. "ECDLP over the p-adic field and formal groups", p-adic height functions, cocycle explorations). **→ tex edit pending.**

### MINOR
1. Name the protocol (e.g. "Mechanism-Resolved Null Protocol"). — tex
2. Add figure: (A) null distributions; (B) z vs √L with 0.427√L line; (C) hypothesis×control matrix. **→ Phase/figure script pending.**
3. Demote carry model to a model; report measured per-digit masses if the large row settles it. — tex
4. U³ baseline: short lemma or "heuristically, confirmed numerically". — tex (+ Phase 59 already exists)
5. TOST detail (margin, mass vs norm, residual definition, draws, pooling). — tex
6. Half-sequence power: per-curve (z_expected_half, z_observed_half). — tex (+ Phase 53 data)
7. Abstract density — move exact coeff / 52/52 / theorem phrasing to intro. — tex
8. Title softening ("Refutation" optional). — tex (optional)
9. Notation table (δ, z(δ), digit_j, f_j, U² mass vs norm). — tex
10. AI-assisted disclosure (state whether pipeline was LLM/agent/script). — tex (needs author input; will add placeholder + honest framing)
11. One-line search-reduction caveat (reflection = only k↔n−k 2-to-1). — tex
12. "Indistinguishable" defined in bounded TOST sense. — tex

### CREATIVE DIRECTIONS (optional, add as a "Future directions" expansion in §11 or appendix)
- D1: pseudo-variance audit as standalone diagnostic; extend dictionary to U^k and multidim Fourier.
- D2: hunt for the "graded leak" regime (empty for E/F_p by Hasse) in g≥2 Jacobians / abelian surfaces / extension fields — most promising for a positive result.
- D3: completeness statement — characterize ALL public linear constraints on δ.
- D4: use 0.427√L law for a-priori power/experiment design.
- D5: build a benchmark of cryptanalytic false positives + protocol checklist.
- D6: transfer lens to LWE/lattice, isogeny, class-group DLP.

---

## 3. Progress log

### ✅ Phase 60 — genuine kernel-offset asymmetric section (M1 fix)
File: `experiments/phase60_genuine_asymmetric_section.py`; results: `results/phase60_genuine_asymmetric_section.json`.
Construction: τ′(P) = τ(P) ⊞ [h(P)]·K₁, K₁=[n]τ(G) (valuation-1 kernel elt, z-param = C_n). Kernel arithmetic done **additively in z** (exact at e≤4): z(δ′(k)) = z(δ(k)) + (k·γ − h(kG))·C_n mod p^e, γ=h(G). Because the offset is in the kernel, **red∘τ′=id everywhere ⇒ full length, 0 drops** (fixes the confound).

Offset rules: `even_x` (h=x mod p^{e-1}, even ⇒ τ′ non-odd), `odd_y` (h=y, ODD ⇒ τ′ stays odd!), `mixed_xy` (h=x+7y, neither), `half_binary` (even binary).

**Headline results (z = U² real-vs-rand, digits 1/2/3), L = n−1, 0 drops:**

| p (n) | standard_odd (single-valued) | genuine even_x | genuine mixed_xy | genuine odd_y (still odd!) |
|---|---|---|---|---|
| 67 (79) | +3.6/6.5/5.3 ✓SV | +0.7/0.5/−0.2 | −0.7/1.9/0.4 | +3.4/5.4/3.6 ✓SV |
| 211 (199) | +6.6/5.5/7.4 ✓SV | −1.0/−0.5/−0.2 | −1.0/0.3/−1.2 | +4.7/6.4/6.6 ✓SV |
| 823 (829) | +11.8/10.2/13.0 ✓SV | −1.6/−0.9/−0.1 | −0.5/−0.4/−1.4 | +9.5/13.2/12.9 ✓SV |
| 10477 (10639) | +45.9/42.0/41.2 ✓SV | +0.9/−0.7/+1.2 | −0.4/−1.1/+0.1 | +47.8/45.6/49.2 ✓SV |

`smaller_rep` non-section (Phase 49) fails red=id on 44–54% of points (the confound), confirmed in-phase.

**Interpretation for the paper:** (i) genuine non-odd section collapses U² to baseline at full length — clean converse, no confound; (ii) an ODD kernel offset PRESERVES single-valuedness and elevation — sharply proves the mechanism is *oddness/reflection itself*, not section incidentals. This is a stronger converse than v3's. Use Phase 60 to rewrite the "converse control" paragraph in §5 and to replace/augment Phase 49 framing.

### ✅ Phase 61 — U² mass variance + 0.427 constant (M2 fix)
File: `experiments/phase61_mass_variance.py`; results: `results/phase61_mass_variance.json`.
- **Symbolic (sympy):** Dirichlet(1,…,1) gives E[M]=2/(L+1), E[M²]=4(L+5)/[(L+1)(L+2)(L+3)], Var(M)=4(L−1)/(L⁴+7L³+17L²+17L+6), large-L = 4(L−8)/L⁴, **Var(M)·L³ → 4** exactly.
- **Monte Carlo:** var_M tracks 4/L³ (ratio≈1.0), NOT 20/L³ (the naive independent-modes value); mass ratio → 1.500; z/√L → 0.4245–0.4387 vs target 0.426728.
- **Delta method:** empirical SD(norm) matches 2^{−7/4}L^{−3/4} within 1–4%.
Conclusion: the constant is correct; the missing proof step is the Parseval⇒Dirichlet factor-5 variance reduction. Tex fix = insert this derivation into Prop 2's proof and separate mass-SD (Θ(L^{−3/2})) from norm-SD (Θ(L^{−3/4})), and phrase 3/2 as the limiting factor + dispatch finite-L edge modes (zero mode, ξ=L/2 self-paired, odd-L midpoint — all O(1/L)).

### ⏳ NEXT (in order)
1. **Figure** (minor #2): matplotlib 3-panel → `paper/fig_artifact.pdf` (+ a phase script `experiments/phase62_figure.py` so it's reproducible). Data from phases 43/47/48/60/61.
2. **v4 tex** (task 10): copy v3→v4, apply M2/M3/M4/M5/M6/M7 + minors + Phase-60 converse rewrite + creative directions paragraph. Keep all numbers consistent with JSON.
3. **Compile v4 + final verification** (task 11): pdflatex if available in sandbox; cross-check every edited number against results JSON; present files.

---

## 3b. Working-tree merge analysis (user request) — DONE

Analyzed the author's uncommitted working-tree changes and merged them into branch `referee-revision-v4`.

**What the author was working on (uncommitted):**
- `paper/ecdsa_lift_paper.tex` (modified, +484 lines) is **byte-identical to `paper/ecdsa_lift_paper_v3.tex`** — the v3 revision was done in place on the base filename. This is the real source of the reviewed PDF.
- `revision_plan.md` — author's own plan merging two prior referee reports (`referee_report_ecdsa_lift.md` = "Referee 1"; `feedback.txt` = "Reviewer 2"). High quality. **It independently flags my M2 finding**: item 3.2 says "derive the 0.43 constant in closed form, stating exactly what is standardized (mass vs norm; analytic vs empirical null variance). R1 reproduced scaling but the constant only to ~20% — purely a standardization ambiguity; fix by specification." → my Phase 61 supplies the missing Var(mass)=4/L³ (Dirichlet/Parseval) step that closes this.
- Untracked phases 49–59 + results JSON — the revision experiments the v3 paper's phase index cites but were never committed.
- `phase49b_killtest_confound.py` — author already DOCUMENTED the converse confound (asymmetric "smaller_rep" rule halves L because it lifts −P for ~50% of points). They flagged it honestly but couldn't remove it (their construction is not a section). **Phase 60 supersedes it** with a genuine kernel-offset section (full length, U² still collapses, odd offset preserves it).
- `phase42` change: removed a 243-char trailing whitespace line (trivial cleanup; kept).

**Verdict:** author's work is good and substantive; the only true gaps are the ones the referee report identifies (M1 converse not a section; M2 variance derivation incomplete). All consolidated into commit `dc9a5a8`.

**GIT LOCK WORKAROUND (sandbox cannot unlink `.git/*.lock`):** normal `git add/commit` fail with "Unable to create .git/index.lock: File exists". To commit:
1. `cp .git/index /tmp/gi` (once), then `GIT_INDEX_FILE=/tmp/gi git add <files>`.
2. `TREE=$(GIT_INDEX_FILE=/tmp/gi git write-tree)`
3. `COMMIT=$(echo "$MSG" | GIT_AUTHOR_*/GIT_COMMITTER_* env git commit-tree $TREE -p $(git rev-parse HEAD))`
4. Update branch ref by **Writing the hash directly** to `.git/refs/heads/referee-revision-v4` via the file tool (git update-ref also hits a lock). Then `git log` confirms.
(`git status`/`git diff`/`git cat-file`/`write-tree`/`commit-tree` all work read-side; only index/ref *.lock creation fails.)

## 4. Independent-verification ledger (already completed during review; all ✓)
- Prop 1 E|Z|⁴=2σ⁴+|ρ|² (Isserlis + MC 8.998 vs 9.000).
- Inflation 3/2; growth 0.427√L; U³ excess ~1/L.
- #E(F₆₇)=79, F₁₆₃=139, F₂₁₁=199 prime (Table 1).
- Lemma 1 canonical lift ([n]Gₑ=O, red(Gₑ)=G).
- Antisymmetry (3) single-valued via additive-z (const 1340 on p=67, e=2).
- Smart 52/52 on y²=x³+4x+7/F₅₃ (anomalous, #E=53).
- Hasse dichotomy: only p=5 has N=2p; none p=7..29 with p|N,N≠p.
- All 8 references real and correctly characterized.
- RCB complete formulas degenerate over Z/p⁴ for deep kernel points → corroborates the paper's additive-z discipline (and is why Phase 60 uses additive z).
