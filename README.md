# ecdsa-lift

A reproducible methodology case study: how an automated *p*-adic ECDLP
lift-error search produced a persuasive **false positive**, and how
mechanism-resolved null controls traced it to a single public identity.

**Paper (canonical):** `paper/ecdsa_lift_paper_v4.pdf`
(source `paper/ecdsa_lift_paper_v4.tex`). Older drafts are in
`paper/archive/`. This is **not** a claim that ECDSA is broken or safe;
the experiments are diagnostic, at toy sizes, and the cryptanalytic
verdict is the expected negative one.

## The one-sentence result
For ordinary curves with `gcd(n,p)=1`, an apparently growing Gowers/Fourier
anomaly in the section lift error `δ(k)` is entirely explained by the
public antisymmetry `δ(k)+δ(n−k)=[n]τ(G)`; a *Mechanism-Resolved Null
Protocol* (mechanism-derived nulls + positive/converse/sham controls)
audits such automated searches better than a generic permutation control.

## Reproduce
```bash
python -m pip install -r requirements.txt        # numpy 2.2, sympy 1.14, matplotlib 3.10
python experiments/phase43_revision_controls.py  # Table 1  (~2 min)
python experiments/phase61_mass_variance.py      # the 0.427√L constant (<1 min)
python experiments/phase62_figure.py             # the figure
bash run_all.sh                                  # every phase, fixed seed 20260609 (<1 hr)
```
Each phase writes a JSON record under `results/`. Seeds are fixed
(`20260609`); the analysis plan was fixed before the clean re-run.

## Layout
- `src/` — curve / projective / formal-group / lift primitives
- `experiments/` — `phaseNN_*.py`, one standalone script per phase (1–63)
- `results/` — committed JSON outputs
- `paper/` — LaTeX source, the compiled PDF, and `archive/` for old drafts
- `notes/` — internal revision notes

## Phase map (high level)
1–37 exploration (all inert); 37 the Gowers anomaly; 40–43 mechanism
controls + clean re-run; 41–42, 51 theorem + anomalous-boundary checks;
47 large-prime confirmatory; 48,55–59,61,63 the closed-form artifact
(3/2 mass, 0.427√L, U³ baseline, pseudo-variance, finite-alphabet);
49/49b first converse attempt (confounded; superseded); 50 CM-independence;
52–54 calibration; 60 genuine-section converse; 62 figure.

## Status
Research artifact accompanying the paper. An archival DOI (Zenodo) will be
minted for the camera-ready; until then cite the `paper-v4r2` release.
License: MIT (see `LICENSE`); citation metadata in `CITATION.cff`.
