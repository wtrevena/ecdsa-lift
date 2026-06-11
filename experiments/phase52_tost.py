"""Phase 52 - TOST / minimum-detectable-effect equivalence bounds (real vs syn).

The paper's decisive claim is that real delta sequences are statistically
INDISTINGUISHABLE from reflection-only synthetics: z(real vs syn) in [-1.4,+1.5].
"Not significant" is not "equivalent". Here we put an upper bound on any residual
"H3" structure (extra U^2 mass beyond the reflection baseline) that the data can
EXCLUDE at 95%, and express it as a fraction of the real-vs-random ELEVATION.

Setup (per clean ordinary j=0 curve, digit j in {1,2,3}):
  * regenerate the REAL z(delta(k)) sequence from the curve (same degeneracy
    filter as Phase 43), encode digit j, U2_real = ||f||_{U^2}.
  * SYN = 300 reflection-only synthetic draws (same RNG conventions as Phase 43):
    mean mu_syn, sd sd_syn.
  * RAND baseline: 300 S^1-uniform draws: mu_rand. Elevation E = U2_real - mu_rand.
  * observed standardized gap z_obs = (U2_real - mu_syn)/sd_syn.

TOST upper bound on the residual (one-sided 95%): the 95% upper confidence limit
on the true mean difference (mu_real - mu_syn) is
    Delta_excl = (z_obs + z_0.95) * sd_syn         [z_0.95 = 1.645]
(if z_obs<0 we use |center|+margin conservatively via max(0,z_obs)+1.645).
We EXCLUDE any residual H3 component larger than Delta_excl. As a fraction of the
elevation:  resid_frac = Delta_excl / E.

Interpretation printed: "we exclude any residual structure larger than X% of the
elevation."

seed 20260609, SYN=RAND=300, digits {1,2,3}; reuses Phase-43 code.
"""
from __future__ import annotations
import json, math, pathlib
import sys
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from experiments.phase43_revision_controls import (
    z_delta_clean, u2_of_digit, synthetic_reflection_u2, ORDINARY, E_PREC)
from experiments.phase37_gowers import random_baseline

SEED = 20260609
SYN = 300
RAND = 300
DIGITS = (1, 2, 3)
Z95 = 1.6448536269514722   # one-sided 95%


def main():
    rng = np.random.default_rng(SEED)
    np.random.seed(SEED)
    out = {"phase": 52, "seed": SEED, "syn_draws": SYN, "rand_draws": RAND,
           "digits": list(DIGITS), "z95_one_sided": Z95,
           "method": ("Delta_excl=(max(0,z_obs)+1.645)*sd_syn is the 95% upper "
                      "confidence limit on residual U2 beyond the reflection "
                      "baseline; reported as a fraction of the elevation "
                      "E=U2_real-mu_rand."),
           "curves": []}
    for p, Nexp in ORDINARY:
        kept, n, M, Cn = z_delta_clean(p, 7, E_PREC)
        ks = [k for k, _ in kept]; zs = [z for _, z in kept]
        L = len(zs)
        rmean, rstd, _ = random_baseline(L, 2, num_samples=RAND)
        rec = {"p": p, "n": n, "L": L, "digits": {}}
        for j in DIGITS:
            u2_real, _ = u2_of_digit(zs, p, j)
            syn, _ = synthetic_reflection_u2(ks, Cn, M, p, j, rng, SYN)
            mu_syn = float(syn.mean()); sd_syn = float(syn.std())
            E = u2_real - rmean                      # elevation over random
            z_obs = (u2_real - mu_syn) / sd_syn
            delta_excl = (max(0.0, z_obs) + Z95) * sd_syn
            resid_frac = delta_excl / E if E > 0 else float("nan")
            rec["digits"][f"d{j}"] = {
                "U2_real": round(u2_real, 5),
                "mu_syn": round(mu_syn, 5), "sd_syn": round(sd_syn, 5),
                "mu_rand": round(rmean, 5),
                "elevation_E": round(E, 5),
                "z_real_vs_syn": round(z_obs, 2),
                "z_real_vs_rand": round((u2_real - rmean) / rstd, 2),
                "Delta_excl_95": round(delta_excl, 5),
                "resid_excluded_frac_of_elevation": round(resid_frac, 3),
                "resid_excluded_pct": round(100 * resid_frac, 1),
            }
        d1 = rec["digits"]["d1"]
        print(f"p={p:>3d} L={L:>3d}  z(real/syn)={d1['z_real_vs_syn']:+.2f}  "
              f"E={d1['elevation_E']:.4f}  exclude residual > "
              f"{d1['resid_excluded_pct']:.1f}% of elevation (d1)")
        out["curves"].append(rec)
    # aggregate: median excluded fraction across curves x digits
    fracs = [d["resid_excluded_frac_of_elevation"]
             for c in out["curves"] for d in c["digits"].values()
             if not math.isnan(d["resid_excluded_frac_of_elevation"])]
    out["summary"] = {
        "median_resid_excluded_pct": round(100 * float(np.median(fracs)), 1),
        "min_resid_excluded_pct": round(100 * float(np.min(fracs)), 1),
        "max_resid_excluded_pct": round(100 * float(np.max(fracs)), 1),
        "interpretation": ("with 300 synthetic draws we exclude, per curve, any "
                           "residual H3 component larger than the listed % of the "
                           "real-vs-random elevation; the reflection model accounts "
                           "for the elevation up to that residual."),
    }
    print(f"SUMMARY: median excluded residual = "
          f"{out['summary']['median_resid_excluded_pct']}% of elevation "
          f"(range {out['summary']['min_resid_excluded_pct']}-"
          f"{out['summary']['max_resid_excluded_pct']}%)")
    res = pathlib.Path(__file__).resolve().parent.parent / "results"
    res.mkdir(exist_ok=True)
    (res / "phase52_tost.json").write_text(json.dumps(out, indent=2, default=str))
    print("wrote results/phase52_tost.json")


if __name__ == "__main__":
    main()
