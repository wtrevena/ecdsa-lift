"""Phase 61 - the U^2 mass variance, and the closed-form 0.427 constant.

Referee objection (v3 review, M2): the proof of Proposition 2 asserts the
i.i.d. null standard deviation SD = 2^{-7/4} L^{-3/4} but (a) attaches it to
the word "mass" when it is the SD of the U^2 NORM, and (b) never derives the
variance of the mass, which is the only non-trivial ingredient behind the
0.427 coefficient.

This phase pins down the variance both ways.

Setup. For f : Z/L -> S^1, the U^2 "mass" is M = sum_xi |fhat(xi)|^4 with
fhat(xi) = (1/L) sum_x f(x) e^{-2pi i x xi / L}. Parseval fixes
sum_xi |fhat(xi)|^2 = (1/L) sum_x |f(x)|^2 = 1 EXACTLY. So the periodogram
values p_xi := |fhat(xi)|^2 sum to 1; for a random-phase f they behave (in
the large-L limit) like a uniform Dirichlet(1,...,1) vector on the
(L-1)-simplex, NOT like independent exponentials.

Two predictions for Var(M):
  * Independent proper-complex-Gaussian modes (IGNORING Parseval):
        Var(|Z|^4) = 20 sigma^8  (|Z|^2 ~ Exp(sigma^2), sigma^2 = 1/L)
        => Var(M) ~ L * 20/L^4 = 20/L^3.
    This would give a growth constant z/sqrt(L) ~ 0.19, contradicting data.
  * Uniform Dirichlet periodogram (RESPECTING Parseval): with the Dirichlet
    moment formula E[prod p_i^{a_i}] = prod(a_i!) * Gamma(L)/Gamma(L+sum a_i),
        E[M]   = L * 2/(L(L+1))                 = 2/(L+1)
        E[M^2] = L*24/D + L(L-1)*4/D, D=L(L+1)(L+2)(L+3)
               = 4(L+5)/[(L+1)(L+2)(L+3)]
        Var(M) = 4(L+5)/[(L+1)(L+2)(L+3)] - 4/(L+1)^2
               = 4/L^3 + O(L^-4).
    The factor-of-5 reduction (20 -> 4) is exactly the Parseval constraint.

Then via the delta method on the NORM = M^{1/4} (mean (2/L)^{1/4}):
    SD(norm) = (1/4) * (2/L)^{-3/4} * SD(M)
             = (1/4) * 2^{-3/4} L^{3/4} * 2 L^{-3/2}
             = 2^{-7/4} L^{-3/4},     using SD(M) = 2 L^{-3/2}.
Mean gap (reflection 3/L vs random 2/L):
    (3/L)^{1/4} - (2/L)^{1/4} = (3^{1/4} - 2^{1/4}) L^{-1/4}.
Standardized:
    z(syn/rand) = (3^{1/4}-2^{1/4}) L^{-1/4} / (2^{-7/4} L^{-3/4})
                = 2^{7/4} (3^{1/4} - 2^{1/4}) sqrt(L)  ~ 0.4267 sqrt(L).

This script (1) does the exact symbolic Dirichlet computation and its
large-L expansion, (2) Monte-Carlo measures Var(M) for random S^1
sequences and checks it tracks 4/L^3 (not 20/L^3), and (3) confirms the
standardized growth constant equals 2^{7/4}(3^{1/4}-2^{1/4}).

Outputs results/phase61_mass_variance.json.
"""
from __future__ import annotations
import sys, pathlib, json, time
import numpy as np
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SEED = 20260609


def symbolic_dirichlet():
    """Exact Var(sum p_i^2) for p ~ Dirichlet(1,...,1) on L parts, and its
    large-L expansion. Returns dict of sympy expressions as strings."""
    L = sp.symbols("L", positive=True)
    # Dirichlet(1,...,1): E[prod p_i^{a_i}] = prod(a_i!) * Gamma(L)/Gamma(L+sum a_i)
    # E[p_i^2]   : a_i=2  -> 2! * Gamma(L)/Gamma(L+2) = 2/(L(L+1))
    # E[p_i^4]   : 4! * Gamma(L)/Gamma(L+4) = 24/(L(L+1)(L+2)(L+3))
    # E[p_i^2 p_j^2] (i!=j): 2!2! * Gamma(L)/Gamma(L+4) = 4/(L(L+1)(L+2)(L+3))
    Em = L * (2 / (L * (L + 1)))                      # E[M]
    D = L * (L + 1) * (L + 2) * (L + 3)
    Em2 = L * (24 / D) + L * (L - 1) * (4 / D)         # E[M^2]
    Var = sp.simplify(Em2 - Em**2)
    Var_series = sp.series(Var, L, sp.oo, 5).removeO()
    return {
        "E[M]": str(sp.simplify(Em)),
        "E[M^2]": str(sp.simplify(Em2)),
        "Var(M)_exact": str(Var),
        "Var(M)_largeL": str(sp.simplify(Var_series)),
        "leading_times_L3": str(sp.limit(Var * L**3, L, sp.oo)),
    }


def monte_carlo_mass_variance(rng, L, ndraw):
    """Return (mean_M, var_M, mean_norm, sd_norm) for random S^1 sequences."""
    f = np.exp(1j * rng.uniform(0, 2 * np.pi, (ndraw, L)))
    fhat = np.fft.fft(f, axis=1) / L
    M = np.sum(np.abs(fhat) ** 4, axis=1)
    norm = M ** 0.25
    return float(M.mean()), float(M.var()), float(norm.mean()), float(norm.std())


def reflection_mass(rng, L, ndraw):
    f = np.exp(1j * rng.uniform(0, 2 * np.pi, (ndraw, L)))
    j = np.arange(L); k = L - 1 - j; m = j < k
    f[:, k[m]] = np.conj(f[:, j[m]])    # index-reversal coupling
    fhat = np.fft.fft(f, axis=1) / L
    M = np.sum(np.abs(fhat) ** 4, axis=1)
    return float(M.mean()), float((M ** 0.25).mean())


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    target = 2 ** (7 / 4) * (3 ** 0.25 - 2 ** 0.25)
    out = {
        "params": {"seed": SEED, "target_constant": round(target, 6)},
        "symbolic_dirichlet": symbolic_dirichlet(),
        "monte_carlo": [],
    }
    print("Symbolic Dirichlet:", out["symbolic_dirichlet"])
    print(f"\nTarget growth constant 2^(7/4)(3^(1/4)-2^(1/4)) = {target:.6f}")
    print(f"{'L':>6} {'mean_M':>10} {'2/L':>10} {'var_M':>11} {'4/L^3':>11} {'20/L^3':>11} "
          f"{'ratio':>7} {'z/sqrtL':>8}")
    for L in [78, 138, 198, 312, 612, 828, 3000]:
        ndraw = 8000
        mM, vM, mN, sN = monte_carlo_mass_variance(rng, L, ndraw)
        mMr, mNr = reflection_mass(rng, L, ndraw)
        ratio = mMr / mM
        z = (mNr - mN) / sN
        rec = {
            "L": L, "mean_M": round(mM, 7), "two_over_L": round(2 / L, 7),
            "var_M": vM, "four_over_L3": 4 / L**3, "twenty_over_L3": 20 / L**3,
            "var_ratio_to_4overL3": round(vM / (4 / L**3), 3),
            "mass_ratio_refl_over_rand": round(ratio, 4),
            "z_per_sqrtL": round(z / np.sqrt(L), 4),
        }
        out["monte_carlo"].append(rec)
        print(f"{L:>6} {mM:>10.6f} {2/L:>10.6f} {vM:>11.3e} {4/L**3:>11.3e} {20/L**3:>11.3e} "
              f"{ratio:>7.4f} {z/np.sqrt(L):>8.4f}")
    # delta-method consistency: predicted SD(norm) = 2^{-7/4} L^{-3/4}
    print("\ndelta-method check: SD(norm) vs 2^{-7/4} L^{-3/4}")
    dm = []
    for rec in out["monte_carlo"]:
        L = rec["L"]
        # recompute sd_norm from this L
        mM, vM, mN, sN = monte_carlo_mass_variance(rng, L, 8000)
        pred = 2 ** (-7 / 4) * L ** (-3 / 4)
        dm.append({"L": L, "sd_norm_emp": round(sN, 6), "sd_norm_pred": round(pred, 6),
                   "ratio": round(sN / pred, 3)})
        print(f"   L={L:>5} emp={sN:.6f} pred(2^-7/4 L^-3/4)={pred:.6f} ratio={sN/pred:.3f}")
    out["delta_method"] = dm
    out["elapsed_s"] = round(time.time() - t0, 1)
    (ROOT / "results").mkdir(exist_ok=True)
    path = ROOT / "results" / "phase61_mass_variance.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
