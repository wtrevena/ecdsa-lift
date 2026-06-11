"""Phase 63 - pseudo-variance mechanism (conjugation) and finite-alphabet baseline.

Two referee requests (R3, 2026-06-10):

(1) The reflection-induced inflation requires the \emph{conjugate} value
    reflection f(L-1-k)=\overline{f(k)} (equivalently the phase image of
    s(n-k)=C-s(k)). We confirm the mechanism is per-mode self-impropriety,
    not a xi<->-xi coupling: measure the per-mode pseudo-variance ratio
    <|rho(xi)|/sigma^2> over nonzero modes. The conjugate reflection drives
    it to 1 (each mode improper, |rho|=sigma^2 -> U^2 mass 3/L); the
    NON-conjugate constraint f(L-1-k)=beta f(k) leaves it ~0 (proper modes,
    no inflation, mass 2/L). This pins the conjugation as essential and
    identifies the exact Fourier-domain constraint for Proposition 1.

(2) Finite-alphabet baseline. The paper's null uses continuous S^1 phases;
    the real digit encoding uses p-th roots of unity (a finite alphabet).
    We show the i.i.d. mass mean (2/L), the reflection inflation (3/2), and
    the per-mode pseudo-variance are alphabet-independent: continuous S^1
    and uniform p-th roots of unity give the same results to finite-size
    noise, for several p. This preempts the objection that the null
    distribution depends on the alphabet.

Outputs results/phase63_pseudovariance_and_alphabet.json.
"""
from __future__ import annotations
import sys, pathlib, json, time
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SEED = 20260609


def u2_mass(F):
    fh = np.fft.fft(F, axis=1) / F.shape[1]
    return np.sum(np.abs(fh) ** 4, axis=1)


def pseudovar_ratio(F):
    """<|E[Z_xi^2]| / E|Z_xi|^2> averaged over nonzero modes."""
    fh = np.fft.fft(F, axis=1) / F.shape[1]
    sig2 = np.mean(np.abs(fh) ** 2, axis=0)
    rho = np.mean(fh ** 2, axis=0)
    return float(np.mean(np.abs(rho[1:]) / sig2[1:]))


def draw_phases(rng, L, nd, alphabet):
    if alphabet == "S1":
        return np.exp(1j * rng.uniform(0, 2 * np.pi, (nd, L)))
    p = int(alphabet)  # p-th roots of unity
    d = rng.integers(0, p, (nd, L))
    return np.exp(2j * np.pi * d / p)


def apply_conj_reflection(F):
    F = F.copy(); L = F.shape[1]; j = np.arange(L); k = L - 1 - j; m = j < k
    F[:, k[m]] = np.conj(F[:, j[m]]); return F


def apply_noconj_reflection(F, beta=1.0):
    F = F.copy(); L = F.shape[1]; j = np.arange(L); k = L - 1 - j; m = j < k
    F[:, k[m]] = beta * F[:, j[m]]; return F


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    out = {"params": {"seed": SEED, "note": "per-mode pseudovariance and finite-alphabet invariance"},
           "conjugation_test": [], "alphabet_invariance": []}

    print("=== (1) conjugation is essential; mechanism is per-mode impropriety ===")
    L, nd = 400, 20000
    base = draw_phases(rng, L, nd, "S1")
    for name, F in [("random", base),
                    ("conj_reflection", apply_conj_reflection(base)),
                    ("nonconj_reflection", apply_noconj_reflection(base))]:
        M = u2_mass(F); ratio = pseudovar_ratio(F)
        rec = {"case": name, "mass_mean_times_L": round(float(M.mean()) * L, 4),
               "pseudovar_ratio_mean": round(ratio, 4)}
        out["conjugation_test"].append(rec)
        print(f"  {name:18s} mass*L={rec['mass_mean_times_L']:.3f}  "
              f"<|rho|/sigma^2>={rec['pseudovar_ratio_mean']:.3f}")

    print("\n=== (2) alphabet invariance: S^1 vs uniform p-th roots of unity ===")
    L, nd = 300, 12000
    for alphabet in ["S1", "67", "211", "823"]:
        base = draw_phases(rng, L, nd, alphabet)
        Mr = u2_mass(base)
        Mc = u2_mass(apply_conj_reflection(base))
        rec = {"alphabet": alphabet,
               "rand_mass_mean_times_L": round(float(Mr.mean()) * L, 4),
               "refl_mass_mean_times_L": round(float(Mc.mean()) * L, 4),
               "inflation_ratio": round(float(Mc.mean() / Mr.mean()), 4),
               "refl_pseudovar_ratio": round(pseudovar_ratio(apply_conj_reflection(base)), 4)}
        out["alphabet_invariance"].append(rec)
        print(f"  alphabet={alphabet:>3}  rand*L={rec['rand_mass_mean_times_L']:.3f}  "
              f"refl*L={rec['refl_mass_mean_times_L']:.3f}  ratio={rec['inflation_ratio']:.4f}  "
              f"<|rho|/s^2>={rec['refl_pseudovar_ratio']:.3f}")

    out["elapsed_s"] = round(time.time() - t0, 1)
    (ROOT / "results").mkdir(exist_ok=True)
    path = ROOT / "results" / "phase63_pseudovariance_and_alphabet.json"
    path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {path}")


if __name__ == "__main__":
    main()
