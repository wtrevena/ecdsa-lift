"""
Phase 59 -- U^3 corollary: the reflection's higher-order inflation vanishes.

The paper reports empirically that U^3 and higher "sat at baseline" while
U^2 was elevated. Here we explain it.

SETUP.  For f: Z/L -> S^1,
   ||f||_{U^2}^4 = sum_xi |hat f(xi)|^4,
   ||f||_{U^3}^8 = (1/L) sum_h ||Delta_h f||_{U^2}^4 ,   Delta_h f(x)=f(x+h) conj f(x).
Baselines (i.i.d. proper f):  E||f||_{U^2}^4 = 2/L,  E||f||_{U^3}^8 = O(1/L^2).

CLAIM.  Under the reflection f(L+1-k)=beta conj f(k):
  * U^2 mass is inflated by the constant factor 3/2 (Phase 57/Prop 1):
        E||f||_{U^2}^4 = 3/L,   ratio -> 3/2.
  * U^3 mass is inflated by a factor 1 + O(1/L): the reflection contributes
    only a measure-zero (O(1/L)) fraction of the (h, derivative) phase space,
    so  E||f||_{U^3}^8 = (1+O(1/L)) * baseline,  ratio -> 1.

WHY (derivation sketch).  ||f||_{U^3}^8 averages ||Delta_h f||_{U^2}^4 over all
shifts h. The reflection symmetry of f induces a reflection symmetry on
Delta_h f ONLY for the single self-paired shift h ~ 0 (more precisely the
derivative Delta_h f inherits an exact pairwise constraint on an O(L)-sized
set of the L^2 (h,x) pairs, i.e. an O(1/L) fraction). For all other h the
two reflection partners of f land on x-locations whose Delta_h values are
NOT conjugate-reflected, so Delta_h f is asymptotically proper and its U^2
mass is at the unconstrained 2/L level. Averaging:
    E||f||_{U^3}^8 = (1/L)[ (#structured h)*O(1/L) + (L - O(1)) * baseline ]
                   = baseline * (1 + O(1/L)).
So the U^3 RATIO -> 1, even as the U^2 ratio -> 3/2. The reflection is a
purely SECOND-order (Fourier/quadratic-mass) phenomenon; it does not survive
the extra derivative that defines U^3. Same argument kills U^k, k>=3.

This file numerically confirms:
  (1) U^2 ratio (syn/rand) -> 3/2 as L grows,
  (2) U^3 ratio (syn/rand) -> 1 as L grows (the excess decays like 1/L),
  (3) tabulates both vs L to exhibit the separation.
"""
from __future__ import annotations
import json, math, pathlib
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
SEED = 20260609


def u2_mass(g):
    L = len(g)
    F = np.fft.fft(g) / L
    return float(np.sum(np.abs(F) ** 4))


def u2_mass_pow4_from_array(g):
    """||g||_{U2}^4 via FFT autocorrelation (same value as u2_mass)."""
    return u2_mass(g)


def u3_pow8(g):
    """||g||_{U^3}^8 = (1/L) sum_h ||Delta_h g||_{U^2}^4, Delta_h g(x)=g(x+h)conj g(x)."""
    L = len(g)
    cg = np.conj(g)
    total = 0.0
    for h in range(L):
        dh = np.roll(g, -h) * cg
        total += u2_mass(dh)
    return total / L


def make_reflection(L, p, rng, beta=None):
    h = L // 2
    if beta is None:
        beta = np.exp(2j*np.pi*rng.random())
    d = rng.integers(0, p, h)
    f = np.empty(L, dtype=complex)
    f[:h] = np.exp(2j*np.pi*d/p)
    f[L-h:] = beta * np.conj(np.exp(2j*np.pi*d[::-1]/p))
    if L % 2 == 1:
        f[h] = np.exp(2j*np.pi*rng.random())
    return f


def make_random(L, rng):
    return np.exp(2j*np.pi*rng.random(L))


def measure(L, rng, nrep, p=10007):
    """Return mean U2^4 and U3^8 for random and reflection ensembles."""
    r2 = np.empty(nrep); r3 = np.empty(nrep)
    s2 = np.empty(nrep); s3 = np.empty(nrep)
    for t in range(nrep):
        g = make_random(L, rng)
        r2[t] = u2_mass(g)
        r3[t] = u3_pow8(g)
        f = make_reflection(L, p, rng)
        s2[t] = u2_mass(f)
        s3[t] = u3_pow8(f)
    return {
        "L": L,
        "U2^4_rand*L": round(float(r2.mean())*L, 4),
        "U2^4_syn*L": round(float(s2.mean())*L, 4),
        "U2_ratio": round(float(s2.mean()/r2.mean()), 4),
        "U3^8_rand*L^2": round(float(r3.mean())*L*L, 4),
        "U3^8_syn*L^2": round(float(s3.mean())*L*L, 4),
        "U3_ratio": round(float(s3.mean()/r3.mean()), 4),
        "U3_excess": round(float(s3.mean()/r3.mean()) - 1.0, 4),
    }


def main():
    rng = np.random.default_rng(SEED)
    out = {"phase": 59, "seed": SEED,
           "claim": "U2 ratio -> 3/2 (const), U3 ratio -> 1 (excess ~ 1/L)",
           "rows": []}
    # U3 is O(L^2) per draw, so keep L modest but spread out to see the trend;
    # nrep large enough for stable ratios.
    grid = [(24, 4000), (40, 3000), (64, 2000), (96, 1500),
            (144, 900), (216, 500), (320, 280)]
    print(f"{'L':>5} {'U2ratio':>9} {'U3ratio':>9} {'U3excess':>9} {'excess*L':>9}")
    for L, nrep in grid:
        r = measure(L, rng, nrep)
        r["U3_excess*L"] = round(r["U3_excess"] * L, 3)
        out["rows"].append(r)
        print(f"{L:5d} {r['U2_ratio']:9.4f} {r['U3_ratio']:9.4f} "
              f"{r['U3_excess']:9.4f} {r['U3_excess*L']:9.3f}")

    # Summaries: U2 ratio should hover at 1.5; U3 excess*L roughly constant
    u2r = [r["U2_ratio"] for r in out["rows"]]
    u3x = [r["U3_excess"] for r in out["rows"]]
    out["U2_ratio_mean"] = round(float(np.mean(u2r)), 4)
    out["U3_excess_decays"] = bool(abs(u3x[-1]) < abs(u3x[0]))
    out["U3_excessL_roughly_const"] = [round(r["U3_excess*L"], 3) for r in out["rows"]]
    (ROOT/"results").mkdir(exist_ok=True)
    (ROOT/"results"/"phase59_u3.json").write_text(json.dumps(out, indent=2))
    print(f"\nU2 ratio mean = {out['U2_ratio_mean']} (target 1.5)")
    print(f"U3 excess decays with L: {out['U3_excess_decays']}")
    print("wrote results/phase59_u3.json")


if __name__ == "__main__":
    main()
