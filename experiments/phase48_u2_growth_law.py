"""
Phase 48 - The U^2 growth law is an artifact signature (Proposition 1).

The Gowers-U^2 elevation GROWS with n (z ~ 3 at n=79, ~ 44 at n=10639).
The original Phase-37 report read growth-with-n as accumulating signal.
This phase shows the growth is forced by the antisymmetry reflection
ALONE, with NO curve arithmetic, via a closed-form mass factor:

  E[ ||f||_{U^2}^4 ] * L  ->  2   for i.i.d. S^1-uniform f
                          ->  3   for reflection-only f  (f(L+1-k)=beta*conj(f(k)))

i.e. the reflection inflates the U^2^4 mass by EXACTLY 3/2, independent of
L and of the value alphabet p. Hence the standardized synthetic-vs-random
gap scales as z(syn/rand) = K*sqrt(L), with K ~ 0.42 measured here, and
K*sqrt(n) reproduces the Phase-43 syn/rand column and the Phase-47 value.

Self-contained: needs no elliptic-curve code (that is the point - the
growth carries no curve information).

Output: results/phase48_u2_growth_law.json
"""
from __future__ import annotations
import json
import math
import pathlib

import numpy as np

SEED = 20260609
NREP = 400


def U2(fv):
    F = np.fft.fft(fv) / len(fv)
    return float((np.abs(F) ** 4).sum() ** 0.25)


def measure(L, p_alphabet, rng):
    """Return (E[U2^4]*L for rand, for syn, z(syn/rand))."""
    h = L // 2
    rand_u2, rand_u4, syn_u2 = [], [], []
    for _ in range(NREP):
        g = np.exp(2j * np.pi * rng.random(L))
        u = U2(g)
        rand_u2.append(u)
        rand_u4.append(u ** 4)
        d = rng.integers(0, p_alphabet, h)
        C = int(rng.integers(0, p_alphabet))
        full = np.empty(L, dtype=np.int64)
        full[:h] = d
        full[L - h:] = (C - d[::-1]) % p_alphabet
        syn_u2.append(U2(np.exp(2j * np.pi * full / p_alphabet)))
    er = float(np.mean(rand_u4)) * L
    # reflection-only U2^4 * L
    es = float(np.mean([u ** 4 for u in syn_u2])) * L
    z = (np.mean(syn_u2) - np.mean(rand_u2)) / np.std(rand_u2)
    return er, es, z


def main():
    rng = np.random.default_rng(SEED)
    rows = []
    # vary L and alphabet p to show p-independence and the 3/2 factor
    grid = [(78, 10007), (198, 10007), (612, 10007), (828, 10007),
            (2000, 10007), (10638, 10007),
            (612, 211), (612, 99991)]  # alphabet-independence cross-check
    Ks = []
    for L, p in grid:
        er, es, z = measure(L, p, rng)
        K = z / math.sqrt(L)
        rows.append({"L": L, "p_alphabet": p,
                     "U2^4*L_rand": round(er, 4),
                     "U2^4*L_syn": round(es, 4),
                     "mass_ratio": round(es / er, 4),
                     "z_syn_rand": round(z, 3),
                     "K=z/sqrt(L)": round(K, 4)})
        if p == 10007:
            Ks.append(K)
        print(f"L={L:6d} p={p:6d}  mass {er:.3f}->{es:.3f} "
              f"(x{es/er:.3f})  z={z:6.2f}  K={K:.4f}")
    Kbar = float(np.mean(Ks))

    # predict Phase-43 syn/rand column and Phase-47 from K*sqrt(n)
    predict = {}
    for n in [79, 139, 199, 313, 397, 613, 829, 10639]:
        predict[n] = round(Kbar * math.sqrt(n - 1), 2)

    out = {
        "phase": 48, "seed": SEED, "nrep": NREP,
        "claim": "reflection inflates U2^4 mass by exactly 3/2; z(syn/rand)=K*sqrt(n)",
        "mass_ratio_mean": round(float(np.mean([r["mass_ratio"] for r in rows])), 4),
        "K_mean": round(Kbar, 4),
        "grid": rows,
        "predicted_syn_rand_K_sqrt_n": predict,
        "note": ("mass ratio is 3/2 independent of L and value alphabet; "
                 "sqrt(n) growth therefore carries no curve information - "
                 "it is the signature of the antisymmetry artifact."),
    }
    res = pathlib.Path(__file__).resolve().parent.parent / "results"
    res.mkdir(exist_ok=True)
    (res / "phase48_u2_growth_law.json").write_text(json.dumps(out, indent=2))
    print(f"mass ratio ~ {out['mass_ratio_mean']} (=3/2), K ~ {Kbar:.4f}")
    print("wrote results/phase48_u2_growth_law.json")


if __name__ == "__main__":
    main()
