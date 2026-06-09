"""
Phase 37d — Permutation control for the Phase 37c residual signal.

If the residual U^2 elevation reflects a true *cocycle/index* structure
(δ(k) is approximately polynomial in k), permuting k uniformly at random
should destroy it: the multiset of digit values is the same, but the
order is scrambled.  If the U^2 stays elevated under permutation, then
the elevation is purely an artifact of the empirical *value distribution*
of digit_j(z(δ(k))) — i.e. the digits aren't uniform on F_p — and has
nothing to do with the group/cocycle structure.  Such a signal would not
leak DLP information.
"""
from __future__ import annotations
import sys
import pathlib
import math
import json
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator
from experiments.phase37_gowers import (
    z_delta_sequence, encode_digit, gowers_norm_complex, random_baseline,
)


def run(p, b, e=4, num_perm=100, num_random=200):
    E_p = Curve(a=0, b=b, N=p)
    G, g_ord = find_generator(E_p)
    seq, _, N = z_delta_sequence(p, b, e, max_k=g_ord - 1)
    n = len(seq)
    if n < 8:
        return {"p": p, "status": "too small"}

    rand_mean, rand_std, _ = random_baseline(n, 2, num_samples=num_random)

    results = {}
    for j in range(1, e):
        f = np.array([encode_digit(z, p, j) for z in seq], dtype=complex)
        u2_real = gowers_norm_complex(f, 2)
        # Permute the values uniformly num_perm times
        perm_u2 = []
        for _ in range(num_perm):
            perm = np.random.permutation(n)
            f_perm = f[perm]
            perm_u2.append(gowers_norm_complex(f_perm, 2))
        perm_arr = np.array(perm_u2)
        results[f"d{j}"] = {
            "U2_real": u2_real,
            "U2_perm_mean": float(perm_arr.mean()),
            "U2_perm_std": float(perm_arr.std()),
            "z_real_vs_perm": round(
                (u2_real - perm_arr.mean()) / perm_arr.std(), 3),
            "U2_random_mean": rand_mean,
            "U2_random_std": rand_std,
            "z_real_vs_random": round((u2_real - rand_mean) / rand_std, 3),
            "z_perm_vs_random": round(
                (perm_arr.mean() - rand_mean) / rand_std, 3),
        }
    return {"p": p, "b": b, "n": n, "g_ord": g_ord, "digits": results}


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    SECP = [(p, 7) for p in (31, 67, 127, 167, 211, 257, 311, 401, 521)]
    reports = []
    for p, b in SECP:
        print(f"[phase37d] p={p} b={b}")
        try:
            r = run(p, b, e=4, num_perm=80, num_random=150)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": repr(exc)}
        reports.append(r)
        if "digits" in r:
            print(f"   n={r['n']}")
            for j in (1, 2, 3):
                d = r["digits"][f"d{j}"]
                print(f"   d{j}: real={d['U2_real']:.4f}  perm_mean={d['U2_perm_mean']:.4f}  "
                      f"rand={d['U2_random_mean']:.4f}  "
                      f"z(real/perm)={d['z_real_vs_perm']:+.2f}  "
                      f"z(perm/rand)={d['z_perm_vs_random']:+.2f}")
        print()
    (out_dir / "phase37d_permutation.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
