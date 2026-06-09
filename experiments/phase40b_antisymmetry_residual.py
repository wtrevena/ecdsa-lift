"""
Phase 40b — Is the Phase 40 ordering-signal just the Phase 21b antisymmetry?

Phase 40 showed the U^2 anomaly in s(k) = z(delta(k)) is in the *ordering*
of the cocycle increments (shuffling increments -> random level; real order
-> 2..12 sigma above).  Before calling that a NEW structure we must rule
out the one universal identity that already constrains the ordering:

    Phase 21b:  delta(k) + delta(n-k) = [n].tau(G)   (in the formal group)
    => z(delta(k)) + z(delta(n-k)) ~ const (mod N)   to leading formal order.

This reflection about k = n/2 is a real, ordered, shuffle-destroyed
structure -- so it could BE the whole Phase 40 signal.  Test: restrict to
the first half k = 1 .. floor((n-1)/2).  Within a single half the
antisymmetry relates each k only to a point OUTSIDE the window, so it
imposes NO constraint on the half-sequence's internal ordering.  Re-run the
exact Phase 40 increment-shuffle control on the half.

  * half real U^2 >> half increment-shuffle null  => the ordering structure
    is INTRINSIC, not merely the global antisymmetry. Genuinely new.
  * half real U^2 ~ half null                      => Phase 40's signal was
    just Phase 21b re-detected. Reduces to a known universal identity.
"""
from __future__ import annotations
import sys
import pathlib
import json
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.phase37_gowers import (
    z_delta_sequence, gowers_norm_complex, random_baseline,
)
from experiments.phase40_increment_shuffle import (
    digits_to_phase, increment_shuffle_null,
)
from src.curves import Curve, find_generator


def run_prime(p, b, e=4, num_shuffles=300, num_random=300, seed=12345):
    rng = np.random.default_rng(seed + p)
    _, g_ord = find_generator(Curve(a=0, b=b, N=p))
    seq, g_ord2, N = z_delta_sequence(p, b, e, max_k=g_ord - 1)
    seq = [z for z in seq if z is not None]
    n_full = len(seq)
    half = (n_full) // 2
    seq_half = seq[:half]                 # k = 1 .. half  (strictly first half)
    n = len(seq_half)
    if n < 24:
        return {"p": p, "status": "too small", "n_half": n}

    rand_mean, rand_std, _ = random_baseline(n, 2, num_samples=num_random)

    digits = {}
    for j in range(1, e):
        f_real = digits_to_phase(seq_half, p, j)
        u2_real = gowers_norm_complex(f_real, 2)
        null = increment_shuffle_null(seq_half, N, p, j, num_shuffles, rng)
        nm, ns = float(null.mean()), float(null.std())
        digits[f"d{j}"] = {
            "U2_real": round(u2_real, 5),
            "U2_incShuffle_mean": round(nm, 5),
            "U2_incShuffle_std": round(ns, 5),
            "z_real_vs_incShuffle": round((u2_real - nm) / ns, 3) if ns > 0 else None,
            "z_real_vs_random": round((u2_real - rand_mean) / rand_std, 3) if rand_std > 0 else None,
            "z_incShuffle_vs_random": round((nm - rand_mean) / rand_std, 3) if rand_std > 0 else None,
        }
    return {"p": p, "b": b, "n_full": n_full, "n_half": n, "g_ord": g_ord,
            "e": e, "digits": digits}


SECP_FAMILY = [(p, 7) for p in (67, 97, 103, 127, 167, 211, 257, 311, 401, 521)]


def main():
    out_dir = ROOT / "results"
    reports = []
    for p, b in SECP_FAMILY:
        print(f"[phase40b] p={p} b={b}")
        try:
            r = run_prime(p, b, e=4, num_shuffles=300, num_random=300)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": repr(exc)}
        reports.append(r)
        if "digits" in r:
            print(f"   n_half={r['n_half']} (full {r['n_full']})")
            for j in (1, 2, 3):
                d = r["digits"][f"d{j}"]
                print(f"   d{j}: U2_real={d['U2_real']:.4f}  "
                      f"incShuffle={d['U2_incShuffle_mean']:.4f}±{d['U2_incShuffle_std']:.4f}  "
                      f"z(real/incShuf)={d['z_real_vs_incShuffle']:+.2f}  "
                      f"z(incShuf/rand)={d['z_incShuffle_vs_random']:+.2f}")
        else:
            print(f"   {r}")
        print()
    (out_dir / "phase40b_antisymmetry_residual.json").write_text(
        json.dumps(reports, indent=2, default=str))
    print("Wrote results/phase40b_antisymmetry_residual.json")


if __name__ == "__main__":
    main()
