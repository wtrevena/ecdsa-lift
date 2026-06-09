"""
Phase 40 — The decisive control for the Phase 37 Fourier / Gowers anomaly.

Phase 37 reported a persistent degree-1 (U^2) Fourier anomaly in the
p-adic digit sequence of the Teichmuller lift error
        s(k) := z(delta(k)) = z([k]tau(G) - tau(kG)),  k = 1 .. g_ord-1.
Phase 37d's permutation control shuffled the *final* digit values and
showed the elevation vanishes -> "the structure is in the index ordering".

But that control cannot distinguish two very different explanations:

  (A) GENUINE cocycle structure: the cocycle is approximately a
      degree-1 polynomial phase in k, a real new property.

  (B) PREFIX-SUM ARTIFACT: s(k) is the running (prefix) sum of its own
      first differences d(k) = s(k+1) - s(k).  Prefix-summing is a
      low-pass filter -- in Fourier space  s_hat(xi) = d_hat(xi)/(1 - e^{2pi i xi/n})
      which blows up like 1/xi near xi = 0.  ANY sequence with the same
      multiset of increments, once prefix-summed, has elevated low-frequency
      mass and hence elevated U^2, regardless of the increments' order.

Phase 37d destroys BOTH the genuine structure and the artifact at once
(shuffling final values kills all ordering), so it cannot tell them apart.

THE DECISIVE CONTROL (this file):
  s is, tautologically, the cumulative sum (mod N) of its first
  differences d.  So:
    * REAL      : encode digit_j(s(k)), measure U^2.
    * NULL (N2) : randomly permute the increments d -> d', rebuild
                  s'(k) = s(1) + sum_{i<k} d'(i)  (mod N), encode
                  digit_j(s'(k)), measure U^2.  Repeat many times.

  N2 has the *identical multiset of increments* as the real sequence,
  only re-ordered, then prefix-summed exactly as the real one is.

Interpretation:
  * real U^2 ~ N2 distribution  (but both >> S^1 random)  => the anomaly
    is FULLY explained by prefix-summing this increment multiset. It is
    NOT a cocycle-specific index structure. Phase 37 is an artifact.
  * real U^2 >> N2 distribution                            => there is
    genuine index structure beyond prefix-summing. Phase 37 is real.

We report, per prime and per digit j in {1,2,3}:
  U2_real, mean/std of U2 over N2 shuffles, z(real vs N2),
  and for context z(real vs S^1 random) (the original Phase 37 number).
"""
from __future__ import annotations
import sys
import pathlib
import json
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from experiments.phase37_gowers import (
    z_delta_sequence, encode_digit, gowers_norm_complex, random_baseline,
)


def digits_to_phase(values, p, j):
    """Encode digit j of each integer in `values` as a phase on S^1."""
    return np.array([encode_digit(int(v), p, j) for v in values], dtype=complex)


def increment_shuffle_null(s, N, p, j, num_shuffles, rng):
    """U^2 distribution for: shuffle first-differences of s, prefix-sum, encode.

    s : list[int]    the real sequence s(k) = z(delta(k)), already mod N.
    Returns array of U^2 values, one per shuffle.
    """
    s = np.asarray(s, dtype=object)
    n = len(s)
    # First differences (mod N). d has length n-1; s is recovered by
    # s(0) + cumulative-sum of d.  Tautology: real prefix sum reproduces s.
    d = [(int(s[i + 1]) - int(s[i])) % N for i in range(n - 1)]
    base = int(s[0])

    out = np.empty(num_shuffles)
    for t in range(num_shuffles):
        d_perm = list(d)
        rng.shuffle(d_perm)
        # Rebuild prefix sums: s'(0)=base, s'(k)=base+sum_{i<k} d_perm[i]  (mod N)
        s_prime = [base]
        acc = base
        for di in d_perm:
            acc = (acc + di) % N
            s_prime.append(acc)
        f = digits_to_phase(s_prime, p, j)
        out[t] = gowers_norm_complex(f, 2)
    return out


def run_prime(p, b, e=4, num_shuffles=300, num_random=300, seed=12345):
    rng = np.random.default_rng(seed + p)
    # z_delta_sequence needs max_k; use g_ord-1 like Phase 37.
    from src.curves import Curve, find_generator
    _, g_ord = find_generator(Curve(a=0, b=b, N=p))
    seq, g_ord2, N = z_delta_sequence(p, b, e, max_k=g_ord - 1)
    assert g_ord == g_ord2
    # Drop any None (non-invertible Y); they break the prefix-sum arithmetic.
    if any(z is None for z in seq):
        seq = [z for z in seq if z is not None]
    n = len(seq)
    if n < 16:
        return {"p": p, "status": "too small", "n": n}

    rand_mean, rand_std, _ = random_baseline(n, 2, num_samples=num_random)

    digits = {}
    for j in range(1, e):
        f_real = digits_to_phase(seq, p, j)
        u2_real = gowers_norm_complex(f_real, 2)
        null = increment_shuffle_null(seq, N, p, j, num_shuffles, rng)
        null_mean, null_std = float(null.mean()), float(null.std())
        digits[f"d{j}"] = {
            "U2_real": round(u2_real, 5),
            "U2_incShuffle_mean": round(null_mean, 5),
            "U2_incShuffle_std": round(null_std, 5),
            "z_real_vs_incShuffle": round((u2_real - null_mean) / null_std, 3)
            if null_std > 0 else None,
            "U2_random_mean": round(rand_mean, 5),
            "U2_random_std": round(rand_std, 5),
            "z_real_vs_random": round((u2_real - rand_mean) / rand_std, 3)
            if rand_std > 0 else None,
            "z_incShuffle_vs_random": round((null_mean - rand_mean) / rand_std, 3)
            if rand_std > 0 else None,
        }
    return {"p": p, "b": b, "n": n, "g_ord": g_ord, "e": e, "digits": digits}


SECP_FAMILY = [(p, 7) for p in
    (31, 43, 67, 79, 97, 103, 127, 167, 211, 257, 311, 401, 521)]


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in SECP_FAMILY:
        print(f"[phase40] p={p} b={b}")
        try:
            r = run_prime(p, b, e=4, num_shuffles=300, num_random=300)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": repr(exc)}
        reports.append(r)
        if "digits" in r:
            print(f"   n={r['n']}  g_ord={r['g_ord']}")
            for j in (1, 2, 3):
                d = r["digits"][f"d{j}"]
                print(f"   d{j}: U2_real={d['U2_real']:.4f}  "
                      f"incShuffle={d['U2_incShuffle_mean']:.4f}"
                      f"±{d['U2_incShuffle_std']:.4f}  "
                      f"z(real/incShuf)={d['z_real_vs_incShuffle']:+.2f}  "
                      f"z(real/rand)={d['z_real_vs_random']:+.2f}  "
                      f"z(incShuf/rand)={d['z_incShuffle_vs_random']:+.2f}")
        else:
            print(f"   {r}")
        print()
    (out_dir / "phase40_increment_shuffle.json").write_text(
        json.dumps(reports, indent=2, default=str))
    print("Wrote results/phase40_increment_shuffle.json")


if __name__ == "__main__":
    main()
