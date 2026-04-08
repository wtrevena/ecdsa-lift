"""
Phase 1 — Reproduce: the Teichmüller-lift error δ(k) is injective at
sufficient p-adic precision.

For a fixed small curve with j = 0 (to stay aligned with secp256k1's
CM structure) and generator G of prime order n, compute

    δ(k) = [k]G̃ - τ([k]G)   in E(Z/p^prec Z)

for every k in [1, n-1] and check whether the map k → δ(k) is
injective. Sweep precision = p, p^2, p^3, p^4, p^5.

Outputs JSON to results/phase1_<p>.json.
"""
from __future__ import annotations
import json
import sys
import pathlib
from collections import defaultdict

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator
from src.lifts import teichmuller_lift


def delta_fingerprint(k_tilde_G, tau_kG) -> tuple:
    """
    Both points reduce to the same point of E(F_p), so their x-
    coordinates agree mod p. We fingerprint their difference as a
    pair of raw coordinate differences mod p^e — this is the natural
    coordinate on the formal neighborhood of the common reduction.
    """
    if k_tilde_G is None and tau_kG is None:
        return ("0", "0")
    if k_tilde_G is None or tau_kG is None:
        return ("inf", None)
    return (k_tilde_G[0] - tau_kG[0], k_tilde_G[1] - tau_kG[1])


def run_prime(p: int, a: int, b: int, max_prec_power: int = 5) -> dict:
    assert a == 0, "we want j = 0 curves to match secp256k1"
    E_p = Curve(a=a, b=b, N=p)
    n = naive_order(E_p)
    G, g_ord = find_generator(E_p)
    report: dict = {
        "prime": p,
        "curve": {"a": a, "b": b},
        "order": n,
        "generator": list(G) if G else None,
        "generator_order": g_ord,
        "precisions": {},
    }
    if G is None:
        return {**report, "error": "no generator found"}

    for e in range(1, max_prec_power + 1):
        N = p ** e
        E_lifted = Curve(a=a, b=b, N=N)
        tau_G = teichmuller_lift(E_p, E_lifted, G)
        seen: dict[tuple, list[int]] = defaultdict(list)
        for k in range(1, g_ord):
            kG = E_p.mul(k, G)
            tau_kG = teichmuller_lift(E_p, E_lifted, kG)
            k_tilde_G = E_lifted.mul(k, tau_G)
            key = delta_fingerprint(k_tilde_G, tau_kG)
            seen[key].append(k)
        collisions = {k: v for k, v in seen.items() if len(v) > 1}
        report["precisions"][f"p^{e}"] = {
            "N": N,
            "distinct_delta_values": len(seen),
            "total_scalars": g_ord - 1,
            "injective": len(collisions) == 0,
            "collision_count": sum(len(v) - 1 for v in collisions.values()),
        }
    return report


def main() -> None:
    # j = 0 curves over primes where -3 is a square (so we have CM).
    # p ≡ 1 mod 3 is the condition.
    # j = 0 curves over primes p ≡ 1 (mod 3). We skip very small primes
    # because they give pathological group structures (e.g. (Z/3)^2) that
    # don't exercise the phenomenon we care about.
    candidates = [
        (13, 0, 3),
        (19, 0, 1),
        (31, 0, 5),
        (37, 0, 2),
        (43, 0, 3),
        (61, 0, 1),
        (67, 0, 2),
        (79, 0, 3),
        (97, 0, 5),
        (103, 0, 2),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, a, b in candidates:
        print(f"[phase1] running p={p}  y^2 = x^3 + {a}x + {b}")
        report = run_prime(p, a, b)
        (out_dir / f"phase1_p{p}.json").write_text(json.dumps(report, indent=2))
        for tag, info in report.get("precisions", {}).items():
            mark = "OK" if info["injective"] else "!!"
            print(f"   {mark} {tag}: injective={info['injective']} "
                  f"({info['distinct_delta_values']}/{info['total_scalars']})")


if __name__ == "__main__":
    main()
