"""
Phase 51 — Smart's attack, ALL k (implementation validation).

The paper reports "12/12 random secrets" on the anomalous curve
y^2 = x^3 + 4x + 7 over F_53 (G = (0, 22), n = p = 53). Here we run
Smart's p-adic-lift attack recovering k for EVERY k in [1, n-1] (52
secrets) and report the success count. Expected: 52/52.

This is a deterministic, exhaustive validation of the smart_dlog
implementation reused from phase42_anomalous_boundary.
"""
from __future__ import annotations
import sys
import pathlib
import json

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, point_order
from experiments.phase42_anomalous_boundary import smart_dlog


def main():
    p, a, b = 53, 4, 7
    G = (0, 22)
    E = Curve(a, b, p)
    n = naive_order(E)
    assert n == p, f"curve is not anomalous: #E={n}, p={p}"
    g_ord = point_order(E, G, n)
    print(f"Curve y^2 = x^3 + {a}x + {b} over F_{p}, #E = {n} (=p: {n==p}), "
          f"G = {G}, ord(G) = {g_ord}")

    results = []
    failures = []
    for k in range(1, n):
        Q = E.mul(k, G)
        krec = smart_dlog(a, b, p, G, Q)
        ok = (krec == k)
        results.append({"k": k, "recovered": krec, "ok": ok})
        if not ok:
            failures.append({"k": k, "recovered": krec})

    success = sum(r["ok"] for r in results)
    total = len(results)
    print(f"Smart attack recovered {success}/{total} secrets "
          f"(k = 1..{n-1}).")
    if failures:
        print(f"FAILURES: {failures}")
    else:
        print("No failures.")

    out = {
        "label": "implementation_validation",
        "curve": {"p": p, "a": a, "b": b, "order": n, "anomalous": n == p},
        "G": list(G),
        "g_ord": g_ord,
        "k_range": [1, n - 1],
        "num_secrets": total,
        "success_count": success,
        "success_rate": success / total,
        "failures": failures,
        "all_results": results,
    }
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "phase51_smart_allk.json").write_text(
        json.dumps(out, indent=2, default=str))
    print("Wrote results/phase51_smart_allk.json")


if __name__ == "__main__":
    main()
