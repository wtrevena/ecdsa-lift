"""
Phase 21 — Symmetry probes: k ↔ n−k  and  k ↔ k^{-1} (mod n).

On E(F_p), [−1]P has the same x-coordinate and negated y. The
Teichmüller lift of [−1]P matches τ(P) with Y negated (because τ
lifts x uniquely via x^p = x and then y via Hensel). So
    τ([−1]P) = [−1]_Z τ(P)    exactly
and consequently δ(n − k) = [n−k]τ(G) − τ((n−k)G). In Ê this
relates to δ(k) by a fixed involution. The question is: what
involution? And does the (k, δ(k)) ↔ (n−k, δ(n−k)) pairing give a
2-to-1 reduction?

Separately: k ↔ k^{-1} (mod n) is the inversion automorphism of
(Z/n)^×. δ(k^{-1}) vs δ(k): any relation?

Phase 5 tested CM λ-action, Phase 2 tested additivity in k. This
phase tests the OTHER natural group actions on Z/n: negation and
inversion.
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def build_delta(E_p, C, G, g_ord):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    deltas = {}
    for k in range(1, g_ord):
        kG = E_p.mul(k, G)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
        X, Y, Z = diff
        try:
            deltas[k] = (-X * modinv(Y % C.N, C.N)) % C.N
        except ZeroDivisionError:
            deltas[k] = 0
    return deltas


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    delta = build_delta(E_p, C, G, g_ord)

    # Negation: δ(n−k) = ?
    neg_plus_hits = 0    # δ(n-k) == -δ(k)
    neg_eq_hits = 0      # δ(n-k) == δ(k)
    neg_zero_hits = 0    # δ(n-k) + δ(k) == 0
    neg_total = 0
    for k in range(1, g_ord):
        nk = (g_ord - k) % g_ord
        if nk == 0:
            continue
        neg_total += 1
        if delta[nk] == delta[k]:
            neg_eq_hits += 1
        if delta[nk] == (-delta[k]) % N:
            neg_plus_hits += 1
        if (delta[nk] + delta[k]) % N == 0:
            neg_zero_hits += 1

    # Inversion k → k^{-1} mod g_ord
    inv_eq_hits = 0
    inv_ratio_hits = 0
    inv_total = 0
    from collections import Counter
    quot_ctr = Counter()
    for k in range(1, g_ord):
        try:
            ki = pow(k, -1, g_ord)
        except ValueError:
            continue
        if ki == 0 or ki >= g_ord:
            continue
        inv_total += 1
        if delta[ki] == delta[k]:
            inv_eq_hits += 1
        try:
            q = (delta[ki] * modinv(delta[k], N)) % N
            quot_ctr[q] += 1
        except ZeroDivisionError:
            pass

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord,
        "negation_eq_hits": f"{neg_eq_hits}/{neg_total}",
        "negation_opp_hits": f"{neg_plus_hits}/{neg_total}",
        "negation_sum_zero": f"{neg_zero_hits}/{neg_total}",
        "inversion_eq_hits": f"{inv_eq_hits}/{inv_total}",
        "inversion_distinct_quotients": len(quot_ctr),
        "inversion_top3_quotients": quot_ctr.most_common(3),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase21] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase21_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
