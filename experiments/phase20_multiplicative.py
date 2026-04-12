"""
Phase 20 — Multiplicative structure probes on δ(k).

δ(jk) vs δ(j)·δ(k), δ(j)+δ(k), δ(j)⊕δ(k). Also test
δ(j)·δ(k) − δ(jk), etc., for constancy-in-jk, linearity-in-(j,k).

Motivation: "pseudorandom under additive probes" doesn't rule out
multiplicative structure. If δ(jk) = δ(j)·δ(k)·f(something), we
have closed-form factorization.
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

    mul_hits = 0
    add_hits = 0
    xor_hits = 0
    total = 0
    mul_over_hits = 0
    rel_residues = []
    for j in range(1, min(g_ord, 15)):
        for k in range(1, min(g_ord, 15)):
            jk = (j * k) % g_ord
            if jk == 0:
                continue
            total += 1
            if (delta[j] * delta[k]) % N == delta[jk]:
                mul_hits += 1
            if (delta[j] + delta[k]) % N == delta[jk]:
                add_hits += 1
            if (delta[j] ^ delta[k]) == delta[jk]:
                xor_hits += 1
            # Test: is delta(jk) / (delta(j)*delta(k)) constant?
            try:
                q = (delta[jk] * modinv(delta[j] * delta[k] % N, N)) % N
                rel_residues.append(q)
            except ZeroDivisionError:
                pass

    # Does the "quotient" q have low entropy / take few values?
    from collections import Counter
    ctr = Counter(rel_residues) if rel_residues else Counter()
    most_common = ctr.most_common(5)

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord,
        "mul_hits": f"{mul_hits}/{total}",
        "add_hits": f"{add_hits}/{total}",
        "xor_hits": f"{xor_hits}/{total}",
        "distinct_quotients": len(ctr),
        "top5_quotients": [(v, c) for v, c in most_common],
        "verdict": ("MULTIPLICATIVE STRUCTURE" if mul_hits > total * 0.1
                    else "no multiplicative structure"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase20] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase20_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
