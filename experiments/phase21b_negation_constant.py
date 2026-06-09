"""
Phase 21b — Verify δ(k) + δ(n-k) = [n]·τ(G) (a fixed constant).

If the relation holds it's a free 2-to-1 reduction but does NOT
invert δ; combined with Phase 2 pseudorandomness it still doesn't
break anything. Testing anyway because any exact algebraic relation
is worth pinning down.
"""
from __future__ import annotations
import sys, pathlib, json
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def run(p, b, e=3):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)

    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    n_tau_G = C.mul(g_ord, tau_G)
    # z-coord of [n]τ(G)
    X, Y, Z = n_tau_G
    try:
        const_z = (-X * modinv(Y % N, N)) % N
    except ZeroDivisionError:
        const_z = None

    deltas = {}
    for k in range(1, g_ord):
        kG = E_p.mul(k, G)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
        X, Y, Z = diff
        try:
            deltas[k] = (-X * modinv(Y % N, N)) % N
        except ZeroDivisionError:
            deltas[k] = None

    # Check: delta(k) + delta(g_ord - k) == const_z for all k?
    hits = 0
    total = 0
    sums = set()
    for k in range(1, g_ord):
        nk = g_ord - k
        if nk == 0 or nk not in deltas:
            continue
        if deltas[k] is None or deltas[nk] is None:
            continue
        total += 1
        s = (deltas[k] + deltas[nk]) % N
        sums.add(s)
        if const_z is not None and s == const_z:
            hits += 1

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "n_times_tau_G_z": const_z,
        "relation_hits": f"{hits}/{total}",
        "distinct_sums": len(sums),
        "sample_sums": list(sums)[:5],
        "verdict": ("δ(k)+δ(n-k)=const VERIFIED" if total > 0 and hits == total
                    else ("sums take few values" if len(sums) < 5
                          else "no constant relation")),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    for p, b in candidates:
        print(f"[phase21b] p={p} b={b}")
        r = run(p, b)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
