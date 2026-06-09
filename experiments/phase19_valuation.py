"""
Phase 19 — p-adic valuation distribution of δ(k) and its differences.

If δ(k) were truly uniform in Z/p^e, then v_p(δ(k)) should be
distributed as Geom(1/p): P(v_p = 0) = 1 − 1/p, P(v_p = 1) = (1 −
1/p)/p, etc. The same applies to finite differences δ(k+1) − δ(k)
and higher-order differences.

A sharp deviation from this distribution — especially mass at high
valuations — would mean δ has structure "hidden in the low bits"
that bypasses the linear/polynomial tests.

We tabulate v_p distributions and compare to theoretical uniform.
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


def vp(x, p, e):
    if x % (p ** e) == 0:
        return e
    v = 0
    while x % p == 0 and v < e:
        x //= p
        v += 1
    return v


def build_delta(E_p, C, G, g_ord):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    deltas = []
    for k in range(1, g_ord):
        kG = E_p.mul(k, G)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
        X, Y, Z = diff
        try:
            d = (-X * modinv(Y % C.N, C.N)) % C.N
        except ZeroDivisionError:
            d = 0
        deltas.append(d)
    return deltas


def distribution(seq, p, e):
    hist = [0] * (e + 1)
    for x in seq:
        hist[vp(x, p, e)] += 1
    return hist


def expected_geometric(n, p, e):
    """Expected counts for uniform Z/p^e."""
    exp = []
    remaining = n
    for v in range(e):
        mass = remaining * (p - 1) / p
        exp.append(mass)
        remaining -= mass
    exp.append(remaining)  # v = e: exactly 0
    return exp


def chi_square(obs, exp):
    s = 0.0
    for o, e in zip(obs, exp):
        if e > 0:
            s += (o - e) ** 2 / e
    return s


def run(p: int, b: int, e: int = 4) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    deltas = build_delta(E_p, C, G, g_ord)
    L = len(deltas)

    d0 = distribution(deltas, p, e)
    ed0 = expected_geometric(L, p, e)
    chi0 = chi_square(d0, ed0)

    # 1st differences
    d1_seq = [(deltas[i + 1] - deltas[i]) % N for i in range(L - 1)]
    d1 = distribution(d1_seq, p, e)
    ed1 = expected_geometric(L - 1, p, e)
    chi1 = chi_square(d1, ed1)

    # 2nd differences
    d2_seq = [(d1_seq[i + 1] - d1_seq[i]) % N for i in range(L - 2)]
    d2 = distribution(d2_seq, p, e)
    ed2 = expected_geometric(L - 2, p, e)
    chi2 = chi_square(d2, ed2)

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "length": L,
        "v_p_distribution_delta": d0,
        "v_p_expected_delta": [round(x, 2) for x in ed0],
        "chi_square_delta": round(chi0, 3),
        "v_p_distribution_diff1": d1,
        "chi_square_diff1": round(chi1, 3),
        "v_p_distribution_diff2": d2,
        "chi_square_diff2": round(chi2, 3),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase19] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase19_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
