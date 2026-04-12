"""
Phase 18 — Berlekamp–Massey linear complexity of δ(k).

If the sequence δ(1), δ(2), ..., δ(n−1) satisfies a short linear
recurrence mod p (or mod p²), it has a closed form that leaks the
scalar. Berlekamp–Massey computes the minimum-length recurrence.

Expected for a truly pseudorandom sequence of length L over F_p:
linear complexity ≈ L/2. Anything significantly smaller is a red
flag; anything ≤ O(log L) is a break.

We test the top p-adic digit (δ(k) div p^{e−1} mod p), the middle
digit, and the bottom digit of δ(k).
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


def bm_over_fp(s, p):
    """Berlekamp–Massey over F_p. Returns the shortest recurrence
    length L such that s_j = Σ_{i=1..L} c_i · s_{j−i} for j ≥ L."""
    n = len(s)
    C = [1] + [0] * n
    B = [1] + [0] * n
    L = 0
    m = 1
    b = 1
    for k in range(n):
        d = s[k] % p
        for i in range(1, L + 1):
            d = (d + C[i] * s[k - i]) % p
        if d == 0:
            m += 1
        elif 2 * L <= k:
            T = C[:]
            try:
                coef = (d * pow(b, -1, p)) % p
            except ValueError:
                return None  # b not invertible
            for i in range(n + 1 - m):
                if i + m <= n:
                    C[i + m] = (C[i + m] - coef * B[i]) % p
            L = k + 1 - L
            B = T
            b = d
            m = 1
        else:
            try:
                coef = (d * pow(b, -1, p)) % p
            except ValueError:
                return None
            for i in range(n + 1 - m):
                if i + m <= n:
                    C[i + m] = (C[i + m] - coef * B[i]) % p
            m += 1
    return L


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


def run(p: int, b: int, e: int = 4) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 12:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    deltas = build_delta(E_p, C, G, g_ord)

    L = len(deltas)
    baseline_pseudorandom = L // 2

    # Digit extraction: δ(k) has e digits in base p
    digits = []
    for d in deltas:
        ds = []
        x = d
        for _ in range(e):
            ds.append(x % p)
            x //= p
        digits.append(ds)

    result = {"prime": p, "curve_b": b, "g_ord": g_ord,
              "sequence_length": L, "baseline": baseline_pseudorandom}
    for digit_idx in range(e):
        seq = [d[digit_idx] for d in digits]
        cplx = bm_over_fp(seq, p)
        result[f"lin_complexity_digit_{digit_idx}"] = cplx
        result[f"ratio_digit_{digit_idx}"] = (
            None if cplx is None else f"{cplx / baseline_pseudorandom:.3f}")

    # Also try sum of digits, δ mod p·p (two-digit window)
    wseq = [(d[0] + p * d[1]) % (p * p) for d in digits]
    # For window of two digits we work mod p (dominant digit) as a sanity pass
    result["lin_complexity_window01_modp"] = bm_over_fp([x % p for x in wseq], p)

    result["verdict"] = ("LOW linear complexity — structural leak"
                         if any(v is not None and v < baseline_pseudorandom * 0.8
                                for v in [result.get(f"lin_complexity_digit_{i}")
                                          for i in range(e)])
                         else "linear complexity near L/2, consistent with pseudorandom")
    return result


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase18] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase18_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
