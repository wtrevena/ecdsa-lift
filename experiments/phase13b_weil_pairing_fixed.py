"""
Phase 13b — Weil pairing Miller's algorithm with auxiliary divisors.

The vanilla Miller calls f_{n,P}(Q) where Q is a single point. If Q
lands on any intermediate kP during doubling, we hit a divisor
collision and the denominator vanishes. The fix: use an auxiliary
point R and compute f_{n,P}(Q + R)/f_{n,P}(R).

Here we:
  (1) Work on the lifted curve E(Z/p^e).
  (2) Compute f_{n,τ(G)}((kG+R)_lift) / f_{n,τ(G)}(R_lift) for a
      chosen R unrelated to G, kG.
  (3) Check if the result (which IS well-defined in (Z/p^e)*)
      depends linearly on k in any p-adic digit.
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv, points
from src.lifts import teichmuller_lift


def line(P1, P2, Q, N):
    x1, y1 = P1
    x2, y2 = P2
    xq, yq = Q
    if (x1 - x2) % N == 0:
        if (y1 + y2) % N == 0:
            return (xq - x1) % N
        num = (3 * x1 * x1) % N
        den = (2 * y1) % N
    else:
        num = (y2 - y1) % N
        den = (x2 - x1) % N
    s = (num * modinv(den, N)) % N
    return (yq - y1 - s * (xq - x1)) % N


def vert(P, Q, N):
    return (Q[0] - P[0]) % N


def miller_f(P, Q, n, E_lifted, N):
    bits = bin(n)[2:]
    T = P
    f_num = 1
    f_den = 1
    for bit in bits[1:]:
        lnum = line(T, T, Q, N)
        T_new = E_lifted.add(T, T)
        if T_new is None:
            return None
        vden = vert(T_new, Q, N)
        f_num = (f_num * f_num * lnum) % N
        f_den = (f_den * f_den * vden) % N
        T = T_new
        if bit == '1':
            lnum = line(T, P, Q, N)
            T_new = E_lifted.add(T, P)
            if T_new is None:
                return None
            vden = vert(T_new, Q, N)
            f_num = (f_num * lnum) % N
            f_den = (f_den * vden) % N
            T = T_new
    # Check denominator is a unit
    try:
        return (f_num * modinv(f_den, N)) % N
    except ZeroDivisionError:
        return None


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n_ord = naive_order(E_p)
    if n_ord % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}
    N = p ** e
    E_lifted = Curve(a=0, b=b, N=N)

    # Pick an auxiliary point R not in <G> span (or just a second
    # random point on E(F_p)). For small p we pick a high-index kR.
    pts = points(E_p)
    R = None
    for pt in pts:
        if pt is None or pt == G or pt == E_p.neg(G):
            continue
        R = pt
        break
    if R is None:
        return {"prime": p, "status": "no auxiliary point"}

    tau_G = teichmuller_lift(E_p, E_lifted, G)
    tau_R = teichmuller_lift(E_p, E_lifted, R)

    miller_vals = {}
    failures = 0
    for k in range(2, min(g_ord, 25)):
        kG = E_p.mul(k, G)
        try:
            tau_kG = teichmuller_lift(E_p, E_lifted, kG)
            # Q' = kG + R
            Q = E_lifted.add(tau_kG, tau_R)
            num = miller_f(tau_G, Q, g_ord, E_lifted, N)
            den = miller_f(tau_G, tau_R, g_ord, E_lifted, N)
            if num is None or den is None:
                failures += 1
                continue
            val = (num * modinv(den, N)) % N
            miller_vals[k] = val
        except (ValueError, ZeroDivisionError):
            failures += 1

    if len(miller_vals) < 4:
        return {"prime": p, "status": f"miller failed {failures}"}

    # Test multiplicativity
    ks = sorted(miller_vals.keys())
    mul_hits = mul_total = 0
    for j in ks:
        for k in ks:
            if j + k in miller_vals:
                mul_total += 1
                if (miller_vals[j] * miller_vals[k]) % N == miller_vals[j + k]:
                    mul_hits += 1

    # Test exponentiation: miller(k) == miller(base)^(k/base) ?
    base = ks[0]
    try:
        inv_base = pow(base, -1, g_ord)
    except ValueError:
        inv_base = None
    exp_hits = 0
    if inv_base is not None:
        for k in ks:
            target = pow(miller_vals[base], (k * inv_base) % g_ord, N)
            if target == miller_vals[k]:
                exp_hits += 1

    # Look for linear-in-k signal in the p-adic digits of miller_vals
    # Express each value modulo p, p², p³ and test sequence linearity
    top_digit_linear = 0
    if len(ks) >= 3:
        top = [miller_vals[k] % p for k in ks]
        if top[1] != top[0]:
            delta_td = (top[1] - top[0]) % p
            base_td = (top[0] - delta_td) % p
            # position: k=ks[i], value = base_td + delta_td * ks[i]? No:
            # Simple linear fit on index
            for i in range(len(ks)):
                pred = (top[0] + (top[1] - top[0]) * (ks[i] - ks[0])) % p
                if pred == top[i]:
                    top_digit_linear += 1

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "miller_failures": failures,
        "miller_computed": len(miller_vals),
        "multiplicativity": f"{mul_hits}/{mul_total}",
        "exponentiation_hits": f"{exp_hits}/{len(miller_vals)}",
        "top_digit_linear_in_k": f"{top_digit_linear}/{len(ks)}",
        "sample_values": [(k, miller_vals[k]) for k in ks[:4]],
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase13b] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase13b_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
