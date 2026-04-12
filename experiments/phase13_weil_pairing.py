"""
Phase 13 — Weil pairing on the lifted curve via Miller's algorithm.

The Weil pairing e_n : E[n] × E[n] → μ_n is bilinear:
    e_n([k]P, Q) = e_n(P, Q)^k.
On F_p this gives the familiar MOV attack only when μ_n sits in a
small extension. For ordinary curves with huge embedding degree, no.

The question: what if we compute the pairing on the LIFT
E(Z/p^e)[n] instead of the base curve? For points of order n in
the lift, the "Weil pairing" e_n lives in (Z/p^e)* / (n-th powers),
which is a larger group than μ_n ⊂ F_p. If the higher p-adic digits
of e_n([k]τ(G), τ(G)) encode k linearly (mod p^{e-1}), that would
be a direct attack.

The subtlety: τ(P) is NOT an n-torsion point in E(Z/p^e) in general
(Phase 2: [n]τ(G) is a nonzero element of Ê). So the "Weil pairing
on τ(G), τ(Q)" isn't literally the n-torsion pairing — it's what
Miller's algorithm outputs when run for n steps on non-n-torsion
input. That output is well-defined in (Z/p^e)* / (n-th powers) only
if we normalize.

We test the simpler object: Miller's function f_{n,P}(Q) computed
for P = τ(G), Q = τ(kG), see if the p-adic expansion of the
resulting value carries linear-in-k signal.
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.lifts import teichmuller_lift


def line_affine(P1, P2, Q, N):
    """Line through P1 and P2 (affine), evaluated at Q. Returns
    a value in Z/N. Assumes the necessary differences are units."""
    x1, y1 = P1
    x2, y2 = P2
    xq, yq = Q
    if (x1 - x2) % N == 0:
        if (y1 + y2) % N == 0:
            # Vertical line at x = x1
            return (xq - x1) % N
        # Tangent: slope = (3 x1²) / (2 y1) for a = 0
        num = (3 * x1 * x1) % N
        den = (2 * y1) % N
        s = (num * modinv(den, N)) % N
    else:
        num = (y2 - y1) % N
        den = (x2 - x1) % N
        s = (num * modinv(den, N)) % N
    # Line: y - y1 = s (x - x1)
    # Value at Q: yq - y1 - s(xq - x1)
    return (yq - y1 - s * (xq - x1)) % N


def vertical_affine(P, Q, N):
    """Vertical line at x = P.x, evaluated at Q."""
    return (Q[0] - P[0]) % N


def miller(P, Q, n, E_lifted, N):
    """Compute f_{n, P}(Q) via Miller's algorithm on the lifted
    curve. Assumes P, Q are affine non-identity points with
    coordinate differences that are units in Z/N at all
    intermediate steps."""
    if n == 1:
        return 1
    # Binary expansion
    bits = bin(n)[2:]
    T = P
    f = 1
    for bit in bits[1:]:
        # Doubling step
        l = line_affine(T, T, Q, N)
        T_new = E_lifted.add(T, T)
        if T_new is None:
            # T_new is identity; can't divide
            raise ValueError("Miller hit identity")
        v = vertical_affine(T_new, Q, N)
        f = (f * f * l * modinv(v, N)) % N
        T = T_new
        if bit == '1':
            l = line_affine(T, P, Q, N)
            T_new = E_lifted.add(T, P)
            if T_new is None:
                raise ValueError("Miller hit identity")
            v = vertical_affine(T_new, Q, N)
            f = (f * l * modinv(v, N)) % N
            T = T_new
    return f


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

    # Lift G and kG via Teichmuller
    tau_G = teichmuller_lift(E_p, E_lifted, G)
    # Pick a "target" Q ≠ G
    # We'll vary kG as "Q" and compute Miller(tau_G, tau_kG, n).
    # Since tau_G is not n-torsion, Miller's output is just f_{n, τG}(τQ).
    # Check if the result varies linearly in k.

    miller_vals = {}
    failures = 0
    for k in range(2, min(g_ord, 30)):
        kG = E_p.mul(k, G)
        try:
            tau_kG = teichmuller_lift(E_p, E_lifted, kG)
            val = miller(tau_G, tau_kG, g_ord, E_lifted, N)
            miller_vals[k] = val
        except (ValueError, ZeroDivisionError):
            failures += 1

    if len(miller_vals) < 5:
        return {"prime": p, "status": f"miller failed ({failures}) — divisor collision"}

    # Test linearity: log-style, f_{n,P}(kG) = f_{n,P}(G)^k in (Z/p^e)*
    # Equivalent to: val(k) / val(k-1) is a constant ratio, OR
    # val(k) == val(1)^k  (if val(1) is known)
    # We computed from k=2 so use val(2) as base
    base_k = min(miller_vals.keys())
    base_val = miller_vals[base_k]
    exp_hits = 0
    for k, v in miller_vals.items():
        expected = pow(base_val, k, N)
        if expected == v:
            exp_hits += 1

    # Multiplicativity: miller_vals[j+k] = miller_vals[j] * miller_vals[k]?
    mul_hits = 0
    mul_total = 0
    ks = sorted(miller_vals.keys())
    for j in ks:
        for k in ks:
            jk = j + k
            if jk not in miller_vals:
                continue
            mul_total += 1
            if (miller_vals[j] * miller_vals[k]) % N == miller_vals[jk]:
                mul_hits += 1

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "miller_failures": failures,
        "miller_computed": len(miller_vals),
        "exponentiation_test_hits": f"{exp_hits}/{len(miller_vals)}",
        "multiplicativity_hits": f"{mul_hits}/{mul_total}",
        "sample_miller_vals": [(k, miller_vals[k]) for k in ks[:5]],
        "verdict": ("LINEAR LEAK IN MILLER OUTPUT"
                    if mul_total > 0 and mul_hits > mul_total * 0.9
                    else "no linearity"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase13] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase13_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
