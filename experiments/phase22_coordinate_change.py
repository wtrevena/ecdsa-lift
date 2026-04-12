"""
Phase 22 — Coordinate change: does rewriting E in a different model
reveal structure the Weierstrass form hides?

The Teichmüller lift τ is a model-dependent construction: it lifts
x via x^p = x, then lifts y via Hensel. A point that is "Teichmüller"
in Weierstrass form is NOT generally Teichmüller in (say) Jacobi
quartic or Montgomery form.

Concretely: for j = 0 the curve y² = x³ + b can be converted to a
Jacobi quartic y² = 1 + d x⁴ or similar via a rational change of
coordinates. The "Teichmüller error δ" depends on WHICH form you
use to define τ. Different forms give different δ, and they cannot
all be pseudorandom — at most one privileged model would have δ
matching the formal logarithm's coordinates linearly.

We test a pragmatic variant: compute δ in the Weierstrass model
AND in a twisted Weierstrass model y² = x³ + b' with b' related
to b by multiplication by a small unit u. This is a trivial
automorphism of the curve (isomorphism via (x, y) ↦ (u²x, u³y))
but τ is not equivariant under it. Measure how δ changes.

If δ₁(k) and δ₂(k) under two different models are related by a
non-trivial affine map in k, we have structure one model couldn't
see.
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


def twist(curve, b_new, g_ord, G_original):
    """Given E: y² = x³ + b and an isomorphic E': y² = x³ + b_new
    where b_new = u^6 · b for some unit u, the isomorphism is
    (x, y) ↦ (u² x, u³ y). Here we assume b_new chosen so the curve
    over F_p has a compatible u."""
    # For simplicity, pick u = 1 (identity) or u = -1 (which gives
    # (x, -y), i.e. the negation map — same curve, same δ).
    # Real coordinate change: use b' = u^6 · b for small u.
    p = curve.N
    for u in range(2, p):
        if (pow(u, 6, p) * curve.b) % p == b_new % p:
            return u
    return None


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}

    # Find another b' in the same isomorphism class: b' = u^6 · b (mod p)
    b_alts = []
    for u in range(2, min(p, 8)):
        b2 = (pow(u, 6, p) * b) % p
        if b2 != b and b2 != 0:
            b_alts.append((u, b2))
    if not b_alts:
        return {"prime": p, "status": "no isomorphism class variants"}

    u, b2 = b_alts[0]
    # Map G to the twisted curve: G' = (u² · G.x, u³ · G.y) mod p
    Gx, Gy = G
    G2 = ((u * u * Gx) % p, (u * u * u * Gy) % p)
    E_p2 = Curve(a=0, b=b2, N=p)
    # Verify G2 is on E_p2
    on_curve = (Gy * Gy * u**6) % p == (Gx**3 * u**6 + b * u**6) % p  # trivially true
    assert (G2[1] * G2[1]) % p == (pow(G2[0], 3, p) + b2) % p

    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    C2 = ProjCurve(a=0, b=b2, N=N)

    d1 = build_delta(E_p, C, G, g_ord)
    d2 = build_delta(E_p2, C2, G2, g_ord)

    # Relationship tests
    eq_hits = 0
    opp_hits = 0
    lin_hits_u2 = 0  # d2 == u² · d1
    lin_hits_u3 = 0  # d2 == u³ · d1
    for k in range(1, g_ord):
        if d1[k] == d2[k]:
            eq_hits += 1
        if (d1[k] + d2[k]) % N == 0:
            opp_hits += 1
        if d2[k] == (pow(u, 2, N) * d1[k]) % N:
            lin_hits_u2 += 1
        if d2[k] == (pow(u, 3, N) * d1[k]) % N:
            lin_hits_u3 += 1

    # Test more general: is d2[k] = α · d1[k] for a fixed α?
    alpha_ctr = {}
    for k in range(1, g_ord):
        if d1[k] == 0:
            continue
        try:
            a = (d2[k] * modinv(d1[k], N)) % N
            alpha_ctr[a] = alpha_ctr.get(a, 0) + 1
        except ZeroDivisionError:
            pass
    dominant_alpha = max(alpha_ctr.items(), key=lambda x: x[1]) if alpha_ctr else (None, 0)

    return {
        "prime": p, "curve_b": b, "curve_b_twist": b2, "isomorphism_u": u,
        "g_ord": g_ord, "N": N,
        "d_equal_hits": f"{eq_hits}/{g_ord - 1}",
        "d_opposite_hits": f"{opp_hits}/{g_ord - 1}",
        "d_scaled_u2_hits": f"{lin_hits_u2}/{g_ord - 1}",
        "d_scaled_u3_hits": f"{lin_hits_u3}/{g_ord - 1}",
        "distinct_alphas": len(alpha_ctr),
        "dominant_alpha_fraction": f"{dominant_alpha[1]}/{g_ord - 1}",
        "dominant_alpha_value": dominant_alpha[0],
        "verdict": ("SCALAR RELATION" if dominant_alpha[1] > (g_ord - 1) * 0.9
                    else "no global scalar; models differ pseudorandomly"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase22] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase22_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
