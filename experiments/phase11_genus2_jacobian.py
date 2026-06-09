"""
Phase 11 — Genus-2 Jacobian cover for j = 0 curves.

For j = 0 elliptic curves E: y² = x³ + b, there is a "classical"
covering map from a genus-2 curve:

    C : y² = x⁶ + b    (genus 2, even model)
    φ : C → E, (x, y) ↦ (x², y)

Jac(C) is 2-dimensional and isogenous to E × E^twist. Under the
cover, a point P ∈ E(F_p) pulls back to a divisor on C with
components above x = ±√(P.x) (if P.x is a square in F_p).

We work with the ODD model because it has a single point at
infinity, which makes Mumford arithmetic simpler:

    y² = f(x) = x⁵ + b·x    ← genus 2 curve, degree-5 model

The map φ: (x, y) ↦ (x, y) to E': y² = x³ + b (substitute x = x)
isn't right here. Let me think again.

Actually the simpler construction: any genus-2 curve y² = f(x) with
f degree 5 has a 2-dimensional Jacobian. We pick f with coefficients
related to our j = 0 curve in a way that makes the embedding E →
Jac(C) explicit. We'll pick a concrete degree-5 polynomial
    f(x) = x^5 + b*x   (this has a 3-torsion structure related to
                        the x -> ζ·x automorphism)
and directly measure: for a divisor D on Jac(C), does the "extra"
coordinate (beyond what the base E captures) carry information
about the scalar k that's not in E alone?

Specifically we:
  (1) Find a divisor D of order n on Jac(C)(F_p) for small p.
  (2) Lift D to Jac(C)(Z/p^e) coefficient-wise via Teichmüller.
  (3) Compute k·D for k = 1..n-1.
  (4) Extract the second (non-E) coordinate of Mumford u(x) and
      check whether it has detectable linear structure in k.

If any Mumford coefficient carries linear-in-k signal, the genus-2
cover IS the structural lever we're looking for.
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import modinv
from src import jac2 as J
from src import poly_modN as P


def curve_order_genus2(f_coeffs, p):
    """Count points on y² = f(x) over F_p, add 1 for ∞.
    Only correct for degree-5 (or 6 with one infinity point)."""
    count = 1  # infinity
    for x in range(p):
        fx = 0
        for c in reversed(f_coeffs):
            fx = (fx * x + c) % p
        # y² = fx
        if fx == 0:
            count += 1
        else:
            # is fx a QR?
            if pow(fx, (p - 1) // 2, p) == 1:
                count += 2
    return count


def find_jac_generator(f_coeffs, p, max_order):
    """Find a divisor D on Jac(C)(F_p) of large order via random
    search. Returns (D, ord). Bound on order from the Hasse-Weil
    bound: |#Jac(C)(F_p) - (p+1)²| ≤ 4·p^{3/2}."""
    # Pick random points on C over F_p
    for x in range(1, p):
        fx = 0
        for c in reversed(f_coeffs):
            fx = (fx * x + c) % p
        if fx == 0:
            continue
        # Find sqrt(fx)
        y = None
        for yt in range(1, p):
            if (yt * yt) % p == fx:
                y = yt
                break
        if y is None:
            continue
        D = J.divisor_from_point(x, y, p)
        # Compute order of D: find smallest k with k·D = identity
        order = None
        cur = D
        for k in range(1, max_order + 1):
            if J.is_identity(cur):
                order = k
                break
            try:
                cur = J.compose(cur, D, f_coeffs, p)
            except (ZeroDivisionError, ValueError):
                cur = None
                break
        if order is not None and order > 3:
            return D, order
    return None, None


def lift_mumford(D, N):
    """Coefficient-wise identity lift (not Teichmüller — for
    Jac(C) the natural lift to Z/N is to keep the coefficients
    the same). The Teichmüller lift of a divisor is a lift of the
    underlying points; for Mumford coordinates this gets messy.
    We use the identity lift and accept the ambiguity."""
    u, v = D
    return (list(u), list(v))


def run(p: int, b: int, e: int = 2) -> dict:
    N = p ** e
    f_p = [0, b, 0, 0, 0, 1]  # x^5 + b*x
    order = curve_order_genus2(f_p, p)
    # Find a generator
    D, D_ord = find_jac_generator(f_p, p, min(order, 200))
    if D is None:
        return {"prime": p, "status": "no usable divisor"}

    f_lifted = [c % N for c in f_p]
    D_lifted = lift_mumford(D, N)

    # Compute k·D for k = 1..D_ord-1 in Jac(C)(Z/N)
    results = []
    try:
        cur = D_lifted
        for k in range(1, min(D_ord, 30)):
            if k > 1:
                cur = J.compose(cur, D_lifted, f_lifted, N)
            u, v = cur
            u_pad = u + [0] * (3 - len(u))
            v_pad = v + [0] * (2 - len(v))
            results.append((k, u_pad, v_pad))
    except (ValueError, ZeroDivisionError) as exc:
        return {"prime": p, "status": f"compose failed: {exc}"}

    # Build feature vectors and check for linearity in k
    # Feature: (u0, u1, u2, v0, v1) at each k
    features = [list(r[1]) + list(r[2]) for r in results]
    L = len(features)
    if L < 4:
        return {"prime": p, "status": "not enough samples"}

    # For each feature dim, test if it's linear in k mod N
    linear_hits = []
    for dim in range(5):
        seq = [features[i][dim] for i in range(L)]
        # Fit: seq[k-1] ≈ A*k + B mod N using (k=1, k=2)
        try:
            A = (seq[1] - seq[0]) % N
            B = (seq[0] - A) % N
            good = sum(1 for i in range(L) if (A * (i + 1) + B - seq[i]) % N == 0)
            linear_hits.append((dim, good, L, A, B))
        except Exception:
            linear_hits.append((dim, None, L, None, None))

    return {
        "prime": p, "curve_b": b, "p_over_F_p_order": order,
        "generator_order": D_ord, "N": N,
        "computed_samples": L,
        "per_dim_linear_fit": linear_hits,
        "first_3_mumford": results[:3],
        "verdict": ("LINEAR LEAK in genus-2 coordinates"
                    if any(lh[1] is not None and lh[1] >= L - 1
                           for lh in linear_hits)
                    else "no linear-in-k Mumford structure"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase11] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase11_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
