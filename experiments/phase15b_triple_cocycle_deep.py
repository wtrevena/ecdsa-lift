"""
Phase 15b — Deeper look at the 3-cochain γ of τ.

Phase 15 showed that γ(iG, jG, kG) has z-coordinate = 0 in ~94% of
tested (i,j,k) triples at precision p^3. That could be either:
  (α) γ ≡ identity (i.e. c IS a coboundary modulo projective
      representation noise), OR
  (β) γ lies in DEEPER p-adic filtration — Ê(p^2 Z/p^e) rather
      than Ê(pZ/p^e) — most of the time.

These are very different claims. (α) would say the explicit
homomorphic retraction from Phase 4 is recoverable at the cochain
level without the [n]^{-1} trick. (β) would say γ has p-adic
valuation ≥ 2 for most inputs, revealing a hidden graded structure
on the cocycle.

We test by measuring v_p of the X/Z coordinates of γ directly,
without the "divide by Y" z_coord projection that can mask the
structure.
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


def tau(E_p, C, P):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    return C.from_affine(teichmuller_lift(E_p, E_aff, P))


def cocycle_c(E_p, C, P, Q):
    tauP = tau(E_p, C, P)
    tauQ = tau(E_p, C, Q)
    PQ = E_p.add(P, Q)
    tauPQ = tau(E_p, C, PQ)
    return C.add(C.add(tauP, tauQ), C.neg(tauPQ))


def gamma(E_p, C, P, Q, R):
    QR = E_p.add(Q, R)
    PQ = E_p.add(P, Q)
    t1 = cocycle_c(E_p, C, P, QR)
    t2 = cocycle_c(E_p, C, PQ, R)
    t3 = cocycle_c(E_p, C, P, Q)
    t4 = cocycle_c(E_p, C, Q, R)
    s = C.add(t1, C.neg(t2))
    s = C.add(s, t3)
    s = C.add(s, C.neg(t4))
    return s


def vp(x, p, cap):
    if x == 0:
        return cap
    v = 0
    while x % p == 0 and v < cap:
        x //= p
        v += 1
    return v


def point_depth(P, p, cap):
    """Depth in formal filtration. A point (X:Y:Z) is in Ê(p^k)
    iff v_p(X/Z) ≥ k (informally, after normalization). For an
    Ê-point we have Y a unit, so the relevant quantity is v_p(X)
    when Z is a unit — or v_p(Z) when Z is not a unit."""
    X, Y, Z = P
    X = X % (p ** cap)
    Y = Y % (p ** cap)
    Z = Z % (p ** cap)
    # If Z is a unit, depth = v_p(X)
    if Z % p != 0:
        return vp(X, p, cap)
    # Otherwise normalize: multiply through by inverse of a unit scale
    # Simplest: look at min(vp(X), vp(Z))
    return min(vp(X, p, cap), vp(Z, p, cap))


def is_identity_strict(C, P):
    X, Y, Z = P
    return X % C.N == 0 and Z % C.N == 0


def run(p: int, b: int, e: int = 4) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 8:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    depth_histo_c = {}
    depth_histo_gamma = {}
    identity_c = 0
    identity_gamma = 0
    c_total = 0
    g_total = 0

    for i in range(1, min(g_ord, 6)):
        for j in range(1, min(g_ord, 6)):
            c_total += 1
            Pi = E_p.mul(i, G)
            Pj = E_p.mul(j, G)
            c = cocycle_c(E_p, C, Pi, Pj)
            if is_identity_strict(C, c):
                identity_c += 1
            d = point_depth(c, p, e)
            depth_histo_c[d] = depth_histo_c.get(d, 0) + 1

            for k in range(1, min(g_ord, 5)):
                g_total += 1
                Pk = E_p.mul(k, G)
                g = gamma(E_p, C, Pi, Pj, Pk)
                if is_identity_strict(C, g):
                    identity_gamma += 1
                d = point_depth(g, p, e)
                depth_histo_gamma[d] = depth_histo_gamma.get(d, 0) + 1

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "e": e, "N": N,
        "c_total": c_total,
        "c_identity_count": identity_c,
        "c_depth_histogram": dict(sorted(depth_histo_c.items())),
        "gamma_total": g_total,
        "gamma_identity_count": identity_gamma,
        "gamma_depth_histogram": dict(sorted(depth_histo_gamma.items())),
        "notes": (
            "c should concentrate at depth 1 (Ê(p)). "
            "If γ concentrates at depth ≥ 2, we have graded structure."
        ),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase15b] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase15b_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
