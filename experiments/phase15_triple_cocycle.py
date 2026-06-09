"""
Phase 15 — Triple cocycle of τ (3-cochain alternating sum).

The Teichmüller lift τ : E(F_p) → E(Z/p^e) is NOT a homomorphism.
Its 2-cocycle is
    c(P, Q) := τ(P) + τ(Q) − τ(P+Q)  ∈  Ê.
Phase 6 looked at c(G, kG) and found it pseudorandom.

Going one level up: the 3-cochain
    γ(P, Q, R) := c(P, Q+R) − c(P+Q, R) + c(P, Q) − c(Q, R)
measures whether c itself is a 2-coboundary. Since H^2 = 0 for our
setting (Phase 3.5), γ must be a 3-coboundary — but the EXPLICIT
formula γ might carry structural information that c did not.

Specifically: γ is trilinear if τ is "quadratic" (i.e., τ's
non-homomorphism error is a bilinear form). We test whether γ is
trilinear or symmetric in (P, Q, R).
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
    """c(P, Q) = τ(P) + τ(Q) − τ(P+Q) in projective (Z/N)."""
    tauP = tau(E_p, C, P)
    tauQ = tau(E_p, C, Q)
    PQ = E_p.add(P, Q)
    tauPQ = tau(E_p, C, PQ)
    return C.add(C.add(tauP, tauQ), C.neg(tauPQ))


def z_coord(C, P):
    X, Y, Z = P
    if Y % C.N == 0:
        return None
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


def gamma(E_p, C, P, Q, R):
    """γ(P,Q,R) = c(P,Q+R) − c(P+Q,R) + c(P,Q) − c(Q,R)."""
    QR = E_p.add(Q, R)
    PQ = E_p.add(P, Q)
    t1 = cocycle_c(E_p, C, P, QR)
    t2 = cocycle_c(E_p, C, PQ, R)
    t3 = cocycle_c(E_p, C, P, Q)
    t4 = cocycle_c(E_p, C, Q, R)
    # γ = t1 − t2 + t3 − t4
    s = C.add(t1, C.neg(t2))
    s = C.add(s, t3)
    s = C.add(s, C.neg(t4))
    return s


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 8:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # 2-cocycle vanishing check (sanity)
    # For a homomorphism τ, c(P,Q) = 0. For us it shouldn't be.
    c_zero = 0
    c_total = 0
    for i in range(1, min(g_ord, 6)):
        for j in range(1, min(g_ord, 6)):
            c_total += 1
            c = cocycle_c(E_p, C, E_p.mul(i, G), E_p.mul(j, G))
            z = z_coord(C, c)
            if z == 0:
                c_zero += 1

    # 3-cochain γ vanishing?
    g_zero = 0
    g_total = 0
    g_values = []
    for i in range(1, min(g_ord, 5)):
        for j in range(1, min(g_ord, 5)):
            for k in range(1, min(g_ord, 5)):
                g_total += 1
                Pi = E_p.mul(i, G)
                Pj = E_p.mul(j, G)
                Pk = E_p.mul(k, G)
                g = gamma(E_p, C, Pi, Pj, Pk)
                z = z_coord(C, g)
                if z == 0:
                    g_zero += 1
                else:
                    g_values.append(((i, j, k), z))

    # Is γ a symmetric trilinear form in (i,j,k)?
    # Test: γ(i,j,k) ≡ i·j·k · γ(1,1,1) mod N ?
    if g_values:
        try:
            γ111 = next(v for (ijk, v) in g_values if ijk == (1, 1, 1))
            tri_hits = 0
            tri_total = 0
            for (i, j, k), v in g_values:
                tri_total += 1
                if (v - i * j * k * γ111) % N == 0:
                    tri_hits += 1
            trilinear = f"{tri_hits}/{tri_total}"
        except StopIteration:
            trilinear = "no γ(1,1,1) available"
    else:
        trilinear = "γ identically zero"

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "c_zero_hits": f"{c_zero}/{c_total}",
        "gamma_zero_hits": f"{g_zero}/{g_total}",
        "trilinear_in_ijk": trilinear,
        "verdict": ("γ ≡ 0 — cocycle trivial"
                    if g_zero == g_total
                    else "γ nonzero — checking structure"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase15] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase15_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
