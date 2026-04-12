"""
Phase 16 — Frobenius canonical-subgroup error.

On E(F_p) the Frobenius Frob_p(x, y) = (x^p, y^p) is the identity
(since everything is in F_p). On the lift E(Z/p^e) it is NOT
generally the identity — the map F: (X, Y) ↦ (X^p, Y^p) gives a
point which IS on the curve but differs from (X, Y) by a formal
group element:
    ε(P)  :=  F(τ(P))  −  τ(P)       in Ê.

For supersingular / anomalous curves, ε has a linear relationship
to [p]. For ordinary curves, ε = 0 mod p but nonzero mod p^2 and
encodes the Hasse–Witt invariant.

Question: is ε(kG) a linear-in-k function of k mod p^{e-1}? If so,
it's extractable from a single sample via formal log.
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


def frob_lift(C: ProjCurve, P, p: int):
    """Naive lift of Frobenius: (X:Y:Z) → (X^p : Y^p : Z^p)."""
    X, Y, Z = P
    N = C.N
    return (pow(X, p, N), pow(Y, p, N), pow(Z, p, N))


def z_coord(C, P):
    X, Y, Z = P
    if Y % C.N == 0:
        return None
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


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
    E_aff = Curve(a=0, b=b, N=N)

    epsilons = {}
    for k in range(1, g_ord):
        kG = E_p.mul(k, G)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        # F(tau) is not generally on E(Z/p^e); it's on the curve
        # y^2 = x^3 + b^p after Frobenius. For b ∈ F_p we have b^p = b
        # so the curve is preserved. Verify.
        F_tau = frob_lift(C, tau_kG, p)
        # ε = F(τ) − τ
        eps = C.add(F_tau, C.neg(tau_kG))
        z = z_coord(C, eps)
        epsilons[k] = z if z is not None else 0

    # Check linearity: ε(kG) ≡ k · ε(G) mod p^e
    linear_hits = 0
    for k in range(1, g_ord):
        if (epsilons[k] - k * epsilons[1]) % N == 0:
            linear_hits += 1

    # Check if all ε = 0 (trivial Frobenius on lift)
    zero_count = sum(1 for v in epsilons.values() if v == 0)

    # Berlekamp-Massey on ε(k) mod p
    eps_top = [(epsilons[k] // (p ** (e - 1))) % p for k in range(1, g_ord)]
    # Simple linear test: is eps_top a linear function of k?
    lin_k_hits = 0
    if epsilons[1] != 0:
        try:
            r = epsilons[1]
            for k in range(1, g_ord):
                if (epsilons[k] - k * r) % p == 0:
                    lin_k_hits += 1
        except Exception:
            pass

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "epsilon_zero_count": f"{zero_count}/{g_ord - 1}",
        "epsilon_linear_in_k_mod_N": f"{linear_hits}/{g_ord - 1}",
        "epsilon_linear_in_k_mod_p": f"{lin_k_hits}/{g_ord - 1}",
        "epsilon_sample_first_5": [epsilons[k] for k in range(1, min(g_ord, 6))],
        "verdict": ("ε identically zero" if zero_count == g_ord - 1
                    else ("ε linear — ATTACK POSSIBLE" if linear_hits > (g_ord - 1) * 0.9
                          else "ε nontrivial but non-linear")),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase16] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase16_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
