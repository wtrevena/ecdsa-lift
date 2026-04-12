"""
Phase 12 — GLV / CRT split in the lift.

For j = 0 curves with p ≡ 1 mod 3 and n = |E(F_p)| prime, the group
E(F_p) has CM by Z[ω]. The scalar ring Z/n sometimes splits under
the action of ω: if n is not inert in Z[ω] (which happens roughly
half the time), then
    Z/n ≅ Z[ω]/π × Z[ω]/π̄
for the two prime factors π, π̄ of n in Z[ω]. This is the GLV
decomposition: any k ∈ Z/n decomposes as (k₁, k₂) with k_i ≈ √n.

The question Phase 12 asks: does the canonical lift distinguish k₁
from k₂? Concretely, is there a projection π : Ê → Z/p^{e-1} such
that π(δ(k)) depends only on k₁ (or only on k₂)?

If so, we have a half-dimension reduction: the DLP on E(F_p)
becomes a joint pair of DLPs on smaller groups, each of size √n.
Combined with Pollard-rho that's a ~√n speed-up (the usual GLV
offensive — still not a break).

But the real prize would be if the lift gave a LINEAR extractor
for k₁ (or k₂), not just a split. Then we'd get k₁ in O(1) time,
reducing the total problem to DLP on a √n group — which is
~2^128 for secp256k1, not a break, but a huge advance.
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


def primitive_cube_root(p):
    if (p - 1) % 3 != 0:
        return None
    for g in range(2, p):
        r = pow(g, (p - 1) // 3, p)
        if r != 1:
            return r
    return None


def glv_eigenvalue(n):
    """Find λ mod n with λ² + λ + 1 ≡ 0 mod n. Requires n to split
    in Z[ω]: exists iff -3 is a square mod n."""
    # Try to find sqrt(-3) mod n; then λ = (−1 + √−3)/2.
    # For small n we brute-force.
    for s in range(2, n):
        if (s * s + 3) % n == 0:
            # s² ≡ -3
            try:
                lam = ((-1 + s) * pow(2, -1, n)) % n
                if (lam * lam + lam + 1) % n == 0:
                    return lam
            except ValueError:
                pass
    return None


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}

    lam = glv_eigenvalue(g_ord)
    if lam is None:
        return {"prime": p, "status": "n is inert in Z[ω] — no CRT split"}

    # Split Z/g_ord: the idempotents for the decomposition
    # Z/g_ord ≅ Z/π × Z/π̄ come from λ. Specifically:
    # e₁ = (λ + 1) · 3^{-1} mod g_ord (if 3 is a unit)
    # e₂ = 1 − e₁
    # and the decomposition is k ↦ (k · e₁ mod g_ord, k · e₂ mod g_ord).
    # Actually this isn't quite the CRT split — the right idempotents
    # come from factoring g_ord in Z[ω]. Simpler: use the GLV
    # decomposition k → (k₁, k₂) with k = k₁ + k₂ · λ, |k_i| ≤ √g_ord.

    # Run GLV decomposition via lattice reduction on [1, λ] vs [0, g_ord]
    # (the standard balanced-representation algorithm)
    # For small g_ord brute-force works.
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    delta = build_delta(E_p, C, G, g_ord)

    # Decompose each k as k = k1 + k2·λ with small k1, k2
    # via rounding: k1 = k − k2·λ, k2 = round(k · v / g_ord) for the
    # right v. For brevity: brute force over small k2.
    import math
    bound = int(math.isqrt(g_ord)) + 1

    decomp = {}
    for k in range(1, g_ord):
        best = None
        for k2 in range(-bound, bound + 1):
            k1 = (k - k2 * lam) % g_ord
            if k1 > g_ord // 2:
                k1 -= g_ord
            if abs(k1) <= bound:
                if best is None or abs(k1) + abs(k2) < abs(best[0]) + abs(best[1]):
                    best = (k1, k2)
        decomp[k] = best

    # Bias check: does δ(k) depend linearly on k1, on k2?
    # Fit δ(k) = A·k1 + B·k2 mod N using first two samples, verify on rest
    ks = [k for k in range(1, g_ord) if decomp[k] is not None]
    if len(ks) < 3:
        return {"prime": p, "status": "decomposition failed"}
    k_a, k_b = ks[0], ks[1]
    ka1, ka2 = decomp[k_a]
    kb1, kb2 = decomp[k_b]
    det = (ka1 * kb2 - ka2 * kb1) % N
    try:
        det_inv = modinv(det % N, N)
        A = ((kb2 * delta[k_a] - ka2 * delta[k_b]) * det_inv) % N
        B = ((ka1 * delta[k_b] - kb1 * delta[k_a]) * det_inv) % N
        good = sum(1 for k in ks
                   if (A * decomp[k][0] + B * decomp[k][1] - delta[k]) % N == 0)
    except ZeroDivisionError:
        good = None

    # Also test: δ(k) linear in just k1 (not k2)?
    try:
        A1 = (delta[k_a] * modinv(ka1 if ka1 != 0 else 1, N)) % N
        good1 = sum(1 for k in ks if decomp[k][0] != 0
                    and (A1 * decomp[k][0] - delta[k]) % N == 0)
    except ZeroDivisionError:
        good1 = 0

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "lambda": lam,
        "bound_sqrt_n": bound,
        "N": N,
        "decomposition_sample": {k: decomp[k] for k in list(decomp.keys())[:6]},
        "joint_linear_fit_hits": (f"{good}/{len(ks)}" if good is not None else "fit failed"),
        "k1_only_fit_hits": f"{good1}/{len(ks)}",
        "verdict": ("GLV LINEAR LEAK" if (good is not None and good > len(ks) * 0.9)
                    else "no linear structure in GLV components"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase12] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase12_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
