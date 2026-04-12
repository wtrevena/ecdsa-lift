"""
Phase 5 — CM linearization via the GLV endomorphism.

For j = 0 curves y^2 = x^3 + b over F_p with p ≡ 1 (mod 3), there is a
non-trivial endomorphism φ(x, y) = (ζ·x, y) where ζ is a primitive
cube root of unity in F_p. On E(F_p) it acts as multiplication by a
scalar λ ∈ Z/n satisfying λ² + λ + 1 ≡ 0 (mod n). This is what the
GLV method exploits defensively.

We ask:  can the same endomorphism, applied to the Teichmüller lift
(NOT the homomorphic section), produce structure that the homomorphic
section erased?

Specifically define the "endomorphism-Teichmüller error":
    ψ(P)  :=  [λ_Z] · τ(P)  -  τ(λ · P)    in E(Z/p^e Z)

where λ_Z is any integer lift of λ. By the same cocycle argument as
Phase 2, ψ(P) lies in Ê(p Z/p^e). The question is whether ψ carries
DLP information — i.e. whether ψ(k·G) is an efficiently computable,
efficiently invertible function of k.

Tests performed:
  (a) Confirm ψ(P) ∈ Ê (reduces to identity mod p).
  (b) Injectivity of k ↦ ψ([k]G) at precision p^e.
  (c) Additivity:  ψ([j+k]G)  ?=  ψ([j]G) + ψ([k]G)   in Ê.
  (d) CM equivariance: [λ_Z] · ψ(P)  ?=  ψ(λ·P).
  (e) If additive, try to recover a random k from ψ([k]G).

If (c) holds we have a linear map Z/n → Ê. That's the cleanest possible
structure — it would make the DLP trivial — so we should expect it to
fail, and the experiment is mainly to quantify HOW it fails.
"""
from __future__ import annotations
import json
import sys
import pathlib
import random

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, points, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def primitive_cube_root(p: int) -> int | None:
    """A primitive cube root of unity mod p, when p ≡ 1 (mod 3)."""
    if (p - 1) % 3 != 0:
        return None
    for g in range(2, p):
        r = pow(g, (p - 1) // 3, p)
        if r != 1 and pow(r, 3, p) == 1:
            return r
    return None


def glv_eigenvalue(E_p: Curve, G, g_ord: int, zeta: int) -> int | None:
    """Find λ ∈ Z/g_ord such that φ(G) = [λ] G, where
    φ(x, y) = (ζ x, y). Brute force for small g_ord."""
    x, y = G
    phi_G = (zeta * x % E_p.N, y % E_p.N)
    Q = phi_G
    # We want the k with [k]G = phi_G; brute force up to g_ord.
    cur = G
    for k in range(1, g_ord + 1):
        if cur == Q:
            return k
        cur = E_p.add(cur, G)
    return None


def psi(E_p: Curve, C_proj: ProjCurve, P, lam_Z: int):
    """ψ(P) := [λ_Z] · τ(P) - τ([λ] · P) in E(Z/p^e Z), projective."""
    E_aff = Curve(a=C_proj.a, b=C_proj.b, N=C_proj.N)
    tau_P = C_proj.from_affine(teichmuller_lift(E_p, E_aff, P))
    lam_P = E_p.mul(lam_Z, P)
    tau_lamP = C_proj.from_affine(teichmuller_lift(E_p, E_aff, lam_P))
    return C_proj.add(C_proj.mul(lam_Z, tau_P), C_proj.neg(tau_lamP))


def in_formal_group(C: ProjCurve, P) -> bool:
    """Is P ∈ Ê(pZ/p^e)?  Equivalent to: P reduces to (0:1:0) mod p."""
    p = None
    N = C.N
    # We don't know p here; infer from N via the smallest prime factor.
    for q in range(2, 1000):
        if N % q == 0:
            p = q
            break
    assert p is not None
    X, Y, Z = P
    return X % p == 0 and Z % p == 0 and Y % p != 0


def proj_to_formal_z(C: ProjCurve, P):
    """Formal group parameter z = -X/Y of a point in Ê. Requires Y
    to be a unit mod N."""
    X, Y, Z = P
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous — skip"}
    G, g_ord = find_generator(E_p)
    if G is None or g_ord < 5:
        return {"prime": p, "status": "tiny group"}

    zeta = primitive_cube_root(p)
    if zeta is None:
        return {"prime": p, "status": "no cube root of unity"}

    lam = glv_eigenvalue(E_p, G, g_ord, zeta)
    if lam is None:
        return {"prime": p, "status": "no GLV eigenvalue"}

    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Compute ψ([k]G) for every k.
    psi_table = {}
    for k in range(0, g_ord):
        kG = E_p.mul(k, G) if k > 0 else None
        if kG is None:
            psi_table[k] = C.identity()
        else:
            psi_table[k] = psi(E_p, C, kG, lam)

    # (a) every ψ([k]G) must be in Ê
    in_A = sum(1 for k in range(1, g_ord) if in_formal_group(C, psi_table[k]))

    # (b) injectivity of k → ψ([k]G)
    seen = {}
    for k in range(0, g_ord):
        key = psi_table[k]
        # use formal-group z parameter as a cleaner fingerprint
        z = proj_to_formal_z(C, key)
        seen.setdefault(z, []).append(k)
    injective = max(len(v) for v in seen.values()) == 1

    # (c) additivity in Ê: ψ(jG) + ψ(kG) ?= ψ((j+k)G)
    add_hits = add_total = 0
    first_add_fail = None
    for j in range(1, min(g_ord, 25)):
        for k in range(1, min(g_ord, 25)):
            add_total += 1
            lhs = psi_table[(j + k) % g_ord]
            rhs = C.add(psi_table[j], psi_table[k])
            if C.equal(lhs, rhs):
                add_hits += 1
            elif first_add_fail is None:
                first_add_fail = (j, k)
    additive = add_hits == add_total

    # (d) CM equivariance: [λ]·ψ(P) ?= ψ(λP)
    cm_hits = cm_total = 0
    for k in range(1, min(g_ord, 40)):
        cm_total += 1
        lhs = C.mul(lam, psi_table[k])
        rhs = psi_table[(lam * k) % g_ord]
        if C.equal(lhs, rhs):
            cm_hits += 1

    # (e) attempt scalar recovery if linear
    attack_ok = None
    if additive:
        random.seed(0)
        k_secret = random.randrange(2, g_ord)
        target = psi_table[k_secret]
        # If ψ is a group homomorphism Z/g_ord → Ê, and ψ(1G) generates
        # a cyclic subgroup of Ê, we can find k by taking the formal
        # log ratio in Ê (Ê ≅ Z/p^(e-1) as abelian group).
        z_target = proj_to_formal_z(C, target)
        z_unit = proj_to_formal_z(C, psi_table[1])
        if z_unit and z_unit % p == 0 and z_target is not None:
            # Ê's formal group parameter z is a uniformiser; division
            # z_target / z_unit in the formal group gives k mod p^(e-1)
            # up to higher-order corrections. For small p and small e
            # this is exact to leading order.
            try:
                k_guess = (z_target * modinv(z_unit // p, p ** (e - 1))) % (p ** (e - 1))
                # Reduce mod g_ord
                attack_ok = (k_guess % g_ord) == k_secret
            except Exception:
                attack_ok = False

    return {
        "prime": p, "curve_b": b, "n": n, "g_ord": g_ord,
        "zeta": zeta, "lambda": lam,
        "in_formal_group_count": f"{in_A}/{g_ord - 1}",
        "injective": injective,
        "additive": f"{add_hits}/{add_total}",
        "additive_first_fail": first_add_fail,
        "cm_equivariant": f"{cm_hits}/{cm_total}",
        "attack_ok": attack_ok,
    }


def main():
    candidates = [
        (13, 2), (19, 2), (31, 3), (43, 7),
        (67, 2), (73, 13), (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    for p, b in candidates:
        print(f"[phase5] p={p} b={b}")
        report = run(p, b)
        (out_dir / f"phase5_p{p}.json").write_text(json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
