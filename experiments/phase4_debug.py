"""
Phase 4 diagnostics. Walk through the cohomological construction on a
single small curve and verify each property we think is true.
"""
from __future__ import annotations
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, points, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def z_proj(C, P):
    """Point's projective Z coordinate mod p (to check if in ker of reduction)."""
    return P[2] % C.N


def reduces_to_identity_modp(C, P, p):
    """Check whether P ∈ Ê(pZ_p). In projective: X ≡ 0 and Z ≡ 0 (mod p)."""
    X, Y, Z = P
    return X % p == 0 and Z % p == 0 and Y % p != 0


def main():
    p, a, b, e = 13, 0, 3, 3
    N = p ** e
    E_p = Curve(a=a, b=b, N=p)
    E_lifted_aff = Curve(a=a, b=b, N=N)
    C = ProjCurve(a=a, b=b, N=N)
    n = naive_order(E_p)
    print(f"p={p}  n=|E(F_p)|={n}  prec=p^{e}={N}")

    pts = points(E_p)
    print(f"E(F_p) has {len(pts)} points")
    # pick a non-identity point
    P0 = next(P for P in pts if P is not None)
    print(f"\nUsing P0 = {P0}")

    tau_P0_aff = teichmuller_lift(E_p, E_lifted_aff, P0)
    print(f"τ(P0) affine in Z/p^{e}: {tau_P0_aff}")
    # Check τ(P0) is on the curve
    x, y = tau_P0_aff
    lhs = (y * y) % N
    rhs = (x * x * x + a * x + b) % N
    print(f"  is on lifted curve? y²-rhs = {(lhs - rhs) % N}")

    tau_P0 = C.from_affine(tau_P0_aff)
    print(f"τ(P0) projective: {tau_P0}")

    # Compute [n]·τ(P0)
    nt = C.mul(n, tau_P0)
    print(f"[n]·τ(P0) = {nt}")
    print(f"  in Ê (ker of reduction)? {reduces_to_identity_modp(C, nt, p)}")
    # Check [n]·τ(P0) is on the curve projectively
    X, Y, Z = nt
    on_curve = (Y * Y * Z - X**3 - a * X * Z * Z - b * Z**3) % N
    print(f"  on projective curve? {on_curve == 0}")

    # Affinize nt to see it in (x,y) form if possible
    try:
        aff_nt = C.to_affine(nt)
        print(f"  affine: {aff_nt}")
    except ZeroDivisionError:
        print("  affine: Z not a unit")

    # Now φ(τ(P0)) = multiply nt by n^{-1} mod p^(e-1)
    exp_A = p ** (e - 1)
    inv_n = modinv(n % exp_A, exp_A)
    print(f"\nμ = mult-by-({inv_n}) (= n^(-1) mod {exp_A})")
    phi = C.mul(inv_n, nt)
    print(f"φ(τ(P0)) = {phi}")
    print(f"  in Ê? {reduces_to_identity_modp(C, phi, p)}")

    # Verify the retraction property: apply φ to a known element of A and
    # check it acts as identity.
    # Take a ∈ A: use the difference τ(P0) - τ(P0) = O ... trivial. Instead,
    # get nt = [n]·τ(P0) which we already know is in A, then φ(nt) should = nt.
    phi_nt = C.mul(inv_n, C.mul(n, nt))
    print(f"φ(nt) equal to nt? {C.equal(phi_nt, nt)}")

    # Retraction check on a nontrivial A-element:
    # Build a_elem = [n]·τ(Q0) - [n]·τ(P0) which is still in A (difference
    # of A-elements).
    Q0 = [P for P in pts if P is not None and P != P0][0]
    tau_Q0 = C.from_affine(teichmuller_lift(E_p, E_lifted_aff, Q0))
    nQ = C.mul(n, tau_Q0)
    a_elem = C.add(nQ, C.neg(nt))
    print(f"\na_elem in Ê? {reduces_to_identity_modp(C, a_elem, p)}")
    phi_a = C.mul(inv_n, C.mul(n, a_elem))
    print(f"φ(a_elem) == a_elem? {C.equal(phi_a, a_elem)}")

    # Test homomorphism of s on a pair
    s_P0 = C.add(tau_P0, C.neg(phi))
    # Pick Q0
    tau_Q0_aff = teichmuller_lift(E_p, E_lifted_aff, Q0)
    tau_Q0 = C.from_affine(tau_Q0_aff)
    phi_Q0 = C.mul(inv_n, C.mul(n, tau_Q0))
    s_Q0 = C.add(tau_Q0, C.neg(phi_Q0))

    PQ = E_p.add(P0, Q0)
    tau_PQ_aff = teichmuller_lift(E_p, E_lifted_aff, PQ)
    tau_PQ = C.from_affine(tau_PQ_aff)
    phi_PQ = C.mul(inv_n, C.mul(n, tau_PQ))
    s_PQ = C.add(tau_PQ, C.neg(phi_PQ))

    lhs = s_PQ
    rhs = C.add(s_P0, s_Q0)
    print(f"\ns(P+Q) == s(P) + s(Q)? {C.equal(lhs, rhs)}")
    # Difference is in A?
    diff = C.add(lhs, C.neg(rhs))
    print(f"difference in Ê? {reduces_to_identity_modp(C, diff, p)}")
    print(f"  diff = {diff}")
    # Maybe the difference is 0 in A but the projective rep is not (0,1,0)?
    # Try reducing modulo p and see
    Xd, Yd, Zd = diff
    print(f"  diff reduced mod p: ({Xd % p}, {Yd % p}, {Zd % p})")

    # Exhaustive sweep: for all pairs (P, Q), does s(P+Q) == s(P)+s(Q)?
    print("\nExhaustive sweep:")
    # Precompute s_map
    s_map = {}
    for P in pts:
        if P is None:
            s_map[None] = C.identity()
            continue
        tP = C.from_affine(teichmuller_lift(E_p, E_lifted_aff, P))
        phi = C.mul(inv_n, C.mul(n, tP))
        s_map[P] = C.add(tP, C.neg(phi))

    fails = 0
    first_fail = None
    for P in pts:
        for Q in pts:
            PQ = E_p.add(P, Q)
            sP = s_map.get(P, C.identity()) if P is not None else C.identity()
            sQ = s_map.get(Q, C.identity()) if Q is not None else C.identity()
            sPQ = s_map.get(PQ, C.identity()) if PQ is not None else C.identity()
            rhs = C.add(sP, sQ)
            if not C.equal(sPQ, rhs):
                fails += 1
                if first_fail is None:
                    first_fail = (P, Q, PQ, sPQ, rhs)
    print(f"  fails: {fails}")
    if first_fail:
        P, Q, PQ, sPQ, rhs = first_fail
        print(f"  first fail: P={P}, Q={Q}, P+Q={PQ}")
        print(f"    s(P+Q)   = {sPQ}")
        print(f"    s(P)+s(Q)= {rhs}")
        # Check if s(P+Q) and s(P)+s(Q) are projectively equal but
        # stored with different scalings.
        X1, Y1, Z1 = sPQ
        X2, Y2, Z2 = rhs
        print(f"    X1 Z2 - X2 Z1 mod N = {(X1*Z2 - X2*Z1) % N}")
        print(f"    Y1 Z2 - Y2 Z1 mod N = {(Y1*Z2 - Y2*Z1) % N}")


if __name__ == "__main__":
    main()
