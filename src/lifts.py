"""
Lifting points from E(F_p) up to E(Z/p^k Z).

We implement two sections:

1. `teichmuller_lift`      — set-theoretic lift via Hensel's lemma.
                              Matches the "Teichmüller" construction:
                              lift x coordinate by repeated p-th powering
                              (the usual limit definition of τ), then
                              Hensel-lift y.
2. `averaged_lift`         — the homomorphic section candidate from
                              Phase 4. Exists in principle whenever
                              gcd(|E(F_p)|, p) = 1 and uses the
                              average-over-group trick to turn the
                              set-theoretic τ into an additive section.

The averaged lift is the primary object this repo exists to study.
"""

from __future__ import annotations
from .curves import Curve, Point, points as enumerate_points


def hensel_lift_y(E_lifted: Curve, x: int, y0: int) -> int:
    """
    Given x mod p^k and y0 mod p with y0^2 = f(x) mod p, return y mod p^k
    with y^2 = f(x) mod p^k. Requires y0 != 0 mod p (non-2-torsion).
    """
    N = E_lifted.N
    p = _extract_prime(N)
    y = y0 % p
    prec = p
    while prec < N:
        prec = min(prec * p, N)
        f = (x * x * x + E_lifted.a * x + E_lifted.b) % prec
        from .curves import modinv
        inv2y = modinv(2 * y, prec)
        y = (y - (y * y - f) * inv2y) % prec
    return y


def hensel_lift_x_root(E_lifted: Curve, x0: int) -> int:
    """
    For a 2-torsion point (y = 0), we need x such that f(x) = x^3+ax+b = 0
    in Z/p^k. Lift the F_p-root x0 via Newton on f. Requires f'(x0) to be
    a unit mod p, i.e. 3x0^2 + a != 0 (mod p), which holds iff (x0, 0)
    is a smooth point of the curve.
    """
    from .curves import modinv
    N = E_lifted.N
    p = _extract_prime(N)
    a, b = E_lifted.a, E_lifted.b
    x = x0 % p
    prec = p
    while prec < N:
        prec = min(prec * p, N)
        f = (x * x * x + a * x + b) % prec
        fp = (3 * x * x + a) % prec
        x = (x - f * modinv(fp, prec)) % prec
    return x


def teichmuller_lift(
    E_p: Curve, E_lifted: Curve, P: Point
) -> Point:
    """
    Teichmüller-style set-theoretic lift of a point in E(F_p) to
    E(Z/p^k Z). Not a group homomorphism in general. Phase 2 showed
    this specific section produces a pseudorandom error δ(k).
    """
    if P is None:
        return None
    x0, y0 = P
    p = E_p.N
    N = E_lifted.N
    if y0 == 0:
        # 2-torsion: Teichmüller of x coordinate doesn't in general lie on
        # the lifted curve. We must lift x as a root of f(x) = 0 instead.
        x = hensel_lift_x_root(E_lifted, x0)
        return (x, 0)
    # Generic case: Teichmüller lift of x (root of X^p - X), Hensel-lift y.
    x = x0 % N
    for _ in range(64):
        x_new = pow(x, p, N)
        if x_new == x:
            break
        x = x_new
    y = hensel_lift_y(E_lifted, x, y0)
    return (x, y)


def splitting_section_proj(E_p: Curve, E_lifted_proj, P: Point, p_power_e: int):
    """
    Phase 4 — homomorphic section from the cohomology-vanishing proof,
    computed in projective coordinates.

    Setup.  We have the exact sequence of abelian groups
        0 → A → B → G → 0
    where G = E(F_p), B = E(Z/p^e Z), A = Ê(p Z/p^e) (the kernel of
    reduction). Writing n := |G|, we assume gcd(n, p) = 1, so n is a
    unit on A (because |A| = p^(e-1) is coprime to n).

    The map φ : B → A defined by
        φ(b) := μ( [n] · b )
    with μ = (mult-by-n)^(-1) on A is a retraction: for any b ∈ A, we
    have [n]·b ∈ A and μ([n]·b) = b, so φ|_A = id_A. Consequently
        s(g) := τ(g) - φ(τ(g))
    is a group homomorphism G → B and a section of reduction.

    Concretely: μ is multiplication by n^(-1) mod p^(e-1).

        s(g)  =  τ(g)  -  (n^(-1) mod p^(e-1)) · ([n] · τ(g))

    This is what we compute. The subtraction `τ(g) - ...` is
    division-free in projective coordinates, so it doesn't care that
    both summands reduce to the same F_p point.
    """
    from .curves import modinv
    if P is None:
        return E_lifted_proj.identity()
    E_lifted_aff = Curve(a=E_lifted_proj.a, b=E_lifted_proj.b,
                         N=E_lifted_proj.N)
    tau_P_aff = teichmuller_lift(E_p, E_lifted_aff, P)
    tau_P = E_lifted_proj.from_affine(tau_P_aff)
    n = len(enumerate_points(E_p))
    # [n] · τ(P) lies in the kernel of reduction because n·P = 0 in E(F_p).
    n_tauP = E_lifted_proj.mul(n, tau_P)
    # On the kernel (order p^(e-1)), multiplication by n is invertible;
    # its inverse is multiplication by n^(-1) mod p^(e-1).
    p = E_p.N
    exp_A = max(1, p ** (p_power_e - 1))
    inv_n = modinv(n % exp_A, exp_A) if exp_A > 1 else 0
    phi_P = E_lifted_proj.mul(inv_n, n_tauP)
    # s(P) = τ(P) - φ(τ(P))
    return E_lifted_proj.add(tau_P, E_lifted_proj.neg(phi_P))


def _extract_prime(N: int) -> int:
    """Given N = p^k with p prime, recover p (its smallest prime factor)."""
    i = 2
    while i * i <= N:
        if N % i == 0:
            return i
        i += 1
    if N > 1:
        return N  # N already prime (k = 1)
    raise ValueError(f"Could not factor {N} to extract the prime")
