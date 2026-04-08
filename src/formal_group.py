"""
Formal group logarithm for short Weierstrass curves, computed to a
given p-adic precision.

For E: y^2 = x^3 + a*x + b, the formal group is obtained by the change
of variables z = -x/y, w = -1/y. The group law on E linearizes in z
up to an explicit power series:

    log_Ê(z) = z + (higher-order terms in z)

with rational coefficients whose denominators are integers coprime
to most primes of interest. For our experiments we only need the
logarithm on points that live in Ê(pZ_p) — i.e. points whose z
coordinate is divisible by p.

This module is deliberately low-precision and Python-only. For
real precision (beyond p^10 or so) we should switch to SageMath's
built-in formal group machinery.
"""

from __future__ import annotations
from .curves import Curve, Point, modinv


def z_of(P: Point, N: int) -> int | None:
    """
    The formal-group parameter z = -x/y for an affine point.
    Returns None if y is not a unit in Z/N (the point is not in Ê).
    """
    if P is None:
        return 0  # identity of the formal group
    x, y = P
    try:
        return (-x * modinv(y, N)) % N
    except ZeroDivisionError:
        return None


def formal_log(E: Curve, z: int, terms: int = 20) -> int:
    """
    The formal group logarithm of a short-Weierstrass curve, truncated
    to `terms` terms, evaluated at z in Z/N. Only correct on inputs
    z that live in pZ_p (otherwise the series doesn't converge p-adically).

    The series is derived from w(z) = z^3 + a*z^7 + b*z^9 + ...
    and integrating dx/(2y + ... ) but for a first pass we use the
    known truncated form:

        log(z) = z + (a/5) z^5 + (b/7) z^7 + (a^2/9) z^9 + ...

    This is the canonical differential's antiderivative; coefficients
    come from expanding dx/y in terms of z. For our small-p
    experiments the low-order part is enough to detect structure.
    """
    N = E.N
    a, b = E.a, E.b
    result = z % N
    # Leading-order correction: integrate w(z) = z^3 + a z^7 + b z^9 + ...
    # term_{2m+1} depends on curve invariants. We'll use a minimal
    # explicit expansion. See Silverman "Arithmetic of Elliptic Curves"
    # Ch IV §1 for the general recursion.
    # z^5 coefficient: a/5
    # z^7 coefficient: b/7
    # z^9 coefficient: a^2/9
    # z^11 coefficient: (ab)/11
    # z^13 coefficient: (a^3 + 2b^2)/13
    coeffs = {
        5: (a, 5),
        7: (b, 7),
        9: (a * a, 9),
        11: (a * b, 11),
        13: (a * a * a + 2 * b * b, 13),
    }
    z_pow = pow(z, 1, N)
    z_sq = (z * z) % N
    # z_pow tracks z^k; we bump it by z^2 twice to step k by 2.
    current = (z * z_sq) % N  # z^3 (we'll never add this — it's the w start)
    current = (current * z_sq) % N  # z^5
    for k in range(5, 2 * terms, 2):
        if k in coeffs:
            num, den = coeffs[k]
            try:
                inv_den = modinv(den, N)
            except ZeroDivisionError:
                # log not defined mod this precision for this curve
                return result
            result = (result + num * inv_den * current) % N
        current = (current * z_sq) % N
    return result
