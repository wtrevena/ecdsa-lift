"""
Genus-2 hyperelliptic Jacobian arithmetic over Z/N.

Curve: C : y^2 = f(x)  with f of degree 5 or 6.
We use the imaginary model: f monic of degree 5 (odd), with one
point at infinity. Mumford representation: a reduced divisor on
Jac(C) is encoded as a pair (u(x), v(x)) with
  • u monic of degree ≤ 2,
  • deg v < deg u,
  • u divides v^2 − f.

The identity element is (u = 1, v = 0) — i.e. the empty divisor.

Composition and reduction follow Cantor's algorithm. Over Z/p^e
we require that non-unit pivots don't appear during polynomial
GCD — this holds for generic inputs on small primes, and we
raise ZeroDivisionError in degenerate cases.

For our experiment we use y^2 = x^5 + b·x (which is a degree-5
model associated to a j = 0 curve via a coordinate change), so
the Jacobian is isogenous to E × E^d for some twist E^d. Points
are pushed from E(Z/N) into Jac(C)(Z/N) via the cover map.
"""
from __future__ import annotations
from typing import Tuple, List

from . import poly_modN as P
from .curves import modinv


Mumford = Tuple[List[int], List[int]]  # (u, v)


def identity() -> Mumford:
    return ([1], [0])


def is_identity(D: Mumford) -> bool:
    u, v = D
    return P.deg(u) == 0 and u[0] == 1 and all(c == 0 for c in v)


def neg(D: Mumford, N: int) -> Mumford:
    u, v = D
    return (list(u), P.sub([0], v, N))


def reduce(u: List[int], v: List[int], f: List[int], N: int,
           genus: int = 2) -> Mumford:
    """Reduce (u, v) until deg u ≤ genus."""
    while P.deg(u) > genus:
        # u' = (f - v^2) / u
        vsq = P.mul(v, v, N)
        num = P.sub(f, vsq, N)
        u_new, rem = P.divmod_poly(num, u, N)
        if rem != [0]:
            raise ValueError(f"reduction non-exact: rem={rem}")
        # v' = -v mod u'
        neg_v = P.sub([0], v, N)
        _, v_new = P.divmod_poly(neg_v, u_new, N)
        # Make u_new monic
        lead = u_new[-1]
        lead_inv = modinv(lead, N)
        u_new = P.scale(u_new, lead_inv, N)
        u, v = u_new, v_new
    return (u, v)


def compose(D1: Mumford, D2: Mumford, f: List[int], N: int) -> Mumford:
    """Cantor composition. Handles the coprime case simply; otherwise
    dispatches to the general formula."""
    u1, v1 = D1
    u2, v2 = D2
    if is_identity(D1):
        return D2
    if is_identity(D2):
        return D1
    # Try coprime case: gcd(u1, u2) = 1
    try:
        d, a, b = P.gcd_poly(u1, u2, N)
    except ZeroDivisionError:
        raise
    # If d is constant (i.e. gcd = 1), we can do the easy case.
    if P.deg(d) == 0:
        # u3 = u1 * u2
        u3 = P.mul(u1, u2, N)
        # v3 = a*u1*v2 + b*u2*v1 (mod u3)
        t1 = P.mul(P.mul(a, u1, N), v2, N)
        t2 = P.mul(P.mul(b, u2, N), v1, N)
        v3_raw = P.add(t1, t2, N)
        _, v3 = P.divmod_poly(v3_raw, u3, N)
    else:
        # General case: compute gcd(d, v1 + v2)
        v_sum = P.add(v1, v2, N)
        d2, s2, t2poly = P.gcd_poly(d, v_sum, N)
        # d2 = s2 * d + t2poly * (v1 + v2)
        # Pull it together: write d2 = e1*u1 + e2*u2 + e3*(v1+v2)
        e1 = P.mul(s2, a, N)
        e2 = P.mul(s2, b, N)
        e3 = t2poly
        # u3 = u1 * u2 / d2^2
        u1u2 = P.mul(u1, u2, N)
        d2sq = P.mul(d2, d2, N)
        u3, rem = P.divmod_poly(u1u2, d2sq, N)
        if rem != [0]:
            raise ValueError("compose: u1*u2 not divisible by d2^2")
        # v3 = (e1*u1*v2 + e2*u2*v1 + e3*(v1*v2+f)) / d2 (mod u3)
        t1 = P.mul(P.mul(e1, u1, N), v2, N)
        t2a = P.mul(P.mul(e2, u2, N), v1, N)
        t3 = P.mul(e3, P.add(P.mul(v1, v2, N), f, N), N)
        v_raw = P.add(P.add(t1, t2a, N), t3, N)
        v_div, rem = P.divmod_poly(v_raw, d2, N)
        if rem != [0]:
            raise ValueError("compose: v divisor rem != 0")
        _, v3 = P.divmod_poly(v_div, u3, N)
    # Monic u3
    lead = u3[-1]
    lead_inv = modinv(lead, N)
    u3 = P.scale(u3, lead_inv, N)
    return reduce(u3, v3, f, N)


def mul_k(k: int, D: Mumford, f: List[int], N: int) -> Mumford:
    """Scalar multiplication via double-and-add."""
    if k < 0:
        return mul_k(-k, neg(D, N), f, N)
    R = identity()
    Q = D
    while k:
        if k & 1:
            R = compose(R, Q, f, N)
        k >>= 1
        if k:
            Q = compose(Q, Q, f, N)
    return R


def equal(D1: Mumford, D2: Mumford) -> bool:
    u1, v1 = D1
    u2, v2 = D2
    return P.trim(u1) == P.trim(u2) and P.trim(v1) == P.trim(v2)


def divisor_from_point(x: int, y: int, N: int) -> Mumford:
    """(P) - (∞) for a point P = (x, y) on C : y^2 = f(x).
    Mumford form: u = x - x_P, v = y_P (constant)."""
    u = [(-x) % N, 1]
    v = [y % N]
    return (u, v)
