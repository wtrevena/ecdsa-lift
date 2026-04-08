"""
Elliptic curve arithmetic over Z/p^k Z, in pure Python.

Deliberately minimal: we need just enough to run lifting experiments
on small primes. Anything that needs real p-adic fields, formal
groups to arbitrary precision, canonical lifts, or isogeny graphs
belongs in the SageMath-backed modules under src/sage/.

Conventions:
- Curves are short Weierstrass: y^2 = x^3 + a*x + b  (mod N).
- The point at infinity is represented as None.
- N is allowed to be a prime p OR a prime power p^k. All inversions
  go through egcd and raise ZeroDivisionError if an element is not
  a unit (caller's problem to decide what that means).
"""

from __future__ import annotations
from dataclasses import dataclass
from math import gcd
from typing import Optional, Tuple

Point = Optional[Tuple[int, int]]  # None == identity (point at infinity)


def egcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return (a, 1, 0)
    g, x1, y1 = egcd(b, a % b)
    return (g, y1, x1 - (a // b) * y1)


def modinv(a: int, n: int) -> int:
    a %= n
    g, x, _ = egcd(a, n)
    if g != 1:
        raise ZeroDivisionError(f"{a} is not a unit mod {n} (gcd = {g})")
    return x % n


@dataclass(frozen=True)
class Curve:
    a: int
    b: int
    N: int  # p or p^k

    def is_on(self, P: Point) -> bool:
        if P is None:
            return True
        x, y = P
        return (y * y - x * x * x - self.a * x - self.b) % self.N == 0

    def neg(self, P: Point) -> Point:
        if P is None:
            return None
        x, y = P
        return (x % self.N, (-y) % self.N)

    def add(self, P: Point, Q: Point) -> Point:
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        N = self.N
        if (x1 - x2) % N == 0:
            if (y1 + y2) % N == 0:
                return None
            # doubling
            m = ((3 * x1 * x1 + self.a) * modinv(2 * y1, N)) % N
        else:
            m = ((y2 - y1) * modinv((x2 - x1) % N, N)) % N
        x3 = (m * m - x1 - x2) % N
        y3 = (m * (x1 - x3) - y1) % N
        return (x3, y3)

    def mul(self, k: int, P: Point) -> Point:
        if P is None or k % self._order_hint() == 0 if False else False:
            return None  # (no-op; the real reduction happens in the loop)
        if k < 0:
            return self.mul(-k, self.neg(P))
        R: Point = None
        Q: Point = P
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.add(Q, Q)
            k >>= 1
        return R

    def _order_hint(self) -> int:
        # placeholder so `mul` is well-defined without needing a known order
        return 1 << 62


def naive_order(E: Curve) -> int:
    """
    Brute-force group order for small primes only. For prime N only —
    counting |E(F_p)|. Uses the Legendre symbol formulation.
    Good for p up to a few thousand.
    """
    p = E.N
    count = 1  # the point at infinity
    for x in range(p):
        rhs = (x * x * x + E.a * x + E.b) % p
        if rhs == 0:
            count += 1
        elif pow(rhs, (p - 1) // 2, p) == 1:
            count += 2
    return count


def points(E: Curve) -> list[Point]:
    """Enumerate every F_p-point. Small primes only."""
    p = E.N
    pts: list[Point] = [None]
    for x in range(p):
        rhs = (x * x * x + E.a * x + E.b) % p
        if rhs == 0:
            pts.append((x, 0))
            continue
        if pow(rhs, (p - 1) // 2, p) != 1:
            continue
        # tonelli-shanks would be proper; for small p we brute force sqrt
        for y in range(p):
            if (y * y) % p == rhs:
                pts.append((x, y))
                pts.append((x, (-y) % p))
                break
    return pts


def point_order(E: Curve, P: Point, bound: int) -> int:
    """Order of P in E(F_p). `bound` should be >= the group order."""
    if P is None:
        return 1
    Q: Point = P
    for k in range(1, bound + 1):
        if Q is None:
            return k
        Q = E.add(Q, P)
    return bound + 1  # unreachable if bound is correct


def find_generator(E: Curve) -> tuple[Point, int]:
    """Point of maximum order in E(F_p). Small primes only."""
    n = naive_order(E)
    best: Point = None
    best_ord = 0
    for P in points(E):
        if P is None:
            continue
        k = point_order(E, P, n)
        if k > best_ord:
            best, best_ord = P, k
    return best, best_ord
