"""
Projective-coordinate arithmetic on short Weierstrass curves over Z/N.

Purpose: allow group operations on lifted points in E(Z/p^k Z) without
requiring coordinate differences to be units. This is essential for any
experiment that manipulates multiple lifts of the same F_p-point (the
averaging construction, the formal group, translations-by-τ(P), etc.).

Representation:  P = (X, Y, Z)  with the curve given by
    Y^2 Z = X^3 + a X Z^2 + b Z^3.

The identity is (0, 1, 0). Addition and doubling formulas below are
division-free — they use only +, -, * in Z/N — so they don't care
whether intermediate quantities are units.

References: Bernstein–Lange EFD / Cohen "A Course in Computational
Algebraic Number Theory" §7.1. Formulas are for general a, b.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Tuple

ProjPoint = Tuple[int, int, int]  # (X, Y, Z)


@dataclass(frozen=True)
class ProjCurve:
    a: int
    b: int
    N: int

    def identity(self) -> ProjPoint:
        return (0, 1, 0)

    def is_identity(self, P: ProjPoint) -> bool:
        """(0:1:0) projectively — so both X and Z must vanish."""
        return P[0] % self.N == 0 and P[2] % self.N == 0

    def from_affine(self, P) -> ProjPoint:
        if P is None:
            return self.identity()
        return (P[0] % self.N, P[1] % self.N, 1)

    def to_affine(self, P: ProjPoint):
        """Project to affine; returns None for the identity, or the pair
        (x, y) with Z-denominator cleared modulo N. Raises if Z is not
        a unit — caller has to decide what to do in that case."""
        from .curves import modinv
        X, Y, Z = P
        if Z % self.N == 0:
            return None
        inv = modinv(Z % self.N, self.N)
        return ((X * inv) % self.N, (Y * inv) % self.N)

    def neg(self, P: ProjPoint) -> ProjPoint:
        X, Y, Z = P
        return (X % self.N, (-Y) % self.N, Z % self.N)

    def equal(self, P: ProjPoint, Q: ProjPoint) -> bool:
        """Equality of projective points over Z/N.

        Over a ring with zero divisors (Z/p^e), two projective reps
        are equal iff ALL THREE 2×2 minors of the matrix [[X1,Y1,Z1],
        [X2,Y2,Z2]] vanish. Using only two minors (the naive check
        over a field) is wrong here: if both points have Z = 0 the
        minors involving Z vanish trivially regardless of X, Y, so
        the naive check reports every point-at-infinity as equal to
        every other.
        """
        N = self.N
        X1, Y1, Z1 = P
        X2, Y2, Z2 = Q
        return ((X1 * Y2 - X2 * Y1) % N == 0
                and (X1 * Z2 - X2 * Z1) % N == 0
                and (Y1 * Z2 - Y2 * Z1) % N == 0)

    def double(self, P: ProjPoint) -> ProjPoint:
        """Doubling on Y^2 Z = X^3 + aXZ^2 + bZ^3, division-free."""
        N = self.N
        if self.is_identity(P):
            return self.identity()
        X, Y, Z = P
        W = (self.a * Z * Z + 3 * X * X) % N
        S = (Y * Z) % N
        B = (X * Y * S) % N
        H = (W * W - 8 * B) % N
        X3 = (2 * H * S) % N
        Y3 = (W * (4 * B - H) - 8 * Y * Y * S * S) % N
        Z3 = (8 * S * S * S) % N
        return (X3, Y3, Z3)

    def add(self, P: ProjPoint, Q: ProjPoint) -> ProjPoint:
        """General projective addition, division-free. Handles P==Q by
        dispatching to `double` first (necessary — the generic formula
        collapses when the inputs are equal)."""
        N = self.N
        if self.is_identity(P):
            return Q
        if self.is_identity(Q):
            return P
        if self.equal(P, Q):
            return self.double(P)
        X1, Y1, Z1 = P
        X2, Y2, Z2 = Q
        U1 = (Y2 * Z1) % N
        U2 = (Y1 * Z2) % N
        V1 = (X2 * Z1) % N
        V2 = (X1 * Z2) % N
        if (V1 - V2) % N == 0:
            # same x projectively but different y → P = -Q → identity
            return self.identity()
        U = (U1 - U2) % N
        V = (V1 - V2) % N
        V2sq = (V * V) % N
        V3 = (V2sq * V) % N
        W = (Z1 * Z2) % N
        A = (U * U * W - V3 - 2 * V2sq * V2) % N
        X3 = (V * A) % N
        Y3 = (U * (V2sq * V2 - A) - V3 * U2) % N
        Z3 = (V3 * W) % N
        return (X3, Y3, Z3)

    def mul(self, k: int, P: ProjPoint) -> ProjPoint:
        if k < 0:
            return self.mul(-k, self.neg(P))
        R = self.identity()
        Q = P
        while k:
            if k & 1:
                R = self.add(R, Q)
            Q = self.double(Q)
            k >>= 1
        return R
