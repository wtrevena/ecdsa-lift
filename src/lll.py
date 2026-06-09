"""
Tiny LLL implementation in pure Python.

Enough to attack small lattices (rank up to maybe 20-40 on the
small-integer problems we care about). Not numerically robust for
lattices with huge condition number — we use Fraction for exact
arithmetic to sidestep that.

Reference: Cohen, "A Course in Computational Algebraic Number
Theory", algorithm 2.6.3.
"""
from __future__ import annotations
from fractions import Fraction
from typing import Iterable


Vec = list          # list of Fractions / ints


def _dot(a: Vec, b: Vec) -> Fraction:
    return sum(x * y for x, y in zip(a, b))


def lll(basis: list[list[int]], delta: Fraction = Fraction(3, 4)) -> list[list[int]]:
    """
    Run LLL on the given integer basis (a list of row vectors). Returns
    a reduced basis (also integer vectors).
    """
    b = [list(row) for row in basis]
    n = len(b)
    if n == 0:
        return []
    mu = [[Fraction(0)] * n for _ in range(n)]
    bstar: list[list[Fraction]] = [[Fraction(0)] * len(b[0]) for _ in range(n)]
    B: list[Fraction] = [Fraction(0)] * n

    def gs_update():
        nonlocal mu, bstar, B
        for i in range(n):
            bstar[i] = [Fraction(x) for x in b[i]]
            for j in range(i):
                if B[j] == 0:
                    mu[i][j] = Fraction(0)
                    continue
                mu[i][j] = _dot([Fraction(x) for x in b[i]], bstar[j]) / B[j]
                bstar[i] = [bstar[i][k] - mu[i][j] * bstar[j][k]
                            for k in range(len(bstar[i]))]
            B[i] = _dot(bstar[i], bstar[i])

    gs_update()

    k = 1
    while k < n:
        # Size reduction
        for j in range(k - 1, -1, -1):
            if abs(mu[k][j]) > Fraction(1, 2):
                r = int(round(mu[k][j]))
                b[k] = [b[k][i] - r * b[j][i] for i in range(len(b[k]))]
                gs_update()
        # Lovász condition
        if B[k] >= (delta - mu[k][k - 1] ** 2) * B[k - 1]:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            gs_update()
            k = max(1, k - 1)

    return b
