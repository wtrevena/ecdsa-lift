"""
Polynomial arithmetic over Z/N.

Polynomials are represented as lists of coefficients in increasing
degree: [a0, a1, a2, ...] means a0 + a1 x + a2 x^2 + ...

Division and GCD require leading coefficients to be UNITS in Z/N.
When they aren't, operations raise ZeroDivisionError — the caller
decides what to do.

Enough to run Cantor's algorithm on a genus-2 Jacobian over Z/p^e.
"""
from __future__ import annotations
from typing import List

from .curves import modinv

Poly = List[int]


def trim(p: Poly) -> Poly:
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(a: Poly, b: Poly, N: int) -> Poly:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % N
    return trim(out)


def sub(a: Poly, b: Poly, N: int) -> Poly:
    n = max(len(a), len(b))
    out = [0] * n
    for i in range(n):
        out[i] = ((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % N
    return trim(out)


def mul(a: Poly, b: Poly, N: int) -> Poly:
    if not a or not b:
        return [0]
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x == 0:
            continue
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % N
    return trim(out)


def scale(a: Poly, s: int, N: int) -> Poly:
    return trim([(x * s) % N for x in a])


def deg(a: Poly) -> int:
    a = trim(a)
    return len(a) - 1 if a != [0] else -1


def divmod_poly(a: Poly, b: Poly, N: int):
    """a = q * b + r with deg(r) < deg(b). Requires b's leading
    coefficient to be a unit in Z/N."""
    a = trim(a)
    b = trim(b)
    if b == [0]:
        raise ZeroDivisionError("divide by zero poly")
    b_lead = b[-1]
    b_lead_inv = modinv(b_lead, N)
    q = [0] * max(1, len(a) - len(b) + 1)
    r = list(a)
    while len(trim(r)) >= len(b) and trim(r) != [0]:
        r = trim(r)
        shift = len(r) - len(b)
        coef = (r[-1] * b_lead_inv) % N
        q[shift] = coef
        # r -= coef * x^shift * b
        for j, bj in enumerate(b):
            r[shift + j] = (r[shift + j] - coef * bj) % N
        r = trim(r)
    return trim(q), trim(r)


def gcd_poly(a: Poly, b: Poly, N: int):
    """Extended gcd: returns (d, s, t) with d = s*a + t*b.
    Requires leading coefficients of intermediate remainders to be
    units. Raises ZeroDivisionError if a non-unit pivot is
    encountered."""
    old_r, r = trim(a), trim(b)
    old_s, s = [1], [0]
    old_t, t = [0], [1]
    while r != [0]:
        q, rem = divmod_poly(old_r, r, N)
        old_r, r = r, rem
        qs = mul(q, s, N)
        new_s = sub(old_s, qs, N)
        old_s, s = s, new_s
        qt = mul(q, t, N)
        new_t = sub(old_t, qt, N)
        old_t, t = t, new_t
    # Normalize old_r to monic
    if old_r == [0]:
        return old_r, old_s, old_t
    lead = old_r[-1]
    lead_inv = modinv(lead, N)
    old_r = scale(old_r, lead_inv, N)
    old_s = scale(old_s, lead_inv, N)
    old_t = scale(old_t, lead_inv, N)
    return old_r, old_s, old_t


def eval_poly(p: Poly, x: int, N: int) -> int:
    v = 0
    for c in reversed(p):
        v = (v * x + c) % N
    return v
