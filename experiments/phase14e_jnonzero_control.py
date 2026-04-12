"""
Phase 14e — Control on j ≠ 0 curves.

Hypothesis (from Phase 14d): the leading-digit rank of M[i][k] =
(δ_{iG}(k) / p) mod p equals the number of non-trivial Aut(E)-orbits
on E(F_p).

For j = 0: Aut(E) = μ₆, rank = (g_ord - 1) / 6.
For j ≠ 0, j ≠ 1728: Aut(E) = {±1}, rank should be (g_ord - 1) / 2.
For j = 1728: Aut(E) = μ₄, rank should be (g_ord - 1) / 4.

We construct j ≠ 0 curves y² = x³ + ax + b with a ≠ 0 and measure.
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


def delta_leading(E_p, C, p, H, k):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kH = E_p.mul(k, H)
    tau_H = C.from_affine(teichmuller_lift(E_p, E_aff, H))
    tau_kH = C.from_affine(teichmuller_lift(E_p, E_aff, kH))
    diff = C.add(C.mul(k, tau_H), C.neg(tau_kH))
    X, Y, Z = diff
    try:
        v = (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return 0
    return (v // p) % p


def rank_over_Fp(M, p):
    M = [row[:] for row in M]
    rows = len(M); cols = len(M[0]) if M else 0
    rank = 0; col = 0
    while rank < rows and col < cols:
        pivot = None
        for r in range(rank, rows):
            if M[r][col] % p != 0:
                pivot = r; break
        if pivot is None:
            col += 1; continue
        M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col] % p, -1, p)
        M[rank] = [(x * inv) % p for x in M[rank]]
        for r in range(rows):
            if r == rank: continue
            f = M[r][col] % p
            if f == 0: continue
            M[r] = [(M[r][c] - f * M[rank][c]) % p for c in range(cols)]
        rank += 1; col += 1
    return rank


def find_good_curve(p, j_type):
    """Find (a, b) over F_p giving a curve with:
       j_type = 'j0'    → a = 0, b ≠ 0  (Aut = μ₆ if p ≡ 1 mod 3)
       j_type = 'j1728' → a ≠ 0, b = 0  (Aut = μ₄ if p ≡ 1 mod 4)
       j_type = 'generic' → a, b ≠ 0, j ≠ 0, 1728  (Aut = ±1)
    with group order not divisible by p, g_ord reasonable."""
    from math import gcd
    for a in range(p):
        for b in range(p):
            if a == 0 and b == 0:
                continue
            disc = (-16 * (4 * a ** 3 + 27 * b ** 2)) % p
            if disc == 0:
                continue
            if j_type == 'j0' and a != 0:
                continue
            if j_type == 'j1728' and b != 0:
                continue
            if j_type == 'generic':
                if a == 0 or b == 0:
                    continue
                # j = 1728 * 4a^3 / (4a^3 + 27b^2)
                num = (1728 * 4 * pow(a, 3, p)) % p
                den = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
                if den == 0:
                    continue
                j = (num * pow(den, -1, p)) % p
                if j == 0 or j == 1728 % p:
                    continue
            E = Curve(a=a, b=b, N=p)
            n = naive_order(E)
            if n % p == 0:
                continue
            G, g_ord = find_generator(E)
            if g_ord < 20:
                continue
            return (a, b, E, G, g_ord)
    return None


def run(p, j_type, e=3, size=None):
    found = find_good_curve(p, j_type)
    if found is None:
        return {"prime": p, "j_type": j_type, "status": "no curve found"}
    a, b, E_p, G, g_ord = found
    if size is None:
        size = min(g_ord - 2, 30)
    N = p ** e
    C = ProjCurve(a=a, b=b, N=N)

    ks = list(range(2, size + 2))
    M = []
    for i in range(1, size + 1):
        H = E_p.mul(i, G)
        row = [delta_leading(E_p, C, p, H, k) for k in ks]
        M.append(row)
    rank = rank_over_Fp(M, p)
    expected_j0 = (g_ord - 1) // 6
    expected_j1728 = (g_ord - 1) // 4
    expected_generic = (g_ord - 1) // 2
    return {
        "prime": p, "j_type": j_type, "a": a, "b": b,
        "g_ord": g_ord, "size": size, "rank": rank,
        "expected_rank_if_aut_mu6": expected_j0,
        "expected_rank_if_aut_mu4": expected_j1728,
        "expected_rank_if_aut_pm1": expected_generic,
        "matches": (
            "μ₆" if rank == expected_j0 else
            "μ₄" if rank == expected_j1728 else
            "±1" if rank == expected_generic else
            "none — other structure?"),
    }


def main():
    cases = [
        (31, 'j0'), (31, 'j1728'), (31, 'generic'),
        (43, 'j0'), (43, 'j1728'), (43, 'generic'),
        (67, 'j0'), (67, 'j1728'), (67, 'generic'),
        (73, 'j0'), (73, 'j1728'), (73, 'generic'),
        (79, 'j0'), (79, 'j1728'), (79, 'generic'),
        (97, 'j0'), (97, 'j1728'), (97, 'generic'),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    all_reports = []
    for p, jt in cases:
        print(f"[phase14e] p={p} j_type={jt}")
        try:
            report = run(p, jt)
        except Exception as exc:
            report = {"prime": p, "j_type": jt, "error": repr(exc)}
        all_reports.append(report)
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase14e_all.json").write_text(
        json.dumps(all_reports, indent=2, default=str))


if __name__ == "__main__":
    main()
