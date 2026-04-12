"""
Phase 39 — Scaling sweep on the Phase 35 LLL cocycle attack.

Phase 35 ran LLL on the upper-triangular 2-cocycle entries for a single
size (8) at e=4 on 7 primes, and found only the trivial cocycle
identities. The natural follow-up: does the shortest non-trivial
relation grow, shrink, or stay constant as we increase n (the matrix
side) and as we increase p?

A genuine structural relation would manifest as: the shortest LLL
vector has max-coefficient that grows much slower than the
Gaussian-heuristic rate p^(e/dim).

We sweep p ∈ secp256k1 family (b=7), size ∈ {6, 8, 10, 12}, and report:
  - dimension dim of the lattice
  - shortest LLL coefficient max |c_i|
  - Gaussian heuristic expected_random_max ≈ p^(e/dim)
  - ratio (shortest_max / expected) — < 1 means LLL did better than random
"""
from __future__ import annotations
import sys
import pathlib
import math
import json

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator
from src.projective import ProjCurve
from experiments.phase31_2cocycle_2var import z_cocycle


def lattice_attack(M_entries, modulus):
    """Build a [I | scaled_modular] lattice and reduce with fpylll.
    Returns shortest non-trivial vector coefficients and stats."""
    from fpylll import IntegerMatrix, LLL
    K = len(M_entries)
    W = modulus
    rows = []
    for i in range(K):
        row = [0] * (K + 1)
        row[i] = 1
        row[K] = M_entries[i] * (W // 1)
        rows.append(row)
    rows.append([0] * (K + 1))
    rows[K][K] = modulus * W
    Mfp = IntegerMatrix.from_matrix(rows)
    LLL.reduction(Mfp)
    best = None
    best_norm2 = None
    for r in range(Mfp.nrows):
        v = [Mfp[r, c] for c in range(Mfp.ncols)]
        coeffs = v[:K]
        if all(c == 0 for c in coeffs):
            continue
        norm2 = sum(c * c for c in coeffs)
        if best is None or norm2 < best_norm2:
            best_norm2 = norm2
            best = coeffs
    if best is None:
        return None
    max_abs = max(abs(c) for c in best)
    nonzero = sum(1 for c in best if c != 0)
    expected = modulus ** (1.0 / (K + 1))
    return {
        "K": K,
        "shortest_norm2": best_norm2,
        "shortest_max_abs": max_abs,
        "shortest_nonzero_count": nonzero,
        "expected_random_max": expected,
        "ratio_max_to_expected": max_abs / expected,
    }


def run(p, b, e, size):
    E_p = Curve(a=0, b=b, N=p)
    if naive_order(E_p) % p == 0:
        return {"status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord <= size + 4:
        return {"status": f"small generator g_ord={g_ord}"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    entries = []
    for i in range(1, size + 1):
        for j in range(i, size + 1):
            v = z_cocycle(E_p, C, G, i, j)
            entries.append(v if v is not None else 0)
    res = lattice_attack(entries, N)
    if res is None:
        return {"status": "lll failed"}
    res.update({"p": p, "b": b, "e": e, "size": size, "g_ord": g_ord})
    return res


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    PRIMES = [31, 43, 67, 79, 97, 103, 167, 211, 257, 311, 401, 521]
    SIZES = [6, 8, 10, 12]
    e = 4
    reports = []
    print(f"=== Phase 39 LLL scaling sweep on y^2 = x^3 + 7, e={e} ===")
    print(f"{'p':>6}{'size':>6}{'dim':>6}{'max|c|':>10}{'expected':>14}{'ratio':>10}{'nonzero':>10}")
    for p in PRIMES:
        for size in SIZES:
            r = run(p, 7, e, size)
            r["p_input"] = p
            r["size_input"] = size
            reports.append(r)
            if "shortest_max_abs" in r:
                print(f"{p:>6}{size:>6}{r['K'] + 1:>6}{r['shortest_max_abs']:>10}"
                      f"{r['expected_random_max']:>14.2f}"
                      f"{r['ratio_max_to_expected']:>10.4f}"
                      f"{r['shortest_nonzero_count']:>10}")
            else:
                print(f"{p:>6}{size:>6}  {r}")
    (out_dir / "phase39_lll_scaling.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
