"""
Phase 35 — LLL / lattice attack on the 2-cocycle matrix M[i,j] = z(c(iG, jG)).

Phase 31 established that M is symmetric, has full leading-digit
F_p rank, and shows no obvious bilinear or polynomial structure.
But "rank deficient" is only one possible signal of structure;
short integer relations among the entries are another.

Setup:
  Treat the upper-triangular entries of M as integers in [0, p^e).
  Build a lattice L spanned by:
        e_{ij} := (0, ..., 0, M[i,j], 0, ..., 0, p^e)
  in Z^{N+1}, where N = (size choose 2 + size). The last coordinate
  is p^e to enforce the modular arithmetic.  A short vector in L
  with small first N coordinates corresponds to an integer linear
  combination Σ c_{ij} M[i,j] ≡ 0 (mod p^e) with small c_{ij}.

  Such a relation, if it exists with much smaller c_{ij} than
  random, indicates hidden algebraic structure: the cocycle entries
  satisfy an unexpected congruence.

If LLL finds nothing short (i.e., shortest vector is of typical
"random lattice" size ≈ p^(e·N/(N+1))), the cocycle is
algebraically generic with respect to short Z-relations.

This is the same trick that breaks weak RSA / NTRU / DSA (many
small biased nonces).
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
from experiments.phase31_2cocycle_2var import z_cocycle


def fpylll_or_olll_reduce(B):
    """Reduce a basis B via Sage if available, then fpylll, then fallback."""
    try:
        from sage.all import Matrix, ZZ
        M = Matrix(ZZ, B)
        L = M.LLL()
        return [list(L.row(i)) for i in range(L.nrows())]
    except ImportError:
        pass
    try:
        from fpylll import IntegerMatrix, LLL
        m = IntegerMatrix.from_matrix(B)
        LLL.reduction(m)
        return [[m[i, j] for j in range(m.ncols)] for i in range(m.nrows)]
    except ImportError:
        pass
    return _olll(B)


def _olll(B):
    """Pedestrian integer LLL with delta=0.75. O(d^4 log) so small only."""
    from fractions import Fraction
    B = [list(row) for row in B]
    n = len(B)
    if n == 0:
        return B

    def dot(u, v):
        return sum(a * b for a, b in zip(u, v))

    def gram_schmidt(B):
        Bs = []
        mu = [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            bi = [Fraction(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = Fraction(
                    sum(Fraction(B[i][k]) * Bs[j][k] for k in range(len(Bs[j])))
                ) / Fraction(sum(bs * bs for bs in Bs[j]))
                bi = [bi[k] - mu[i][j] * Bs[j][k] for k in range(len(bi))]
            Bs.append(bi)
        return Bs, mu

    delta = Fraction(3, 4)
    k = 1
    while k < n:
        Bs, mu = gram_schmidt(B)
        # size reduce
        for j in range(k - 1, -1, -1):
            q = round(mu[k][j])
            if q != 0:
                B[k] = [B[k][i] - q * B[j][i] for i in range(len(B[k]))]
        Bs, mu = gram_schmidt(B)
        # Lovasz check
        norm_k = sum(b * b for b in Bs[k])
        norm_km1 = sum(b * b for b in Bs[k - 1])
        if norm_k >= (delta - mu[k][k - 1] * mu[k][k - 1]) * norm_km1:
            k += 1
        else:
            B[k], B[k - 1] = B[k - 1], B[k]
            k = max(k - 1, 1)
    return B


def lattice_attack(M_entries, modulus):
    """Build the lattice [I | scaled_entries] / modulus row, run LLL,
    return the shortest vector and its 'expected' size (Gaussian
    heuristic) for comparison."""
    N = len(M_entries)
    # Basis: N+1 rows of length N+1
    # rows 0..N-1: e_i augmented by M_entries[i]
    # row N:        zeros augmented by modulus
    # We scale the modular column by a large weight W to make LLL
    # prefer small linear combos that vanish modularly.
    W = modulus  # weight on modular coordinate
    B = []
    for i in range(N):
        row = [0] * (N + 1)
        row[i] = 1
        row[N] = M_entries[i] * W // 1
        B.append(row)
    B.append([0] * (N + 1))
    B[N][N] = modulus * W

    reduced = fpylll_or_olll_reduce(B)
    # Find shortest vector among reduced basis
    def norm2(v):
        return sum(x * x for x in v[:N])  # only the c_i part matters
    sorted_b = sorted(reduced, key=norm2)
    shortest = sorted_b[0]
    # Drop the modular column for reporting
    coeffs = shortest[:N]
    relation_residue = sum(c * m for c, m in zip(coeffs, M_entries)) % modulus
    return {
        "N": N,
        "shortest_norm2": norm2(shortest),
        "shortest_max_abs": max(abs(x) for x in coeffs),
        "shortest_coeffs": coeffs[:10],
        "relation_residue": relation_residue,
        "modulus": modulus,
        "expected_random_max": modulus ** (1.0 / N),  # rough heuristic
    }


def run(p, b, e=4, size=10):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"p": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < size + 4:
        return {"p": p, "status": "small generator"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Build the upper-triangular cocycle entries (i ≤ j)
    entries = []
    coords = []
    for i in range(1, size + 1):
        for j in range(i, size + 1):
            v = z_cocycle(E_p, C, G, i, j)
            entries.append(v if v is not None else 0)
            coords.append((i, j))

    res = lattice_attack(entries, N)
    res["p"] = p
    res["b"] = b
    res["g_ord"] = g_ord
    res["size"] = size
    res["e"] = e
    res["coords"] = coords[:10]
    return res


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3),
                 (97, 5), (103, 5)]:
        print(f"[phase35] p={p} b={b}")
        try:
            r = run(p, b, e=4, size=8)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": repr(exc)}
        reports.append(r)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase35_lll_cocycle.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
