"""
Phase 31 — The 2-cocycle as a function of TWO variables.

All previous phases studied δ(k) := [k]τ(G) − τ(kG) along the single
direction k → [k]G.  But the underlying object is the symmetric
2-cocycle

    c(P, Q) := τ(P + Q) − τ(P) − τ(Q)  ∈ Ê(p Z/p^e)

a function E(F_p) × E(F_p) → Ê.  We have NEVER probed it as a
2-variable object.  Open questions:

  (Q1) Is `z(c(P, Q))` (the formal-group parameter of the 2-cocycle)
       BILINEAR in (P, Q), or merely symmetric and continuous?
       Bilinear here means: there exists a Z_p-bilinear pairing
            β : E(F_p) × E(F_p) → Z_p
       with `z(c(P, Q)) ≡ β(P, Q) (mod p^?)`.
  (Q2) Equivalently, viewing `f(i, j) := z(c(iG, jG))` as a function
       Z/n × Z/n → Z/p^e, what is the rank of the matrix (f(i, j))?
       A bilinear form has rank ≤ 1 over Q_p; an arbitrary symmetric
       function has rank up to n.
  (Q3) Are there partial-derivative-like recurrences:
            f(i+1, j) − f(i, j) = g(i, j)
       where g has lower complexity?  This is the discrete analogue
       of "f is locally polynomial".
  (Q4) Does the matrix have any eigenstructure under the
       Aut(E) = μ₆ action (P → ζ·P)?

Run with the standard Python venv (no Sage required for the basic
2-variable test).
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


def z_cocycle(E_p, C, G, i, j):
    """Compute z(c(iG, jG)) = z( τ(iG+jG) − τ(iG) − τ(jG) ) in Z/N."""
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    iG = E_p.mul(i, G)
    jG = E_p.mul(j, G)
    ijG = E_p.mul((i + j) % naive_order(E_p), G)  # = (i+j)G
    tau_iG = C.from_affine(teichmuller_lift(E_p, E_aff, iG))
    tau_jG = C.from_affine(teichmuller_lift(E_p, E_aff, jG))
    tau_ijG = C.from_affine(teichmuller_lift(E_p, E_aff, ijG))
    diff = C.add(tau_ijG, C.add(C.neg(tau_iG), C.neg(tau_jG)))
    X, Y, Z = diff
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


def matrix_rank_mod_N(M, N, p, e):
    """Compute the F_p-rank of M after stripping common p-factors
    on each entry. Returns the leading-digit rank."""
    rows = len(M)
    cols = len(M[0]) if M else 0
    # Reduce each entry: divide by p once and take mod p
    M2 = []
    for r in M:
        row = []
        for x in r:
            if x is None:
                row.append(0)
            else:
                row.append((x // p) % p)
        M2.append(row)
    # F_p Gaussian elimination
    rk = 0
    col = 0
    while rk < rows and col < cols:
        piv = None
        for r in range(rk, rows):
            if M2[r][col] % p != 0:
                piv = r
                break
        if piv is None:
            col += 1
            continue
        M2[rk], M2[piv] = M2[piv], M2[rk]
        inv = pow(M2[rk][col] % p, -1, p)
        M2[rk] = [(x * inv) % p for x in M2[rk]]
        for r in range(rows):
            if r == rk:
                continue
            f = M2[r][col] % p
            if f == 0:
                continue
            M2[r] = [(M2[r][c] - f * M2[rk][c]) % p for c in range(cols)]
        rk += 1
        col += 1
    return rk


def diagonal_difference_test(M, p):
    """Test whether the discrete partial derivative f(i+1, j) - f(i, j)
    has lower rank than f itself. If so, f is "polynomial-like" in i."""
    rows = len(M)
    cols = len(M[0]) if M else 0
    if rows < 2:
        return None
    Dx = [[((M[i+1][j] - M[i][j]) // p) % p if M[i+1][j] is not None else 0
           for j in range(cols)]
          for i in range(rows - 1)]
    rkD = 0
    M2 = [r[:] for r in Dx]
    rk = 0
    col = 0
    while rk < rows - 1 and col < cols:
        piv = None
        for r in range(rk, rows - 1):
            if M2[r][col] % p != 0:
                piv = r
                break
        if piv is None:
            col += 1
            continue
        M2[rk], M2[piv] = M2[piv], M2[rk]
        inv = pow(M2[rk][col] % p, -1, p)
        M2[rk] = [(x * inv) % p for x in M2[rk]]
        for r in range(rows - 1):
            if r == rk:
                continue
            f = M2[r][col] % p
            if f == 0:
                continue
            M2[r] = [(M2[r][c] - f * M2[rk][c]) % p for c in range(cols)]
        rk += 1
        col += 1
    return rk


def symmetry_check(M, N):
    """Verify M[i][j] == M[j][i] mod N (the cocycle is symmetric)."""
    rows = len(M)
    cols = len(M[0]) if M else 0
    if rows != cols:
        return False
    for i in range(rows):
        for j in range(cols):
            if M[i][j] is None or M[j][i] is None:
                continue
            if (M[i][j] - M[j][i]) % N != 0:
                return False
    return True


def run(p, b, e=3, size=18):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < size + 4:
        return {"prime": p, "status": "small generator"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Build the matrix M[i][j] = z(c(iG, jG)) for i, j = 1..size
    M = []
    for i in range(1, size + 1):
        row = []
        for j in range(1, size + 1):
            row.append(z_cocycle(E_p, C, G, i, j))
        M.append(row)

    # Symmetry sanity check
    sym = symmetry_check(M, N)

    # Leading-digit F_p rank of (M / p) mod p
    rank = matrix_rank_mod_N(M, N, p, e)

    # Discrete partial derivative rank
    deriv_rank = diagonal_difference_test(M, p)

    # Bilinearity test: a Z_p-bilinear function would have leading-digit
    # rank at most 1.  We test whether the leading digit f(i, j)/p mod p
    # has rank 1 (= bilinear) or higher (= more structure).
    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N, "size": size,
        "symmetric": sym,
        "leading_digit_rank": rank,
        "is_rank_1_bilinear": rank == 1,
        "discrete_partial_rank": deriv_rank,
        "rank_drops_under_partial": (
            deriv_rank is not None and deriv_rank < rank),
        "first_3x3_leading": [[(M[i][j] // p) % p if M[i][j] is not None else None
                               for j in range(3)] for i in range(3)],
    }


def main():
    candidates = [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5)]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in candidates:
        print(f"[phase31] p={p} b={b}")
        try:
            r = run(p, b)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"prime": p, "error": repr(exc)}
        reports.append(r)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase31_2cocycle_2var.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
