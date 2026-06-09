"""
Phase 14 — Does the Teichmüller error δ_H(k) factor as f(H)·g(k)?

Motivation. We've only tested δ(k) with a fixed base point G. If the
matrix M[H][k] := δ_H(k) = z-coord of ([k]τ(H) − τ(kH)) has rank 1
over Z/p^e, then δ_H(k) = f(H)·g(k) and g(k) is a linear-in-k
invariant independent of H — recoverable from a single column.

Rank 1 means every row is a scalar multiple of every other row.
Equivalently: δ_H(k)·δ_{H'}(k') ≡ δ_H(k')·δ_{H'}(k) mod p^e.

More generally, low rank r means g is r-dimensional and extractable
via SVD/Smith normal form. Even rank 2 would give a 2D lattice
attack.
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv, points
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def delta_H_k(E_p: Curve, C: ProjCurve, H, k: int):
    """z-coord of ([k]·τ(H) − τ(kH)) in Ê."""
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kH = E_p.mul(k, H)
    tau_H = C.from_affine(teichmuller_lift(E_p, E_aff, H))
    tau_kH = C.from_affine(teichmuller_lift(E_p, E_aff, kH))
    diff = C.add(C.mul(k, tau_H), C.neg(tau_kH))
    X, Y, Z = diff
    if Y % C.N == 0:
        return None
    return (-X * modinv(Y % C.N, C.N)) % C.N


def matrix_rank_mod_N(M, N):
    """Gaussian elimination mod N with unit-pivot check. Returns
    the number of rows we could reduce."""
    M = [row[:] for row in M]
    rows = len(M)
    cols = len(M[0]) if M else 0
    rank = 0
    col = 0
    while rank < rows and col < cols:
        pivot = None
        for r in range(rank, rows):
            v = M[r][col] % N
            if v == 0:
                continue
            try:
                modinv(v, N)
                pivot = r
                break
            except ZeroDivisionError:
                continue
        if pivot is None:
            col += 1
            continue
        M[rank], M[pivot] = M[pivot], M[rank]
        inv = modinv(M[rank][col] % N, N)
        M[rank] = [(x * inv) % N for x in M[rank]]
        for r in range(rows):
            if r == rank:
                continue
            f = M[r][col] % N
            if f == 0:
                continue
            M[r] = [(M[r][c] - f * M[rank][c]) % N for c in range(cols)]
        rank += 1
        col += 1
    return rank


def run(p: int, b: int, e: int = 3, max_H: int = 10, max_k: int = 15) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < max_k + 2:
        return {"prime": p, "status": "too small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Pick max_H distinct non-identity base points.
    # Use multiples of G so we stay in the cyclic subgroup.
    Hs = [E_p.mul(i, G) for i in range(1, min(max_H, g_ord - 1) + 1)]
    ks = list(range(1, min(max_k, g_ord - 1) + 1))

    M = []
    for H in Hs:
        row = []
        for k in ks:
            v = delta_H_k(E_p, C, H, k)
            row.append(0 if v is None else v)
        M.append(row)

    rank = matrix_rank_mod_N(M, N)
    rank_modp = matrix_rank_mod_N(M, p)
    rank_modp2 = matrix_rank_mod_N(M, p * p)

    # Is the matrix rank-1? Test via 2x2 minors vanishing.
    rank1_hits = rank1_total = 0
    for i in range(len(Hs)):
        for j in range(i + 1, len(Hs)):
            for a in range(len(ks)):
                for c in range(a + 1, len(ks)):
                    rank1_total += 1
                    det = (M[i][a] * M[j][c] - M[i][c] * M[j][a]) % N
                    if det == 0:
                        rank1_hits += 1

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "matrix_shape": [len(Hs), len(ks)],
        "rank_mod_N": rank,
        "rank_mod_p": rank_modp,
        "rank_mod_p2": rank_modp2,
        "rank1_minor_hits": f"{rank1_hits}/{rank1_total}",
        "verdict": ("rank-1 factorization — ATTACK POSSIBLE"
                    if rank == 1
                    else f"generic rank {rank}, no low-rank structure"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase14] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase14_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        print(f"   rank_mod_N: {report.get('rank_mod_N')}  "
              f"rank_mod_p: {report.get('rank_mod_p')}  "
              f"rank1_minors: {report.get('rank1_minor_hits')}")
        print(f"   {report.get('verdict', report.get('status'))}")


if __name__ == "__main__":
    main()
