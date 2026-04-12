"""
Phase 14b — Multi-base δ factorization, FIXED.

The Phase 14 rank was reported as 0 because every δ value is
divisible by p (δ ∈ Ê means δ ≡ 0 mod p). The correct object is
δ/p ∈ Z/p^{e-1}. We redo the rank computation on the normalized
matrix M'[H][k] = δ_H(k) / p mod p^{e-1}.

A rank-1 M' would factor δ_H(k) = p · f(H) · g(k), giving a
one-parameter family of linear functions in k per H.
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


def delta_H_k(E_p, C, H, k):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kH = E_p.mul(k, H)
    tau_H = C.from_affine(teichmuller_lift(E_p, E_aff, H))
    tau_kH = C.from_affine(teichmuller_lift(E_p, E_aff, kH))
    diff = C.add(C.mul(k, tau_H), C.neg(tau_kH))
    X, Y, Z = diff
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return 0


def rank_over_Fq(M, q):
    """Gaussian elimination over F_q (q prime)."""
    M = [row[:] for row in M]
    rows = len(M)
    cols = len(M[0]) if M else 0
    rank = 0
    col = 0
    while rank < rows and col < cols:
        pivot = None
        for r in range(rank, rows):
            if M[r][col] % q != 0:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col] % q, -1, q)
        M[rank] = [(x * inv) % q for x in M[rank]]
        for r in range(rows):
            if r == rank:
                continue
            f = M[r][col] % q
            if f == 0:
                continue
            M[r] = [(M[r][c] - f * M[rank][c]) % q for c in range(cols)]
        rank += 1
        col += 1
    return rank


def run(p, b, e=4, max_H=12, max_k=20):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < max_k + 2:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    Hs = [E_p.mul(i, G) for i in range(1, min(max_H, g_ord - 1) + 1)]
    # Skip k = 1 (trivially 0) and k = g_ord - 1 (related to k=1 by Phase 21b)
    ks = list(range(2, min(max_k + 2, g_ord // 2 + 1)))

    M_raw = []
    for H in Hs:
        row = []
        for k in ks:
            v = delta_H_k(E_p, C, H, k)
            row.append(v)
        M_raw.append(row)

    # Check every entry is divisible by p (sanity: they are in Ê)
    all_p = all(e % p == 0 for row in M_raw for e in row)

    # Normalize: divide by p
    M_norm = [[(e // p) % (p ** (e_dummy - 1)) for e in row] for row in M_raw
              for e_dummy in [e]]  # hack to use e in list comp
    M_norm = [[(e // p) % (p ** (4 - 1)) for e in row] for row in M_raw]

    # Rank over F_p (reduce entries mod p)
    M_mod_p = [[e % p for e in row] for row in M_norm]
    rank_p = rank_over_Fq(M_mod_p, p)

    # Also rank over Z/(p^(e-1)) — requires handling of non-unit pivots
    # Simpler: compute a sequence of ranks over F_p of (M_norm / p^j) for j=0..e-2
    ranks = [rank_p]
    for j in range(1, e - 1):
        Mj = [[(e // (p ** j)) % p for e in row] for row in M_norm]
        ranks.append(rank_over_Fq(Mj, p))

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "matrix_shape": [len(Hs), len(ks)],
        "all_entries_in_E_hat": all_p,
        "ranks_over_Fp_at_layers": ranks,
        "max_rank": max(ranks),
        "verdict": ("low-rank leak" if max(ranks) < min(len(Hs), len(ks)) // 2
                    else "full rank / generic"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase14b] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase14b_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
