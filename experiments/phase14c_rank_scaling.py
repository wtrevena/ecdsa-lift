"""
Phase 14c — Rank scaling of the leading-digit δ matrix.

Phase 14b found that M[i][k] = (δ_{iG}(k) / p) mod p is rank-deficient
on every test curve (e.g. 6/12 for p=31, 10/12 for p=97), with the
deficit confined to layer 0 (the leading p-adic coefficient).

This could be either:
  (A) Genuine low-rank structure — rank saturates below min(rows, cols)
      as size grows, implying a factorization of δ/p mod p.
  (B) Small-matrix artifact — rank grows proportionally with size until
      it hits the bound.

14c settles it: for each curve, build matrices at several sizes
s = 4, 6, 8, 10, ..., up to ≈ g_ord/2, and record the rank. If rank
saturates, report the saturation point and the row-kernel basis;
otherwise record the linear growth.
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
    """Return δ_H(k) = z-coord of [k]τ(H) − τ(kH) as an element of Z/N."""
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


def rank_and_kernel_over_Fp(M, p):
    """Gaussian elimination over F_p. Returns (rank, row_kernel_basis).
    row_kernel_basis is a list of vectors v such that v · M = 0."""
    M = [row[:] for row in M]
    rows = len(M)
    cols = len(M[0]) if M else 0
    # Track row operations to recover the kernel
    aug = [row[:] + [0] * i + [1] + [0] * (rows - i - 1) for i, row in enumerate(M)]
    rank = 0
    col = 0
    pivot_rows = []
    while rank < rows and col < cols:
        pivot = None
        for r in range(rank, rows):
            if aug[r][col] % p != 0:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        aug[rank], aug[pivot] = aug[pivot], aug[rank]
        inv = pow(aug[rank][col] % p, -1, p)
        aug[rank] = [(x * inv) % p for x in aug[rank]]
        for r in range(rows):
            if r == rank:
                continue
            f = aug[r][col] % p
            if f == 0:
                continue
            aug[r] = [(aug[r][c] - f * aug[rank][c]) % p for c in range(len(aug[r]))]
        pivot_rows.append(rank)
        rank += 1
        col += 1
    # Rows of aug with all-zero first `cols` entries are in the kernel.
    # Their second half (last `rows` columns) gives the kernel combination.
    kernel = []
    for r in range(rows):
        if all(aug[r][c] % p == 0 for c in range(cols)):
            kernel.append(aug[r][cols:])
    return rank, kernel


def measure_rank_curve(p, b, e=3, max_size=None):
    E_p = Curve(a=0, b=b, N=p)
    n_ord = naive_order(E_p)
    if n_ord % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 8:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    if max_size is None:
        max_size = min(g_ord - 2, 32)

    sizes = list(range(4, max_size + 1, 2))
    results = []

    # Build the full max_size × max_size matrix once, then take submatrices
    full_H = [E_p.mul(i, G) for i in range(1, max_size + 1)]
    # Skip k = 1 (δ = 0 trivially). Use k = 2 .. max_size + 1.
    full_k = list(range(2, max_size + 2))
    full_M = []
    for H in full_H:
        row = []
        for k in full_k:
            v = delta_H_k(E_p, C, H, k)
            # Normalize: divide by p and reduce mod p (layer 0)
            row.append((v // p) % p)
        full_M.append(row)

    for s in sizes:
        M = [row[:s] for row in full_M[:s]]
        rank, kernel = rank_and_kernel_over_Fp(M, p)
        results.append({
            "size": s,
            "rank": rank,
            "kernel_dim": s - rank,
            "saturated_below_size": rank < s,
        })

    # Also extract kernel at max size
    M_max = [row[:max_size] for row in full_M[:max_size]]
    rank_max, kernel_max = rank_and_kernel_over_Fp(M_max, p)

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "max_size": max_size,
        "rank_progression": results,
        "final_rank": rank_max,
        "final_kernel_dim": max_size - rank_max,
        "kernel_basis_sample": kernel_max[:3] if kernel_max else [],
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase14c] p={p} b={b}")
        try:
            report = measure_rank_curve(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase14c_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            if k == "rank_progression":
                print(f"   {k}:")
                for r in v:
                    print(f"      size={r['size']:3d} rank={r['rank']:3d} "
                          f"ker_dim={r['kernel_dim']}")
            else:
                print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
