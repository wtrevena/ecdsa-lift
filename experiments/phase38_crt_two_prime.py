"""
Phase 38 — CRT two-prime lift of the cocycle (Tier B1 from meta_reflection.md).

The proposal: lift y^2 = x^3 + b to Z/(p^e · q^f) via Hensel + CRT and
look at the joint structure of the cocycle.  By the Chinese Remainder
Theorem, E(Z/(p^e · q^f)) ≅ E(Z/p^e) × E(Z/q^f), so any *function* of
points decomposes as a direct sum. The 2-cocycle c(P,Q) = τ(P+Q) -
τ(P) - τ(Q) therefore satisfies
        c_{pq}(P, Q) = (c_p(P_p, Q_p), c_q(P_q, Q_q))
componentwise.  This means there is NO additional structure to be
gained by gluing — the components are independent.

But there is one non-tautological question hidden inside this:
given two primes (p, b=7) and (q, b=7) — i.e. the *same curve*
y^2 = x^3 + 7 reduced at two different primes — is there a *correlation*
between the F_p cocycle and the F_q cocycle that is invisible at either
single prime?  Specifically:
  - Pick generators G_p ∈ E(F_p) and G_q ∈ E(F_q).
  - Build the 2-cocycle matrices M_p[i,j] = z_p(c_p(iG_p, jG_p))
    and M_q[i,j] = z_q(c_q(iG_q, jG_q)) for small i,j.
  - Compute correlation coefficients between vec(M_p) and vec(M_q).
  - If non-zero, there is a curve-intrinsic (b-dependent) feature
    that survives across primes.
  - Run LLL on the joint vector [M_p / p^e ‖ M_q / q^f] to look for
    a short integer relation across the two primes.

If everything is at noise level, this confirms the CRT direct-sum
decomposition empirically and closes B1.
"""
from __future__ import annotations
import sys
import pathlib
import math
import json
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift
from experiments.phase31_2cocycle_2var import z_cocycle


def cocycle_matrix(p, b, e=4, size=10):
    """Return M[i,j] = z(c(iG, jG)) mod p^e for i,j = 1..size."""
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return None, None
    G, g_ord = find_generator(E_p)
    if g_ord <= size + 2:
        return None, None
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    M = np.zeros((size, size), dtype=np.int64)
    for i in range(1, size + 1):
        for j in range(1, size + 1):
            v = z_cocycle(E_p, C, G, i, j)
            M[i - 1, j - 1] = v if v is not None else 0
    return M, N


def normalize_to_unit_box(M, N):
    """Map entries from [0, N) to [-0.5, 0.5)."""
    return ((M.astype(np.float64) / N) + 0.5) % 1.0 - 0.5


def cross_prime_correlation(p_list, b=7, e=4, size=8):
    """Compute pairwise correlation of normalized cocycle matrices."""
    Ms = {}
    Ns = {}
    for p in p_list:
        M, N = cocycle_matrix(p, b, e, size)
        if M is None:
            print(f"   skip p={p} (small generator or anomalous)")
            continue
        Ms[p] = M
        Ns[p] = N
    keys = sorted(Ms)
    out = {"primes": keys, "correlations": {}}
    # Pairwise pearson correlation between normalized vec(M)
    for i, p in enumerate(keys):
        for q in keys[i + 1:]:
            v1 = normalize_to_unit_box(Ms[p], Ns[p]).flatten()
            v2 = normalize_to_unit_box(Ms[q], Ns[q]).flatten()
            v1 -= v1.mean()
            v2 -= v2.mean()
            denom = math.sqrt(float((v1 * v1).sum() * (v2 * v2).sum()))
            corr = float((v1 * v2).sum() / denom) if denom > 0 else 0.0
            out["correlations"][f"{p},{q}"] = round(corr, 6)
    return out, Ms, Ns


def cross_prime_lll(Ms, Ns):
    """Stack the cocycle entries from multiple primes into a single vector
    and run LLL on the lattice [I ‖ scaled_entries] to find short integer
    relations across primes.

    For each entry index k, the lattice column is:
      (e_k, M_p[k]·W/N_p, M_q[k]·W/N_q, ...)
    where W is a weight; the 'modular' rows are p, q, ... ensuring the
    relation holds mod each prime power. A short combination of M-rows
    that vanishes simultaneously mod every N corresponds to an integer
    relation across primes."""
    try:
        from fpylll import IntegerMatrix, LLL
        use_fpylll = True
    except ImportError:
        use_fpylll = False
    if not use_fpylll:
        return {"error": "fpylll not available; skipping cross-prime LLL"}
    primes = sorted(Ms)
    if len(primes) < 2:
        return {"error": "need ≥ 2 primes"}
    # Flatten matrices
    flats = {p: Ms[p].flatten().tolist() for p in primes}
    K = len(flats[primes[0]])
    W = max(Ns.values())  # large weight
    # Lattice rows: K + len(primes) rows of dimension K + len(primes)
    rows = []
    for k in range(K):
        row = [0] * (K + len(primes))
        row[k] = 1
        for i, p in enumerate(primes):
            row[K + i] = flats[p][k] * (W // Ns[p])
        rows.append(row)
    for i, p in enumerate(primes):
        row = [0] * (K + len(primes))
        row[K + i] = Ns[p] * (W // Ns[p])
        rows.append(row)
    Mfp = IntegerMatrix.from_matrix(rows)
    LLL.reduction(Mfp)
    shortest = None
    shortest_norm = None
    for r in range(Mfp.nrows):
        v = [Mfp[r, c] for c in range(Mfp.ncols)]
        coeffs = v[:K]
        if all(c == 0 for c in coeffs):
            continue
        norm2 = sum(c * c for c in coeffs)
        max_abs = max(abs(c) for c in coeffs)
        if shortest is None or norm2 < shortest_norm:
            shortest_norm = norm2
            shortest = (coeffs, max_abs, norm2)
    return {
        "primes": primes,
        "K": K,
        "shortest_norm2": shortest[2] if shortest else None,
        "shortest_max_abs": shortest[1] if shortest else None,
        "shortest_coeffs_first10": shortest[0][:10] if shortest else None,
        "expected_random_max": int(round(W ** (1.0 / (K + len(primes))))),
    }


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    SECP_PRIMES = [31, 43, 67, 79, 97, 103, 127, 167, 211, 257]
    print("=== Phase 38: cross-prime cocycle correlations on y^2 = x^3 + 7 ===")
    corr, Ms, Ns = cross_prime_correlation(SECP_PRIMES, b=7, e=4, size=8)
    print("Pairwise Pearson correlation of normalized cocycle vec(M):")
    print(f"  primes used: {corr['primes']}")
    vals = list(corr["correlations"].values())
    print(f"  min/max/mean correlation: {min(vals):+.4f} / {max(vals):+.4f} / {sum(vals)/len(vals):+.4f}")
    print("  per-pair (first 10):")
    for k, v in list(corr["correlations"].items())[:10]:
        print(f"    {k}: {v:+.4f}")
    print()

    print("=== Cross-prime LLL on stacked cocycle entries ===")
    lll_result = cross_prime_lll(Ms, Ns)
    for k, v in lll_result.items():
        print(f"   {k}: {v}")

    (out_dir / "phase38_crt_two_prime.json").write_text(
        json.dumps({"correlations": corr, "lll": lll_result},
                   indent=2, default=str))


if __name__ == "__main__":
    main()
