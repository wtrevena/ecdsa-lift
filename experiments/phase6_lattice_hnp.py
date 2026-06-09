"""
Phase 6 — Hidden Number Problem on the p-adic lifting data.

Setup. Given a j = 0 curve E/F_p with generator G of prime order n,
and a point Q = [k]G with k secret, we want to recover k from
"noisy" p-adic information extracted from the Teichmüller / canonical
lifts.

HNP framing. The Boneh–Venkatesan Hidden Number Problem says: if we
are given
    t_i · k  +  a_i  ≡  b_i  (mod n)
with k unknown and (a_i, b_i) satisfying |a_i| < n / 2^ℓ for some ℓ
bits, then LLL on the right lattice recovers k once we have enough
samples (ℓ ≥ log n / m for m samples).

Phase 6 experiment. We need to produce (t_i, b_i) pairs where b_i
leaks the top bits of t_i · k mod n, from the p-adic lift data.

Naive attempt. For each trial k ∈ [0, n), compute
    z_k  :=  zcoord(τ(G + [k]G)  -  τ(G)  -  τ([k]G))   ∈ Ê
i.e. the formal-group parameter of the Teichmüller cocycle. In
principle z_k is a function of k — Phase 5 showed it's injective
but neither additive nor CM-equivariant.

If z_k / z_1 (mod p^(e-1)) agrees with k mod p^(e-1) to ANY leading
number of bits, we have exploitable bias.

This experiment:
  (a) Tabulates z_k for k = 1..n-1 across several curves.
  (b) For each curve, checks the correlation between z_k and k.
  (c) Attempts an HNP-style LLL attack using (k, z_k) pairs for a
      small sample of k values (pretending we don't know the full
      correspondence).

A negative result here is still interesting: it rules out a large
family of hidden-number lattice attacks on the p-adic lift data.
"""
from __future__ import annotations
import json
import random
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift
from src.lll import lll


def cocycle_z(E_p, C, G, k: int):
    """Compute z(c) where c = τ(G) + τ([k]G) - τ(G + [k]G) ∈ Ê."""
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kG = E_p.mul(k, G)
    sum_p = E_p.add(G, kG)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
    tau_sum = C.from_affine(teichmuller_lift(E_p, E_aff, sum_p))
    cocycle = C.add(C.add(tau_G, tau_kG), C.neg(tau_sum))
    X, Y, Z = cocycle
    if Y % C.N == 0:
        return None
    return (-X * modinv(Y % C.N, C.N)) % C.N


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "tiny"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Tabulate z_k for k = 1..g_ord-1
    zs = {}
    for k in range(1, g_ord):
        zs[k] = cocycle_z(E_p, C, G, k)

    # Linearity check: does z_k = k * z_1 mod N / N' for some N'?
    z1 = zs[1]
    linear_hits = 0
    for k in range(2, min(g_ord, 40)):
        if zs[k] is None or z1 is None:
            continue
        if zs[k] == (k * z1) % N:
            linear_hits += 1

    # Bias check: pseudo-random or biased?
    # Compute the "residue": z_k - k*z1 mod N, see if it has small top bits
    residues = []
    for k in range(1, min(g_ord, 80)):
        if zs[k] is None or z1 is None:
            continue
        r = (zs[k] - k * z1) % N
        residues.append((k, r))
    max_res = max(r for _, r in residues)
    avg_res = sum(r for _, r in residues) / max(len(residues), 1)

    # HNP attack attempt. Build a lattice with rows
    #   [ N   0    0   ...  0 ]
    #   [ t_1 1    0   ...  0 ]   with t_i = z_{k_i} for sampled k_i
    #   [ t_2 0    1   ...  0 ]
    #   [ ... ]
    #   [ B   0    0   ...  1 ]   B = some target
    # and see if LLL finds a short vector related to the secret.
    random.seed(p)
    m = 8  # number of samples
    sample_ks = random.sample(range(2, min(g_ord, 100)), min(m, g_ord - 3))
    ts = [zs[k] for k in sample_ks]
    if any(t is None for t in ts):
        return {"prime": p, "status": "z computation failure"}

    # This is a naive lattice; LLL is not expected to find the secret.
    rows = []
    rows.append([N] + [0] * m)
    for i in range(m):
        row = [ts[i]] + [0] * m
        row[i + 1] = 1
        rows.append(row)
    reduced = lll(rows)
    shortest_norm_sq = min(sum(x * x for x in row) for row in reduced)

    return {
        "prime": p, "curve_b": b, "n": n, "g_ord": g_ord,
        "linear_hits": f"{linear_hits}/{min(g_ord, 40) - 2}",
        "residues_max": max_res,
        "residues_avg_ratio_to_N": f"{avg_res / N:.3f}",
        "lattice_shortest_norm_sq": shortest_norm_sq,
        "notes": "non-trivial short vector would indicate HNP-style bias",
    }


def main():
    candidates = [
        (13, 2), (19, 2), (31, 3), (43, 7),
        (67, 2), (73, 13), (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    for p, b in candidates:
        print(f"[phase6] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase6_p{p}.json").write_text(json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
