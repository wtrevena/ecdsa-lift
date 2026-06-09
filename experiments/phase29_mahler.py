"""
Phase 29 — Mahler expansion of δ(k).

Any continuous f: Z_p → Z_p has a unique Mahler expansion
        f(x) = Σ_{n≥0} a_n · C(x, n),        |a_n|_p → 0,
where C(x, n) = x(x-1)...(x-n+1)/n! is the binomial polynomial.
The coefficients are recovered from finite differences:
        a_n = Σ_{k=0}^n (-1)^{n-k} C(n,k) · f(k) = (Δ^n f)(0).

Question: does k ↦ z(δ(k)), viewed as a function Z_p → Z_p (extended
from k = 0..g_ord−1), have a Mahler expansion with any exploitable
structure?

Specific probes:
  (1) Distribution of p-adic valuations v_p(a_n): does it grow fast
      (coefficients die, meaning the function is "smooth" p-adically)
      or slow (meaning the function is p-adically "rough")?
  (2) Support pattern: are the nonzero a_n concentrated at indices
      with some arithmetic property (e.g., n ≡ 0 mod 6, echoing the
      μ₆ symmetry of j=0)?
  (3) Do the Mahler coefficients satisfy any obvious recurrence?
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


def z_delta(E_p, C, G, k):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kG = E_p.mul(k, G)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
    diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
    X, Y, Z = diff
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


def finite_differences(seq, N):
    """Return the Newton forward-difference triangle. table[n][0] are
    the Mahler coefficients a_n = (Δ^n f)(0)."""
    rows = [list(seq)]
    while len(rows[-1]) > 1:
        prev = rows[-1]
        rows.append([(prev[i+1] - prev[i]) % N for i in range(len(prev) - 1)])
    return [row[0] for row in rows]


def vp(x, p, cap):
    if x % (p ** cap) == 0:
        return cap
    v = 0
    while x % p == 0:
        v += 1
        x //= p
    return v


def run(p, b, e=5):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 12:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Sample δ(k) for k = 0..g_ord-1 (δ(0) = 0)
    f = [0]
    for k in range(1, g_ord):
        z = z_delta(E_p, C, G, k)
        f.append(z if z is not None else 0)

    # Mahler coefficients via finite differences on the FULL sequence
    mahler = finite_differences(f, N)
    # Trim trailing zeros modulo p^e
    # Valuations histogram
    vals = [vp(a, p, e) for a in mahler]
    # Growth: does v_p(a_n) grow linearly with n? (It should if
    # f is locally analytic; bounded if f is merely continuous.)
    # We compute the minimum index with v_p(a_n) >= t for various t.
    first_hit = {t: None for t in range(1, e + 1)}
    for i, v in enumerate(vals):
        for t in range(1, e + 1):
            if first_hit[t] is None and v >= t:
                first_hit[t] = i

    # Support pattern: for each residue class mod 6, count nonzero a_n
    mod6_nonzero = {r: 0 for r in range(6)}
    for i, a in enumerate(mahler):
        if a != 0:
            mod6_nonzero[i % 6] += 1

    # Nonzero count
    nonzero = sum(1 for a in mahler if a != 0)

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "num_coeffs": len(mahler),
        "nonzero_coeffs": nonzero,
        "valuations_histogram": {v: vals.count(v) for v in sorted(set(vals))},
        "first_index_with_vp_ge": first_hit,
        "nonzero_by_residue_mod_6": mod6_nonzero,
        "first_8_valuations": vals[:8],
        "last_8_valuations": vals[-8:],
    }


def main():
    candidates = [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5)]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in candidates:
        print(f"[phase29] p={p} b={b}")
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
    (out_dir / "phase29_mahler.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
