"""
Phase 8 — Can the Teichmüller error be fit by a low-degree polynomial?

Phase 2 showed that finite differences of δ(k) don't vanish at order
up to 4. That's a necessary condition for polynomial behavior, not a
sufficient one — δ could still be *well approximated* by a polynomial,
with enough slack that the finite differences don't look zero modulo
p^e but do look small compared to p^e. If so, Lagrange interpolation
on (k, δ(k)) pairs would give a polynomial whose evaluation at an
unknown k leaks bits of k.

We test three increasingly generous models:

  (i)   Z-linear:     δ(k)  ≡  A·k  (mod p^e)
  (ii)  Z-quadratic:  δ(k)  ≡  A·k + B·k²  (mod p^e)
  (iii) Z[ω]-linear:  δ(k)  ≡  A·k + B·(λ·k)  (mod p^e),
        where λ is the GLV eigenvalue. This is what j=0 curves might
        give via the CM endomorphism.

Each model is a linear system in the unknowns (A, B, ...). We solve
it using the first few (k, δ(k)) points and then test the fit on the
remaining ones. "Fits" means the predicted δ matches the actual δ
modulo p^e for every test point.

If any model fits globally, that model gives a closed-form attack.
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


def build_delta_series(E_p: Curve, C: ProjCurve, G, g_ord: int):
    """δ(k) = z-coordinate of ([k]·τ(G) - τ([k]G)) ∈ Ê, as a Z/p^e
    integer. This IS the information Phase 1 showed is injective in k."""
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    out = {}
    for k in range(1, g_ord):
        kG = E_p.mul(k, G)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
        X, Y, Z = diff
        if Y % C.N == 0:
            out[k] = None
        else:
            out[k] = (-X * modinv(Y % C.N, C.N)) % C.N
    return out


def fit_linear_mod_N(pairs, N):
    """Solve A·k ≡ δ(k) (mod N) using the first pair and test on all."""
    if not pairs:
        return None, 0
    k0, d0 = pairs[0]
    try:
        A = (d0 * modinv(k0 % N, N)) % N
    except ZeroDivisionError:
        return None, 0
    good = sum(1 for k, d in pairs if (A * k - d) % N == 0)
    return {"A": A}, good


def fit_quadratic_mod_N(pairs, N):
    """Solve A·k + B·k² ≡ δ(k) using two pairs; test all."""
    if len(pairs) < 2:
        return None, 0
    (k1, d1), (k2, d2) = pairs[0], pairs[1]
    # | k1  k1^2 | |A|   |d1|
    # | k2  k2^2 | |B| = |d2|
    det = (k1 * k2 * k2 - k2 * k1 * k1) % N
    try:
        det_inv = modinv(det, N)
    except ZeroDivisionError:
        return None, 0
    A = ((k2 * k2 * d1 - k1 * k1 * d2) * det_inv) % N
    B = ((k1 * d2 - k2 * d1) * det_inv) % N
    good = sum(1 for k, d in pairs
               if (A * k + B * k * k - d) % N == 0)
    return {"A": A, "B": B}, good


def fit_zomega_linear_mod_N(pairs, N, lam):
    """A·k + B·(λ k) ≡ δ(k) — degenerate since both terms are
    proportional to k. Returns only if consistent."""
    # This is actually equivalent to (A + B·λ)·k ≡ δ(k) which is
    # just the linear model in disguise. Include for completeness.
    return fit_linear_mod_N(pairs, N)


def fit_polynomial(pairs, N, degree):
    """Fit δ(k) = c_0 + c_1 k + c_2 k² + ... + c_d k^d mod N using
    the first (degree+1) pairs; test all."""
    if len(pairs) < degree + 1:
        return None, 0
    # Vandermonde system over Z/N
    rows = []
    rhs = []
    for i in range(degree + 1):
        k, d = pairs[i]
        rows.append([pow(k, j, N) for j in range(degree + 1)])
        rhs.append(d % N)
    # Gaussian elimination with explicit modular inverses — fails
    # cleanly on non-unit pivots.
    M = [row[:] + [rhs[i]] for i, row in enumerate(rows)]
    size = degree + 1
    for col in range(size):
        pivot_row = None
        for r in range(col, size):
            if M[r][col] % N != 0:
                try:
                    modinv(M[r][col] % N, N)
                    pivot_row = r
                    break
                except ZeroDivisionError:
                    continue
        if pivot_row is None:
            return None, 0
        M[col], M[pivot_row] = M[pivot_row], M[col]
        inv = modinv(M[col][col] % N, N)
        M[col] = [(x * inv) % N for x in M[col]]
        for r in range(size):
            if r == col:
                continue
            factor = M[r][col] % N
            if factor == 0:
                continue
            M[r] = [(M[r][c] - factor * M[col][c]) % N for c in range(size + 1)]
    coeffs = [M[i][size] for i in range(size)]
    # Test
    good = 0
    for k, d in pairs:
        v = sum(coeffs[i] * pow(k, i, N) for i in range(size)) % N
        if v == d % N:
            good += 1
    return {"coeffs": coeffs}, good


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 15:
        return {"prime": p, "status": "too small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    delta = build_delta_series(E_p, C, G, g_ord)
    pairs = [(k, d) for k, d in delta.items() if d is not None]
    total = len(pairs)

    # Degrees to try
    results = {}
    for degree in range(1, min(total, 8)):
        params, good = fit_polynomial(pairs, N, degree)
        results[f"degree_{degree}"] = {
            "fit_count": good,
            "total": total,
            "ratio": good / total if total else 0,
        }

    best = max(results.values(), key=lambda r: r["fit_count"])
    return {
        "prime": p,
        "curve_b": b,
        "g_ord": g_ord,
        "N": N,
        "total_scalars": total,
        "fits": results,
        "best_ratio": best["ratio"],
        "status": "LINEAR ATTACK POSSIBLE" if best["ratio"] > 0.95
                  else "pseudorandom under polynomial fit",
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    for p, b in candidates:
        print(f"[phase8] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase8_p{p}.json").write_text(json.dumps(report, indent=2, default=str))
        print(f"   {report.get('status')}")
        if "fits" in report:
            for deg, info in report["fits"].items():
                print(f"   {deg}: {info['fit_count']}/{info['total']}")


if __name__ == "__main__":
    main()
