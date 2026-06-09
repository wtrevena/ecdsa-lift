"""
Phase 26 — Formal logarithm of δ.

The most important new experiment. Taking log_Ê(δ(k)) linearizes the
formal group, so if the non-linearity of δ is "cohomological", log
should reveal it.

Strategy:
  1. Work on E(Z/p^e) with e large enough that higher-order log
     terms matter (for j = 0 curves the log has terms at z, z^7,
     z^13, ..., so we need e ≥ 8).
  2. For k = 1..g_ord-1, compute the z-coordinate of δ(k) ∈ Ê.
  3. Evaluate log_Ê(z) via the truncated power series (coefficients
     from PARI).
  4. Test linear-in-k fit, polynomial fits, and look at successive
     differences of log(δ(k)).

Also computes log(c(jG, G)) where c is the Teichmüller 2-cocycle —
this is the "step function" whose partial sums give log(δ).
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

import cypari2
pari = cypari2.Pari()


def get_formal_log_coeffs(a, b, terms=60):
    """Return list of (exponent, rational) pairs for the formal log of
    y^2 = x^3 + a*x + b, up to x^{terms}. Rationals are (num, den)."""
    E = pari.ellinit([int(a), int(b)])
    flog = pari.ellformallog(E, terms + 1)
    # Parse the power series. flog is a PARI t_SER.
    # Convert to polynomial and iterate monomials.
    polytype = pari("t_POL")
    # Use the coefficient extraction via polcoef-like:
    # ellformallog returns a power series in x. We can get coefficients
    # with polcoef(truncate(...), i).
    series = pari.truncate(flog)
    out = []
    for i in range(terms + 1):
        c = pari.polcoef(series, i)
        if c == 0:
            continue
        # Each c is a rational (t_FRAC) or int
        num = int(pari.numerator(c))
        den = int(pari.denominator(c))
        out.append((i, num, den))
    return out


def log_in_ring(coeffs, z, N):
    """Evaluate log at z ∈ Z/N using the rational-coefficient power
    series. Requires every denominator to be a unit in Z/N (true when
    gcd(den, p) = 1, which holds for coeffs of the formal log of a
    short Weierstrass curve since all denominators are coprime to p
    for our test primes)."""
    result = 0
    for k, num, den in coeffs:
        # z^k mod N
        zk = pow(z, k, N)
        if zk == 0:
            continue
        try:
            inv_den = modinv(den, N)
        except ZeroDivisionError:
            continue
        result = (result + num * inv_den * zk) % N
    return result


def z_coord_of_delta(E_p, C, G, k):
    """Compute z = -X/Y of δ(k) = [k]τ(G) − τ(kG) in projective E(Z/p^e).
    Returns the z-coordinate as an integer in Z/N, or None."""
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kH = E_p.mul(k, G)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    tau_kH = C.from_affine(teichmuller_lift(E_p, E_aff, kH))
    diff = C.add(C.mul(k, tau_G), C.neg(tau_kH))
    X, Y, Z = diff
    try:
        # z = -X / Y  (in the formal group parameterization)
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


def linear_fit_mod(seq, N):
    """Fit seq[k] ≈ A*k + B mod N. Returns (A, B, hits) where hits
    is the count of exact matches."""
    if len(seq) < 2 or seq[0] is None or seq[1] is None:
        return None, None, 0
    A = (seq[1] - seq[0]) % N
    B = (seq[0] - A) % N
    hits = sum(1 for i, v in enumerate(seq)
               if v is not None and (A * i + B - v) % N == 0)
    return A, B, hits


def poly_fit_mod(seq, N, deg):
    """Fit a degree-deg polynomial to seq mod N via finite differences.
    Returns coefficient list and hit count."""
    # Use first deg+1 non-None values to fit via Newton's forward
    # difference. Then check full sequence.
    good = [(i, v) for i, v in enumerate(seq) if v is not None]
    if len(good) < deg + 1:
        return None, 0
    # Newton forward difference on first deg+1 points
    xs = [good[i][0] for i in range(deg + 1)]
    ys = [good[i][1] for i in range(deg + 1)]
    # Build finite differences assuming uniform spacing (not quite true
    # if some are None, but good enough for detection)
    if any(xs[i+1] - xs[i] != 1 for i in range(deg)):
        return None, 0
    table = [ys[:]]
    for d in range(deg):
        table.append([(table[d][i+1] - table[d][i]) % N for i in range(len(table[d]) - 1)])
    # Coefficients of p(x) expressed in the basis {(x-x0 choose k)}
    newton = [row[0] for row in table]
    x0 = xs[0]
    def eval_newton(x):
        acc = 0
        for k, c in enumerate(newton):
            # binomial(x - x0, k)
            b = 1
            for j in range(k):
                b = (b * (x - x0 - j)) // (j + 1)
            acc = (acc + c * b) % N
        return acc
    hits = sum(1 for i, v in enumerate(seq)
               if v is not None and eval_newton(i) == v)
    return newton, hits


def run(p, b, e=10, max_k=None):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 8:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    if max_k is None:
        max_k = min(g_ord - 1, 40)

    # Formal log series for this curve
    log_coeffs = get_formal_log_coeffs(0, b, terms=60)

    # Compute δ(k) z-coord and log(δ(k)) for k = 1..max_k
    z_vals = [None]  # index 0 unused; δ(0) = identity
    log_vals = [None]
    for k in range(1, max_k + 1):
        z = z_coord_of_delta(E_p, C, G, k)
        z_vals.append(z)
        if z is None or z == 0:
            log_vals.append(None)
        else:
            log_vals.append(log_in_ring(log_coeffs, z, N))

    # Analysis: linear fit
    valid = [(k, v) for k, v in enumerate(log_vals) if v is not None and k > 0]
    if len(valid) < 4:
        return {"prime": p, "status": "not enough samples"}

    # Linear fit (ignoring index 0)
    seq = [v for _, v in valid]
    A_lin, B_lin, hits_lin = linear_fit_mod(seq, N)

    # Polynomial fits up to degree 4
    poly_hits = {}
    for d in range(2, 5):
        _, h = poly_fit_mod(seq, N, d)
        poly_hits[d] = h

    # Successive differences of log(δ(k))
    diffs = [(seq[i+1] - seq[i]) % N for i in range(len(seq) - 1)]
    # Are the differences constant?
    const_diff = all(d == diffs[0] for d in diffs) if diffs else False
    # Linear fit of differences
    A_d, B_d, hits_d = linear_fit_mod(diffs, N) if len(diffs) >= 2 else (None, None, 0)

    # Valuations: v_p(log(δ(k)))
    def vp(x):
        if x == 0:
            return e  # saturate at e
        v = 0
        while x % p == 0:
            v += 1
            x //= p
        return v
    valuations = [vp(v) for v in seq]

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "max_k": max_k, "samples": len(valid),
        "log_coeffs_degrees": [c[0] for c in log_coeffs[:8]],
        "first_5_z": z_vals[1:6],
        "first_5_log": seq[:5],
        "valuations_histogram": {
            v: valuations.count(v) for v in sorted(set(valuations))},
        "linear_fit_hits": f"{hits_lin}/{len(seq)}",
        "poly_fit_hits": {d: f"{h}/{len(seq)}" for d, h in poly_hits.items()},
        "const_successive_diff": const_diff,
        "linear_fit_of_diffs_hits": f"{hits_d}/{len(diffs)}" if diffs else None,
        "verdict": (
            "LINEAR IN k — BREAK" if hits_lin == len(seq) else
            "QUADRATIC fit — strong structure" if poly_hits.get(2, 0) == len(seq) else
            "CUBIC fit — moderate structure" if poly_hits.get(3, 0) == len(seq) else
            "no low-degree polynomial structure"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase26] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase26_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
