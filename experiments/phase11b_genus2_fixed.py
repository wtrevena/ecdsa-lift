"""
Phase 11b — Fixed genus-2 Jacobian experiment.

Fixes to Phase 11:
  (1) Use Tonelli–Shanks for square roots (not linear scan).
  (2) Start from a literal C(F_p) point, Hensel-lift both coordinates
      to C(Z/p^e), and BUILD the Mumford (u, v) from the lifted point.
      This guarantees v² ≡ f (mod u) exactly over Z/p^e, so Cantor
      reduction never hits a non-exact remainder.
  (3) Stop after 20 successful compositions (enough to test linearity).

The question is still: does any Mumford coordinate of k·D carry a
linear-in-k signal?
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import modinv
from src import jac2 as J
from src import poly_modN as P


def tonelli_shanks(n, p):
    if n % p == 0:
        return 0
    if pow(n, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        r = pow(n, (p + 1) // 4, p)
        return r
    # Tonelli–Shanks for p ≡ 1 mod 4
    q = p - 1
    s = 0
    while q % 2 == 0:
        q //= 2
        s += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    while True:
        if t == 1:
            return r
        i = 0
        temp = t
        while temp != 1:
            temp = (temp * temp) % p
            i += 1
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p


def eval_poly(coeffs, x, N):
    r = 0
    for c in reversed(coeffs):
        r = (r * x + c) % N
    return r


def hensel_lift_y(x, y0, f_coeffs, p, e):
    """Given y0² ≡ f(x) (mod p), lift y to Z/p^e via Newton."""
    N = p ** e
    x = x % N
    target = eval_poly(f_coeffs, x, N)
    y = y0 % p
    pk = p
    while pk < N:
        pk = min(pk * pk, N)
        # y ← y − (y² − target)/(2y)  mod pk
        num = (y * y - target) % pk
        inv = modinv((2 * y) % pk, pk)
        y = (y - num * inv) % pk
    return y


def make_weight1_divisor(x_fp, y_fp, f_coeffs, p, e):
    """Build Mumford (u, v) for the divisor (P) − ∞ with P = (x, y)
    ∈ C(Z/p^e). u = x − x_lift, v = y_lift (constant)."""
    y_lift = hensel_lift_y(x_fp, y_fp, f_coeffs, p, e)
    N = p ** e
    u = [(-x_fp) % N, 1]
    v = [y_lift]
    return (u, v)


def find_point_on_C(f_coeffs, p):
    """Find an F_p point (x, y) on y² = f(x)."""
    for x in range(1, p):
        fx = eval_poly(f_coeffs, x, p)
        if fx == 0:
            continue
        y = tonelli_shanks(fx, p)
        if y is not None and y != 0:
            return (x, y)
    return None


def estimate_order(D, f_lifted, N, max_k=120):
    """Compute k·D until identity or max_k."""
    cur = D
    for k in range(1, max_k + 1):
        if J.is_identity(cur):
            return k
        try:
            cur = J.compose(cur, D, f_lifted, N)
        except Exception:
            return None
    return None


def run(p, b, e=2):
    N = p ** e
    f_p = [0, b, 0, 0, 0, 1]  # x^5 + b*x over F_p
    f_lifted = [c % N for c in f_p]

    pt = find_point_on_C(f_p, p)
    if pt is None:
        return {"prime": p, "status": "no point on C"}
    x, y = pt
    D = make_weight1_divisor(x, y, f_p, p, e)

    # Verify v^2 ≡ f (mod u) over Z/N
    u, v = D
    vsq = P.mul(v, v, N)
    diff = P.sub(f_lifted, vsq, N)
    _, rem = P.divmod_poly(diff, u, N)
    if any(c != 0 for c in rem):
        return {"prime": p, "status": f"lift inconsistent: rem={rem}"}

    # Compute k·D for k = 1..L
    L = 24
    results = []
    cur = D
    try:
        for k in range(1, L + 1):
            if k > 1:
                cur = J.compose(cur, D, f_lifted, N)
            if J.is_identity(cur):
                results.append((k, "identity"))
                break
            u_k, v_k = cur
            u_pad = list(u_k) + [0] * (3 - len(u_k))
            v_pad = list(v_k) + [0] * (2 - len(v_k))
            results.append((k, u_pad, v_pad))
    except Exception as exc:
        return {"prime": p, "status": f"compose at k={k} failed: {exc}",
                "partial": results}

    valid = [r for r in results if len(r) == 3]
    if len(valid) < 4:
        return {"prime": p, "status": "not enough valid samples",
                "results": results}

    features = [list(r[1]) + list(r[2]) for r in valid]
    Lf = len(features)
    # Test each feature dim for linear-in-k fit mod N
    linear_report = []
    for dim in range(5):
        seq = [features[i][dim] for i in range(Lf)]
        A = (seq[1] - seq[0]) % N
        B = (seq[0] - A) % N
        hits = sum(1 for i in range(Lf)
                   if (A * (i + 1) + B - seq[i]) % N == 0)
        linear_report.append({"dim": dim, "hits": f"{hits}/{Lf}",
                              "A": A, "B": B})

    return {
        "prime": p, "curve_b": b, "N": N,
        "initial_point": (x, y),
        "samples": Lf,
        "linear_fit_per_dim": linear_report,
        "first_4_mumford": valid[:4],
        "verdict": ("LINEAR LEAK found"
                    if any(int(r["hits"].split("/")[0]) >= Lf - 1
                           and r["A"] != 0 for r in linear_report)
                    else "no linear-in-k structure"),
    }


def main():
    candidates = [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5)]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in candidates:
        print(f"[phase11b] p={p} b={b}")
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
    (out_dir / "phase11b_all.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
