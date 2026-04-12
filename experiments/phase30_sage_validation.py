"""
Phase 30 — Sage validation of the three load-bearing claims.

Re-runs the most important negative results using SageMath's
battle-tested implementations:

  (A) Phase 21b z-additive antisymmetry on E(Z_p) using Sage's
      EllipticCurve over Zp(p, prec).
  (B) Phase 26 formal log of δ using Sage's E.formal_group().log()
      (which gives the series symbolically — no manual coefficient
      extraction).
  (C) The leading-digit Aut-orbit rank rule on j = 1728 curves
      (which we never verified in earlier phases because we couldn't
      easily find prime-order curves).

Run with:  ~/miniforge3/bin/conda run -n sage sage experiments/phase30_sage_validation.py
"""
from sage.all import (EllipticCurve, GF, Zp, ZZ, QQ, Integer, Mod,
                      next_prime, is_prime, randint, set_random_seed)
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


# ----------------------------------------------------------------------
# (A) Sage Phase 21b — z-additive antisymmetry on E(Z_p)
# ----------------------------------------------------------------------
def sage_phase21b(p, b, prec=6):
    """Verify z(δ(k)) + z(δ(n−k)) ≡ z([n]τ(G))  in Z_p / p^prec, where
    z(P) = -X/Y is the formal-group parameter on Ê. Uses Sage's lift
    to Zp."""
    Fp = GF(p)
    E_fp = EllipticCurve(Fp, [0, b])
    n = E_fp.cardinality()
    # find a generator
    G_fp = None
    for P in E_fp.points():
        if P.order() == n:
            G_fp = P
            break
    if G_fp is None:
        # no cyclic generator → take any maximal-order point
        G_fp = E_fp.gens()[0]
        n = G_fp.order()

    R = Zp(p, prec=prec, type='capped-rel', print_mode='terse')
    E_zp = EllipticCurve(R, [0, b])

    def lift_pt(Q):
        if Q.is_zero():
            return E_zp(0)
        x_fp, y_fp = Q.xy()
        # Hensel-lift y from y_fp ∈ F_p
        x_lift = R(ZZ(x_fp))
        # solve y² = x³ + b in Z_p starting from y_fp mod p
        y_lift = R(ZZ(y_fp))
        # Newton refinement
        for _ in range(prec + 2):
            f = y_lift * y_lift - (x_lift ** 3 + R(b))
            df = 2 * y_lift
            y_lift = y_lift - f / df
        return E_zp(x_lift, y_lift)

    tau_G = lift_pt(G_fp)
    nTauG = Integer(n) * tau_G
    # z = -X/Y on the formal group
    def z_of(P):
        if P.is_zero():
            return R(0)
        return -P[0] / P[1]

    const_z = z_of(nTauG)

    samples = []
    hits = 0
    total = 0
    for k in [2, 3, 5, 7, 11, 13, 17, 19, 23]:
        if k >= n - 2:
            continue
        kG_fp = Integer(k) * G_fp
        nkG_fp = Integer(n - k) * G_fp
        d_k = Integer(k) * tau_G - lift_pt(kG_fp)
        d_nk = Integer(n - k) * tau_G - lift_pt(nkG_fp)
        s = z_of(d_k) + z_of(d_nk)
        ok = (s - const_z).valuation() >= prec - 1
        total += 1
        if ok:
            hits += 1
        samples.append((k, bool(ok)))
    return {
        "prime": p, "b": b, "n": int(n), "prec": prec,
        "const_z_first_digits": str(const_z),
        "hits": f"{hits}/{total}",
        "samples": samples,
    }


# ----------------------------------------------------------------------
# (B) Sage Phase 26 — formal log of δ via E.formal_group().log()
# ----------------------------------------------------------------------
def sage_phase26(p, b, prec=8):
    """Compute z(δ(k)) for k = 1..g_ord-1 and apply Sage's symbolic
    formal-group logarithm. Test linear and polynomial fits in k."""
    Fp = GF(p)
    E_fp = EllipticCurve(Fp, [0, b])
    n = E_fp.cardinality()
    G_fp = None
    for P in E_fp.points():
        if P.order() == n:
            G_fp = P
            break
    if G_fp is None:
        return {"prime": p, "status": "no cyclic gen"}

    R = Zp(p, prec=prec, type='capped-rel', print_mode='terse')
    E_zp = EllipticCurve(R, [0, b])

    def lift_pt(Q):
        if Q.is_zero():
            return E_zp(0)
        x_fp, y_fp = Q.xy()
        x_lift = R(ZZ(x_fp))
        y_lift = R(ZZ(y_fp))
        for _ in range(prec + 2):
            f = y_lift * y_lift - (x_lift ** 3 + R(b))
            df = 2 * y_lift
            y_lift = y_lift - f / df
        return E_zp(x_lift, y_lift)

    tau_G = lift_pt(G_fp)

    # Sage's formal group logarithm as a power series
    fg = E_zp.formal_group()
    log_series = fg.log(prec * 2)  # truncated to degree 2*prec

    def z_of(P):
        if P.is_zero():
            return R(0)
        return -P[0] / P[1]

    log_vals = [None]
    for k in range(1, min(n, 30)):
        kG_fp = Integer(k) * G_fp
        d_k = Integer(k) * tau_G - lift_pt(kG_fp)
        z = z_of(d_k)
        if z == 0:
            log_vals.append(R(0))
        else:
            log_vals.append(log_series(z))

    # Linear fit on log_vals[1..]
    seq = log_vals[1:]
    valid = [(i + 1, v) for i, v in enumerate(seq) if v is not None]
    if len(valid) < 4:
        return {"prime": p, "b": b, "status": "few samples"}

    # Try a*k + b fit using the first two
    k1, v1 = valid[0]
    k2, v2 = valid[1]
    try:
        slope = (v2 - v1) / (k2 - k1)
        intercept = v1 - slope * k1
        hits = 0
        for k, v in valid:
            pred = slope * k + intercept
            if (v - pred).valuation() >= prec - 1:
                hits += 1
    except Exception as e:
        slope = None
        hits = 0

    # Quadratic fit using first three
    quad_hits = 0
    if len(valid) >= 3:
        k1, v1 = valid[0]; k2, v2 = valid[1]; k3, v3 = valid[2]
        # solve [a*k² + b*k + c = v] for three points
        try:
            from sage.all import matrix, vector
            M = matrix(R, [[k1*k1, k1, 1], [k2*k2, k2, 1], [k3*k3, k3, 1]])
            rhs = vector(R, [v1, v2, v3])
            coeffs = M.solve_right(rhs)
            for k, v in valid:
                pred = coeffs[0]*k*k + coeffs[1]*k + coeffs[2]
                if (v - pred).valuation() >= prec - 1:
                    quad_hits += 1
        except Exception:
            pass

    return {
        "prime": p, "b": b, "n": int(n), "prec": prec,
        "samples": len(valid),
        "linear_fit_hits": f"{hits}/{len(valid)}",
        "quadratic_fit_hits": f"{quad_hits}/{len(valid)}",
        "first_5_log_vals": [str(v) for v in seq[:5]],
    }


# ----------------------------------------------------------------------
# (C) j = 1728 rank verification
# ----------------------------------------------------------------------
def find_prime_order_j1728(p_min, p_max):
    """Find a prime p and curve y² = x³ + a*x with prime |E(F_p)|."""
    p = next_prime(p_min)
    while p < p_max:
        for a in range(1, p):
            E = EllipticCurve(GF(p), [a, 0])
            n = E.cardinality()
            if is_prime(n) and n > 30:
                return (p, int(a), int(n))
        p = next_prime(p)
    return None


def sage_phase14_j1728(p, a, n, e=3, size=14):
    """Compute the leading-digit rank matrix for a j=1728 curve and
    verify rank = (n-1)/4 (μ₄ Aut equivariance)."""
    Fp = GF(p)
    E_fp = EllipticCurve(Fp, [a, 0])
    G_fp = None
    for P in E_fp.points():
        if P.order() == n:
            G_fp = P
            break
    R = Zp(p, prec=e + 2, type='capped-rel', print_mode='terse')
    E_zp = EllipticCurve(R, [a, 0])

    def lift_pt(Q):
        if Q.is_zero():
            return E_zp(0)
        xf, yf = Q.xy()
        xl = R(ZZ(xf))
        yl = R(ZZ(yf))
        for _ in range(e + 4):
            f = yl * yl - (xl ** 3 + R(a) * xl)
            df = 2 * yl
            yl = yl - f / df
        return E_zp(xl, yl)

    def delta_leading(H_fp, k):
        kH_fp = Integer(k) * H_fp
        tau_H = lift_pt(H_fp)
        tau_kH = lift_pt(kH_fp)
        d = Integer(k) * tau_H - tau_kH
        if d.is_zero():
            return 0
        z = -d[0] / d[1]
        # leading digit (after stripping one factor of p)
        zl = z.lift()
        return (zl // p) % p

    M = []
    for i in range(1, size + 1):
        H = Integer(i) * G_fp
        row = []
        for k in range(2, size + 2):
            row.append(delta_leading(H, k))
        M.append(row)

    # F_p rank
    M2 = [r[:] for r in M]
    rk = 0
    col = 0
    while rk < size and col < size:
        piv = None
        for r in range(rk, size):
            if M2[r][col] % p != 0:
                piv = r
                break
        if piv is None:
            col += 1
            continue
        M2[rk], M2[piv] = M2[piv], M2[rk]
        inv = pow(M2[rk][col] % p, -1, p)
        M2[rk] = [(x * inv) % p for x in M2[rk]]
        for r in range(size):
            if r == rk:
                continue
            f = M2[r][col] % p
            if f == 0:
                continue
            M2[r] = [(M2[r][c] - f * M2[rk][c]) % p for c in range(size)]
        rk += 1
        col += 1

    return {
        "prime": p, "a": a, "n": n, "size": size,
        "rank": rk,
        "expected_mu4": min(size, (n - 1) // 4),
        "expected_pm1": min(size, (n - 1) // 2),
        "matches_mu4": rk == min(size, (n - 1) // 4),
    }


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    report = {}

    print("=== (A) Sage Phase 21b ===")
    p21b = []
    for p, b in [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3)]:
        try:
            r = sage_phase21b(p, b, prec=6)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"prime": p, "b": b, "error": str(exc)}
        print(f"  p={p} b={b}:  {r.get('hits', r.get('error'))}")
        p21b.append(r)
    report["phase21b_sage"] = p21b

    print("\n=== (B) Sage Phase 26 ===")
    p26 = []
    for p, b in [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3)]:
        try:
            r = sage_phase26(p, b, prec=8)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"prime": p, "b": b, "error": str(exc)}
        print(f"  p={p} b={b}:  linear={r.get('linear_fit_hits')}  quad={r.get('quadratic_fit_hits')}")
        p26.append(r)
    report["phase26_sage"] = p26

    print("\n=== (C) j=1728 rank check ===")
    j1728 = []
    found = find_prime_order_j1728(30, 5000)
    if found is None:
        print("  no prime-order j=1728 curve found in range")
    else:
        p, a, n = found
        print(f"  found: p={p} a={a} n={n}")
        try:
            r = sage_phase14_j1728(p, a, n)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"error": str(exc)}
        print(f"  result: {r}")
        j1728.append(r)
    report["j1728_check"] = j1728

    (out_dir / "phase30_sage_validation.json").write_text(
        json.dumps(report, indent=2, default=str))
    print(f"\nWritten {out_dir / 'phase30_sage_validation.json'}")


if __name__ == "__main__":
    main()
