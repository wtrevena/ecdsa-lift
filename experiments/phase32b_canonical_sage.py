"""
Phase 32b — Serre-Tate canonical lift in SAGE.

Phase 32 (pure Python) hit precision overflow in projective coordinates
when scalar-multiplying kernel-of-reduction points. Sage's
EllipticCurve(Zp(p, prec)) handles p-adic point arithmetic correctly,
so we redo it there.

Approach: Newton iteration on f(Q) = π·Q − Q, where π = a + b·ω is the
Frobenius element of Z[ω] and ω lifts to Z_p via Hensel.

  f'(Q)|_{tangent} = (π_tangent − 1)  where π_t = a + b·ω_lift ∈ Z_p
  Newton step:  Q ← Q − [(π_t − 1)^{-1}]_int · (π·Q − Q)

The integer (π_t − 1)^{-1} mod p^prec is treated as an integer scalar.
This is correct because the formal-group action of u ∈ Z_p on Ê
agrees with integer scalar multiplication when u is reduced to its
representative in [0, p^prec).

Run:  ~/miniforge3/bin/conda run -n sage sage experiments/phase32b_canonical_sage.py
"""
from sage.all import EllipticCurve, GF, Zp, ZZ, Integer
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def find_pi(p, ap):
    disc = 12 * p - 3 * ap * ap
    if disc < 0:
        return None
    s = ZZ(disc).isqrt()
    if s * s != disc:
        return None
    for sign in (1, -1):
        num = 3 * ap + sign * s
        if num % 6 == 0:
            a = num // 6
            b = 2 * a - ap
            if a * a - a * b + b * b == p and 2 * a - b == ap:
                return (int(a), int(b))
    return None


def lift_point_naive(E_zp, P_fp, R, b_curve):
    if P_fp.is_zero():
        return E_zp(0)
    xf, yf = P_fp.xy()
    x_lift = R(ZZ(xf))
    y_lift = R(ZZ(yf))
    for _ in range(10):
        f = y_lift * y_lift - (x_lift ** 3 + R(b_curve))
        df = 2 * y_lift
        y_lift = y_lift - f / df
    return E_zp(x_lift, y_lift)


def lift_omega(p, R):
    omega_fp = None
    for w in range(2, p):
        if pow(w, 3, p) == 1 and w != 1:
            omega_fp = w
            break
    if omega_fp is None:
        return None
    w = R(ZZ(omega_fp))
    for _ in range(20):
        f = w * w + w + 1
        df = 2 * w + 1
        w = w - f / df
    return w


def cm_action_omega(P, omega_lift):
    """ω·(x, y) = (ω·x, y). Bypass strict validation since
    omega^3 = 1 only to precision."""
    if P.is_zero():
        return P
    x, y = P[0], P[1]
    E = P.curve()
    return E.point([omega_lift * x, y, x.parent().one()], check=False)


def pi_action(P, a_pi, b_pi, omega_lift):
    return Integer(a_pi) * P + Integer(b_pi) * cm_action_omega(P, omega_lift)


def canonical_lift_newton(P_fp, E_zp, R, a_pi, b_pi, omega_lift, b_curve,
                          max_iter=20):
    """Newton iteration on f(Q) = π·Q − Q with tangent-space
    inverse (π_t − 1)^{-1} as the Newton scaling."""
    Q = lift_point_naive(E_zp, P_fp, R, b_curve)
    pi_t = R(a_pi) + R(b_pi) * omega_lift
    pim1_t = pi_t - R(1)
    inv_pim1 = R(1) / pim1_t
    # Convert to integer scalar
    inv_int = ZZ(inv_pim1.lift())
    for _ in range(max_iter):
        piQ = pi_action(Q, a_pi, b_pi, omega_lift)
        Rpt = piQ - Q
        if Rpt.is_zero():
            break
        # If Rpt has very high precision-level zero, also break
        try:
            v = Rpt[2].valuation()
            if v >= R.precision_cap() - 1:
                break
        except Exception:
            pass
        correction = (-inv_int) * Rpt
        Q_new = Q + correction
        if Q_new == Q:
            break
        Q = Q_new
    return Q


def z_of(P):
    if P.is_zero():
        return None
    return -P[0] / P[1]


def run(p, b_curve, prec=10, max_k=12):
    Fp = GF(p)
    E_fp = EllipticCurve(Fp, [0, b_curve])
    n = E_fp.cardinality()
    if n % p == 0:
        return {"p": p, "status": "anomalous"}
    if p % 3 != 1:
        return {"p": p, "status": "not split"}
    G_fp = None
    for P in E_fp.points():
        if P.order() == n:
            G_fp = P
            break
    if G_fp is None:
        return {"p": p, "status": "no cyclic gen"}
    ap = p + 1 - int(n)
    pi = find_pi(p, ap)
    if pi is None:
        return {"p": p, "status": f"no pi for ap={ap}"}
    a_pi, b_pi = pi

    R = Zp(p, prec=prec, type='capped-rel', print_mode='terse')
    E_zp = EllipticCurve(R, [0, b_curve])
    omega_lift = lift_omega(p, R)
    if omega_lift is None:
        return {"p": p, "status": "no omega"}

    def naive(P):
        return lift_point_naive(E_zp, P, R, b_curve)

    def canon(P):
        return canonical_lift_newton(P, E_zp, R, a_pi, b_pi, omega_lift,
                                     b_curve)

    # Verify canonical lift is Frobenius-fixed
    tau_G = canon(G_fp)
    pi_tauG = pi_action(tau_G, a_pi, b_pi, omega_lift)
    R_check = pi_tauG - tau_G
    if R_check.is_zero():
        pi_fix_val = "exact"
    else:
        try:
            v = R_check[2].valuation()
            pi_fix_val = f"v_p={v}/{prec}"
        except Exception:
            pi_fix_val = "nontrivial"

    delta_can_z = []
    delta_naive_z = []
    for k in range(1, max_k + 1):
        kG = Integer(k) * G_fp
        d_can = Integer(k) * canon(G_fp) - canon(kG)
        d_naive = Integer(k) * naive(G_fp) - naive(kG)
        delta_can_z.append(z_of(d_can))
        delta_naive_z.append(z_of(d_naive))

    can_zero = sum(1 for z in delta_can_z if z is None or z == 0
                   or (hasattr(z, 'valuation') and z.valuation() >= prec - 2))
    naive_zero = sum(1 for z in delta_naive_z if z is None or z == 0
                     or (hasattr(z, 'valuation') and z.valuation() >= prec - 2))

    return {
        "p": p, "b": b_curve, "n": int(n), "ap": ap, "pi": pi, "prec": prec,
        "max_k": max_k,
        "pi_fixes_tau_G": pi_fix_val,
        "canonical_zeros": f"{can_zero}/{max_k}",
        "naive_zeros": f"{naive_zero}/{max_k}",
        "first_3_canonical_z": [str(z) for z in delta_can_z[:3]],
        "first_3_naive_z": [str(z) for z in delta_naive_z[:3]],
    }


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5),
                 (103, 5)]:
        print(f"[phase32b] p={p} b={b}")
        try:
            r = run(p, b, prec=10, max_k=12)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": str(exc)}
        reports.append(r)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase32b_canonical_sage.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
