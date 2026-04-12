"""
Phase 13c — Miller's algorithm tracked as (numerator, denominator)
fractions to avoid inversion failures in Z/p^e.

The vanilla Miller / phase13b aborted whenever a line or vertical
denominator was a non-unit mod p^e. We fix this by:

  (1) Keeping num and den as SEPARATE ring elements throughout —
      never invert until the end.
  (2) At the very end, compute v_p(num), v_p(den) and divide out the
      common p-power factor. In the Weil pairing f_{n,P}(Q+R)/f_{n,P}(R)
      these p-factors should cancel.
  (3) Representing line slopes as (slope_num, slope_den) pairs too, so
      the line evaluations `y - y1 - slope*(x - x1)` become
          (slope_den * (y - y1) - slope_num * (x - x1)) / slope_den.
      And we pull the extra slope_den factor into the denominator
      accumulator.

With these tracked fractions, Miller runs to completion in Z/p^e with
no inversions.
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv, points
from src.lifts import teichmuller_lift


def vp(x, p):
    if x == 0:
        return None  # infinite
    v = 0
    while x % p == 0:
        v += 1
        x //= p
    return v


def line_num_den(P1, P2, Q, a, N):
    """Return (numerator, denominator) for the line through P1, P2
    evaluated at Q, with NO inversions. The geometric value is
       num / den
    and at the end of Miller we divide out shared p-powers with the
    accumulator denominator."""
    x1, y1 = P1
    x2, y2 = P2
    xq, yq = Q
    if (x1 - x2) % N == 0 and (y1 + y2) % N == 0:
        # vertical line x = x1
        return ((xq - x1) % N, 1)
    if (x1 - x2) % N == 0:
        # tangent at P1 = P2
        s_num = (3 * x1 * x1 + a) % N
        s_den = (2 * y1) % N
    else:
        s_num = (y2 - y1) % N
        s_den = (x2 - x1) % N
    # line: y - y1 = (s_num/s_den)(x - x1)
    # value at Q: (s_den*(yq - y1) - s_num*(xq - x1)) / s_den
    num = (s_den * (yq - y1) - s_num * (xq - x1)) % N
    den = s_den % N
    return (num, den)


def miller(P, Q, n, E_lifted, N, a):
    bits = bin(n)[2:]
    T = P
    num = 1
    den = 1
    for bit in bits[1:]:
        # Doubling step
        ln_num, ln_den = line_num_den(T, T, Q, a, N)
        T2 = E_lifted.add(T, T)
        if T2 is None:
            # [2]T = ∞: vertical line = x - T.x; contributes as number 1
            # (we're essentially at a 2-torsion point of P, which for
            # prime n shouldn't happen until the last step)
            return None
        # f ← f² · ln / vert(T2)
        vnum = (Q[0] - T2[0]) % N
        vden = 1
        num = (num * num * ln_num * vden) % N
        den = (den * den * ln_den * vnum) % N
        T = T2
        if bit == '1':
            ln_num, ln_den = line_num_den(T, P, Q, a, N)
            T_new = E_lifted.add(T, P)
            if T_new is None:
                return None
            vnum = (Q[0] - T_new[0]) % N
            vden = 1
            num = (num * ln_num * vden) % N
            den = (den * ln_den * vnum) % N
            T = T_new
    return (num, den)


def reduce_fraction(num, den, p, e):
    """Strip common p-powers from (num, den) in Z/p^e.  Returns
    (num', den', shift) where shift = v_p(num) − v_p(den) after
    cancellation; num'/den' is the 'regularized' value."""
    N = p ** e
    vn = vp(num, p) if num != 0 else e
    vd = vp(den, p) if den != 0 else e
    s = min(vn, vd)
    if s == 0:
        return num, den
    ps = p ** s
    return (num // ps) % N, (den // ps) % N


def run(p, b, e=3):
    E_p = Curve(a=0, b=b, N=p)
    n_ord = naive_order(E_p)
    if n_ord % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}
    N = p ** e
    E_lifted = Curve(a=0, b=b, N=N)

    # Auxiliary R
    pts = [pt for pt in points(E_p)
           if pt is not None and pt != G and pt != E_p.neg(G)]
    if not pts:
        return {"prime": p, "status": "no aux point"}
    R = pts[len(pts) // 2]

    tau_G = teichmuller_lift(E_p, E_lifted, G)
    tau_R = teichmuller_lift(E_p, E_lifted, R)

    miller_vals = {}
    failures = 0
    nonunit_events = 0
    for k in range(2, min(g_ord, 25)):
        kG = E_p.mul(k, G)
        try:
            tau_kG = teichmuller_lift(E_p, E_lifted, kG)
            Q = E_lifted.add(tau_kG, tau_R)
            num_Q, den_Q = miller(tau_G, Q, g_ord, E_lifted, N, 0)
            num_R, den_R = miller(tau_G, tau_R, g_ord, E_lifted, N, 0)
        except Exception:
            failures += 1
            continue
        # f(Q)/f(R) = (num_Q * den_R) / (den_Q * num_R)
        top = (num_Q * den_R) % N
        bot = (den_Q * num_R) % N
        # Cancel common p-powers
        tn, td = reduce_fraction(top, bot, p, e)
        try:
            val = (tn * modinv(td, N)) % N
        except ZeroDivisionError:
            nonunit_events += 1
            # Record the "p-adic" value as (tn, td) debris
            continue
        miller_vals[k] = val

    if len(miller_vals) < 4:
        return {"prime": p, "status": f"miller usable {len(miller_vals)}",
                "failures": failures, "nonunit": nonunit_events}

    ks = sorted(miller_vals.keys())
    # Multiplicativity test: miller(k1)*miller(k2) vs miller(k1+k2)?
    # (The Weil/Miller pairing is bilinear, so for fixed G and varying Q=kG
    # the value should satisfy e(G, kG) = e(G, G)^k. Test that.)
    # Raise miller_vals[ks[0]] to (k / ks[0]) mod g_ord.
    base_k = ks[0]
    base_val = miller_vals[base_k]
    try:
        inv_base = pow(base_k, -1, g_ord)
    except ValueError:
        inv_base = None
    exp_hits = 0
    total = 0
    exp_samples = []
    if inv_base is not None:
        for k in ks:
            exp = (k * inv_base) % g_ord
            target = pow(base_val, exp, N)
            total += 1
            ok = (target == miller_vals[k])
            if ok:
                exp_hits += 1
            if len(exp_samples) < 4:
                exp_samples.append((k, ok))

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "miller_computed": len(miller_vals),
        "miller_failures": failures,
        "nonunit_after_cancel": nonunit_events,
        "exponential_rule_hits": f"{exp_hits}/{total}",
        "exp_samples": exp_samples,
        "first_4_values": [(k, miller_vals[k]) for k in ks[:4]],
    }


def main():
    candidates = [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5)]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in candidates:
        print(f"[phase13c] p={p} b={b}")
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
    (out_dir / "phase13c_all.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
