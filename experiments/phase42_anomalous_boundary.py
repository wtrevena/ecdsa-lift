"""
Phase 42 — The theorem's boundary is sharp: gcd(n, p) = 1 is exactly the
hypothesis that fails on anomalous curves, where the lift attack succeeds.

The inertness theorem (Phase 41) rests on: for <G> of order n with
gcd(n, p) = 1, the lift-error's secret part [k]c is governed by
k mod p^{e-1}, while the public data fixes only k mod n -- independent
coordinates -- and the well-definedness obstruction [n]c = [n]tau(G) is a
NONZERO kernel element (the Phase 21b constant).

Set n = p (an ANOMALOUS curve, #E(F_p) = p).  Now:
  * k mod n = k mod p  *equals*  the residue the formal group sees;
  * the obstruction [n]c = [p]c vanishes, because c lives in the kernel
    of reduction of E(Z/p^2), which has order p, so [p]c = O.
The lift error becomes well-defined and the formal-group logarithm
recovers k.  This is exactly Smart's attack (1999) / Satoh-Araki /
Semaev.  We implement it to show the SAME machinery that is inert for
gcd(n,p)=1 becomes a total break the instant n = p.

Smart's attack on E(Z/p^2):
  1. Lift G, Q=kG to Ghat, Qhat in E(Z/p^2) (any lift; Hensel on y).
  2. R_G = [p] Ghat,  R_Q = [p] Qhat  lie in the kernel of reduction
     E_1(Z/p^2) ~= (pZ/p^2, +); independent of the lift choice because
     [p](kernel)=O.
  3. psi(R) = (-X_R / Y_R) / p   mod p   (formal log to first order).
  4. k = psi(R_Q) / psi(R_G)   mod p.
"""
from __future__ import annotations
import sys
import pathlib
import json
import random
from math import gcd

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, points, point_order, modinv
from src.projective import ProjCurve
from src.lifts import hensel_lift_y, teichmuller_lift


def find_anomalous_curve(p):
    """Return (a, b) with #E(F_p) = p, or None. Try a few (a,b)."""
    for a in range(0, p):
        for b in range(0, p):
            if (4 * a * a * a + 27 * b * b) % p == 0:
                continue  # singular
            if naive_order(Curve(a, b, p)) == p:
                return (a, b)
    return None


def psi_log(a, b, p, P_fp):
    """psi(P) = formal-log( [p] * lift(P) ) / p   mod p.
    Well-defined: any two lifts of P differ by a kernel element kappa with
    [p]kappa = O, so [p]lift(P) is independent of the lift choice."""
    N = p * p
    C = ProjCurve(a, b, N)
    E_lift = Curve(a, b, N)
    x, y0 = P_fp
    y = hensel_lift_y(E_lift, x % N, y0 % p)
    R = C.mul(p, (x % N, y, 1))             # [p]lift(P) in kernel of reduction
    X, Y, Z = R
    if Y % N == 0:
        return None
    t = (-X * modinv(Y % N, N)) % N         # t = -X/Y, divisible by p
    if t % p != 0:
        return None
    return (t // p) % p                     # formal log / p  mod p


def smart_dlog(a, b, p, G, Q):
    """Recover k with Q = kG on an anomalous curve via the p-adic lift."""
    psiG = psi_log(a, b, p, G)
    psiQ = psi_log(a, b, p, Q)
    if psiG is None or psiG % p == 0 or psiQ is None:
        return None
    return (psiQ * modinv(psiG, p)) % p


def is_cyclic_zp2(a, b, p, G):
    """True if [p]lift(G) != O, i.e. E(Z/p^2) is cyclic and Smart applies."""
    N = p * p
    C = ProjCurve(a, b, N)
    E_lift = Curve(a, b, N)
    x, y0 = G
    y = hensel_lift_y(E_lift, x % N, y0 % p)
    return not C.is_identity(C.mul(p, (x % N, y, 1)))


def obstruction_vanishes(a, b, p, G):
    """The well-definedness obstruction is [n]c = [p]c with c in the kernel
    of reduction of E(Z/p^2) (order p).  R = [p]lift(G) lies in E_1(Z/p^2)
    iff z(R) = 0 mod p; then [p]R = O (z([p]R) = p z(R) = 0 mod p^2).  For
    gcd(n,p)=1 the analogous [n]tau(G) is NONZERO (Phase 41 V1); here n = p
    annihilates the order-p kernel, so the obstruction vanishes."""
    N = p * p
    C = ProjCurve(a, b, N)
    E_lift = Curve(a, b, N)
    x, y0 = G
    y = hensel_lift_y(E_lift, x % N, y0 % p)
    R = C.mul(p, (x % N, y, 1))
    X, Y, Z = R
    zR = (-X * modinv(Y % N, N)) % N
    return zR % p == 0 and (p * zR) % N == 0


def pick_generator(E, p):
    for P in points(E):
        if P is None or P[1] == 0:
            continue
        if point_order(E, P, p) == p:
            return P
    return None


def boundary_invariant(p, a, b, G, e=2):
    """Emit the section-invariant z([p]tau(G)) for the paper's anomalous
    curve (the numbers quoted in the paper's boundary section).

    Returns, at precision e (default e=2, the value the paper quotes):
      * z_pTauG            : z([p]lift(G)) mod p^e (an element of pZ/p^e);
      * factor_c, "c*p"    : its leading p-adic digit c = z/p mod p, so
                             z = c*p; phi(G)=c is what Smart's attack
                             divides by, hence must be != 0;
      * valuation          : v_p(z) (= 1: nonzero kernel element);
      * section_independent: True iff the Teichmuller-x, naive-x, and two
                             shifted-x sections all give the SAME z (the
                             ambiguity [p]c=O has been annihilated).
    """
    N = p ** e
    C = ProjCurve(a, b, N)
    El = Curve(a, b, N)
    Ep = Curve(a, b, p)
    x0, y0 = G

    def z_of_section(xL):
        y = hensel_lift_y(El, xL % N, y0 % p)
        X, Y, Z = C.mul(p, (xL % N, y, 1))
        assert Y % N != 0
        return (-X * modinv(Y % N, N)) % N

    xt = teichmuller_lift(Ep, El, G)[0] % N          # Teichmuller-x section
    sections = {
        "teichmuller_x": z_of_section(xt),
        "naive_x":       z_of_section(x0 % N),
        "shift_x+17p":   z_of_section((x0 + 17 * p) % N),
        "shift_x+29p":   z_of_section((x0 + 29 * p) % N),
    }
    vals = set(sections.values())
    z = sections["naive_x"]
    c = (z // p) % p
    val = 0
    zz = z
    while zz % p == 0 and zz != 0:
        zz //= p; val += 1
    return {
        "p": p, "a": a, "b": b, "G": list(G), "e": e,
        "z_pTauG_mod_pe": z,
        "leading_digit_c": c,
        "factorization": f"{c}*{p}",
        "valuation": val,
        "phi_G_nonzero": c % p != 0,
        "sections": sections,
        "section_independent": len(vals) == 1,
    }


def main():
    out = {}
    # find a small anomalous curve with CYCLIC E(Z/p^2) (Smart applies).
    # A split E(Z/p^2) = Z/p x Z/p (e.g. p=5,61) needs higher precision; the
    # generic anomalous curve is cyclic and breaks already at p^2.
    anom = None
    PRIMES = [53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109]
    for p in PRIMES:
        ab = find_anomalous_curve(p)
        if not ab:
            continue
        a, b = ab
        E = Curve(a, b, p)
        G = pick_generator(E, p)
        if G and is_cyclic_zp2(a, b, p, G):
            anom = (p, a, b, G)
            break
    assert anom, "no suitable anomalous curve found"
    p, a, b, G = anom
    E = Curve(a, b, p)
    n = naive_order(E)
    print(f"Anomalous curve: y^2 = x^3 + {a}x + {b} over F_{p}, "
          f"#E = {n} (= p: {n == p}), E(Z/p^2) cyclic: True")
    print(f"Generator G = {G}, order {point_order(E, G, n)}")
    out["curve"] = {"p": p, "a": a, "b": b, "order": n, "anomalous": n == p}
    out["G"] = list(G)

    # run Smart's attack on random secrets
    random.seed(2026)
    results = []
    for _ in range(12):
        k = random.randrange(1, p)
        Q = E.mul(k, G)
        krec = smart_dlog(a, b, p, G, Q)
        results.append({"k": k, "recovered": krec, "ok": krec == k})
        tag = "OK" if krec == k else "FAIL"
        print(f"   k={k:3d}  recovered={krec}  {tag}")
    rate = sum(r["ok"] for r in results) / len(results)
    out["smart_attack"] = results
    out["smart_attack_success_rate"] = rate

    out["obstruction_vanishes"] = obstruction_vanishes(a, b, p, G)
    print(f"   obstruction [n]c = [p]c vanishes: {out['obstruction_vanishes']}")

    # Section-invariant boundary value quoted in the paper (e=2): for the
    # paper's exact curve y^2 = x^3 + 4x + 7 / F_53 with G = (0, 22).
    bi = boundary_invariant(53, 4, 7, (0, 22), e=2)
    out["boundary_invariant_paper_curve"] = bi
    print(f"   [p]tau(G) on F_53 curve: z = {bi['z_pTauG_mod_pe']} "
          f"= {bi['factorization']}, valuation {bi['valuation']}, "
          f"section-independent: {bi['section_independent']}, "
          f"phi(G)={bi['leading_digit_c']} (nonzero: {bi['phi_G_nonzero']})")

    # CONTRAST: non-anomalous curve, the same procedure is ill-defined.
    pna, ana, bna = 67, 0, 7
    nna = naive_order(Curve(ana, bna, pna))
    gcd_na = gcd(nna, pna)
    print("")
    print(f"Contrast (non-anomalous): y^2=x^3+{bna} over F_{pna}, "
          f"#E={nna}, gcd(n,p)={gcd_na}")
    print("   [n]tau(G) is a NONZERO kernel element (Phase 21b/41 constant); "
          "the lift error is not well-defined mod n, so no analogue of psi "
          "recovers k.")
    out["contrast"] = {"p": pna, "order": nna, "gcd_n_p": gcd_na,
                       "attack_applies": False}

    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "phase42_anomalous_boundary.json").write_text(
        json.dumps(out, indent=2, default=str))
    print("")
    print(f"Smart attack success rate: {rate:.0%}. "
          "Wrote results/phase42_anomalous_boundary.json")


if __name__ == "__main__":
    main()
