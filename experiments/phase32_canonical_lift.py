"""
Phase 32 — δ under the TRUE Serre–Tate canonical lift (j = 0 curves).

Setup. For an ordinary elliptic curve E/F_p with j = 0 (so CM by
Z[ζ₃]) and p ≡ 1 mod 3, the prime p splits in Z[ω] = Z[ζ₃] as
p = π · π̄ for some Frobenius element π = a + bω with norm
N(π) = a² − ab + b² = p.  The integral model y² = x³ + b is already
the canonical lift of E to Z_p (because the CM locus is rigid under
deformation).  For each F_p point P, the canonical-lift section
τ_can(P) ∈ E(Z_p) is the UNIQUE point reducing to P that is fixed by
the lifted Frobenius:

        π · τ_can(P)  =  τ_can(P).

This is a *different* section than our naive coordinate-fixed Hensel
lift τ_naive (which fixes x but doesn't enforce Frobenius equivariance).

Algorithm.
  1. Lift ω from F_p to Z_p via Hensel: solve x² + x + 1 = 0.
     The CM action is then ω · (x, y) = (ω · x, y) on the j = 0 model.
  2. Solve N(π) = p in Z[ω] for π = a + bω; pick the right sign by
     matching the trace a_p of Frobenius (Tr(π) = 2a − b = a_p).
  3. For each F_p point P:
        Q ← naive Hensel lift of P
        repeat several times:
            R = π · Q − Q                  (an element of Ê)
            Q = Q − formal_inv_pi_minus_1(R)
     where π − 1 acts on Ê via the formal group law and is invertible
     because (π − 1)(π̄ − 1) = N(π − 1) = p + 1 − Tr(π) = #E(F_p),
     which is coprime to p.
  4. Compute δ_can(k) = [k] τ_can(G) − τ_can(kG) and rerun all the
     structural probes (linear fit, polynomial fit, DFT, BM,
     Mahler) on the canonical-lift sequence.

If δ_can has any structure that δ_naive does not, we have finally
found the discrepancy that hides the cryptanalytic lever.  If they
agree (or both are random), we definitively close the canonical-lift
angle.
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


def hensel_lift_root(poly_coeffs_mod_p, root_fp, p, e):
    """Hensel-lift a simple root of a polynomial from F_p to Z/p^e."""
    N = p ** e
    coeffs = [c % N for c in poly_coeffs_mod_p]
    deriv = [(i * c) % N for i, c in enumerate(coeffs)][1:]

    def evalp(cs, x):
        r = 0
        for c in reversed(cs):
            r = (r * x + c) % N
        return r

    x = root_fp % p
    pk = p
    while pk < N:
        pk = min(pk * pk, N)
        f = evalp(coeffs, x) % pk
        df = evalp(deriv, x) % pk
        try:
            dfi = modinv(df, pk)
        except ZeroDivisionError:
            break
        x = (x - f * dfi) % pk
    return x % N


def find_pi(p, ap):
    """Find π = a + bω in Z[ω] with N(π) = p and Tr(π) = ap.
    N(π) = a² − ab + b², Tr(π) = 2a − b = ap (since ω + ω̄ = -1)."""
    # 2a - b = ap  →  b = 2a - ap
    # a² - ab + b² = p
    # Substitute b: a² - a(2a-ap) + (2a-ap)² = p
    #             = a² - 2a² + a·ap + 4a² - 4a·ap + ap²
    #             = 3a² - 3a·ap + ap² = p
    # 3a² - 3a·ap + (ap² - p) = 0
    # a = (3ap ± √(9ap² - 12(ap² - p))) / 6
    #   = (3ap ± √(12p - 3ap²)) / 6
    disc = 12 * p - 3 * ap * ap
    if disc < 0:
        return None
    s = int(round(disc ** 0.5))
    if s * s != disc:
        # Not exact integer sqrt — try nearby
        for d in range(-2, 3):
            cand = s + d
            if cand >= 0 and cand * cand == disc:
                s = cand
                break
        else:
            return None
    if (3 * ap + s) % 6 == 0:
        a = (3 * ap + s) // 6
    elif (3 * ap - s) % 6 == 0:
        a = (3 * ap - s) // 6
    else:
        return None
    b = 2 * a - ap
    # Verify
    if a * a - a * b + b * b != p:
        return None
    if 2 * a - b != ap:
        return None
    return (a, b)


def cm_action_omega(C, P, omega_lift):
    """Apply ω·(x, y) = (ω·x, y) on the j=0 curve y² = x³ + b.
    Operates in projective coordinates: (X:Y:Z) → (ω·X:Y:Z).
    Note this is an automorphism of the curve."""
    X, Y, Z = P
    return ((omega_lift * X) % C.N, Y % C.N, Z % C.N)


def pi_action(C, P, a, b, omega_lift):
    """Apply π = a + b·ω on a projective point P."""
    aP = C.mul(a, P)
    omega_P = cm_action_omega(C, P, omega_lift)
    bw_P = C.mul(b, omega_P)
    return C.add(aP, bw_P)


def formal_neg_kernel(C, R, p):
    """Project a kernel-of-reduction projective point R to its
    formal-group parameter z = -X/Y. Returns z mod N (an element of
    pZ/p^eZ) or None if Y is non-unit."""
    X, Y, Z = R
    try:
        return (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return None


def canonical_lift_point(E_p, C, P_fp, a_pi, b_pi, omega_lift,
                         max_iter=30):
    """Compute the canonical (Frobenius-fixed) lift of P_fp ∈ E(F_p)
    in E(Z/p^e) via Newton iteration on f(Q) = π·Q − Q.

    Tangent-space action: on the formal group Ê (whose tangent space
    at e is one-dimensional Z_p), the endomorphism π acts as the
    scalar π_t = a_pi + b_pi · ω_lift ∈ Z_p (using the embedding
    Z[ω] → Z_p sending ω → ω_lift). So (π − 1)_t = (a_pi − 1) +
    b_pi · ω_lift, which is a UNIT in Z_p iff gcd(#E, p) = 1
    (i.e., ordinary), exactly our setting.

    Newton step:  Q ← Q − f(Q) / f'(Q)
                    = Q − (1/(π−1)_t) · (π·Q − Q)
    where 1/(π−1)_t is computed in (Z/p^e)^* and the multiplication
    by it on the formal group is mul-by-scalar in projective E
    coordinates (which equals formal-group scaling on Ê).
    """
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    Q = C.from_affine(teichmuller_lift(E_p, E_aff, P_fp))
    pi_tangent = (a_pi + b_pi * omega_lift) % C.N
    pim1_tangent = (pi_tangent - 1) % C.N
    try:
        inv_pim1 = modinv(pim1_tangent, C.N)
    except ZeroDivisionError:
        return Q  # singular — give up
    for _ in range(max_iter):
        piQ = pi_action(C, Q, a_pi, b_pi, omega_lift)
        R = C.add(piQ, C.neg(Q))
        if C.is_identity(R):
            return Q
        # Correction: -(1/(π−1)_t) · R   in the formal group
        neg_R = C.neg(R)
        correction = C.mul(inv_pim1, neg_R)
        Q_new = C.add(Q, correction)
        if C.equal(Q_new, Q):
            break
        Q = Q_new
    return Q


def z_delta(E_p, C, G_fp, k, lift_fn):
    """δ(k) = [k]·lift_fn(G) − lift_fn(kG); return its z = -X/Y."""
    kG_fp = E_p.mul(k, G_fp)
    tau_G = lift_fn(G_fp)
    tau_kG = lift_fn(kG_fp)
    diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
    return formal_neg_kernel(C, diff, None)


def linear_fit_test(seq, N):
    if len(seq) < 2 or seq[0] is None or seq[1] is None:
        return None, None, 0
    A = (seq[1] - seq[0]) % N
    B = (seq[0] - A) % N
    hits = sum(1 for i, v in enumerate(seq)
               if v is not None and (A * i + B - v) % N == 0)
    return A, B, hits


def run(p, b, e=4):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G_fp, g_ord = find_generator(E_p)
    if g_ord < 12:
        return {"prime": p, "status": "small"}
    if p % 3 != 1:
        return {"prime": p, "status": f"p % 3 = {p % 3} (need 1 for split)"}
    ap = p + 1 - n
    pi = find_pi(p, ap)
    if pi is None:
        return {"prime": p, "status": f"could not find π: ap={ap}"}
    a_pi, b_pi = pi

    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)

    # Lift ω: root of x² + x + 1 in Z/p^e
    # over F_p, ω = (-1 + √-3)/2 = ((p-1)/2)·(1 + √-3)
    # find a primitive cube root mod p
    omega_fp = None
    for w in range(2, p):
        if pow(w, 3, p) == 1 and w != 1:
            omega_fp = w
            break
    if omega_fp is None:
        return {"prime": p, "status": "no cube root mod p"}
    omega_lift = hensel_lift_root([1, 1, 1], omega_fp, p, e)
    # Verify ω² + ω + 1 ≡ 0 mod p^e
    check = (omega_lift * omega_lift + omega_lift + 1) % N
    if check != 0:
        return {"prime": p, "status": f"ω lift bad: check={check}"}

    # Two lift functions: naive and canonical
    E_aff = Curve(a=C.a, b=C.b, N=C.N)

    def naive_lift(P):
        return C.from_affine(teichmuller_lift(E_p, E_aff, P))

    def canonical_lift(P):
        return canonical_lift_point(E_p, C, P, a_pi, b_pi, omega_lift)

    # Compute z(δ_naive(k)) and z(δ_can(k)) for k = 1..max_k
    max_k = min(g_ord - 1, 30)
    delta_naive = [None]
    delta_can = [None]
    for k in range(1, max_k + 1):
        delta_naive.append(z_delta(E_p, C, G_fp, k, naive_lift))
        delta_can.append(z_delta(E_p, C, G_fp, k, canonical_lift))

    # Sanity: do the two sections AGREE on the generator G?
    nG = naive_lift(G_fp)
    cG = canonical_lift(G_fp)
    same_G = C.equal(nG, cG)

    # If they differ, summarize the difference at z-level
    diff_G_X, diff_G_Y, diff_G_Z = C.add(cG, C.neg(nG))
    try:
        diff_G_z = (-diff_G_X * modinv(diff_G_Y % N, N)) % N
    except ZeroDivisionError:
        diff_G_z = None

    # Linear-fit tests on both sequences
    seq_n = [v for v in delta_naive[1:] if v is not None]
    seq_c = [v for v in delta_can[1:] if v is not None]
    A_n, B_n, h_n = linear_fit_test(seq_n, N)
    A_c, B_c, h_c = linear_fit_test(seq_c, N)

    # Are δ_naive and δ_can ever equal? Compute their pairwise differences
    diffs = []
    for k in range(1, max_k + 1):
        if delta_naive[k] is None or delta_can[k] is None:
            continue
        d = (delta_can[k] - delta_naive[k]) % N
        diffs.append(d)
    nz_diffs = sum(1 for d in diffs if d != 0)

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "ap": ap,
        "pi": (a_pi, b_pi),
        "norm_pi_check": a_pi * a_pi - a_pi * b_pi + b_pi * b_pi,
        "trace_pi_check": 2 * a_pi - b_pi,
        "omega_fp": omega_fp,
        "omega_lift_first_digits": omega_lift % (p ** 2),
        "naive_canonical_agree_on_G": same_G,
        "diff_G_z_first_digits": (diff_G_z % (p ** 2)) if diff_G_z is not None else None,
        "naive_linear_fit": f"{h_n}/{len(seq_n)}",
        "canonical_linear_fit": f"{h_c}/{len(seq_c)}",
        "samples": max_k,
        "nonzero_diffs_naive_vs_canonical": f"{nz_diffs}/{len(diffs)}",
        "first_5_naive": seq_n[:5],
        "first_5_canonical": seq_c[:5],
    }


def main():
    candidates = [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5), (103, 5)]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in candidates:
        print(f"[phase32] p={p} b={b}")
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
    (out_dir / "phase32_canonical_lift.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
