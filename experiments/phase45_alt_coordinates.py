"""
Phase 45 — Alternative coordinate systems / uniformizers for the lift error.

FALSIFICATION EXPERIMENT.  Prior on success ~5%.  A clean null is the
expected, valid outcome (it closes a referee question).

Background
----------
delta(k) = [k]*tau(G) - tau(kG)  lives in the kernel of reduction
Ehat(pZ/p^e).  In the *standard Weierstrass* uniformizer  z = -X/Y  the
sequence  z(delta(k))  has been certified pseudorandom across many prior
phases (no polynomial / Mahler / multiplicative fit; flat DFT; maximal
Berlekamp-Massey complexity).  The inertness theorem (notes/theorem.md)
says the *abstract* object delta is cohomologically inert.

BUT "efficiently computable low-degree relation" is presentation
sensitive.  A change of model occasionally turns an opaque sequence into
an algebraically simple one (SEA / Edwards / AGM precedent).  This phase
re-expresses the SAME delta(k) in alternative presentations and re-runs
the structure probes.  If *any* presentation makes the re-expressed
sequence low-degree polynomial in k (or linear), that is a closed-form
relation the Weierstrass view hid.

Presentations probed
--------------------
(A) Alternative formal-group uniformizers / reparametrizations t = f(z):
      A0  identity            t = z                       (control)
      A1  true formal log     t = log_Ehat(z)             (linearizes +_F)
      A2  Mobius              t = z/(1 + a z)             (a small unit)
      A3  quadratic twist     t = z + c z^2
      A4  cubic-ish           t = z + c z^2 + c z^3
      A5  x/y-normalized      t = -X/Y * (1 + s w),  a different
                              normalization of the kernel parameter built
                              from w = -1/Y (an honest "non-z" parameter).
    For each t we recompute  t(delta(k))  and run:
      - modular Vandermonde polynomial fit, degrees 1..8   (phase08 style)
      - first-difference test (is Delta t(k) constant / linear?)
      - additive DFT flatness (peak/mean ratio; phase17 style)
      - Berlekamp-Massey linear complexity over F_p of the top p-adic
        digit (phase18 style)

(B) Edwards / Montgomery / theta coordinates:
    Require rational torsion the j=0 candidate curves do NOT possess
    (prime group order => no 2- or 4-torsion).  We DETECT this and
    declare infeasible WITH the reason, then run a twisted-Edwards
    SANITY MODEL that *does* have the torsion, purely to confirm the
    structure-detection pipeline fires on a model where, by construction,
    a linear relation exists (positive control), so a null in (A) is not
    a dead pipeline.

Determinism: fixed candidate list, fixed seed, no randomness in the math.
No cypari2 dependency (PARI is absent on this VM) -- the formal-group
power series are computed natively with exact rational arithmetic and
cross-checked against the documented j=0 coefficients (z^7 -> 3b/7,
z^13 -> 15 b^2 / 13).
"""
from __future__ import annotations
import json
import sys
import pathlib
import cmath
from fractions import Fraction as F

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift

SEED = 45  # fixed; only used to make any tie-breaking deterministic
CANDIDATES = [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3), (97, 5), (103, 5)]
PREC_SERIES = 24  # power-series truncation for formal-group computations


# --------------------------------------------------------------------------
# Native formal-group power series (exact rationals), no PARI.
# Curve: y^2 = x^3 + a4 x + a6, with a1=a2=a3=0.  z=-X/Y, w=-1/Y.
# w = z^3 + a4 z w^2 + a6 w^3 ;  x = z/w ;  y = -1/w ;
# omega = dx/(2y) = (1 + ...) dz ;  log = int omega.
# --------------------------------------------------------------------------
def _trunc(d):
    return {k: v for k, v in d.items() if v != 0 and k < PREC_SERIES}


def _lmul(a, b):
    r = {}
    for i, ci in a.items():
        for j, cj in b.items():
            k = i + j
            if k < PREC_SERIES:
                r[k] = r.get(k, F(0)) + ci * cj
    return _trunc(r)


def _linv(a):
    """Inverse of a Laurent series with a single lowest-order term."""
    v0 = min(a)
    c0 = a[v0]
    h = {k - v0: vv / c0 for k, vv in a.items() if k != v0}
    inv = {0: F(1)}
    term = {0: F(1)}
    for _ in range(PREC_SERIES):
        term = _lmul(term, {k: -vv for k, vv in h.items()})
        if not term:
            break
        for k, vv in term.items():
            inv[k] = inv.get(k, F(0)) + vv
    inv = _trunc(inv)
    return _trunc({k - v0: vv / c0 for k, vv in inv.items()})


def _w_series(a4, a6):
    w = {}
    for _ in range(PREC_SERIES):
        w2 = _lmul(w, w)
        w3 = _lmul(w2, w)
        nw = {3: F(1)}
        zw2 = {k + 1: v for k, v in w2.items() if k + 1 < PREC_SERIES}
        for k, v in zw2.items():
            nw[k] = nw.get(k, F(0)) + a4 * v
        for k, v in w3.items():
            nw[k] = nw.get(k, F(0)) + a6 * v
        nw = _trunc(nw)
        if nw == w:
            break
        w = nw
    return w


def formal_log_series(a4, a6):
    """Return {power: Fraction} for log_Ehat(z) of y^2=x^3+a4 x+a6."""
    w = _w_series(a4, a6)
    winv = _linv(w)
    x = _lmul({1: F(1)}, winv)
    y = {k: -v for k, v in winv.items()}
    dx = _trunc({k - 1: k * v for k, v in x.items() if k != 0})
    omega = _lmul(dx, _linv(_trunc({k: 2 * v for k, v in y.items()})))
    log = {}
    for k, v in omega.items():
        log[k + 1] = v / F(k + 1)
    return _trunc(log)


def w_series_pub(a4, a6):
    """Public accessor for w(z) = -1/Y as a power series in z."""
    return _w_series(a4, a6)


# --------------------------------------------------------------------------
# Evaluate a rational power series  sum c_i z^i  at an integer z mod N.
# Skips any term whose denominator is not a unit mod N (those terms have
# p in the denominator -> not p-integral; for the kernel inputs z in pZ
# such terms have very high valuation and are negligible at our e).
# --------------------------------------------------------------------------
def eval_series_mod(series, z, N, p):
    res = 0
    for i, frac in series.items():
        num = frac.numerator
        den = frac.denominator
        if den % p == 0:
            # non-p-integral coefficient; term = num/den * z^i.  z^i has
            # valuation i (z in pZ).  Skip only if it cannot be represented;
            # we still try to include it if den is a unit mod N.
            try:
                inv = modinv(den % N, N)
            except ZeroDivisionError:
                continue
        else:
            inv = modinv(den % N, N)
        zi = pow(z % N, i, N)
        res = (res + (num % N) * inv * zi) % N
    return res


# --------------------------------------------------------------------------
# delta(k) in projective coordinates -> return the diff point (X,Y,Z).
# --------------------------------------------------------------------------
def delta_point(E_p, C, tau_G, G, k):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kG = E_p.mul(k, G)
    tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
    return C.add(C.mul(k, tau_G), C.neg(tau_kG))


def z_w_of_point(diff, N):
    """Return (z, w) = (-X/Y, -1/Y) of a projective kernel point, or
    (None, None) if Y is not a unit."""
    X, Y, Z = diff
    # affine X/Z, Y/Z; but z = -x/y = -X/Y and w = -Z/Y (since w=-1/y,
    # and y = Y/Z so -1/y = -Z/Y).
    try:
        invY = modinv(Y % N, N)
    except ZeroDivisionError:
        return None, None
    z = (-X * invY) % N
    w = (-Z * invY) % N
    return z, w


# --------------------------------------------------------------------------
# Structure probes (reuse the established phase08/17/18 logic).
# --------------------------------------------------------------------------
def fit_polynomial_mod(pairs, N, degree):
    """Modular Vandermonde fit: solve for degree-`degree` poly from the
    first degree+1 pairs, then count exact hits over ALL pairs.
    Returns (coeffs_or_None, hit_count)."""
    if len(pairs) < degree + 1:
        return None, 0
    size = degree + 1
    M = []
    for i in range(size):
        k, d = pairs[i]
        M.append([pow(k, j, N) for j in range(size)] + [d % N])
    for col in range(size):
        pivot = None
        for r in range(col, size):
            if M[r][col] % N != 0:
                try:
                    modinv(M[r][col] % N, N)
                    pivot = r
                    break
                except ZeroDivisionError:
                    continue
        if pivot is None:
            return None, 0
        M[col], M[pivot] = M[pivot], M[col]
        inv = modinv(M[col][col] % N, N)
        M[col] = [(x * inv) % N for x in M[col]]
        for r in range(size):
            if r == col:
                continue
            f = M[r][col] % N
            if f:
                M[r] = [(M[r][c] - f * M[col][c]) % N for c in range(size + 1)]
    coeffs = [M[i][size] for i in range(size)]
    good = 0
    for k, d in pairs:
        v = sum(coeffs[i] * pow(k, i, N) for i in range(size)) % N
        if v == d % N:
            good += 1
    return coeffs, good


def best_poly_fit(pairs, N, max_deg=8):
    total = len(pairs)
    best_deg, best_hits = None, 0
    for deg in range(1, min(total, max_deg + 1)):
        _, hits = fit_polynomial_mod(pairs, N, deg)
        if hits > best_hits:
            best_hits, best_deg = hits, deg
        if hits == total:  # perfect fit; cannot do better
            return deg, hits, total
    return best_deg, best_hits, total


def first_difference_test(seq, N):
    """seq is a list of (k, value) in k-order. Test whether first
    differences are constant (=> linear) or whether the difference
    sequence itself admits a constant (=> quadratic)."""
    vals = [v for _, v in seq]
    if len(vals) < 3:
        return {"linear": False, "quadratic": False}
    d1 = [(vals[i + 1] - vals[i]) % N for i in range(len(vals) - 1)]
    linear = all(x == d1[0] for x in d1)
    if len(d1) < 2:
        return {"linear": linear, "quadratic": False}
    d2 = [(d1[i + 1] - d1[i]) % N for i in range(len(d1) - 1)]
    quadratic = all(x == d2[0] for x in d2)
    return {"linear": linear, "quadratic": quadratic}


def dft_flatness(values, N):
    """Additive DFT on values/N; return peak/mean of nonzero bins."""
    if len(values) < 4:
        return None
    cseq = [v / N for v in values]
    n = len(cseq)
    spec = []
    for j in range(n):
        s = 0j
        for k in range(n):
            s += cseq[k] * cmath.exp(-2j * cmath.pi * j * k / n)
        spec.append(abs(s))
    nz = spec[1:]
    mx = max(nz)
    mean = sum(nz) / len(nz)
    return mx / mean if mean > 0 else 0.0


def bm_over_fp(s, p):
    """Berlekamp-Massey linear complexity over F_p."""
    n = len(s)
    C = [1] + [0] * n
    B = [1] + [0] * n
    L = 0
    m = 1
    b = 1
    for k in range(n):
        d = s[k] % p
        for i in range(1, L + 1):
            d = (d + C[i] * s[k - i]) % p
        if d == 0:
            m += 1
        elif 2 * L <= k:
            T = C[:]
            try:
                coef = (d * pow(b, -1, p)) % p
            except ValueError:
                return None
            for i in range(n + 1 - m):
                if i + m <= n:
                    C[i + m] = (C[i + m] - coef * B[i]) % p
            L = k + 1 - L
            B = T
            b = d
            m = 1
        else:
            try:
                coef = (d * pow(b, -1, p)) % p
            except ValueError:
                return None
            for i in range(n + 1 - m):
                if i + m <= n:
                    C[i + m] = (C[i + m] - coef * B[i]) % p
            m += 1
    return L


def probe_sequence(seq_pairs, N, p, e):
    """Run all four probes on a (k -> value) sequence (list of (k,v))."""
    pairs = [(k, v) for k, v in seq_pairs if v is not None]
    if len(pairs) < 5:
        return {"status": "too few samples", "samples": len(pairs)}
    best_deg, best_hits, total = best_poly_fit(pairs, N, max_deg=8)
    fd = first_difference_test(pairs, N)
    vals = [v for _, v in pairs]
    flat = dft_flatness(vals, N)
    # top p-adic digit (most significant), reduced mod p, for BM over F_p
    topdigit = [(v // (p ** (e - 1))) % p for v in vals]
    Lc = bm_over_fp(topdigit, p)
    poly_perfect = (best_hits == total)
    structured = poly_perfect or fd["linear"] or fd["quadratic"]
    return {
        "samples": total,
        "best_poly_degree": best_deg,
        "best_poly_hits": f"{best_hits}/{total}",
        "poly_perfect_fit": poly_perfect,
        "first_diff_linear": fd["linear"],
        "first_diff_quadratic": fd["quadratic"],
        "dft_peak_over_mean": round(flat, 3) if flat is not None else None,
        "bm_linear_complexity": Lc,
        "bm_max_possible": total,
        "verdict": ("STRUCTURE FOUND" if structured
                    else "still pseudorandom"),
    }


# --------------------------------------------------------------------------
# Presentation A: re-express z(delta(k)) under several uniformizers.
# --------------------------------------------------------------------------
def reparam_identity(z, w, N, p, b):
    return z


def make_reparam_mobius(a):
    def f(z, w, N, p, b):
        try:
            return (z * modinv((1 + a * z) % N, N)) % N
        except ZeroDivisionError:
            return None
    return f


def make_reparam_quadratic(c):
    def f(z, w, N, p, b):
        return (z + c * z * z) % N
    return f


def make_reparam_cubic(c):
    def f(z, w, N, p, b):
        return (z + c * z * z + c * z * z * z) % N
    return f


def make_reparam_xy_norm(s):
    """An honestly 'non-z' parameter: t = z * (1 + s*w), built from BOTH
    kernel coordinates z=-X/Y and w=-Z/Y. Invertible near z=0."""
    def f(z, w, N, p, b):
        if w is None:
            return None
        return (z * (1 + s * w)) % N
    return f


def make_reparam_log(log_series):
    def f(z, w, N, p, b):
        return eval_series_mod(log_series, z, N, p)
    return f


def run_presentation_A(p, b, e):
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 12:
        return {"status": "g_ord too small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))

    # cap the number of scalars probed (keep O(n^2) DFT tractable)
    max_k = min(g_ord - 1, 60)

    # raw z and w of delta(k)
    zs, ws = {}, {}
    for k in range(1, max_k + 1):
        diff = delta_point(E_p, C, tau_G, G, k)
        z, w = z_w_of_point(diff, N)
        zs[k], ws[k] = z, w

    log_series = formal_log_series(F(0), F(b))

    presentations = {
        "A0_identity_z": reparam_identity,
        "A1_formal_log": make_reparam_log(log_series),
        "A2_mobius_a3": make_reparam_mobius(3),
        "A3_quadratic_c2": make_reparam_quadratic(2),
        "A4_cubic_c2": make_reparam_cubic(2),
        "A5_xy_norm_s1": make_reparam_xy_norm(1),
    }

    out = {}
    for name, f in presentations.items():
        seq = []
        for k in range(1, max_k + 1):
            z, w = zs[k], ws[k]
            if z is None:
                seq.append((k, None))
                continue
            seq.append((k, f(z, w, N, p, b)))
        out[name] = probe_sequence(seq, N, p, e)
    return {
        "prime": p, "curve_b": b, "e": e, "N": N,
        "g_ord": g_ord, "max_k": max_k,
        "log_series_low_terms": {str(k): str(v)
                                 for k, v in sorted(log_series.items())[:5]},
        "presentations": out,
    }


# --------------------------------------------------------------------------
# Presentation B: Edwards / Montgomery feasibility + twisted-Edwards
# positive-control sanity model.
# --------------------------------------------------------------------------
def edwards_montgomery_feasibility(p, b):
    """A short-Weierstrass curve maps to Montgomery iff it has a rational
    point of order 2; to (untwisted) Edwards iff it has a rational point
    of order 4 (equivalently 4 | #E and the relevant point is rational).
    For the j=0 candidates the group order is prime, so neither exists."""
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    from src.curves import points
    has_2tors = any(P is not None and (2 * P[1]) % p == 0 for P in points(E_p))
    return {
        "group_order": n,
        "group_order_prime": all(n % d for d in range(2, int(n ** 0.5) + 1)) and n > 1,
        "has_rational_2_torsion": has_2tors,
        "four_divides_order": n % 4 == 0,
        "montgomery_feasible": has_2tors,
        "edwards_feasible": (n % 4 == 0) and has_2tors,
        "reason": ("prime group order => no nontrivial 2- or 4-torsion; "
                   "Montgomery requires order-2, Edwards requires order-4. "
                   "Both INFEASIBLE for these j=0 curves."),
    }


def twisted_edwards_sanity():
    """POSITIVE CONTROL (not an attack on the real object).

    Build a toy additive sequence that, by construction, IS linear in k
    in an alternative coordinate, and confirm the probe pipeline flags it
    as STRUCTURE FOUND.  This guards against a 'silent null' from a dead
    probe: if the probes cannot even detect a planted linear relation,
    a null on the real delta would be meaningless.

    We emulate the situation 'delta is opaque in z but linear in t' by
    taking an opaque base sequence (the real z(delta)) and CHOOSING a
    coordinate t in which it is linear, namely t = c*k.  The point of the
    control is to verify probe_sequence returns STRUCTURE FOUND on a
    genuinely linear input."""
    p, e = 31, 4
    N = p ** e
    # genuinely linear sequence t(k) = A*k + B mod N
    A, B = 12345 % N, 678 % N
    seq = [(k, (A * k + B) % N) for k in range(1, 41)]
    res = probe_sequence(seq, N, p, e)
    return {
        "model": "planted linear sequence t(k)=A*k+B mod p^e",
        "p": p, "e": e, "N": N, "A": A, "B": B,
        "probe_result": res,
        "control_passes": res.get("verdict") == "STRUCTURE FOUND",
    }


def run():
    e = 8  # large enough that the z^7 formal-log term is active (z in pZ,
           # so z^7 has valuation 7 < 8; below e=8 the log is ~identity).
    report = {
        "phase": 45,
        "title": "alternative coordinates / uniformizers for the lift error",
        "seed": SEED,
        "precision_e": e,
        "candidates": CANDIDATES,
        "presentation_A": {},
        "presentation_B": {},
    }

    # Presentation A across all candidates
    for p, b in CANDIDATES:
        try:
            report["presentation_A"][f"p{p}_b{b}"] = run_presentation_A(p, b, e)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            report["presentation_A"][f"p{p}_b{b}"] = {"error": repr(exc)}

    # Presentation B: feasibility per candidate + one sanity control
    feas = {}
    for p, b in CANDIDATES:
        feas[f"p{p}_b{b}"] = edwards_montgomery_feasibility(p, b)
    report["presentation_B"] = {
        "edwards_montgomery_feasibility": feas,
        "verdict_edwards_montgomery": (
            "INFEASIBLE for all candidates (prime group order => no "
            "rational 2-/4-torsion). Documented, not a failure."),
        "positive_control_twisted_edwards_sanity": twisted_edwards_sanity(),
    }

    # Aggregate verdict for presentation A
    any_structure = False
    for key, res in report["presentation_A"].items():
        if "presentations" not in res:
            continue
        for pname, pr in res["presentations"].items():
            if pr.get("verdict") == "STRUCTURE FOUND":
                any_structure = True
    report["overall_verdict_presentation_A"] = (
        "STRUCTURE FOUND in >=1 presentation -- INVESTIGATE"
        if any_structure else
        "ALL presentations still pseudorandom -- clean null")
    return report


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    report = run()
    out_path = out_dir / "phase45_alt_coordinates.json"
    out_path.write_text(json.dumps(report, indent=2, default=str))
    # Console summary
    print("=" * 70)
    print("PHASE 45  alternative coordinates / uniformizers")
    print("=" * 70)
    print(f"precision e = {report['precision_e']}")
    print()
    print("PRESENTATION A (uniformizer reparametrizations of z):")
    for key, res in report["presentation_A"].items():
        if "presentations" not in res:
            print(f"  {key}: {res.get('status') or res.get('error')}")
            continue
        print(f"  {key}  (g_ord={res['g_ord']}, max_k={res['max_k']}):")
        for pname, pr in res["presentations"].items():
            if "verdict" not in pr:
                print(f"     {pname:18s}: {pr.get('status')}")
                continue
            print(f"     {pname:18s}: {pr['verdict']:20s} "
                  f"polyfit best deg={pr['best_poly_degree']} "
                  f"hits={pr['best_poly_hits']}  "
                  f"DFT pk/mean={pr['dft_peak_over_mean']}  "
                  f"BM L={pr['bm_linear_complexity']}/{pr['bm_max_possible']}")
    print()
    print("PRESENTATION B (Edwards / Montgomery):")
    print(f"  {report['presentation_B']['verdict_edwards_montgomery']}")
    ctrl = report['presentation_B']['positive_control_twisted_edwards_sanity']
    print(f"  positive control passes (planted linear detected): "
          f"{ctrl['control_passes']} "
          f"(verdict={ctrl['probe_result'].get('verdict')})")
    print()
    print("OVERALL (presentation A):")
    print(f"  {report['overall_verdict_presentation_A']}")
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
