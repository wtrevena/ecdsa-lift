"""
Phase 46 — "Close the lid" check: does the lift-error sequence z(delta(k))
match any known closed form / OEIS sequence?

Motivation
----------
The inertness theorem (Phase 41, notes/theorem.md) says the section's
lift error delta(k) = [k] tau(G) - tau(kG) lives in the kernel of
reduction and, read through the formal-group parameter z = -X/Y, behaves
like a pseudorandom-looking element of pZ/p^e (subject to the exact
relations V1-V5 of Phase 41: antisymmetry around Cn, etc.). Before the
paper claims "pseudorandom / no exploitable structure", we run one cheap,
falsifiable check: is z(delta(k)) actually a *named* integer sequence in
disguise (OEIS), or expressible as a simple closed form mod p?

Prior on a hit: low (~2%). Expected outcome: NULL. A null is a one-line
footnote, not a result. This script produces:

  1. The integer sequences, for a few small primes:
       - z(delta(k))               for k = 1..n-1   (raw, in Z/p^e)
       - z(delta(k)) mod p         (identically 0: delta is a kernel
                                    element, z in pZ -- kept as a sanity
                                    check, confirmed all-zero)
       - d1(k) = (z(delta(k))//p) mod p   the FIRST NON-TRIVIAL p-adic
                                    digit of z -- the genuinely varying
                                    F_p-valued sequence; this is what we
                                    actually probe / fit
       - leading decimal digit of z(delta(k))  (a deliberately silly
                                    "looks like a sequence" probe -- the
                                    kind of short sequence that matches
                                    OEIS spuriously, included to calibrate
                                    skepticism)
     printed to stdout and saved to results/phase46_oeis.json, formatted
     ready to paste into the OEIS search box.

  2. Offline symbolic-regression / closed-form detection on d1(k) as a
     function of k, with simple modular models:
       - constant            f(k) = c
       - affine / linear     f(k) = a*k + c
       - quadratic           f(k) = a*k^2 + b*k + c
       - cubic               f(k) = a*k^3 + ... + c
       - pure multiplicative f(k) = a*k
       - geometric           f(k) = a * c^k
     Each is fit exactly over Z/p (least-squares has no meaning mod p;
     we instead solve the linear system on a training half and measure
     exact-match accuracy on a held-out half), and compared against the
     baseline accuracy of guessing the most common residue (the
     "random sequence" null). A genuine closed form would generalize
     (~100% held-out); a pseudorandom sequence will sit at baseline.

The web (OEIS) part is NOT done here -- per the task, web queries are run
with the WebSearch/web_fetch tools and the exact URLs + verdicts are
recorded in notes/phase46_oeis.md. This script emits the search strings.

Reuses: src.curves, src.projective, src.lifts.teichmuller_lift, and the
exact z_delta construction from experiments/phase41_theorem_verification.py
(single curve op then read z; degeneracy-flagged) so the sequences are the
same object the theorem is about.
"""
from __future__ import annotations
import sys
import json
import pathlib
from math import gcd
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


# --------------------------------------------------------------------------
# Sequence generation  (z(delta(k)), identical method to Phase 41)
# --------------------------------------------------------------------------

def build_z_delta(p: int, b: int, e: int = 4):
    """Return (n, z_delta_fn, p, e, N) for E: y^2 = x^3 + b over F_p,
    or None if gcd(n, p) != 1 (theorem hypothesis fails -> skip).

    z_delta_fn(k) -> (z, ok): z = z(delta(k)) in Z/p^e and a boolean ok
    flagging whether the pure-Python projective arithmetic stayed
    non-degenerate (ok == False => value not trustworthy; same wall as
    Phase 41 / Phase 32)."""
    N = p ** e
    E_p = Curve(a=0, b=b, N=p)
    G_fp, n = find_generator(E_p)
    if gcd(n, p) != 1:
        return None
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)

    def tau(P_fp):
        return C.from_affine(teichmuller_lift(E_p, E_aff, P_fp))

    tau_G = tau(G_fp)

    def z_delta(k: int):
        kG = E_p.mul(k % n, G_fp) if k % n != 0 else None
        tkG = C.identity() if kG is None else tau(kG)
        d = C.add(C.mul(k, tau_G), C.neg(tkG))
        X, Y, Z = d
        if gcd(Y % N, N) != 1:
            return (None, False)
        z = (-X * modinv(Y % N, N)) % N
        return (z, z % p == 0)  # kernel elements satisfy z = 0 mod p

    return n, z_delta, p, e, N


def sequences_for_prime(p: int, b: int, e: int = 4):
    """Compute the sequences for one prime.

    Returns a dict with the raw / mod-p / first-p-adic-digit / leading-digit
    sequences plus bookkeeping (length, how many entries were degenerate, an
    'arithmetic_reliable' overall flag = no degenerate entries)."""
    built = build_z_delta(p, b, e)
    if built is None:
        return {"p": p, "b": b, "e": e, "status": "skip: gcd(n,p) != 1"}
    n, z_delta, p, e, N = built

    raw, modp, digit1, leading = [], [], [], []
    degenerate = 0
    for k in range(1, n):
        z, ok = z_delta(k)
        if not ok or z is None:
            degenerate += 1
            raw.append(None)
            modp.append(None)
            digit1.append(None)
            leading.append(None)
            continue
        raw.append(z)
        modp.append(z % p)               # ALWAYS 0: delta is a kernel elt (z in pZ)
        digit1.append((z // p) % p)      # first NON-trivial p-adic digit of z
        leading.append(int(str(z)[0]) if z != 0 else 0)

    raw_clean = [x for x in raw if x is not None]
    modp_clean = [x for x in modp if x is not None]
    digit1_clean = [x for x in digit1 if x is not None]
    leading_clean = [x for x in leading if x is not None]

    return {
        "p": p, "b": b, "e": e, "n": n, "N": N,
        "len_full": n - 1,
        "degenerate_entries": degenerate,
        "arithmetic_reliable": degenerate == 0,
        "seq_mod_p_all_zero": all(x == 0 for x in modp_clean),
        "seq_raw": raw,
        "seq_mod_p": modp,
        "seq_first_padic_digit": digit1,
        "seq_leading_digit": leading,
        "seq_raw_clean": raw_clean,
        "seq_mod_p_clean": modp_clean,
        "seq_first_padic_digit_clean": digit1_clean,
        "seq_leading_digit_clean": leading_clean,
        "oeis_query_raw": ",".join(map(str, raw_clean)),
        "oeis_query_mod_p": ",".join(map(str, modp_clean)),
        "oeis_query_first_padic_digit": ",".join(map(str, digit1_clean)),
        "oeis_query_leading_digit": ",".join(map(str, leading_clean)),
    }


# --------------------------------------------------------------------------
# Offline closed-form detection on the first non-trivial p-adic digit of
# z(delta(k)),  d1(k) := (z(delta(k)) // p) mod p.
# --------------------------------------------------------------------------
#
# NOTE: z(delta(k)) mod p is identically 0 (delta is a kernel-of-reduction
# element, z in pZ), so fitting THAT is vacuous -- every model scores 1.0 on
# the constant-zero sequence and the test proves nothing. The genuinely
# varying F_p-valued sequence is the *first non-trivial p-adic digit* d1(k);
# that is what we subject to closed-form detection below.
#
# We fit f(k) = d1(k) with simple *modular* models. There is no meaningful
# least-squares over Z/p, so instead we:
#   - take a training prefix (first half of valid k),
#   - solve the model's parameters exactly over Z/p if invertible,
#   - measure EXACT-match accuracy on the held-out second half.
#
# Baseline (the "random sequence" null): predict the single most common
# residue seen in training. For a uniform-random sequence mod p the
# expected held-out accuracy is ~1/p. A real closed form generalizes to
# ~1.0. We report each model's held-out accuracy vs this baseline.

def _solve_linear_modp(A_rows, y_vals, p):
    """Solve A x = y over Z/p by Gaussian elimination. Returns x or None
    if singular / inconsistent / underdetermined."""
    ncol = len(A_rows[0])
    M = [list(map(lambda v: v % p, row)) + [y % p] for row, y in zip(A_rows, y_vals)]
    nrow = len(M)
    piv_cols = []
    r = 0
    for c in range(ncol):
        piv = None
        for i in range(r, nrow):
            if M[i][c] % p != 0:
                piv = i
                break
        if piv is None:
            continue
        M[r], M[piv] = M[piv], M[r]
        inv = modinv(M[r][c], p)
        M[r] = [(v * inv) % p for v in M[r]]
        for i in range(nrow):
            if i != r and M[i][c] % p != 0:
                f = M[i][c]
                M[i] = [(a - f * b) % p for a, b in zip(M[i], M[r])]
        piv_cols.append(c)
        r += 1
        if r == ncol:
            break
    if len(piv_cols) < ncol:
        return None
    for i in range(r, nrow):
        if M[i][ncol] % p != 0 and all(M[i][c] % p == 0 for c in range(ncol)):
            return None
    x = [0] * ncol
    for idx, c in enumerate(piv_cols):
        x[c] = M[idx][ncol] % p
    return x


def fit_models_mod_p(seq, p, n):
    """seq is indexed 0..n-2 corresponding to k=1..n-1 (None entries are
    degenerate and dropped here together with their k). It should be the
    informative F_p-valued sequence d1(k) = (z//p) mod p."""
    ks = [k for k in range(1, n) if seq[k - 1] is not None]
    ys = [seq[k - 1] for k in ks]
    m = len(ks)
    if m < 6:
        return {"status": f"too few valid points ({m}) to fit/test"}

    split = m // 2
    tr_k, tr_y = ks[:split], ys[:split]
    te_k, te_y = ks[split:], ys[split:]

    results = {}

    def held_out_acc(predict):
        correct = sum(1 for k, y in zip(te_k, te_y) if predict(k) == y % p)
        return correct / len(te_k)

    common = Counter(tr_y).most_common(1)[0][0]
    results["baseline_mode"] = {
        "held_out_accuracy": held_out_acc(lambda k: common),
        "params": {"const": common},
        "n_train": len(tr_k), "n_test": len(te_k),
    }
    results["baseline_uniform_random_expected"] = {
        "held_out_accuracy": round(1.0 / p, 6),
        "note": "expected accuracy guessing for a uniform sequence mod p",
    }

    def make_pred(coef):
        return lambda k: sum(coef[j] * pow(k, j, p)
                             for j in range(len(coef))) % p

    for degree, name in [(0, "constant"), (1, "affine_linear"),
                         (2, "quadratic"), (3, "cubic")]:
        ncoef = degree + 1
        if len(tr_k) < ncoef:
            continue
        # First ask the strong question: is the WHOLE training half exactly a
        # degree<=d polynomial mod p? (overdetermined system, solved exactly).
        A_all = [[pow(k, j, p) for j in range(ncoef)] for k in tr_k]
        coef_exact = _solve_linear_modp(A_all, tr_y, p)
        if coef_exact is not None:
            pred = make_pred(coef_exact)
            tr_acc = sum(1 for k, y in zip(tr_k, tr_y) if pred(k) == y % p) / len(tr_k)
            results[name] = {
                "held_out_accuracy": held_out_acc(pred),
                "train_accuracy": tr_acc,
                "params": {"coeffs_low_to_high": coef_exact},
                "fit_mode": "exact_over_all_training",
                "n_train": len(tr_k), "n_test": len(te_k),
            }
            continue
        # Not globally polynomial on training. Fall back to INTERPOLATION on the
        # first ncoef training points (always solvable; Vandermonde invertible
        # for distinct k mod p), then test how that degree-d hypothesis
        # generalizes. A true closed form -> ~1.0 held-out; noise -> ~1/p.
        A_min = [[pow(tr_k[i], j, p) for j in range(ncoef)] for i in range(ncoef)]
        coef = _solve_linear_modp(A_min, tr_y[:ncoef], p)
        if coef is None:
            results[name] = {"held_out_accuracy": None, "params": None,
                             "note": "singular interpolation (repeated k mod p)"}
            continue
        pred = make_pred(coef)
        tr_acc = sum(1 for k, y in zip(tr_k, tr_y) if pred(k) == y % p) / len(tr_k)
        results[name] = {
            "held_out_accuracy": held_out_acc(pred),
            "train_accuracy": tr_acc,
            "params": {"coeffs_low_to_high": coef},
            "fit_mode": "interpolation_first_%d_points" % ncoef,
            "note": "no degree-%d polynomial fits ALL training points; this is "
                    "the interpolant through the first %d -> generalization is the "
                    "honest test" % (degree, ncoef),
            "n_train": len(tr_k), "n_test": len(te_k),
        }

    a_votes = Counter()
    for k, y in zip(tr_k, tr_y):
        if k % p != 0:
            a_votes[(y * modinv(k % p, p)) % p] += 1
    if a_votes:
        a_best = a_votes.most_common(1)[0][0]
        pred = lambda k: (a_best * k) % p
        tr_acc = sum(1 for k, y in zip(tr_k, tr_y) if pred(k) == y % p) / len(tr_k)
        results["multiplicative_a_k"] = {
            "held_out_accuracy": held_out_acc(pred),
            "train_accuracy": tr_acc,
            "params": {"a": a_best},
            "n_train": len(tr_k), "n_test": len(te_k),
        }

    best = None
    for c in range(1, p):
        a_votes = Counter()
        for k, y in zip(tr_k, tr_y):
            ck = pow(c, k, p)
            if ck != 0:
                a_votes[(y * modinv(ck, p)) % p] += 1
        if not a_votes:
            continue
        a_best = a_votes.most_common(1)[0][0]
        predc = lambda k, a=a_best, cc=c: (a * pow(cc, k, p)) % p
        tr_acc = sum(1 for k, y in zip(tr_k, tr_y) if predc(k) == y % p) / len(tr_k)
        if best is None or tr_acc > best[0]:
            best = (tr_acc, c, a_best)
    if best is not None:
        tr_acc, c, a = best
        pred = lambda k: (a * pow(c, k, p)) % p
        results["geometric_a_cpow_k"] = {
            "held_out_accuracy": held_out_acc(pred),
            "train_accuracy": tr_acc,
            "params": {"a": a, "c": c},
            "n_train": len(tr_k), "n_test": len(te_k),
            "note": "c chosen by best TRAINING fit (search over F_p^*), "
                    "so train_accuracy is optimistic; held-out is the honest number",
        }

    return results


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

# secp256k1 j=0 family (b=7), as in Phase 41. The task asks for
# p in {31, 43, 67, 97}; we include those plus a couple more that Phase 41
# showed are arithmetically reliable, for a longer (less spurious) sequence.
PRIMES = [(31, 7), (43, 7), (67, 7), (97, 7), (109, 7), (151, 7), (163, 7)]


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    report = {"phase": 46, "curve_family": "y^2 = x^3 + 7 (j=0)",
              "e": 4, "primes": [], "fits": {}}

    for p, b in PRIMES:
        print(f"\n[phase46] p={p} b={b}")
        s = sequences_for_prime(p, b, e=4)
        if "status" in s:
            print(f"   {s['status']}")
            report["primes"].append(s)
            continue
        flag = "" if s["arithmetic_reliable"] else \
            f"  *** {s['degenerate_entries']} DEGENERATE entries (untrustworthy) ***"
        print(f"   n={s['n']}  seq length={s['len_full']}{flag}")
        print(f"   z(delta(k))             : {s['seq_raw_clean']}")
        print(f"   z(delta(k)) mod p       : {s['seq_mod_p_clean']}"
              f"   [all-zero: {s['seq_mod_p_all_zero']}  (kernel elt, expected)]")
        print(f"   first p-adic digit d1(k): {s['seq_first_padic_digit_clean']}")
        print(f"   leading decimal digit   : {s['seq_leading_digit_clean']}")
        print(f"   OEIS query (d1, mod p)  : {s['oeis_query_first_padic_digit']}")

        fits = fit_models_mod_p(s["seq_first_padic_digit"], p, s["n"])
        report["fits"][str(p)] = fits
        if "status" in fits:
            print(f"   fit: {fits['status']}")
        else:
            base = fits["baseline_mode"]["held_out_accuracy"]
            rnd = fits["baseline_uniform_random_expected"]["held_out_accuracy"]
            print(f"   --- modular closed-form fit of d1(k) (held-out accuracy) ---")
            print(f"       baseline (mode)        : {base:.3f}")
            print(f"       baseline (uniform 1/p) : {rnd:.3f}")
            for name in ("constant", "affine_linear", "quadratic", "cubic",
                         "multiplicative_a_k", "geometric_a_cpow_k"):
                if name in fits and fits[name].get("held_out_accuracy") is not None:
                    print(f"       {name:22s}: {fits[name]['held_out_accuracy']:.3f}")

        report["primes"].append(s)

    out_path = out_dir / "phase46_oeis.json"
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[phase46] wrote {out_path}")

    print("\n[phase46] OEIS search URLs (for the memo / manual web queries):")
    for s in report["primes"]:
        if "status" in s:
            continue
        for label, q in (("raw", s["oeis_query_raw"]),
                         ("d1_modp", s["oeis_query_first_padic_digit"]),
                         ("leading", s["oeis_query_leading_digit"])):
            url = f"https://oeis.org/search?q={q}&fmt=text"
            print(f"   p={s['p']:>4} [{label:>8}] {url}")


if __name__ == "__main__":
    main()
