"""
Phase 36 — Iwasawa-style p-adic L-function probe.

Idea: take an elliptic curve E/Q whose reduction mod p matches one of
our test curves (same a_p, same j-invariant), compute its anticyclotomic
or cyclotomic p-adic L-function via Sage's `padic_lseries`, and read
off the Iwasawa invariants μ_p(E), λ_p(E).

Why this might matter:
  • The p-adic L-function is a power series L_p(T) ∈ Z_p[[T]] (or its
    field of fractions) whose coefficients encode arithmetic of E
    over the cyclotomic Z_p-extension of Q.
  • Its Newton polygon controls μ and λ.
  • If μ ≠ 0 OR λ has a structured pattern, this is a non-trivial
    p-adic invariant of E that "sees" something the formal-group /
    canonical-lift section does not.
  • Translating this back to the F_p group structure of E(F_p) is
    non-trivial, but at minimum it gives us a new computable
    invariant whose p-adic expansion can be compared against
    δ(k) and the cocycle entries.

If λ = 0 and μ = 0 (which is the "boring" generic case), we close
this angle; if not, we have something new to study.

For our test curves we need to find a global E/Q with the same
mod-p reduction. A simple way: search the Cremona / LMFDB tables
for curves with conductor coprime to p and matching a_p.

Run:  ~/miniforge3/bin/conda run -n sage sage experiments/phase36_iwasawa.py
"""
from sage.all import (EllipticCurve, GF, ZZ, Integer, cremona_curves,
                       prime_range)
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def find_global_curve_with_ap(p, target_ap, target_j=0, max_cond=200):
    """Find a curve E/Q with good reduction at p, same a_p, and same
    j-invariant mod p as our F_p curve."""
    Fp = GF(p)
    for cond in range(11, max_cond):
        try:
            curves = cremona_curves([cond])
        except Exception:
            continue
        for E in curves:
            if cond % p == 0:
                continue
            if E.discriminant() % p == 0:
                continue
            ap = E.ap(p)
            if int(ap) == target_ap:
                jE = E.j_invariant()
                # Check j mod p
                try:
                    jE_fp = Fp(jE)
                    if jE_fp == Fp(target_j):
                        return E
                except Exception:
                    pass
    return None


def run_iwasawa(p, b_curve):
    Fp = GF(p)
    E_fp = EllipticCurve(Fp, [0, b_curve])
    n = E_fp.cardinality()
    ap = p + 1 - int(n)
    j_target = int(E_fp.j_invariant())  # 0 for our curves

    # Try to find a global E/Q with matching reduction
    E_global = find_global_curve_with_ap(p, ap, target_j=j_target,
                                         max_cond=400)
    if E_global is None:
        # Fall back: any curve with good red. and matching a_p
        E_global = find_global_curve_with_ap(p, ap, target_j=None,
                                             max_cond=400)
    if E_global is None:
        return {"p": p, "status": f"no global E/Q with ap={ap}"}

    # Compute the p-adic L-series
    try:
        if ap % p == 0:
            # supersingular case — different machinery
            return {"p": p, "global_curve": str(E_global),
                    "status": "supersingular at p"}
        Lp = E_global.padic_lseries(p)
        # Series expansion to a few digits
        try:
            series = Lp.series(n=5, prec=10)
        except Exception as e:
            series = f"series failed: {e}"
        # Iwasawa invariants
        try:
            mu = Lp.mu_invariant()
        except Exception as e:
            mu = f"mu failed: {e}"
        try:
            lam = Lp.lambda_invariant()
        except Exception as e:
            lam = f"lambda failed: {e}"
        return {
            "p": p, "b_curve": b_curve, "ap": ap,
            "global_curve": str(E_global),
            "global_label": (E_global.cremona_label()
                             if hasattr(E_global, 'cremona_label') else ''),
            "j_invariant_global": str(E_global.j_invariant()),
            "mu_invariant": str(mu),
            "lambda_invariant": str(lam),
            "lp_series": str(series),
        }
    except Exception as exc:
        import traceback
        traceback.print_exc()
        return {"p": p, "global_curve": str(E_global),
                "error": str(exc)}


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3),
                 (97, 5), (103, 5)]:
        print(f"[phase36] p={p} b={b}")
        try:
            r = run_iwasawa(p, b)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": str(exc)}
        reports.append(r)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase36_iwasawa.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
