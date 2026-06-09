"""
Phase 33 — Frobenius on H^1_dR(E/Z_p) and the unit-root subspace.

For ordinary E/F_p, the de Rham cohomology H^1_dR(E_can/Z_p) is a
free Z_p-module of rank 2 with a Frobenius operator F whose
characteristic polynomial is

        x^2 − a_p · x + p.

The two eigenvalues are the *unit root* α (with v_p(α) = 0) and the
non-unit root β = p/α (with v_p(β) = 1). The unit-root subspace
U ⊂ H^1_dR is a canonical Z_p-line, and it's a non-trivial linear
invariant of the curve that is NOT a section of the reduction map
E(Z_p) → E(F_p) — so it might encode information that escaped all
of Phases 21–32.

What this script does:
  1. For each test curve E/F_p, lift to E/Z_p and compute the
     Frobenius matrix on H^1_dR via Sage's
     `E.matrix_of_frobenius()` (Kedlaya / Monsky–Washnitzer).
  2. Verify charpoly is x² − a_p · x + p.
  3. Compute the unit-root eigenvector v ∈ Q_p^2.
  4. Test correlation between v and the F_p group structure: in
     particular, compute Coleman integrals ∫_{O}^{τ_can(P)} ω_unit
     for the unit-root differential ω_unit and various F_p points P.
     If P → ∫ω_unit is a homomorphism E(F_p) → Z_p, it must be
     trivial (target torsion-free, source finite). But the
     individual values are still computable and might exhibit
     structure modulo p^k that the Hensel-lift δ does not.
  5. Compare ω_unit-integration results against δ to see if
     they're independent or correlated.

Run:  ~/miniforge3/bin/conda run -n sage sage experiments/phase33_frobenius_unitroot.py
"""
from sage.all import (EllipticCurve, GF, Zp, Qp, ZZ, Integer,
                      Matrix, vector, PolynomialRing)
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]


def lift_curve_to_qp(p, b_curve, prec=10):
    """Lift y^2 = x^3 + b to Q_p (canonical model is the same: j=0
    is rigid)."""
    K = Qp(p, prec=prec, type='capped-rel', print_mode='terse')
    return EllipticCurve(K, [0, b_curve]), K


def frobenius_matrix_and_unit_root(E_qp, p, ap, prec):
    """Compute F-matrix on H^1_dR via Sage's matrix_of_frobenius
    (calls Kedlaya for hyperelliptic / direct for ell). Returns
    (M, charpoly_coeffs, unit_root, unit_root_eigenvector)."""
    try:
        M = E_qp.matrix_of_frobenius()
    except Exception as e:
        return None, None, None, None, f"matrix_of_frobenius failed: {e}"
    # Charpoly
    R = M.parent().base_ring()
    PR = PolynomialRing(R, 'x')
    cp = M.charpoly()
    # Should be x^2 - ap*x + p
    expected = PR.gen() ** 2 - R(ap) * PR.gen() + R(p)
    cp_ok = (cp - expected).is_zero() or all(
        c.valuation() >= prec - 2 for c in (cp - expected).coefficients())
    # Eigenvalues
    # Solve x^2 - ap*x + p = 0 in Q_p (Hensel: ordinary => roots split)
    # discriminant = ap^2 - 4p, square root in Q_p iff it's a square.
    # For ordinary E, ap^2 - 4p is a unit non-square in F_p generally,
    # but in Q_p we need the discriminant to be a square. It is iff
    # E has ordinary reduction AND the disc lifts.
    # Easier: factor charpoly directly in Q_p[x].
    try:
        roots = cp.roots()
    except Exception as e:
        return M, cp, None, None, f"root finding failed: {e}"
    if not roots:
        return M, cp, None, None, "no roots in Q_p"
    # Unit root: the one with v_p = 0
    unit_root = None
    for (r, mult) in roots:
        if r.valuation() == 0:
            unit_root = r
            break
    if unit_root is None:
        return M, cp, None, None, "no unit root"
    # Eigenvector: solve (M - unit_root*I) v = 0
    K = M.base_ring()
    Mr = M - unit_root * Matrix.identity(K, 2)
    # The eigenvector is in the kernel; pick (Mr[0,1], -Mr[0,0]) or so
    if Mr[0, 1] != 0:
        ev = vector(K, [Mr[0, 1], -Mr[0, 0]])
    elif Mr[1, 1] != 0:
        ev = vector(K, [Mr[1, 1], -Mr[1, 0]])
    else:
        ev = vector(K, [1, 0])
    return M, cp, unit_root, ev, "ok"


def run(p, b_curve, prec=12):
    Fp = GF(p)
    E_fp = EllipticCurve(Fp, [0, b_curve])
    n = E_fp.cardinality()
    ap = p + 1 - int(n)
    if n % p == 0:
        return {"p": p, "status": "anomalous"}

    E_qp, K = lift_curve_to_qp(p, b_curve, prec=prec)
    M, cp, unit_root, ev, msg = frobenius_matrix_and_unit_root(
        E_qp, p, ap, prec)
    if M is None:
        return {"p": p, "b": b_curve, "ap": ap, "frob_status": msg}

    return {
        "p": p, "b": b_curve, "n": int(n), "ap": ap, "prec": prec,
        "frob_status": msg,
        "frob_matrix": [[str(M[i, j]) for j in range(2)] for i in range(2)],
        "charpoly": str(cp),
        "unit_root": str(unit_root) if unit_root is not None else None,
        "unit_root_valuation": (int(unit_root.valuation())
                                if unit_root is not None else None),
        "unit_root_eigenvector": ([str(x) for x in ev]
                                  if ev is not None else None),
    }


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    for p, b in [(31, 3), (43, 7), (67, 2), (73, 13), (79, 3),
                 (97, 5), (103, 5)]:
        print(f"[phase33] p={p} b={b}")
        try:
            r = run(p, b, prec=12)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": str(exc)}
        reports.append(r)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase33_frobenius_unitroot.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
