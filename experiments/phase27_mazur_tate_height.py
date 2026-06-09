"""
Phase 27 — Mazur–Tate cyclotomic p-adic height pairing.

The global Mazur–Tate height h_p : E(Q) × E(Q) → Q_p is bilinear:
    h_p(kP, kP) = k² · h_p(P, P)

If we could compute h_p(G, G) and h_p(Q, Q) for Q = kG with G, Q ∈ E(F_p)
via some F_p-computable proxy, the ratio would give k² and thus k.

Reality check: h_p is defined on E(Q), not E(F_p). The local
contribution at p (on points in Ê) is `-log_p(z(P))` — exactly what
Phase 26 computed. Phase 26 found NO polynomial-in-k structure in
that local height of δ(k).

This phase makes the connection explicit by computing h_p on actual
curves over Q with rank ≥ 1 and checking:
  (1) Does h_p(kP, kP) = k² h_p(P, P) hold numerically? (sanity)
  (2) Can any "p-adic discriminant" of the scalar be read off from
      h_p of a single point, i.e. is there a k-dependent quantity
      that reveals k mod p^i?
"""
from __future__ import annotations
import json
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cypari2
pari = cypari2.Pari()


def find_rank1_curve():
    """Return (E, P, label) for a small curve of rank ≥ 1 with a
    known generator P."""
    # E: y^2 = x^3 + 3, point (1, 2) has infinite order (easy check)
    E = pari.ellinit([0, 3])
    P = pari("[1, 2]")
    # Verify
    assert pari.ellisoncurve(E, P) == 1
    return E, P, "y^2 = x^3 + 3"


def bilinear_test(E, P, p, precision=8):
    """Verify h_p(kP, kP) = k² h_p(P, P) for small k."""
    try:
        hP = pari.ellpadicheight(E, p, precision, P)
    except Exception as exc:
        return {"status": f"height failed: {exc}"}
    results = []
    for k in range(1, 10):
        kP = pari.ellmul(E, P, k)
        if kP == 0:
            results.append((k, "identity"))
            continue
        try:
            hkP = pari.ellpadicheight(E, p, precision, kP)
        except Exception as exc:
            results.append((k, f"err: {exc}"))
            continue
        # The pairing height returns [h_cyc, h_anti_cyc] — two values
        # For bilinearity test, use first component
        # Compute k^2 * h(P) and compare
        try:
            h_first = hP[0] if pari.type(hP) == 't_VEC' else hP
            hk_first = hkP[0] if pari.type(hkP) == 't_VEC' else hkP
            expected = k * k * h_first
            diff = hk_first - expected
            results.append((k, "match" if diff == 0 else f"diff={diff}"))
        except Exception as exc:
            results.append((k, f"compare err: {exc}"))
    return {"bilinearity": results, "h_P": str(hP)}


def run(p):
    E, P, label = find_rank1_curve()
    # Skip primes of bad reduction or anomalous
    try:
        disc = pari.ellglobalred(E)
    except Exception:
        pass
    good = pari.ellap(E, p) % p != 0  # ordinary iff ap ≠ 0 mod p
    if not good:
        return {"prime": p, "status": "supersingular or bad", "ap": int(pari.ellap(E, p))}
    report = {"prime": p, "curve": label, "point": str(P),
              "ap_mod_p": int(pari.ellap(E, p) % p)}
    report.update(bilinear_test(E, P, p))
    return report


def main():
    primes = [5, 7, 11, 13, 17, 19, 23]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p in primes:
        print(f"[phase27] p={p}")
        try:
            report = run(p)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase27_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
