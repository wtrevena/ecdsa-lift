"""
Phase 4 — The whole ballgame.

Group-cohomology argument:  H^2(E(F_p), Ê(pZ_p)) = 0 whenever
gcd(|E(F_p)|, p) = 1, because the coefficient group is a Z_p-module
and |G| is a unit. Therefore the short exact sequence

    0  ->  Ê(pZ_p)  ->  E(Z_p)  ->  E(F_p)  ->  0

splits. A homomorphic section s : E(F_p) -> E(Z_p) provably exists.

This experiment: build a candidate s via the averaging construction
applied to the Teichmüller lift τ, then TEST whether it is actually a
homomorphism. If yes, attack. If no, record exactly how it fails so
we know whether the averaging itself is the problem or the
Teichmüller choice of τ is.

We work in projective coordinates throughout so that the affine group
law's "must be units" requirement doesn't break intermediate sums of
points that reduce to the same F_p-point.
"""
from __future__ import annotations
import json
import random
import sys
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, points, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import splitting_section_proj
from src.formal_group import formal_log


def _proj_equal(C: ProjCurve, A, B) -> bool:
    return C.equal(A, B)


def check_homomorphism(E_p: Curve, C: ProjCurve, s_map):
    pts = points(E_p)
    total = ok = 0
    first_fail = None
    for P in pts:
        for Q in pts:
            total += 1
            PQ = E_p.add(P, Q)
            lhs = s_map[_key(PQ)]
            rhs = C.add(s_map[_key(P)], s_map[_key(Q)])
            if _proj_equal(C, lhs, rhs):
                ok += 1
            elif first_fail is None:
                first_fail = {
                    "P": list(P) if P else None,
                    "Q": list(Q) if Q else None,
                }
    return total, ok, first_fail


def _key(P):
    return None if P is None else (P[0], P[1])


def run(p: int, a: int, b: int, prec_e: int = 3) -> dict:
    E_p = Curve(a=a, b=b, N=p)
    N = p ** prec_e
    C = ProjCurve(a=a, b=b, N=N)

    n_group = naive_order(E_p)
    if n_group % p == 0:
        return {
            "prime": p,
            "status": "anomalous — Smart–SSSA already handles this",
        }

    G, g_ord = find_generator(E_p)
    if G is None:
        return {"prime": p, "error": "no generator"}

    # Build s on every point of E(F_p) using the averaging construction
    # in projective coordinates.
    pts = points(E_p)
    s_map: dict = {}
    for P in pts:
        s_map[_key(P)] = splitting_section_proj(E_p, C, P, prec_e)

    total, ok, first_fail = check_homomorphism(E_p, C, s_map)
    report = {
        "prime": p,
        "curve": {"a": a, "b": b},
        "n_group": n_group,
        "generator": list(G),
        "generator_order": g_ord,
        "precision": f"p^{prec_e}",
        "homomorphism_check": {
            "total": total,
            "ok": ok,
            "pass_rate": ok / total if total else None,
            "first_failure": first_fail,
        },
    }

    if total and ok == total:
        report["status"] = "HOMOMORPHIC — running attack"
        report["attack"] = _run_attack(E_p, C, s_map, G, g_ord, n_group, N)
    else:
        report["status"] = "not a homomorphism"
        # Sanity: confirm that s_map[P] at least REDUCES to P in E(F_p).
        reductions_ok = 0
        for P in pts:
            sP = s_map[_key(P)]
            try:
                aff = C.to_affine(sP)
            except ZeroDivisionError:
                aff = "z-not-unit"
            if aff is None and P is None:
                reductions_ok += 1
            elif aff == "z-not-unit":
                pass
            elif aff is not None and P is not None and (aff[0] - P[0]) % p == 0 and (aff[1] - P[1]) % p == 0:
                reductions_ok += 1
        report["reductions_correct"] = f"{reductions_ok}/{len(pts)}"

    return report


def _run_attack(E_p, C: ProjCurve, s_map, G, g_ord, n_group, N):
    sG = s_map[_key(G)]
    n_sG = C.mul(n_group, sG)
    try:
        aff = C.to_affine(n_sG)
    except ZeroDivisionError:
        return {"error": "[n]·s(G) not affinisable"}
    if aff is None:
        return {"error": "[n]·s(G) is the point at infinity"}
    # z-coordinate of a point in the formal group neighbourhood of O
    # is z = -x/y. For our purposes this formal-group element should
    # satisfy p | z.
    from src.curves import modinv as _modinv
    try:
        zG = (-aff[0] * _modinv(aff[1], N)) % N
    except ZeroDivisionError:
        return {"error": "y(n·s(G)) not a unit"}
    if zG % E_p.N != 0:
        return {"error": f"zG = {zG} not divisible by p"}
    logG = formal_log(C, zG)

    random.seed(0)
    k_secret = random.randrange(2, g_ord)
    Q = E_p.mul(k_secret, G)
    sQ = s_map[_key(Q)]
    n_sQ = C.mul(n_group, sQ)
    affQ = C.to_affine(n_sQ)
    if affQ is None:
        return {"error": "[n]·s(Q) is infinity"}
    try:
        zQ = (-affQ[0] * _modinv(affQ[1], N)) % N
    except ZeroDivisionError:
        return {"error": "y(n·s(Q)) not a unit"}
    logQ = formal_log(C, zQ)
    try:
        k_recovered = (logQ * _modinv(logG, N)) % N
    except ZeroDivisionError:
        return {"error": "log(G) not invertible"}
    return {
        "k_secret": k_secret,
        "k_recovered_mod_N": k_recovered,
        "k_recovered_mod_group_order": k_recovered % g_ord,
        "match": (k_recovered % g_ord) == k_secret,
    }


def main():
    # j=0 curves over primes p ≡ 1 (mod 3) whose (b) gives cyclic,
    # prime, non-anomalous group order. The pathological (Z/3)² cases
    # (where every Teichmüller lift is trivially 3-torsion) are
    # intentionally excluded.
    candidates = [
        (13, 0, 2),   # n = 19
        (19, 0, 2),   # n = 13
        (31, 0, 3),   # n = 43
        (43, 0, 7),   # n = 31
        (67, 0, 2),   # n = 73
        (73, 0, 13),  # n = 67
        (79, 0, 3),   # n = 97
        (97, 0, 5),   # n = 79
        (103, 0, 5),  # n = 97
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, a, b in candidates:
        print(f"[phase4] p={p}")
        try:
            report = run(p, a, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
            print(f"   ERROR: {exc}")
        (out_dir / f"phase4_p{p}.json").write_text(json.dumps(report, indent=2))
        print(f"   {report.get('status', report.get('error', '?'))}")
        hc = report.get("homomorphism_check")
        if hc:
            print(f"   hom: {hc['ok']}/{hc['total']}")
        if "attack" in report:
            print(f"   attack: {report['attack']}")


if __name__ == "__main__":
    main()
