"""
Phase 14f — Clean control: curves with PRIME g_ord across j-types.

Phase 14d verified rank = (g-1)/6 on j=0 curves with prime g_ord.
Phase 14e muddied the waters by using curves with non-prime g_ord.

This phase finds curves with prime g_ord for each j-type:
  • j = 0     (Aut = μ₆, expected rank = (g-1)/6)
  • j = 1728  (Aut = μ₄, expected rank = (g-1)/4)
  • j generic (Aut = ±1, expected rank = (g-1)/2)

If the prediction `rank × |Aut(E)| = g-1` holds across j-types, we've
fully characterized the leading-digit structure and closed the lead.
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


def is_prime(n):
    if n < 2: return False
    if n < 4: return True
    if n % 2 == 0: return False
    i = 3
    while i * i <= n:
        if n % i == 0: return False
        i += 2
    return True


def delta_leading(E_p, C, p, H, k):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    kH = E_p.mul(k, H)
    tau_H = C.from_affine(teichmuller_lift(E_p, E_aff, H))
    tau_kH = C.from_affine(teichmuller_lift(E_p, E_aff, kH))
    diff = C.add(C.mul(k, tau_H), C.neg(tau_kH))
    X, Y, Z = diff
    try:
        v = (-X * modinv(Y % C.N, C.N)) % C.N
    except ZeroDivisionError:
        return 0
    return (v // p) % p


def rank_over_Fp(M, p):
    M = [row[:] for row in M]
    rows = len(M); cols = len(M[0]) if M else 0
    rank = 0; col = 0
    while rank < rows and col < cols:
        pivot = None
        for r in range(rank, rows):
            if M[r][col] % p != 0:
                pivot = r; break
        if pivot is None:
            col += 1; continue
        M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col] % p, -1, p)
        M[rank] = [(x * inv) % p for x in M[rank]]
        for r in range(rows):
            if r == rank: continue
            f = M[r][col] % p
            if f == 0: continue
            M[r] = [(M[r][c] - f * M[rank][c]) % p for c in range(cols)]
        rank += 1; col += 1
    return rank


def search_curve_prime_ord(p, j_type, min_ord=20, max_ord=200):
    """Find (a, b) on F_p with specified j_type and prime g_ord."""
    for a in range(p):
        for b in range(p):
            if a == 0 and b == 0: continue
            disc = (-16 * (4 * a**3 + 27 * b**2)) % p
            if disc == 0: continue
            if j_type == 'j0' and a != 0: continue
            if j_type == 'j1728' and b != 0: continue
            if j_type == 'generic':
                if a == 0 or b == 0: continue
                num = (1728 * 4 * pow(a, 3, p)) % p
                den = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
                if den == 0: continue
                j = (num * pow(den, -1, p)) % p
                if j == 0 or j == 1728 % p: continue
            E = Curve(a=a, b=b, N=p)
            if E.N != p: continue
            try:
                n = naive_order(E)
            except Exception:
                continue
            if n % p == 0: continue
            if not (min_ord <= n <= max_ord): continue
            if not is_prime(n): continue
            G, g_ord = find_generator(E)
            if g_ord != n: continue
            return (a, b, E, G, g_ord)
    return None


def measure(p, j_type, e=3):
    found = search_curve_prime_ord(p, j_type)
    if found is None:
        return {"prime": p, "j_type": j_type, "status": "no prime-order curve"}
    a, b, E_p, G, g_ord = found
    size = min(g_ord - 2, 32)
    N = p ** e
    C = ProjCurve(a=a, b=b, N=N)
    ks = list(range(2, size + 2))
    M = [[delta_leading(E_p, C, p, E_p.mul(i, G), k) for k in ks]
         for i in range(1, size + 1)]
    rank = rank_over_Fp(M, p)
    exp6 = (g_ord - 1) // 6
    exp4 = (g_ord - 1) // 4
    exp2 = (g_ord - 1) // 2
    matches = []
    if rank == exp6: matches.append("μ₆")
    if rank == exp4: matches.append("μ₄")
    if rank == exp2: matches.append("±1")
    return {
        "prime": p, "j_type": j_type, "a": a, "b": b,
        "g_ord": g_ord, "size": size, "rank": rank,
        "exp_rank_mu6": exp6, "exp_rank_mu4": exp4, "exp_rank_pm1": exp2,
        "matches": matches or ["none"],
        "g_over_rank": round(g_ord / rank, 3) if rank else None,
    }


def main():
    # Several primes — search each for all three j-types
    primes = [37, 43, 47, 53, 61, 67, 73, 79, 83, 89]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    all_reports = []
    for p in primes:
        for jt in ['j0', 'j1728', 'generic']:
            print(f"[phase14f] p={p} j_type={jt}")
            try:
                report = measure(p, jt)
            except Exception as exc:
                report = {"prime": p, "j_type": jt, "error": repr(exc)}
            all_reports.append(report)
            for k, v in report.items():
                print(f"   {k}: {v}")
            print()
    (out_dir / "phase14f_all.json").write_text(
        json.dumps(all_reports, indent=2, default=str))
    # Summary table
    print("=" * 70)
    print(f"{'prime':>6} {'j_type':>10} {'g_ord':>6} {'rank':>6} {'g/rank':>8} {'matches':>15}")
    for r in all_reports:
        if 'rank' in r:
            print(f"{r['prime']:>6} {r['j_type']:>10} {r['g_ord']:>6} "
                  f"{r['rank']:>6} {r['g_over_rank']:>8} "
                  f"{','.join(r['matches']):>15}")


if __name__ == "__main__":
    main()
