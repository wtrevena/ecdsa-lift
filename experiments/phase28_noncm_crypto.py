"""
Phase 28 — Non-CM ordinary curves at (near-)cryptographic sizes via PARI.

All previous phases used small j=0 curves (CM by Z[ζ₃]). This phase tests
whether the structural features we found (rank deficit = (g-1)/|Aut|,
Phase 21b antisymmetry δ(k)+δ(n-k) = [n]τ(G), Phase 17 DFT pattern,
Phase 18 BM linear-recurrence absence) change when:

  (a) j ≠ 0, 1728  (non-CM, Aut = ±1)
  (b) p is large (32-bit → 128-bit)
  (c) g_ord is prime

The hypothesis to falsify: maybe at cryptographic scale, some new
structure emerges that small-prime experiments missed. If none appears,
the "attack is dead" verdict from Phase 14 extends to realistic sizes.

We use PARI for:
  • ellcard (Schoof's algorithm)
  • ellmul, elladd on E(F_p) and on E(Z/p^e) via direct formulas
  • ellformallog for the formal log series
"""
from __future__ import annotations
import json
import sys
import pathlib
import random
import time

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import cypari2
pari = cypari2.Pari()

from src.projective import ProjCurve


def modinv(a, n):
    g, x, _ = ext_gcd(a % n, n)
    if g != 1:
        raise ZeroDivisionError(f"no inverse of {a} mod {n}")
    return x % n


def ext_gcd(a, b):
    if b == 0:
        return a, 1, 0
    g, x1, y1 = ext_gcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def ec_add(P, Q, a, N):
    """Affine add on y^2 = x^3 + ax + b over Z/N. Points are (x,y) or None."""
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if (x1 - x2) % N == 0:
        if (y1 + y2) % N == 0:
            return None
        # double
        num = (3 * x1 * x1 + a) % N
        den = (2 * y1) % N
        lam = (num * modinv(den, N)) % N
    else:
        num = (y2 - y1) % N
        den = (x2 - x1) % N
        lam = (num * modinv(den, N)) % N
    x3 = (lam * lam - x1 - x2) % N
    y3 = (lam * (x1 - x3) - y1) % N
    return (x3, y3)


def ec_mul(k, P, a, N):
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R = ec_add(R, Q, a, N)
        Q = ec_add(Q, Q, a, N)
        k >>= 1
    return R


def teichmuller_point(P_fp, a_fp, b_fp, p, e):
    """Hensel-lift an F_p point to E(Z/p^e): fix x, solve for y to higher
    precision by Newton (y² = x³+ax+b).  Then iterate x as well via the
    canonical Teichmüller construction: we want [p]τ(P) = τ([p]P) in
    the limit — but a single-step Hensel lift of coordinates suffices
    for our purposes since p·Ê lives in higher filtration."""
    # Here we just use the naive lift: keep (x mod p, y mod p) and
    # Hensel-lift y. This is NOT the true Teichmüller — but δ is
    # measured as [k]τ(G) − τ(kG) with the SAME lift on both sides, so
    # it still quantifies the lifting obstruction for the chosen section.
    x0, y0 = P_fp
    N = p ** e
    x = x0 % N
    # solve y² = x³+ax+b in Z/p^e starting from y0 mod p
    target = (pow(x, 3, N) + a_fp * x + b_fp) % N
    y = y0 % p
    pk = p
    while pk < N:
        pk *= p
        # y' = y - (y²-target)/(2y)  mod pk
        f = (y * y - target) % pk
        df = (2 * y) % pk
        try:
            y = (y - f * modinv(df, pk)) % pk
        except ZeroDivisionError:
            break
    return (x % N, y % N)


def delta_z(G_fp, k, a_fp, b_fp, p, e, C=None):
    """z-coordinate of δ(k) = [k]τ(G) − τ(kG) in E(Z/p^e).
    Uses projective arithmetic so the kernel-of-reduction point is safe."""
    N = p ** e
    if C is None:
        C = ProjCurve(a=a_fp, b=b_fp, N=N)
    kG_fp = pari_ec_mul(k, G_fp, a_fp, b_fp, p)
    if kG_fp is None:
        return None
    tau_G = teichmuller_point(G_fp, a_fp, b_fp, p, e)
    tau_kG = teichmuller_point(kG_fp, a_fp, b_fp, p, e)
    PG = C.from_affine(tau_G)
    PkG = C.from_affine(tau_kG)
    kPG = C.mul(k, PG)
    diff = C.add(kPG, C.neg(PkG))
    X, Y, Z = diff
    # z = -X/Y in the formal group parameterization
    try:
        return (-X * modinv(Y % N, N)) % N
    except ZeroDivisionError:
        return None


def pari_ec_mul(k, P, a, b, p):
    E = pari.ellinit([int(a), int(b)], int(p))
    Q = pari([int(P[0]), int(P[1])])
    R = pari.ellmul(E, Q, int(k))
    if R == pari(0) or len(R) < 2:
        return None
    return (int(R[0]), int(R[1]))


def find_noncm_curve(p):
    """Find (a,b) over F_p with j≠0,1728, prime group order, and a
    generator G. Returns (a, b, G, n) or None."""
    E_global = None
    for attempt in range(200):
        a = random.randrange(1, p)
        b = random.randrange(1, p)
        disc = (-16 * (4 * pow(a, 3, p) + 27 * pow(b, 2, p))) % p
        if disc == 0:
            continue
        # j = 1728 * 4a^3 / (4a^3 + 27b^2)
        num = (1728 * 4 * pow(a, 3, p)) % p
        den = (4 * pow(a, 3, p) + 27 * pow(b, 2, p)) % p
        j = (num * modinv(den, p)) % p
        if j == 0 or j == 1728 % p:
            continue
        E = pari.ellinit([int(a), int(b)], int(p))
        n = int(pari.ellcard(E))
        if n % p == 0:
            continue  # anomalous
        # prime?
        if not pari.isprime(n):
            continue
        # find a point
        for _ in range(50):
            xg = random.randrange(p)
            rhs = (pow(xg, 3, p) + a * xg + b) % p
            ys = pari.Mod(rhs, p)
            try:
                yg = int(pari.sqrt(ys))
            except Exception:
                continue
            # verify
            if (yg * yg) % p != rhs:
                continue
            G = (xg, yg)
            return (a, b, G, n)
    return None


def phase21b_check(a_fp, b_fp, G, n, p, e=4):
    """Check the Phase 21b z-coordinate identity
        z(δ(k)) + z(δ(n−k)) ≡ z([n]τ(G))  mod p^e
    Returns hits/total and the constant."""
    N = p ** e
    C = ProjCurve(a=a_fp, b=b_fp, N=N)
    tauG = teichmuller_point(G, a_fp, b_fp, p, e)
    PG = C.from_affine(tauG)
    nTauG = C.mul(n, PG)
    try:
        const_z = (-nTauG[0] * modinv(nTauG[1] % N, N)) % N
    except ZeroDivisionError:
        const_z = None

    def z_delta(k):
        kG_fp = pari_ec_mul(k, G, a_fp, b_fp, p)
        if kG_fp is None:
            return None
        tau_kG = C.from_affine(teichmuller_point(kG_fp, a_fp, b_fp, p, e))
        diff = C.add(C.mul(k, PG), C.neg(tau_kG))
        X, Y, Z = diff
        try:
            return (-X * modinv(Y % N, N)) % N
        except ZeroDivisionError:
            return None

    ks = [2, 3, 5, 7, 11, 13, 17, 19]
    hits = 0
    total = 0
    samples = []
    for k in ks:
        if k >= n - 2:
            continue
        zk = z_delta(k)
        znk = z_delta(n - k)
        if zk is None or znk is None or const_z is None:
            continue
        total += 1
        s = (zk + znk) % N
        ok = (s == const_z)
        if ok:
            hits += 1
        samples.append((k, bool(ok)))
    return {"const_z_mod_p2": const_z % (p * p) if const_z is not None else None,
            "hits": f"{hits}/{total}",
            "samples": samples}


def run_one(bits, seed=None):
    if seed is not None:
        random.seed(seed)
    # Pick a random prime of given bit size
    lo = 1 << (bits - 1)
    hi = (1 << bits) - 1
    while True:
        p = random.randrange(lo | 1, hi, 2)
        if pari.isprime(p):
            p = int(p)
            break
    print(f"  p = {p}  ({bits} bits)")
    t0 = time.time()
    found = find_noncm_curve(p)
    if found is None:
        return {"bits": bits, "p": p, "status": "no curve found"}
    a, b, G, n = found
    print(f"  curve y²=x³+{a}x+{b}  order n={n}  prime={pari.isprime(n)}")
    print(f"  (search took {time.time()-t0:.1f}s)")

    # Phase 21b antisymmetry
    pb = phase21b_check(a, b, G, n, p, e=4)

    # Phase 14-style: leading-digit rank on a small block
    # M[i][k] = (δ_{iG}(k)/p) mod p  — do this only if p not too big
    # and we can afford size×size matrix
    size = 16
    N = p ** 3
    M = []
    for i in range(1, size + 1):
        iG = pari_ec_mul(i, G, a, b, p)
        if iG is None:
            M.append([0] * size)
            continue
        row = []
        for k in range(2, size + 2):
            z = delta_z(iG, k, a, b, p, e=3)
            if z is None:
                row.append(0)
            else:
                row.append((z // p) % p)
        M.append(row)

    # F_p rank
    def rank_fp(M, p):
        M = [row[:] for row in M]
        rows = len(M)
        cols = len(M[0]) if M else 0
        r = 0
        c = 0
        while r < rows and c < cols:
            piv = None
            for i in range(r, rows):
                if M[i][c] % p != 0:
                    piv = i
                    break
            if piv is None:
                c += 1
                continue
            M[r], M[piv] = M[piv], M[r]
            inv = pow(M[r][c] % p, -1, p)
            M[r] = [(x * inv) % p for x in M[r]]
            for i in range(rows):
                if i == r:
                    continue
                f = M[i][c] % p
                if f == 0:
                    continue
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[r])]
            r += 1
            c += 1
        return r

    rk = rank_fp(M, p)
    expected_generic = min(size, (n - 1) // 2)

    return {
        "bits": bits, "p": p, "a": a, "b": b, "n": n,
        "n_is_prime": bool(pari.isprime(n)),
        "size": size, "leading_rank": rk,
        "expected_pm1": expected_generic,
        "rank_matches_pm1": rk == expected_generic,
        "phase21b": pb,
    }


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    # 20, 24, 32, 40 bits — keep small enough for matrix + Hensel arithmetic
    for bits in [20, 24, 28, 32, 40]:
        print(f"[phase28] bits={bits}")
        try:
            r = run_one(bits, seed=bits * 1000 + 7)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"bits": bits, "error": repr(exc)}
        reports.append(r)
        for k, v in r.items():
            print(f"   {k}: {v}")
        print()
    (out_dir / "phase28_noncm_crypto.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
