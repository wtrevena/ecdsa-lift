"""Phase 49 - Asymmetric-section KILL TEST. See report."""
from __future__ import annotations
import sys, pathlib, json, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.curves import Curve, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift, hensel_lift_y, hensel_lift_x_root
from experiments.phase37_gowers import gowers_norm_complex, encode_digit, random_baseline
SEED = 20260609; E_PREC = 4; DIGITS = (1, 2, 3); RAND = 300
TARGETS = [(67, 79), (211, 199), (823, 829), (10477, 10639)]

def _tonelli(c, p):
    if pow(c, (p - 1) // 2, p) != 1:
        return None
    if p % 4 == 3:
        return pow(c, (p + 1) // 4, p)
    q = p - 1; ss = 0
    while q % 2 == 0:
        q //= 2; ss += 1
    z = 2
    while pow(z, (p - 1) // 2, p) != p - 1:
        z += 1
    m = ss; cc = pow(z, q, p); t = pow(c, q, p); r = pow(c, (q + 1) // 2, p)
    while t != 1:
        i = 0; tt = t
        while tt != 1:
            tt = tt * tt % p; i += 1
        b = pow(cc, 1 << (m - i - 1), p)
        m = i; cc = b * b % p; t = t * cc % p; r = r * b % p
    return r


def fast_generator(E_p, n):
    """For prime group order n, any non-identity, non-2-torsion point of
    order n generates. Returns (G_fp, n). Assumes n prime."""
    p = E_p.N; a = E_p.a; b = E_p.b
    for x in range(p):
        c = (x * x * x + a * x + b) % p
        if c == 0:
            continue
        if pow(c, (p - 1) // 2, p) != 1:
            continue
        y = _tonelli(c, p)
        if y is None:
            continue
        G = (x, y)
        # order divides n (prime) and G != O, so order is exactly n
        return G, n
    raise RuntimeError("no generator found")



def teichmuller_x(E_p, E_lifted, x0):
    p = E_p.N
    N = E_lifted.N
    x = x0 % N
    for _ in range(64):
        xn = pow(x, p, N)
        if xn == x:
            break
        x = xn
    return x


def asymmetric_lift(E_p, E_lifted, P, canon="smaller_rep"):
    if P is None:
        return None
    x0, y0 = P
    p = E_p.N
    N = E_lifted.N
    if y0 == 0:
        return (hensel_lift_x_root(E_lifted, x0), 0)
    x = teichmuller_x(E_p, E_lifted, x0)
    y = hensel_lift_y(E_lifted, x, y0)
    y_alt = (-y) % N
    cands = [y, y_alt]
    if canon == "smaller_rep":
        ychosen = min(cands)
    elif canon == "even_lsd":
        evens = [c for c in cands if (c % p) % 2 == 0]
        ychosen = evens[0] if len(evens) == 1 else min(cands)
    else:
        raise ValueError(canon)
    return (x, ychosen)


def z_delta_sequence(p, b, e, section, canon=None, n_known=None):
    from math import gcd
    N = p ** e
    E_p = Curve(a=0, b=b, N=p)
    if n_known is not None:
        G_fp, n = fast_generator(E_p, n_known)
    else:
        G_fp, n = find_generator(E_p)
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)
    def lift(P):
        if section == "standard":
            return teichmuller_lift(E_p, E_aff, P)
        return asymmetric_lift(E_p, E_aff, P, canon=canon)
    tauG = C.from_affine(lift(G_fp))
    nT = C.mul(n, tauG)
    Cn = None
    if gcd(nT[1] % N, N) == 1:
        Cn = (-nT[0] * modinv(nT[1] % N, N)) % N
    kept = []
    zmap = {}
    for k in range(1, n):
        kG = E_p.mul(k, G_fp)
        tkG = C.from_affine(lift(kG))
        d = C.add(C.mul(k, tauG), C.neg(tkG))
        X, Y, Z = d
        if gcd(Y % N, N) != 1:
            continue
        z = (-X * modinv(Y % N, N)) % N
        if z % p != 0:
            continue
        kept.append((k, z))
        zmap[k] = z
    return kept, n, N, Cn, zmap


def antisymmetry_single_valued(zmap, n, N):
    sums = set()
    npairs = 0
    sample = []
    for k in list(zmap.keys()):
        partner = n - k
        if partner in zmap and k < partner:
            s = (zmap[k] + zmap[partner]) % N
            sums.add(s)
            npairs += 1
            if len(sample) < 8:
                sample.append(s)
    return (len(sums) == 1, len(sums), npairs, sorted(sample))


def u2_of_digit(values, p, j):
    f = np.array([encode_digit(int(v), p, j) for v in values], dtype=complex)
    return gowers_norm_complex(f, 2)


def main():
    t0 = time.time()
    np.random.seed(SEED)
    out = {"params": {"e": E_PREC, "digits": list(DIGITS), "random_baselines": RAND, "seed": SEED, "canon_rules": ["smaller_rep", "even_lsd"]}, "curves": []}
    for p, Nexp in TARGETS:
        print("[phase49] p=%d (expect n=%d)" % (p, Nexp))
        rec = {"p": p, "expected_n": Nexp, "n": Nexp, "sections": {}}
        for sec_label, canon in [("standard", None), ("asym_smaller_rep", "smaller_rep"), ("asym_even_lsd", "even_lsd")]:
            section = "standard" if sec_label == "standard" else "asym"
            kept, n, N, Cn, zmap = z_delta_sequence(p, 7, E_PREC, section, canon, n_known=(Nexp if p >= 5000 else None))
            assert n == Nexp, "order mismatch"
            zs = [z for _, z in kept]
            L = len(zs)
            single, ndist, npairs, sample = antisymmetry_single_valued(zmap, n, N)
            rmean, rstd, _ = random_baseline(L, 2, num_samples=RAND)
            digits = {}
            for j in DIGITS:
                u2 = u2_of_digit(zs, p, j)
                digits["d%d" % j] = {"U2_real": round(u2, 5), "U2_rand_mean": round(rmean, 5), "U2_rand_std": round(rstd, 5), "z_real_vs_rand": round((u2 - rmean) / rstd, 2)}
            rec["sections"][sec_label] = {"L": L, "dropped": (n - 1) - L, "Cn": Cn, "antisymmetry_single_valued": single, "antisymmetry_num_distinct_sums": ndist, "antisymmetry_num_pairs": npairs, "antisymmetry_sample_sums": sample, "digits": digits}
            z1 = digits["d1"]["z_real_vs_rand"]
            z2 = digits["d2"]["z_real_vs_rand"]
            z3 = digits["d3"]["z_real_vs_rand"]
            print("   %-18s L=%5d single=%-5s distinct_sums=%5d  z(real/rand) d1=%+.2f d2=%+.2f d3=%+.2f" % (sec_label, L, str(single), ndist, z1, z2, z3))
        out["curves"].append(rec)
    out["elapsed_s"] = round(time.time() - t0, 1)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "phase49_asymmetric_section.json").write_text(json.dumps(out, indent=2, default=str))
    print("Wrote results/phase49_asymmetric_section.json")


if __name__ == "__main__":
    main()
