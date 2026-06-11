"""Phase 50 - Non-CM and j!=0 curves under the clean-set protocol.

Table 1 of the paper is all j=0 (y^2=x^3+7). Here we add curves with
other j-invariants and run the identical pre-registered clean-set
protocol (e=4, digits j in {1,2,3}, 300-draw nulls, seed 20260609):
  * j=1728 (y^2=x^3+ax). NOTE: y^2=x^3+ax always carries the 2-torsion
    point (0,0), so #E is ALWAYS even and prime full-order is impossible.
    We therefore run on the generator of the largest prime-order
    (odd) subgroup, giving a cyclic prime-order group with gcd(n,p)=1.
  * 3 generic non-CM curves (random a,b, j not in {0,1728}), prime
    full order n, p in 100-1000.

PREDICTION: the 3/2 inflation (real/inc ~ +(growing)) and real~=syn
equivalence hold identically, independent of CM.

Reports a Table-1-style row per curve: real/inc, inc/rand, half/inc,
syn/rand, real/syn z over digits 1,2,3.
"""
from __future__ import annotations
import sys, pathlib, json, time
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from src.curves import Curve, find_generator, modinv, naive_order, points, point_order
from src.projective import ProjCurve
from src.lifts import teichmuller_lift
from experiments.phase37_gowers import gowers_norm_complex, encode_digit, random_baseline
SEED = 20260609; E_PREC = 4; DIGITS = (1, 2, 3); SHUF = 300; RAND = 300; SYN = 300

# (p, a, b, n_sub, label) ; n_sub = order of the generator we use (prime).
# j=1728: full order even, use largest odd prime subgroup (order 101).
# generic: full order is prime, n_sub = full order.
CURVES = [
    (181, 7, 0, 101, "j=1728 (y^2=x^3+7x), sub-order 101 of #E=202"),
    (653, 390, 144, 613, "generic non-CM, j=446"),
    (163, 100, 84, 173, "generic non-CM, j=6"),
    (521, 51, 471, 479, "generic non-CM, j=316"),
]


def j_invariant(a, b, p):
    disc = (4 * a ** 3 + 27 * b * b) % p
    if disc % p == 0:
        return None
    return (1728 * 4 * a ** 3 * modinv(disc, p)) % p


def find_gen_of_order(E_p, m):
    """Return a point G of E(F_p) with exact order m (m | #E). Brute force
    over points; small p only."""
    n = naive_order(E_p)
    for P in points(E_p):
        if P is None:
            continue
        if point_order(E_p, P, n) == m:
            return P
    raise RuntimeError("no point of order %d" % m)

def z_delta_clean(p, a, b, e, n_sub):
    """z(delta(k)) for k=1..n_sub-1 on a generator of order n_sub, with the
    pre-registered degeneracy filter. Returns (kept (k,z), n_sub, N, Cn)."""
    from math import gcd
    N = p ** e
    E_p = Curve(a, b, p)
    G_fp = find_gen_of_order(E_p, n_sub)
    C = ProjCurve(a, b, N)
    E_aff = Curve(a, b, N)
    tauG = C.from_affine(teichmuller_lift(E_p, E_aff, G_fp))
    nT = C.mul(n_sub, tauG)
    Cn = (-nT[0] * modinv(nT[1] % N, N)) % N
    kept = []
    for k in range(1, n_sub):
        kG = E_p.mul(k, G_fp)
        tkG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        d = C.add(C.mul(k, tauG), C.neg(tkG))
        X, Y, Z = d
        if gcd(Y % N, N) != 1:
            continue
        z = (-X * modinv(Y % N, N)) % N
        if z % p != 0:
            continue
        kept.append((k, z))
    return kept, n_sub, N, Cn

def u2_of_digit(values, p, j):
    f = np.array([encode_digit(int(v), p, j) for v in values], dtype=complex)
    return gowers_norm_complex(f, 2), f


def peak_fourier(values, p, j):
    f = np.array([encode_digit(int(v), p, j) for v in values], dtype=complex)
    F = np.abs(np.fft.fft(f)) / len(f)
    return float(np.sort(F)[::-1][0])


def increment_shuffle_u2(seq, N, p, j, rng, num):
    s = [int(v) for v in seq]
    d = [(s[i + 1] - s[i]) % N for i in range(len(s) - 1)]
    base = s[0]
    out = np.empty(num)
    for t in range(num):
        dp = d[:]
        rng.shuffle(dp)
        acc = base
        sp = [base]
        for di in dp:
            acc = (acc + di) % N
            sp.append(acc)
        out[t], _ = u2_of_digit(sp, p, j)
    return out

def synthetic_reflection_u2(n_idx, Cn, N, p, j, rng, num):
    """U^2 of digit-j of synthetic sequences imposing ONLY the reflection
    s(k)+s(n-k)=Cn on i.i.d.-uniform values. Mirrors phase43 exactly."""
    n = max(n_idx) + 1
    idx = sorted(n_idx)
    idxset = set(idx)
    out = np.empty(num)
    pf = np.empty(num)
    for t in range(num):
        s = {}
        for k in idx:
            if k in s:
                continue
            partner = n - k
            r = int(rng.integers(0, N))
            s[k] = r
            if partner in idxset:
                s[partner] = (Cn - r) % N
        vals = [s[k] for k in idx]
        out[t], _ = u2_of_digit(vals, p, j)
        pf[t] = peak_fourier(vals, p, j)
    return out, pf

def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    np.random.seed(SEED)
    out = {"params": {"e": E_PREC, "digits": list(DIGITS), "shuffles": SHUF, "random_baselines": RAND, "synthetic": SYN, "seed": SEED}, "curves": []}
    for p, a, b, n_sub, label in CURVES:
        Efull = naive_order(Curve(a, b, p))
        jv = j_invariant(a, b, p)
        print("[phase50] p=%d a=%d b=%d  #E=%d  n_sub=%d  j=%s  %s" % (p, a, b, Efull, n_sub, jv, label))
        kept, n, M, Cn = z_delta_clean(p, a, b, E_PREC, n_sub)
        ks = [k for k, _ in kept]
        zs = [z for _, z in kept]
        L = len(zs)
        from math import gcd as _g
        rec = {"p": p, "a": a, "b": b, "full_order": Efull, "n": n, "j_invariant": jv, "label": label, "gcd_n_p": _g(n, p), "L": L, "dropped": (n - 1) - L, "Cn": Cn, "digits": {}}
        rmean, rstd, _ = random_baseline(L, 2, num_samples=RAND)
        half = zs[:L // 2]
        for j in DIGITS:
            u2_real, _ = u2_of_digit(zs, p, j)
            inc = increment_shuffle_u2(zs, M, p, j, rng, SHUF)
            u2_half, _ = u2_of_digit(half, p, j)
            inc_half = increment_shuffle_u2(half, M, p, j, rng, SHUF)
            syn, syn_pf = synthetic_reflection_u2(ks, Cn, M, p, j, rng, SYN)
            pf_real = peak_fourier(zs, p, j)
            rec["digits"]["d" + str(j)] = _digit_rec(u2_real, inc, u2_half, inc_half, syn, syn_pf, rmean, rstd, pf_real)
        d1 = rec["digits"]["d1"]
        print("   #E=%d n=%d L=%d drop=%d j=%s gcd=%d  real/inc=%+.1f inc/rand=%+.1f half/inc=%+.1f syn/rand=%+.1f real/syn=%+.1f" % (Efull, n, L, rec["dropped"], jv, rec["gcd_n_p"], d1["z_real_vs_inc"], d1["z_inc_vs_rand"], d1["z_half_vs_inc"], d1["z_syn_vs_rand"], d1["z_real_vs_syn"]))
        out["curves"].append(rec)
    out["elapsed_s"] = round(time.time() - t0, 1)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "phase50_noncm.json").write_text(json.dumps(out, indent=2, default=str))
    print("Wrote results/phase50_noncm.json")

def _digit_rec(u2_real, inc, u2_half, inc_half, syn, syn_pf, rmean, rstd, pf_real):
    return {
        "U2_real": round(u2_real, 5),
        "U2_inc_mean": round(float(inc.mean()), 5),
        "U2_inc_std": round(float(inc.std()), 5),
        "z_real_vs_inc": round((u2_real - inc.mean()) / inc.std(), 2),
        "z_real_vs_rand": round((u2_real - rmean) / rstd, 2),
        "z_inc_vs_rand": round((inc.mean() - rmean) / rstd, 2),
        "U2_half": round(u2_half, 5),
        "z_half_vs_inc": round((u2_half - inc_half.mean()) / inc_half.std(), 2),
        "U2_syn_mean": round(float(syn.mean()), 5),
        "U2_syn_std": round(float(syn.std()), 5),
        "z_syn_vs_rand": round((syn.mean() - rmean) / rstd, 2),
        "z_real_vs_syn": round((u2_real - syn.mean()) / syn.std(), 2),
        "peak_fourier_real": round(pf_real, 5),
        "peak_fourier_syn": round(float(syn_pf.mean()), 5),
    }


if __name__ == "__main__":
    main()
