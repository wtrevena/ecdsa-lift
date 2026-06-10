"""
Phase 43 — Referee-revision controls on a clean ORDINARY j=0 set.

Fixes raised in review:
  * Table 1 of the paper mixed regimes: some primes p=2 mod 3 give
    SUPERSINGULAR y^2=x^3+7 (N=p+1), and p=127 is ANOMALOUS (N=p). Here
    we use only ORDINARY j=0 curves: p = 1 mod 3 with N=#E(F_p) PRIME and
    N not in {p, p+1}. Each row reports its exact (p, n=N).
  * Pre-specified parameters (no post-hoc tuning):
        e = 4, digits j in {1,2,3}, SHUF = 300 increment shuffles,
        RAND = 300 S^1-uniform baselines, SYN = 300 synthetic draws,
        seed = 20260609. Exclusion rule, fixed in advance: a per-k value
        z(delta(k)) is DROPPED iff its projective lift is degenerate
        (Y not a unit, or z not in pZ — a kernel element must have
        z = 0 mod p); curves are reported regardless of size.
  * Positive control for the antisymmetry hypothesis H2: synthetic
    sequences that impose ONLY the reflection s(k)+s(n-k)=C on otherwise
    i.i.d.-uniform values. If reflection alone reproduces the U^2
    elevation and its growth with n, H2 is confirmed constructively.

Outputs results/phase43_revision_controls.json with, per curve:
  real U^2, increment-shuffle null, half-sequence U^2, synthetic-reflection
  U^2, S^1-random baseline, and peak Fourier magnitude (real vs synthetic)
  for the beta fit.
"""
from __future__ import annotations
import sys
import pathlib
import json
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, find_generator, modinv, naive_order
from src.projective import ProjCurve
from src.lifts import teichmuller_lift
from experiments.phase37_gowers import gowers_norm_complex, encode_digit, random_baseline

SEED = 20260609
E_PREC = 4
DIGITS = (1, 2, 3)
SHUF = 300
RAND = 300
SYN = 300

# clean ordinary j=0 curves (p = 1 mod 3, N prime, N != p, p+1)
ORDINARY = [(67, 79), (163, 139), (211, 199), (349, 313),
            (433, 397), (577, 613), (823, 829)]


def z_delta_clean(p, b, e):
    """z(delta(k)) for k=1..n-1 with the pre-specified degeneracy filter.
    Returns (list of (k, z) kept, n, N=p^e)."""
    N = p ** e
    E_p = Curve(a=0, b=b, N=p)
    G_fp, n = find_generator(E_p)
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)
    tauG = C.from_affine(teichmuller_lift(E_p, E_aff, G_fp))
    nT = C.mul(n, tauG)
    Cn = (-nT[0] * modinv(nT[1] % N, N)) % N
    kept = []
    for k in range(1, n):
        kG = E_p.mul(k, G_fp)
        tkG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        d = C.add(C.mul(k, tauG), C.neg(tkG))
        X, Y, Z = d
        from math import gcd
        if gcd(Y % N, N) != 1:
            continue
        z = (-X * modinv(Y % N, N)) % N
        if z % p != 0:                      # not a kernel element -> degenerate
            continue
        kept.append((k, z))
    return kept, n, N, Cn


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
    """U^2 of digit-j of synthetic sequences that impose ONLY the reflection
    s(k)+s(n-k)=Cn on i.i.d.-uniform values. n_idx = list of indices present
    (k = 1..n-1). We pair k with n-k.

    Returns (u2_array, peak_fourier_array): for each synthetic draw we also
    record the peak Fourier magnitude of its digit-j encoding (same definition
    as peak_fourier() on the real sequence), so the caller can report the mean
    synthetic peak magnitude for the beta fit. The RNG is consumed in exactly
    the same order as before (one draw per index per t), so all previously
    written U^2/z statistics reproduce byte-for-byte."""
    n = max(n_idx) + 1  # n-k pairing uses group order n
    idx = sorted(n_idx)
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
            if partner in set(idx):
                s[partner] = (Cn - r) % N
        vals = [s[k] for k in idx]
        out[t], _ = u2_of_digit(vals, p, j)
        pf[t] = peak_fourier(vals, p, j)
    return out, pf


def main():
    rng = np.random.default_rng(SEED)
    # random_baseline() in phase37 uses the legacy global numpy RNG; seed it
    # with SEED for reproducibility. This is independent of the default_rng(SEED)
    # stream above, so the U^2/synthetic statistics are unaffected.
    np.random.seed(SEED)
    out = {"params": {"e": E_PREC, "digits": list(DIGITS), "shuffles": SHUF,
                      "random_baselines": RAND, "synthetic": SYN, "seed": SEED},
           "curves": []}
    for p, Nexp in ORDINARY:
        print(f"[phase43] p={p} (expect n={Nexp})")
        kept, n, M, Cn = z_delta_clean(p, 7, E_PREC)
        assert n == Nexp, f"order mismatch {n}!={Nexp}"
        ks = [k for k, _ in kept]
        zs = [z for _, z in kept]
        L = len(zs)
        rec = {"p": p, "n": n, "L": L, "dropped": (n - 1) - L,
               "ordinary": naive_order(Curve(0, 7, p)) not in (p, p + 1),
               "digits": {}}
        rmean, rstd, _ = random_baseline(L, 2, num_samples=RAND)
        # half indices (40b): first half of the kept k's
        half = zs[:L // 2]
        for j in DIGITS:
            u2_real, _ = u2_of_digit(zs, p, j)
            inc = increment_shuffle_u2(zs, M, p, j, rng, SHUF)
            u2_half, _ = u2_of_digit(half, p, j)
            inc_half = increment_shuffle_u2(half, M, p, j, rng, SHUF)
            syn, syn_pf = synthetic_reflection_u2(ks, Cn, M, p, j, rng, SYN)
            pf_real = peak_fourier(zs, p, j)
            rec["digits"][f"d{j}"] = {
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
        d1 = rec["digits"]["d1"]
        print(f"   n={n} L={L} drop={rec['dropped']}  "
              f"real/inc={d1['z_real_vs_inc']:+.1f} inc/rand={d1['z_inc_vs_rand']:+.1f} "
              f"half/inc={d1['z_half_vs_inc']:+.1f} syn/rand={d1['z_syn_vs_rand']:+.1f} "
              f"real/syn={d1['z_real_vs_syn']:+.1f}")
        out["curves"].append(rec)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "phase43_revision_controls.json").write_text(
        json.dumps(out, indent=2, default=str))
    print("Wrote results/phase43_revision_controls.json")


if __name__ == "__main__":
    main()
