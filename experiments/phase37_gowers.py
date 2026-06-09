"""
Phase 37 — Gowers U^k norms on the lift-error cocycle.

Phase 26 tested EXACT polynomial fits of z(δ(k)) in k. Negative.
Phase 29 tested Mahler expansion. Negative.
Both are *strict* notions of polynomial structure.

Gowers U^k norms test APPROXIMATE polynomial structure: a function
f: Z/n → C has large U^k norm iff it CORRELATES with a polynomial
phase of degree < k. For random f: Z/n → S^1, |f|_{U^k} ≈ n^{-1/2}
or smaller. For structured f, |f|_{U^k} is bounded below by a
constant.

Inverse theorem (Green-Tao-Ziegler): if |f|_{U^k} > δ then f
correlates with a (k-1)-step nilsequence. For k=2 this is a Fourier
coefficient; for k=3 it's a quadratic phase e^{2πi α k²}; for k=4
it's a cubic phase or a 3-step nilsequence.

Our test:
  Compute z(δ(k)) for k = 1..n-1 on a generator G.
  Form three encodings:
    f_top(k) = exp(2πi · (z(δ(k)) // p^(e-1)) / p)   # top digit
    f_bot(k) = exp(2πi · (z(δ(k)) % p) / p)            # bottom digit
    f_mid(k) = exp(2πi · ((z(δ(k)) // p) % p) / p)     # middle digit
  Compute |f|_{U^2}, |f|_{U^3}, |f|_{U^4} for each, plus the
  Fourier-spectrum (largest Fourier coefficient) for U^2 confirmation.
  Compare against the empirical distribution of |g|_{U^k} for random
  g: Z/n → S^1.

If any |f|_{U^k} is anomalously large (e.g. > random_mean + 3*std),
this is a signal of approximate polynomial structure that all
previous phases missed.
"""
from __future__ import annotations
import json
import sys
import pathlib
import math
import cmath
import random
import itertools
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def z_of_delta(E_p, C, G_proj, k):
    """Compute z(δ(k)) = z([k]·τ(G) − τ(kG)) ∈ Z/N."""
    G_fp = C.to_affine(G_proj) if G_proj[2] != 0 else None
    # We need the F_p generator; recompute by reduction
    N = C.N
    p = next(p for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43,
                          47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
                          103, 107] if N % p == 0 and (N // p) % p == 0)
    # Easier: assume we know p
    raise NotImplementedError("use the cleaner z_delta_seq below")


def z_delta_sequence(p, b, e, max_k):
    """Compute z(δ(k)) for k = 1..max_k on the canonical generator."""
    N = p ** e
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    G_fp, g_ord = find_generator(E_p)
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)

    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G_fp))
    seq = []
    for k in range(1, max_k + 1):
        kG_fp = E_p.mul(k, G_fp)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG_fp))
        diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
        X, Y, Z = diff
        try:
            z = (-X * modinv(Y % N, N)) % N
        except ZeroDivisionError:
            z = None
        seq.append(z)
    return seq, g_ord, N


def _uk_pow(g, k):
    """Return |g|_{U^k}^{2^k} via the recursive identity
        |g|_{U^k}^{2^k} = E_h |Δ_h g|_{U^{k-1}}^{2^{k-1}}
    with Δ_h g(x) = g(x+h) · conj(g(x)).
    Base case k=1:  |g|_{U^1}^2 = |E_x g(x)|^2.
    Complexity: O(n^k)."""
    n = len(g)
    if k == 1:
        m = g.mean()
        return float((m * m.conjugate()).real)
    if k == 2:
        # FFT-based: |g|_{U^2}^4 = (1/n) Σ_h |a(h)|^2 where
        # a(h) = (1/n) Σ_x g(x+h) conj(g(x))
        F = np.fft.fft(g)
        a = np.fft.ifft(F * np.conj(F)) / n  # length-n correlation, normalized
        return float(np.sum(np.abs(a) ** 2).real / n)
    # k >= 3: recurse over shifts h
    cg = np.conj(g)
    total = 0.0
    for h in range(n):
        dh = np.roll(g, -h) * cg  # Δ_h g(x) = g(x+h) · conj(g(x))
        total += _uk_pow(dh, k - 1)
    return total / n


def gowers_norm_complex(f, k):
    """Compute |f|_{U^k} for f: Z/n → C, given as a list of complex numbers.
    Uses the recursive derivative identity (FFT for k=2, direct for k≥3).
    Total complexity: O(n^{k-1} log n) for k=2, O(n^{k-1}) for k≥3."""
    g = np.asarray(f, dtype=complex)
    n = len(g)
    if k == 0:
        return float(abs(g.sum()) / n)
    val = _uk_pow(g, k)
    if val < 0:
        val = abs(val)
    return val ** (1.0 / (2 ** k))


def encode_digit(z, p, j):
    """Encode the j-th p-adic digit of z as a phase in S^1.
    Note: z(δ(k)) is divisible by p (δ ∈ kernel of reduction), so digit 0
    is identically 0 — useless. Use j ≥ 1.
    """
    if z is None:
        return 1.0 + 0j
    d = (z // (p ** j)) % p
    return cmath.exp(2j * math.pi * d / p)


def encode_mod_n(z, p, e, n):
    if z is None:
        return 1.0 + 0j
    return cmath.exp(2j * math.pi * (z % n) / n)


def random_baseline(n, k, num_samples=200):
    """Empirical distribution of |g|_{U^k} for random g: Z/n → S^1."""
    samples = np.empty(num_samples)
    for i in range(num_samples):
        phases = np.random.random(n) * (2 * np.pi)
        g = np.exp(1j * phases)
        samples[i] = gowers_norm_complex(g, k)
    return float(samples.mean()), float(samples.std()), float(samples.max())


def fourier_spectrum(f):
    """Return the magnitudes of the discrete Fourier transform of f."""
    n = len(f)
    spec = []
    for k in range(n):
        s = sum(f[j] * cmath.exp(-2j * math.pi * j * k / n) for j in range(n))
        spec.append(abs(s) / n)
    return spec


def run(p, b, e=4, k_max=4, num_random=200):
    E_p = Curve(a=0, b=b, N=p)
    G_fp, g_ord = find_generator(E_p)
    seq, _, N = z_delta_sequence(p, b, e, max_k=g_ord - 1)
    n = len(seq)
    if n < 8:
        return {"p": p, "status": "too small", "g_ord": g_ord}

    # Use digits 1..e-1 of z(δ(k)). Digit 0 is identically zero (artifact:
    # δ ∈ kernel of reduction so z(δ) ≡ 0 mod p), so we exclude it.
    encodings = [(f"d{j}", [encode_digit(z, p, j) for z in seq])
                 for j in range(1, e)]

    result = {"p": p, "b": b, "g_ord": g_ord, "n": n, "e": e}
    for name, f in encodings:
        norms = {}
        for k in range(2, k_max + 1):
            try:
                norms[f"U{k}"] = gowers_norm_complex(f, k)
            except Exception as exc:
                norms[f"U{k}"] = f"err: {exc}"
        F = np.fft.fft(np.asarray(f, dtype=complex)) / n
        mags = np.abs(F)
        sorted_mags = np.sort(mags)[::-1]
        norms["max_fourier"] = float(sorted_mags[0])
        norms["second_fourier"] = float(sorted_mags[1]) if len(sorted_mags) > 1 else 0.0
        result[f"f_{name}"] = norms

    # Random baselines (one per (n, k))
    baselines = {}
    for k in range(2, k_max + 1):
        mean, std, mx = random_baseline(n, k, num_samples=num_random)
        baselines[f"U{k}_random_mean"] = mean
        baselines[f"U{k}_random_std"] = std
        baselines[f"U{k}_random_max"] = mx
    # Random baseline for max Fourier coefficient
    fourier_max_samples = []
    for _ in range(num_random):
        g = np.exp(1j * np.random.random(n) * 2 * np.pi)
        F = np.fft.fft(g) / n
        fourier_max_samples.append(float(np.max(np.abs(F))))
    baselines["fourier_max_random_mean"] = float(np.mean(fourier_max_samples))
    baselines["fourier_max_random_std"] = float(np.std(fourier_max_samples))
    result["baselines"] = baselines

    # Z-scores
    zscores = {}
    for name, _ in encodings:
        for k in range(2, k_max + 1):
            v = result[f"f_{name}"].get(f"U{k}")
            if isinstance(v, (int, float)):
                z = (v - baselines[f"U{k}_random_mean"]) / baselines[f"U{k}_random_std"]
                zscores[f"f_{name}_U{k}_z"] = round(z, 3)
        v = result[f"f_{name}"]["max_fourier"]
        z = (v - baselines["fourier_max_random_mean"]) / baselines["fourier_max_random_std"]
        zscores[f"f_{name}_fourier_z"] = round(z, 3)
    result["zscores"] = zscores

    return result


SECP_FAMILY = [(p, 7) for p in
    (31, 43, 67, 79, 97, 103, 127, 167, 211, 257, 311, 401, 521, 631, 751, 1009)]


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    reports = []
    # Scaling sweep on the secp256k1 family y^2 = x^3 + 7
    for p, b in SECP_FAMILY:
        print(f"[phase37] p={p} b={b}")
        try:
            r = run(p, b, e=4, k_max=3, num_random=200)
        except Exception as exc:
            import traceback
            traceback.print_exc()
            r = {"p": p, "error": repr(exc)}
        reports.append(r)
        if "f_d1" in r:
            print(f"   n={r['n']}  g_ord={r['g_ord']}")
            for j in (1, 2, 3):
                fk = r[f"f_d{j}"]
                print(f"   d{j}: U2={fk['U2']:.4f}  U3={fk['U3']:.4f}  "
                      f"max_F={fk['max_fourier']:.4f}")
            print(f"   rand U2={r['baselines']['U2_random_mean']:.4f}"
                  f"±{r['baselines']['U2_random_std']:.4f}  "
                  f"max_F={r['baselines']['fourier_max_random_mean']:.4f}"
                  f"±{r['baselines']['fourier_max_random_std']:.4f}")
            print(f"   z: {r['zscores']}")
        else:
            print(f"   {r}")
        print()
    (out_dir / "phase37_gowers.json").write_text(
        json.dumps(reports, indent=2, default=str))


if __name__ == "__main__":
    main()
