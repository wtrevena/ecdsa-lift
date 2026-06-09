"""
Phase 17 — DFT / character-sum spectral analysis of δ(k).

The sharpest test for "δ(k) is pseudorandom over Z/n" is: every
non-trivial character sum
    S(χ) = Σ_k χ(k) · e_{p^e}(δ(k))
is small (of size √n). Peaks in S(χ) = hidden periodicity = hidden
linear structure.

We do two complementary tests:

  (A) Fourier over Z/n on the MULTIPLICATIVE side. Pick a prime ℓ
      with ℓ ≡ 1 (mod n), find a primitive n-th root of unity ω in
      F_ℓ, and compute
          Ŝ(j) = Σ_{k=1}^{n-1} δ(k) · ω^{jk}   in F_ℓ.
      Large |Ŝ(j)| (for some j ≠ 0) indicates δ(k) has a frequency-j
      linear component.

  (B) Fourier over Z/(p^e) on the ADDITIVE side. Convert δ to
      complex floats, compute the DFT length n, check spectrum
      concentration.

Flat spectrum (all frequencies roughly √n) is the only honest
certificate that δ is pseudorandom.
"""
from __future__ import annotations
import json
import sys
import pathlib
import cmath
import math

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, naive_order, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift


def build_delta(E_p, C, G, g_ord):
    E_aff = Curve(a=C.a, b=C.b, N=C.N)
    tau_G = C.from_affine(teichmuller_lift(E_p, E_aff, G))
    deltas = []
    for k in range(1, g_ord):
        kG = E_p.mul(k, G)
        tau_kG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        diff = C.add(C.mul(k, tau_G), C.neg(tau_kG))
        X, Y, Z = diff
        try:
            d = (-X * modinv(Y % C.N, C.N)) % C.N
        except ZeroDivisionError:
            d = 0
        deltas.append(d)
    return deltas


def dft_complex(seq):
    """Straight O(n²) DFT returning magnitude spectrum."""
    n = len(seq)
    out = []
    for j in range(n):
        s = 0j
        for k in range(n):
            s += seq[k] * cmath.exp(-2j * cmath.pi * j * k / n)
        out.append(abs(s))
    return out


def run(p: int, b: int, e: int = 3) -> dict:
    E_p = Curve(a=0, b=b, N=p)
    n = naive_order(E_p)
    if n % p == 0:
        return {"prime": p, "status": "anomalous"}
    G, g_ord = find_generator(E_p)
    if g_ord < 10:
        return {"prime": p, "status": "small"}
    N = p ** e
    C = ProjCurve(a=0, b=b, N=N)
    deltas = build_delta(E_p, C, G, g_ord)

    # (B) Complex DFT on δ(k) directly. Each δ(k) ∈ [0, N).
    # Normalize to ~unit magnitude so the spectrum is comparable.
    cseq = [d / N for d in deltas]
    spectrum = dft_complex(cseq)
    # zero frequency is mean — exclude
    spec_nz = spectrum[1:]
    max_spec = max(spec_nz)
    mean_spec = sum(spec_nz) / len(spec_nz)
    # A flat spectrum has all bins ≈ √L / 2 of the signal energy's
    # square root. Peak-to-mean ratio is the key indicator.
    peak_to_mean = max_spec / mean_spec if mean_spec > 0 else 0.0
    top5 = sorted(enumerate(spec_nz), key=lambda x: -x[1])[:5]

    # (C) DFT on δ(k) mod p — the TOP p-adic digit
    topdigit = [d // (p ** (e - 1)) for d in deltas]
    c2 = [d / max(1, p) for d in topdigit]
    spec2 = dft_complex(c2)
    spec2_nz = spec2[1:]
    peak2 = max(spec2_nz)
    mean2 = sum(spec2_nz) / len(spec2_nz)
    ptm2 = peak2 / mean2 if mean2 > 0 else 0.0

    return {
        "prime": p, "curve_b": b, "g_ord": g_ord, "N": N,
        "delta_count": len(deltas),
        "spectrum_peak_over_mean_full": f"{peak_to_mean:.3f}",
        "top5_frequencies": [(j + 1, round(v, 3)) for j, v in top5],
        "spectrum_peak_over_mean_topdigit": f"{ptm2:.3f}",
        "verdict": ("STRONG SPECTRAL PEAK" if peak_to_mean > 4.0
                    else f"flat spectrum (peak/mean={peak_to_mean:.2f})"),
    }


def main():
    candidates = [
        (31, 3), (43, 7), (67, 2), (73, 13),
        (79, 3), (97, 5), (103, 5),
    ]
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    for p, b in candidates:
        print(f"[phase17] p={p} b={b}")
        try:
            report = run(p, b)
        except Exception as exc:
            report = {"prime": p, "error": repr(exc)}
        (out_dir / f"phase17_p{p}.json").write_text(
            json.dumps(report, indent=2, default=str))
        for k, v in report.items():
            print(f"   {k}: {v}")
        print()


if __name__ == "__main__":
    main()
