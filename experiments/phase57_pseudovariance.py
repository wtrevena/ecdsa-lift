"""
Phase 57 -- Pseudo-variance classification of the U^2 mass-inflation factor.

THEOREM (per-mode pseudo-variance => mass inflation).
Each nonzero DFT mode Z_xi = hat f(xi) is, in the large-L Gaussian limit, a
mean-zero complex Gaussian with variance sigma^2(xi)=E|Z_xi|^2 and
pseudo-variance rho(xi)=E[Z_xi^2]. Such a Gaussian has
    E|Z|^4 = 2 sigma^4 + |rho|^2,
so M(f)=sum_xi E|Z_xi|^4 and the inflation factor over the proper baseline is
    INFLATION = 1 + <|rho(xi)|^2> / (2 sigma^4).
Since |rho|<=sigma^2, the factor lies in [1, 3/2]: 3/2 iff |rho|=sigma^2 a.e.
(reflection), 1 iff rho==0 (proper / shift).
"""
from __future__ import annotations
import json, pathlib
import numpy as np
import sympy as sp

ROOT = pathlib.Path(__file__).resolve().parents[1]
SEED = 20260609


def symbolic_fourth_moment():
    a, b, c = sp.symbols('a b c', real=True)
    EX4, EY4, EX2Y2 = 3*a**2, 3*b**2, a*b + 2*c**2
    E_Z4 = EX4 + EY4 + 2*EX2Y2
    sigma2 = a + b
    target = 2*sigma2**2 + (a-b)**2 + 4*c**2
    diff = sp.simplify(E_Z4 - target)
    return {"E_Z4": str(sp.expand(E_Z4)),
            "2sigma4_plus_rho2": str(sp.expand(target)),
            "difference": str(diff),
            "identity_holds": bool(diff == 0)}


def U2_mass(fv):
    L = len(fv)
    F = np.fft.fft(fv) / L
    return float(np.sum(np.abs(F) ** 4))


def _dphase(d, p):
    return np.exp(2j*np.pi*np.asarray(d)/p)


def _reflect(L, rng, beta, p):
    h = L // 2
    d = rng.integers(0, p, h)
    f = np.empty(L, dtype=complex)
    f[:h] = _dphase(d, p)
    f[L-h:] = beta * np.conj(_dphase(d[::-1], p))
    if L % 2 == 1:
        f[h] = np.exp(2j*np.pi*rng.random())
    return f


def _shift(L, rng, p):
    s = max(1, L // 3)
    omega = np.exp(2j*np.pi*rng.random())
    g = np.exp(2j*np.pi*rng.random(L))
    return omega * np.roll(g, -s)


def _partial_reflect(L, rng, q, p):
    h = L // 2
    d = rng.integers(0, p, h)
    f = np.empty(L, dtype=complex)
    f[:h] = _dphase(d, p)
    partner = np.conj(_dphase(d[::-1], p))
    indep = _dphase(rng.integers(0, p, h), p)
    coupled = rng.random(h) > q
    f[L-h:] = np.where(coupled, partner, indep)
    if L % 2 == 1:
        f[h] = np.exp(2j*np.pi*rng.random())
    return f


def measure_constraint(kind, L, rng, nrep=3000, p=10007, q=0.5, frac=0.5):
    masses = np.empty(nrep)
    abs2 = np.zeros(L)
    sq = np.zeros(L, dtype=complex)
    for t in range(nrep):
        if kind == "C0_none":
            f = np.exp(2j*np.pi*rng.random(L))
        elif kind == "C1_reflection":
            f = _reflect(L, rng, 1.0, p)
        elif kind == "C1p_reflection_phase":
            f = _reflect(L, rng, np.exp(2j*np.pi*rng.random()), p)
        elif kind == "C3_real":
            f = _reflect(L, rng, 1.0, p)
        elif kind == "C2_shift":
            f = _shift(L, rng, p)
        elif kind == "C4_partial":
            f = _partial_reflect(L, rng, q, p)
        elif kind == "C5_block":
            f = _partial_reflect(L, rng, 1.0 - frac, p)
        else:
            raise ValueError(kind)
        masses[t] = U2_mass(f)
        F = np.fft.fft(f) / L
        abs2 += np.abs(F)**2
        sq += F**2
    sigma2 = abs2 / nrep
    rho = sq / nrep
    s2 = sigma2[1:]; rh = np.abs(rho[1:])
    mean_s4 = float(np.mean(s2**2))
    mean_rho2 = float(np.mean(rh**2))
    infl_rho = 1 + mean_rho2/(2*mean_s4) if mean_s4 > 0 else float('nan')
    return {"kind": kind, "L": L,
            "ratio_rho_over_sigma2": round(float(np.mean(rh)/np.mean(s2)), 4),
            "inflation_from_rho": round(infl_rho, 4),
            "inflation_from_mass": round(float(np.mean(masses))*L/2.0, 4)}


def main():
    rng = np.random.default_rng(SEED)
    out = {"phase": 57, "seed": SEED,
           "symbolic_fourth_moment": symbolic_fourth_moment()}
    L = 400
    catalogue = [
        ("C0_none", {}), ("C1_reflection", {}), ("C1p_reflection_phase", {}),
        ("C3_real", {}), ("C2_shift", {}),
        ("C4_partial", {"q": 0.0}), ("C4_partial", {"q": 0.25}),
        ("C4_partial", {"q": 0.5}), ("C4_partial", {"q": 0.75}),
        ("C4_partial", {"q": 1.0}), ("C5_block", {"frac": 0.5}),
    ]
    rows = []
    for kind, kw in catalogue:
        r = measure_constraint(kind, L, rng, nrep=3000, **kw)
        r.update(kw)
        if kind in ("C1_reflection", "C1p_reflection_phase", "C3_real"):
            r["inflation_theory"] = 1.5
        elif kind in ("C0_none", "C2_shift"):
            r["inflation_theory"] = 1.0
        elif kind == "C4_partial":
            r["inflation_theory"] = round(1 + (1 - kw["q"])**2/2, 4)
        elif kind == "C5_block":
            r["inflation_theory"] = round(1 + kw["frac"]**2/2, 4)
        rows.append(r)
        print("%-22s %-14s rho/s2=%.3f infl(mass)=%.3f infl(rho)=%.3f theory=%s" %
              (kind, str(kw), r["ratio_rho_over_sigma2"], r["inflation_from_mass"],
               r["inflation_from_rho"], r["inflation_theory"]))
    out["catalogue"] = rows
    (ROOT/"results").mkdir(exist_ok=True)
    (ROOT/"results"/"phase57_pseudovariance.json").write_text(json.dumps(out, indent=2))
    print("\nsymbolic identity holds:", out["symbolic_fourth_moment"]["identity_holds"])
    print("wrote results/phase57_pseudovariance.json")


if __name__ == "__main__":
    main()
