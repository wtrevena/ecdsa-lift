"""Phase 55 - Fixed-shift discriminator (Prop 1 control)."""
from __future__ import annotations
import json, math, pathlib
import numpy as np

SEED = 20260609
NREP = 400


def mass(fv):
    F = np.fft.fft(fv) / len(fv)
    return float((np.abs(F) ** 4).sum())


def rand_seq(L, p, rng):
    return np.exp(2j * np.pi * rng.integers(0, p, L) / p)


def involution_seq(L, p, partner, rng):
    C = int(rng.integers(0, p)); inv2 = pow(2, -1, p)
    full = np.empty(L, dtype=np.int64); seen = np.zeros(L, dtype=bool)
    for k in range(L):
        if seen[k]:
            continue
        j = partner(k)
        if j == k:
            full[k] = (C * inv2) % p; seen[k] = True
        else:
            r = int(rng.integers(0, p)); full[k] = r; full[j] = (C - r) % p
            seen[k] = seen[j] = True
    return np.exp(2j * np.pi * full / p)


def chained_shift_seq(L, p, m, rng):
    g = math.gcd(L, m); ol = L // g
    full = np.empty(L, dtype=np.int64); seen = np.zeros(L, dtype=bool)
    C = int(rng.integers(0, p)); inv2 = pow(2, -1, p)
    for start in range(L):
        if seen[start]:
            continue
        if ol % 2 == 0:
            k = start; val = int(rng.integers(0, p))
            for step in range(ol):
                full[k] = val if step % 2 == 0 else (C - val) % p
                seen[k] = True; k = (k + m) % L
        else:
            k = start
            for _ in range(ol):
                full[k] = (C * inv2) % p; seen[k] = True; k = (k + m) % L
    return np.exp(2j * np.pi * full / p)


def measure(seqfun, L, p, rng):
    return float(np.mean([mass(seqfun(L, p, rng)) for _ in range(NREP)])) * L


def main():
    rng = np.random.default_rng(SEED); p = 10007; rows = []
    for L in [200, 612, 1000]:
        refl = lambda k, L=L: (L - k) % L
        trans = lambda k, L=L: (k + L // 2) % L
        mr = measure(rand_seq, L, p, rng)
        mref = measure(lambda L, p, r: involution_seq(L, p, refl, r), L, p, rng)
        mfix = measure(lambda L, p, r: involution_seq(L, p, trans, r), L, p, rng)
        m_long = 3 if math.gcd(L, 3) == 1 else 7
        m_art = measure(lambda L, p, r: chained_shift_seq(L, p, m_long, r), L, p, rng)
        rows.append({
            "L": L, "M*L_rand": round(mr, 4),
            "M*L_reflection(xi->-xi)": round(mref, 4),
            "ratio_reflection": round(mref / mr, 4),
            "M*L_fixedshift_involution(k<->k+L/2)": round(mfix, 4),
            "ratio_fixedshift": round(mfix / mr, 4),
            "caveat_chained_shift_m": m_long,
            "M*L_chained_shift_longorbit": round(m_art, 3),
            "ratio_chained_shift": round(m_art / mr, 3)})
        print(f"L={L}: rand {mr:.3f} | reflection {mref:.3f} (x{mref/mr:.3f}) | "
              f"fixedshift k<->k+L/2 {mfix:.3f} (x{mfix/mr:.3f}) | "
              f"[caveat chained m={m_long}: x{m_art/mr:.1f}]")
    out = {"phase": 55, "seed": SEED, "nrep": NREP, "p_alphabet": p,
           "claim": ("reflection (xi->-xi involution) inflates U^2 mass to 3/L "
                     "(ratio 3/2); fixed-shift involution k<->k+L/2 (xi->xi) leaves "
                     "it at 2/L (ratio ~1)."),
           "rows": rows,
           "caveat": ("chained_shift_longorbit is NOT the discriminator: chaining "
                      "s(k+m)=C-s(k) over a long single orbit forces a period-2 "
                      "pattern giving M*L ~ L/3 (separate periodicity artifact). The "
                      "clean fixed-shift involution is k<->k+L/2.")}
    res = pathlib.Path(__file__).resolve().parent.parent / "results"
    res.mkdir(exist_ok=True)
    (res / "phase55_fixedshift.json").write_text(json.dumps(out, indent=2))
    print("wrote results/phase55_fixedshift.json")


if __name__ == "__main__":
    main()
