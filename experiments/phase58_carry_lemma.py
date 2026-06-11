"""
Phase 58 -- Noisy-reflection / carry lemma.

Proposition 1 idealises the digit reflection as EXACT:
    f(L+1-k) = beta * conj f(k)      (carry-free).
The real digit sequence digit_j(z(delta(k))) obeys an exact reflection on
the INTEGER z (z(delta(k)) + z(delta(n-k)) = Cn mod N, Phase 21b/40b), but
base-p SUBTRACTION to recover a single digit introduces CARRIES: digit_j of
a sum/difference is not the (mod p) sum/difference of digits. So the phase-
level reflection holds only for a fraction (1-q_j) of pairs, and acts like
independent noise on the remaining fraction q_j.

LEMMA (noisy reflection).  Model f on Z/L: for a fraction (1-q) of the
reflection pairs (k, L+1-k) the exact relation f(L+1-k)=beta conj f(k)
holds; for the remaining fraction q the two entries are independent
S^1-uniform. Then in the Gaussian limit each mode has pseudo-variance
|rho(xi)| = (1-q) sigma^2, so by the Phase-57 master formula

    E[ ||f||_{U^2}^4 ] = (2 + (1-q)^2) / L  + o(1/L),

i.e. the mass ratio vs the i.i.d. baseline (2/L) is (2 + (1-q)^2)/2,
interpolating between 1 (q=1, no reflection) and 3/2 (q=0, exact Prop 1).

This file:
  (A) DERIVES q_j from REAL delta data: for each clean ordinary curve and
      digit j, q_j = fraction of pairs (k,n-k) for which the phase-level
      reflection digit_j(z(delta(n-k))) = (Cj - digit_j(z(delta(k)))) mod p
      FAILS (Cj fitted as the modal residual). Also gives a first-principles
      carry estimate q_j ~ (p-1)/p * (1 - p^{-1}) ... (see carry_prob).
  (B) VERIFIES the (2+(1-q)^2)/L law synthetically across q.
  (C) Predicts the real U^2 mass from the measured q_j and compares to the
      actually-measured real U^2 mass (Phase-43 numbers / recomputed here),
      converting "residual = finite-size noise" into "residual = carries".
"""
from __future__ import annotations
import sys, json, math, pathlib
from collections import Counter
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, find_generator, modinv, naive_order
from src.projective import ProjCurve
from src.lifts import teichmuller_lift

SEED = 20260609
E_PREC = 4
DIGITS = (1, 2, 3)
# clean ordinary j=0 curves (p = 1 mod 3, N prime), same set as Phase 43
ORDINARY = [(67, 79), (163, 139), (211, 199), (349, 313),
            (433, 397), (577, 613), (823, 829)]


# ---------------- real delta data ----------------
def z_delta_clean(p, b, e):
    """z(delta(k)) for k=1..n-1 with the Phase-43 degeneracy filter.
    Returns kept=[(k,z)], n, N, Cn (the integer reflection constant)."""
    N = p ** e
    E_p = Curve(a=0, b=b, N=p)
    G_fp, n = find_generator(E_p)
    C = ProjCurve(a=0, b=b, N=N)
    E_aff = Curve(a=0, b=b, N=N)
    tauG = C.from_affine(teichmuller_lift(E_p, E_aff, G_fp))
    nT = C.mul(n, tauG)
    Cn = (-nT[0] * modinv(nT[1] % N, N)) % N
    kept = []
    from math import gcd
    for k in range(1, n):
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
    return kept, n, N, Cn


def digit(z, p, j):
    return (z // (p ** j)) % p


def measure_q(p, b, e):
    """For each digit j, measure q_j = fraction of reflection pairs where the
    PHASE-LEVEL reflection digit_j(z(n-k)) == (Cj - digit_j(z(k))) mod p FAILS.
    Cj is fitted as the modal value of (digit_j(z(k)) + digit_j(z(n-k))) mod p.
    Only pairs with BOTH partners present (after the degeneracy filter) count.
    """
    kept, n, N, Cn = z_delta_clean(p, b, e)
    zmap = {k: z for k, z in kept}
    present = set(zmap)
    out = {}
    for j in DIGITS:
        sums = []
        pairs = []
        for k in present:
            partner = n - k
            if partner in present and partner > k:
                dk = digit(zmap[k], p, j)
                dp = digit(zmap[partner], p, j)
                sums.append((dk + dp) % p)
                pairs.append((dk, dp))
        if not pairs:
            out[f"d{j}"] = None
            continue
        Cj = Counter(sums).most_common(1)[0][0]
        fail = sum(1 for (dk, dp) in pairs if (dk + dp) % p != Cj)
        out[f"d{j}"] = {
            "n_pairs": len(pairs),
            "Cj_modal": int(Cj),
            "q_measured": round(fail / len(pairs), 4),
            "matched_frac_(1-q)": round(1 - fail / len(pairs), 4),
        }
    return out, n, len(kept), N


def carry_prob_first_principles(p, j, e):
    """Heuristic first-principles q for digit j of a base-p difference.

    The exact integer reflection is z(k)+z(n-k) = Cn (mod p^e). Recovering
    digit_j of z(n-k) from digit_j of z(k) requires Cn - z(k) with no carry
    crossing position j. For digits of two 'random' base-p integers, a borrow
    propagates into position j unless every lower position already determines
    it. The probability that the digit-j relation is *exactly* the mod-p
    complement (no net borrow disturbance at position j) is, to leading order,
    the probability that the borrow into position j is the SAME on both sides.
    A borrow enters position j with probability ~ 1/2 for uniform digits (the
    classic 'half the time there is a carry' heuristic), and the two sides'
    borrows agree with probability ~1/2 + 1/(2p) (they share the constant Cn).
    Net mismatch probability:
        q_j ~ 1/2 * (1 - 1/p)             for the interior digits 1..e-2,
        q_j ~ 0                            for the top digit j=e-1 (no higher
                                            position to borrow from; the
                                            reflection is exact there).
    This is heuristic; the measured q (above) is the ground truth.
    """
    if j >= e - 1:
        return 0.0
    return 0.5 * (1 - 1.0 / p)


# ---------------- synthetic verification of (2+(1-q)^2)/L ----------------
def U2_mass(fv):
    L = len(fv)
    F = np.fft.fft(fv) / L
    return float(np.sum(np.abs(F) ** 4))


def synth_noisy_reflection_mass(L, q, p_alphabet, rng, nrep=3000):
    h = L // 2
    masses = np.empty(nrep)
    for t in range(nrep):
        d = rng.integers(0, p_alphabet, h)
        f = np.empty(L, dtype=complex)
        f[:h] = np.exp(2j*np.pi*d/p_alphabet)
        partner = np.conj(np.exp(2j*np.pi*d[::-1]/p_alphabet))
        indep = np.exp(2j*np.pi*rng.integers(0, p_alphabet, h)/p_alphabet)
        coupled = rng.random(h) > q
        f[L-h:] = np.where(coupled, partner, indep)
        if L % 2 == 1:
            f[h] = np.exp(2j*np.pi*rng.random())
        masses[t] = U2_mass(f)
    return float(np.mean(masses)) * L


def real_u2_mass(p, b, e, j):
    """Mass M(f)=||f||_{U2}^4 of the real digit-j sequence (no synthetics)."""
    kept, n, N, Cn = z_delta_clean(p, b, e)
    zs = [z for _, z in kept]
    f = np.array([np.exp(2j*np.pi*digit(z, p, j)/p) for z in zs], dtype=complex)
    return U2_mass(f), len(zs)


def main():
    rng = np.random.default_rng(SEED)
    out = {"phase": 58, "seed": SEED,
           "lemma": "E[U2^4] = (2+(1-q)^2)/L ; ratio (2+(1-q)^2)/2",
           "curves": [], "synthetic_law": [], "first_principles": []}

    # (B) synthetic law across q
    L = 400
    print("== synthetic (2+(1-q)^2)/L law ==")
    for q in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
        m = synth_noisy_reflection_mass(L, q, 10007, rng, nrep=3000)
        pred = 2 + (1 - q) ** 2
        out["synthetic_law"].append(
            {"q": q, "E[M]*L_mc": round(m, 4), "predicted_2+(1-q)^2": round(pred, 4),
             "ratio_mc": round(m / 2, 4), "ratio_pred": round(pred / 2, 4)})
        print(f"  q={q:.1f}  E[M]*L={m:.3f}  pred={pred:.3f}  ratio={m/2:.3f}")

    # (A)+(C) real q per curve, and predicted vs measured real mass
    print("\n== real q_j and mass prediction ==")
    for p, Nexp in ORDINARY:
        qrec, n, L, N = measure_q(p, 7, E_PREC)
        rec = {"p": p, "n": n, "L": L, "q": qrec, "digits": {}}
        for j in DIGITS:
            qj = qrec[f"d{j}"]
            mass_real, Lr = real_u2_mass(p, 7, E_PREC, j)
            if qj is not None:
                q = qj["q_measured"]
                pred_mass_L = 2 + (1 - q) ** 2          # E[M]*L predicted
                # synthetic mass at this exact (L, q) for a calibrated compare
                syn_mass_L = synth_noisy_reflection_mass(
                    L, q, 10007, rng, nrep=1500)
                rec["digits"][f"d{j}"] = {
                    "q_measured": q,
                    "mass_real*L": round(mass_real * L, 4),
                    "pred_mass*L_(2+(1-q)^2)": round(pred_mass_L, 4),
                    "syn_mass*L_at_q": round(syn_mass_L, 4),
                    "first_principles_q": round(
                        carry_prob_first_principles(p, j, E_PREC), 4),
                }
        rec["fp"] = {f"d{j}": round(carry_prob_first_principles(p, j, E_PREC), 4)
                     for j in DIGITS}
        out["curves"].append(rec)
        d = rec["digits"]
        print(f"  p={p:4d} n={n:4d} L={L:4d}  "
              + "  ".join(f"d{j}:q={d[f'd{j}']['q_measured']:.2f} "
                          f"M*L real={d[f'd{j}']['mass_real*L']:.2f}/"
                          f"pred={d[f'd{j}']['pred_mass*L_(2+(1-q)^2)']:.2f}"
                          for j in DIGITS if f"d{j}" in d))

    (ROOT/"results").mkdir(exist_ok=True)
    (ROOT/"results"/"phase58_carry_lemma.json").write_text(json.dumps(out, indent=2))
    print("\nwrote results/phase58_carry_lemma.json")


if __name__ == "__main__":
    main()
