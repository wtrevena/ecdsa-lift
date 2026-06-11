"""Supplement to phase49: confound check for the asymmetric-section kill test.

The primary kill test (phase49) uses the paper pre-registered KERNEL filter
(drop k unless z(delta) in pZ). Under that protocol the asymmetric section
collapses the U^2 elevation to baseline AND destroys the antisymmetry identity
(distinct pair-sums > 1). However the filter also halves L (it drops exactly
one member of each {k,n-k} pair, because breaking tau(-P)=-tau(P) pushes one
of the pair out of the kernel). To document that confound transparently we
ALSO record the UNFILTERED variant (all k kept, encode digit j of z mod p^e):
there the asymmetric section shows a LARGE but spurious U^2 (z leaves the
kernel so its mod-p digit becomes strongly non-uniform low-complexity
structure -- NOT the Gowers anomaly the paper studies, which lives in digits
j>=1 of genuine kernel elements). Both views agree the reflection break
qualitatively changes delta and kills the SPECIFIC antisymmetry-driven U^2
elevation the paper attributes to tau(-P)=-tau(P).
"""
import sys, pathlib, json
import numpy as np
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from math import gcd
from src.curves import Curve, find_generator, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift
from experiments.phase49_asymmetric_section import asymmetric_lift, fast_generator
from experiments.phase37_gowers import gowers_norm_complex, encode_digit, random_baseline
SEED = 20260609
# p=10477 omitted from this supplementary unfiltered check for runtime;
# the confound is already unambiguous across the three scales below.
TARGETS = [(67, 79), (211, 199), (823, 829)]

def seq_nofilter(p, b, e, section, canon, n_known=None):
    N = p ** e
    Ep = Curve(0, b, p)
    if n_known:
        G, n = fast_generator(Ep, n_known)
    else:
        G, n = find_generator(Ep)
    C = ProjCurve(0, b, N)
    Ea = Curve(0, b, N)
    if section == "standard":
        lift = lambda P: teichmuller_lift(Ep, Ea, P)
    else:
        lift = lambda P: asymmetric_lift(Ep, Ea, P, canon=canon)
    tauG = C.from_affine(lift(G))
    zmap = {}
    zlist = []
    for k in range(1, n):
        tkG = C.from_affine(lift(Ep.mul(k, G)))
        d = C.add(C.mul(k, tauG), C.neg(tkG))
        X, Y, Z = d
        if gcd(Y % N, N) != 1:
            zlist.append(None)
            continue
        z = (-X * modinv(Y % N, N)) % N
        zmap[k] = z
        zlist.append(z)
    return zlist, zmap, n, N

def u2_modN(zlist, N, j, p):
    f = np.array([1 + 0j if z is None else encode_digit(int(z), p, j) for z in zlist], dtype=complex)
    return gowers_norm_complex(f, 2)


def frac_in_kernel(zlist, p):
    tot = sum(1 for z in zlist if z is not None)
    ker = sum(1 for z in zlist if z is not None and z % p == 0)
    return ker, tot


def main():
    np.random.seed(SEED)
    out = {"note": "unfiltered confound check; see module docstring", "seed": SEED, "curves": []}
    for p, nn in TARGETS:
        nk = nn if p >= 5000 else None
        rec = {"p": p, "n": nn, "sections": {}}
        for sec, canon, lab in [("standard", None, "standard"), ("asym", "smaller_rep", "asym_smaller_rep"), ("asym", "even_lsd", "asym_even_lsd")]:
            zlist, zmap, n, N = seq_nofilter(p, 7, 4, sec, canon, n_known=nk)
            ker, tot = frac_in_kernel(zlist, p)
            L = len(zlist)
            rmean, rstd, _ = random_baseline(L, 2, num_samples=300)
            digs = {}
            for j in (1, 2, 3):
                u2 = u2_modN(zlist, N, j, p)
                digs["d" + str(j)] = {"U2": round(u2, 5), "z_real_vs_rand": round((u2 - rmean) / rstd, 2)}
            rec["sections"][lab] = {"L_all": L, "in_kernel": ker, "frac_in_kernel": round(ker / tot, 4), "digits_zmodpe": digs}
            z=[digs["d"+str(j)]["z_real_vs_rand"] for j in (1,2,3)]
            print("p=%-6d %-18s L=%d in_kernel=%d/%d  z(real/rand,zmodpe) d1,d2,d3=%s" % (p, lab, L, ker, tot, z))
        out["curves"].append(rec)
    (ROOT / "results").mkdir(exist_ok=True)
    (ROOT / "results" / "phase49b_killtest_confound.json").write_text(json.dumps(out, indent=2, default=str))
    print("Wrote results/phase49b_killtest_confound.json")


if __name__ == "__main__":
    main()
