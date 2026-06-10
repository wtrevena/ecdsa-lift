"""
Phase 44 -- Multi-instance / cross-curve statistical aggregation.

This is the ONE attack lane the inertness theorem (notes/theorem.md) does
NOT directly cover. The theorem proves that for gcd(n,p)=1 the *per-
instance* lift error delta(k) = [k]tau(G) - tau(kG) is inert: the only
public structure is

  (H1) an Aut(E)-equivariance of the leading base-p digit of z(delta),
  (H2) the antisymmetry pairing  z(delta(k)) + z(delta(n-k)) = Cn = z([n]tau(G)).

Phase 40/40b already showed the *entire* Phase-37 Fourier "anomaly" was the
Fourier shadow of (H2). So the per-instance lane is well understood and
null beyond (H1),(H2).

HYPOTHESIS (this phase): even though each delta is per-instance inert, the
JOINT distribution of delta-derived features across a LARGE ENSEMBLE of
curves might carry a secret-correlated signal exploitable only in aggregate
-- the differential/linear-cryptanalysis analogy, where single traces look
random but pooled statistics leak. Concretely: does the mutual information
I(secret ; features), POOLED across many curves, exceed its label-shuffled
null floor by more than per-curve MI does, AFTER we quotient out the public
identities (H1),(H2)?

FALSIFIER: if pooled MI (beyond the public identities) sits at the shuffle
null floor (|z| <~ 3 with multiple-comparison context), the lane is dead.

DESIGN (pre-specified, no post-hoc tuning):
  Precision           e = 4   (z is additive at this precision; see phase41).
  Seed                20260609.    Per-curve null SHUF=300; pooled-kNN SHUF_SK=120.
  Ensembles:
    (A) j=0 CM family: y^2 = x^3 + 7 over F_p, p = 1 mod 3, p in 31..200,
        excluding anomalous (N=p) and supersingular (N=p+1).
    (B) generic non-CM: random (a,b), j != 0,1728 (trivial extra Aut),
        p in 31..200, drawn deterministically from the seed.
        Ensemble (B) is the clean cross-curve test because the Aut-
        equivariance identity (H1) is ABSENT there (Aut = {+-1} only).
  Features of z(delta(k)) (kernel element => z = 0 mod p):
    f_digit1,f_digit2,f_digit3 : base-p digits 1,2,3 (digit 0 is identically
        0; digit 1 is the leading digit carrying (H1)).
    f_modp2 : z(delta(k)) mod p^2.   f_val : p-adic valuation capped at e.
  Labels: k (full secret); k mod 2,3,5 (coarse residues); knorm = floor(8k/n)
    (curve-independent normalized secret bucket, used for pooling).
  PUBLIC-IDENTITY QUOTIENT (the critical control):
    * (H2): restrict to FIRST HALF k <= floor((n-1)/2). Within one half the
      reflection k<->n-k maps out of the window, imposing no internal
      constraint (the phase40b control). Full range also reported.
    * (H1): two feature sets -- FULL (incl f_digit1) and BEYOND_H1 (drop
      f_digit1). Headline uses FIRST-HALF + BEYOND_H1, corroborated on the
      non-CM ensemble where (H1) is vacuous.
  MI ESTIMATION (two independent estimators):
    * histogram / plug-in MI (exact discrete MI, nats) -- the workhorse,
      used per-curve and pooled, null = 300 label permutations.
    * sklearn mutual_info_classif CONTINUOUS kNN (Kraskov; tiny jitter) --
      a genuinely different estimator, run on the POOLED data as a cross-
      check (null = 120 permutations); skipped when a class has
      < n_neighbors+1 members. (For discrete_features=True sklearn returns
      exactly the plug-in MI, so the continuous path is the informative one.)
  z-score = (MI_real - mean(MI_null)) / std(MI_null), plus one-sided p_emp.

Outputs results/phase44_aggregation.json and a notes/phase44 memo.
"""
from __future__ import annotations
import sys
import pathlib
import json
import random
import warnings
from math import gcd

import numpy as np

warnings.filterwarnings("ignore")

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.curves import Curve, find_generator, naive_order, modinv
from src.projective import ProjCurve
from src.lifts import teichmuller_lift

try:
    from sklearn.feature_selection import mutual_info_classif
    HAVE_SKLEARN = True
except Exception:  # pragma: no cover
    HAVE_SKLEARN = False

SEED = 20260609
E_PREC = 4
SHUF = 300        # per-feature histogram null (fast, vectorized)
SHUF_SK = 120     # pooled kNN cross-check null (slower)
P_LO, P_HI = 31, 200
N_NONCM = 8


# ----------------------------------------------------------------------
def is_prime(m: int) -> bool:
    if m < 2:
        return False
    i = 2
    while i * i <= m:
        if m % i == 0:
            return False
        i += 1
    return True


def p_adic_val(x: int, p: int, e: int) -> int:
    if x % (p ** e) == 0:
        return e
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


# ----------------------------------------------------------------------
def z_delta_sequence(p: int, a: int, b: int, e: int = E_PREC):
    """(records, n, N, Cn). Pre-specified degeneracy filter (phase41/43):
    drop k whose projective lift is degenerate (Y not a unit) or whose z is
    not a kernel element (z != 0 mod p)."""
    N = p ** e
    E_p = Curve(a=a, b=b, N=p)
    G_fp, n = find_generator(E_p)
    if gcd(n, p) != 1:
        return None
    C = ProjCurve(a=a, b=b, N=N)
    E_aff = Curve(a=a, b=b, N=N)
    tauG = C.from_affine(teichmuller_lift(E_p, E_aff, G_fp))
    nT = C.mul(n, tauG)
    try:
        Cn = (-nT[0] * modinv(nT[1] % N, N)) % N
    except ZeroDivisionError:
        Cn = None
    records = []
    # k=1 is excluded by the public-identity quotient: delta(1) = [1]tau(G) -
    # tau(G) = 0 identically (valuation = e) on EVERY curve, carrying zero
    # secret information. Pooling this fixed point across the ensemble is a
    # pure public identity (the most trivial one), so we start at k=2.
    for k in range(2, n):
        kG = E_p.mul(k, G_fp)
        tkG = C.from_affine(teichmuller_lift(E_p, E_aff, kG))
        d = C.add(C.mul(k, tauG), C.neg(tkG))
        X, Y, Z = d
        if gcd(Y % N, N) != 1:
            continue
        z = (-X * modinv(Y % N, N)) % N
        if z % p != 0:
            continue
        dig = []
        zz = z
        for _ in range(e):
            dig.append(zz % p)
            zz //= p
        records.append({
            "k": k,
            "f_digit1": dig[1] if len(dig) > 1 else 0,
            "f_digit2": dig[2] if len(dig) > 2 else 0,
            "f_digit3": dig[3] if len(dig) > 3 else 0,
            "f_modp2": z % (p * p),
            "f_val": p_adic_val(z, p, e),
        })
    return records, n, N, Cn


# ----------------------------------------------------------------------
# MI estimators
# ----------------------------------------------------------------------
FEATURE_KEYS_FULL = ["f_digit1", "f_digit2", "f_digit3", "f_modp2", "f_val"]
FEATURE_KEYS_BEYOND_H1 = ["f_digit2", "f_digit3", "f_modp2", "f_val"]


def _design_matrix(records, feature_keys):
    return np.array([[r[k] for k in feature_keys] for r in records], dtype=float)


def _codes(col):
    """Map an integer feature column to dense 0..K-1 codes; return (codes, K)."""
    vals, inv = np.unique(col, return_inverse=True)
    return inv.astype(np.int64), len(vals)


def histogram_mi_codes(fcodes, Kf, ycodes, Ky):
    """Plug-in MI (nats) between two code arrays via a 2D bincount."""
    n = len(fcodes)
    if n == 0 or Kf < 2 or Ky < 2:
        return 0.0
    flat = fcodes * Ky + ycodes
    joint = np.bincount(flat, minlength=Kf * Ky).reshape(Kf, Ky).astype(float)
    joint /= n
    pf = joint.sum(axis=1, keepdims=True)
    py = joint.sum(axis=0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = joint * np.log(joint / (pf * py))
    terms[~np.isfinite(terms)] = 0.0
    return float(terms.sum())


def hist_mi_total_and_null(X, labels, rng, shuf):
    """Sum over features of plug-in MI(feature; label), with a shuffle null.
    Returns (real, null_array). Feature codes are precomputed once; the null
    permutes labels only -> fast."""
    ycodes, Ky = _codes(np.asarray(labels))
    fc = [_codes(X[:, j].astype(np.int64)) for j in range(X.shape[1])]
    real = float(sum(histogram_mi_codes(c, K, ycodes, Ky) for c, K in fc))
    null = np.empty(shuf)
    yc = ycodes.copy()
    for s in range(shuf):
        ys = yc[rng.permutation(len(yc))]
        null[s] = float(sum(histogram_mi_codes(c, K, ys, Ky) for c, K in fc))
    return real, null


def sklearn_mi(X, labels, rng_seed):
    """Total continuous-kNN mutual_info_classif (nats). Returns None when the
    smallest class has < n_neighbors+1 members (kNN classif cannot build the
    per-class tree there); the histogram estimator carries those cases."""
    if not HAVE_SKLEARN:
        return None
    n_neighbors = 3
    _, counts = np.unique(np.asarray(labels), return_counts=True)
    if counts.min() < n_neighbors + 1:
        return None
    jr = np.random.default_rng(rng_seed)
    Xc = X + jr.normal(0, 1e-6, size=X.shape)
    mi = mutual_info_classif(Xc, labels, discrete_features=False,
                             random_state=rng_seed, n_neighbors=n_neighbors)
    return float(np.sum(mi))


def _zscore(real, null):
    mu = float(np.mean(null))
    sd = float(np.std(null))
    z = (real - mu) / sd if sd > 0 else 0.0
    p_emp = float((np.sum(null >= real) + 1) / (len(null) + 1))
    return dict(real=real, null_mean=mu, null_std=sd, z=z, p_emp=p_emp)


def mi_block(X, labels, rng, with_sklearn=False, sk_seed=0):
    """Histogram MI (+null) always; kNN MI (+null) only if with_sklearn."""
    real_h, null_h = hist_mi_total_and_null(X, labels, rng, SHUF)
    out = {"hist": _zscore(real_h, null_h)}
    if with_sklearn:
        real_sk = sklearn_mi(X, labels, sk_seed)
        if real_sk is not None:
            null_sk = np.empty(SHUF_SK)
            lab = np.asarray(labels)
            for s in range(SHUF_SK):
                v = sklearn_mi(X, lab[rng.permutation(len(lab))], sk_seed + 1 + s)
                null_sk[s] = v if v is not None else 0.0
            out["sklearn"] = _zscore(real_sk, null_sk)
    return out


# ----------------------------------------------------------------------
def make_labels(records, n, kind):
    ks = np.array([r["k"] for r in records])
    if kind == "k":
        return ks
    if kind == "knorm":
        return np.minimum(np.floor((ks / n) * 8).astype(int), 7)
    if kind.startswith("kmod"):
        return ks % int(kind[4:])
    raise ValueError(kind)


def first_half(records, n):
    half = (n - 1) // 2
    return [r for r in records if r["k"] <= half]


# ----------------------------------------------------------------------
def build_j0_ensemble():
    fam = []
    for p in range(P_LO, P_HI + 1):
        if not is_prime(p) or p % 3 != 1:
            continue
        E = Curve(a=0, b=7, N=p)
        N = naive_order(E)
        if N == p or N == p + 1:
            continue
        G, n = find_generator(E)
        if n < 8 or gcd(n, p) != 1:
            continue
        fam.append({"p": p, "a": 0, "b": 7, "n": n, "N_group": N, "family": "j0"})
    return fam


def build_noncm_ensemble(rng):
    primes = [p for p in range(P_LO, P_HI + 1) if is_prime(p)]
    fam, seen, attempts = [], set(), 0
    while len(fam) < N_NONCM and attempts < 5000:
        attempts += 1
        p = rng.choice(primes)
        a = rng.randrange(1, p)
        b = rng.randrange(1, p)
        if (4 * a * a * a + 27 * b * b) % p == 0:
            continue
        if a % p == 0 or b % p == 0:
            continue
        E = Curve(a=a, b=b, N=p)
        N = naive_order(E)
        if N == p or N == p + 1:
            continue
        G, n = find_generator(E)
        if n < 10 or gcd(n, p) != 1:
            continue
        key = (p, a, b)
        if key in seen:
            continue
        seen.add(key)
        fam.append({"p": p, "a": a, "b": b, "n": n, "N_group": N, "family": "noncm"})
    return fam


# ----------------------------------------------------------------------
def analyze_curve(spec, rng):
    """Per-curve MI (histogram estimator only -- fast). kNN is reserved for
    the pooled headline."""
    res = z_delta_sequence(spec["p"], spec["a"], spec["b"], E_PREC)
    if res is None:
        return None
    records, n, N, Cn = res
    out = {**spec, "kept": len(records), "total": n - 1, "Cn": Cn, "per_curve": {}}
    recs_half = first_half(records, n)
    for rng_label, recs in (("full", records), ("half", recs_half)):
        if len(recs) < 6:
            continue
        block = {}
        for fset_name, fkeys in (("full_feats", FEATURE_KEYS_FULL),
                                 ("beyond_H1", FEATURE_KEYS_BEYOND_H1)):
            X = _design_matrix(recs, fkeys)
            for lab_kind in ("k", "kmod2", "kmod3", "kmod5"):
                labels = make_labels(recs, n, lab_kind)
                if len(set(labels.tolist())) < 2:
                    continue
                block[f"{fset_name}|{lab_kind}"] = mi_block(
                    X, labels, rng, with_sklearn=False)
        out["per_curve"][rng_label] = block
    return out, records, recs_half, n


def pooled_analysis(curve_payloads, rng, family_name):
    """Pool FIRST-HALF rows across the ensemble (H2 quotient). Coarse/
    normalized labels are curve-independent so pooling is legitimate. Both
    estimators (histogram + kNN cross-check) reported here -- this is the
    headline lane."""
    pooled = {}
    for fset_name, fkeys in (("full_feats", FEATURE_KEYS_FULL),
                             ("beyond_H1", FEATURE_KEYS_BEYOND_H1)):
        for lab_kind in ("kmod2", "kmod3", "kmod5", "knorm"):
            X_rows, y_rows = [], []
            for (records, recs_half, n) in curve_payloads:
                if len(recs_half) < 4:
                    continue
                X_rows.append(_design_matrix(recs_half, fkeys))
                y_rows.append(make_labels(recs_half, n, lab_kind))
            if not X_rows:
                continue
            X = np.vstack(X_rows)
            y = np.concatenate(y_rows)
            if len(set(y.tolist())) < 2 or X.shape[0] < 20:
                continue
            pooled[f"{fset_name}|{lab_kind}"] = {
                "n_rows": int(X.shape[0]),
                **mi_block(X, y, rng, with_sklearn=True,
                           sk_seed=hash((family_name, lab_kind)) % 9973),
            }
    return pooled


def summarize_z(block, estimator):
    return {k: round(v[estimator]["z"], 3) for k, v in block.items()
            if k.startswith("beyond_H1|") and estimator in v}


def main():
    out_dir = ROOT / "results"
    out_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    pyrng = random.Random(SEED)

    print("=== Phase 44: multi-instance / cross-curve aggregation ===")
    print(f"sklearn: {HAVE_SKLEARN}; e={E_PREC}; SHUF={SHUF}; SHUF_SK={SHUF_SK}; seed={SEED}")
    j0 = build_j0_ensemble()
    noncm = build_noncm_ensemble(pyrng)
    print(f"j0 curves:    {[c['p'] for c in j0]}  ({len(j0)})")
    print(f"noncm curves: {[(c['p'], c['a'], c['b']) for c in noncm]}  ({len(noncm)})")

    report = {
        "meta": {
            "phase": 44, "seed": SEED, "e": E_PREC, "shuf": SHUF, "shuf_sk": SHUF_SK,
            "p_range": [P_LO, P_HI], "sklearn": HAVE_SKLEARN,
            "estimators": ["histogram_plugin (per-curve + pooled)",
                           "sklearn_continuous_kNN (pooled cross-check)"],
            "feature_keys_full": FEATURE_KEYS_FULL,
            "feature_keys_beyond_H1": FEATURE_KEYS_BEYOND_H1,
            "public_identities_quotiented": {
                "H0_trivial_k1": "k=1 excluded (delta(1)=0 identically, no secret content)",
                "H2_antisymmetry": "restricted to first half k<=floor((n-1)/2)",
                "H1_aut_leading_digit": "dropped f_digit1; corroborated on non-CM ensemble",
            },
        },
        "ensembles": {},
    }

    for fam_name, fam in (("j0", j0), ("noncm", noncm)):
        print(f"\n--- ensemble {fam_name} ({len(fam)} curves) ---", flush=True)
        per_curve_out, payloads = [], []
        for spec in fam:
            r = analyze_curve(spec, rng)
            if r is None:
                print(f"  p={spec['p']} skipped", flush=True)
                continue
            cinfo, records, recs_half, n = r
            per_curve_out.append(cinfo)
            payloads.append((records, recs_half, n))
            hb = cinfo["per_curve"].get("half", {})
            print(f"  p={spec['p']:>3} n={n:>3} kept={cinfo['kept']:>3}/{cinfo['total']:>3}"
                  f"  half beyondH1 z(hist)={summarize_z(hb,'hist')}", flush=True)

        pooled = pooled_analysis(payloads, rng, fam_name)
        zs_pool_h = summarize_z(pooled, "hist")
        zs_pool_sk = summarize_z(pooled, "sklearn")
        print(f"  POOLED beyond_H1 z(hist):    {zs_pool_h}", flush=True)
        print(f"  POOLED beyond_H1 z(sklearn): {zs_pool_sk}", flush=True)

        agg = {}
        for lab in ("kmod2", "kmod3", "kmod5", "k"):
            key = f"beyond_H1|{lab}"
            vals = [c["per_curve"]["half"][key]["hist"]["z"]
                    for c in per_curve_out
                    if "half" in c["per_curve"] and key in c["per_curve"]["half"]]
            if vals:
                agg[f"hist|{lab}"] = {
                    "per_curve_mean_z": round(float(np.mean(vals)), 3),
                    "per_curve_max_abs_z": round(float(np.max(np.abs(vals))), 3),
                    "n_curves": len(vals),
                }

        report["ensembles"][fam_name] = {
            "curves": [{"p": c["p"], "a": c["a"], "b": c["b"], "n": c["n"],
                        "kept": c["kept"], "total": c["total"]} for c in per_curve_out],
            "per_curve_aggregate_headline": agg,
            "pooled": pooled,
            "pooled_headline_z_hist": zs_pool_h,
            "pooled_headline_z_sklearn": zs_pool_sk,
            "per_curve_full_detail": per_curve_out,
        }

    # VERDICT (pre-specified). Headline = first-half + beyond_H1, both
    # estimators, both ensembles. With many tests: |z|>4 pooled OR |z|>5
    # per-curve = hit; <3 = dead; 3-4 = weak.
    pooled_max = 0.0
    percurve_max = 0.0
    for fam_name in ("j0", "noncm"):
        ens = report["ensembles"][fam_name]
        for est in ("hist", "sklearn"):
            for v in ens[f"pooled_headline_z_{est}"].values():
                pooled_max = max(pooled_max, abs(v))
        for v in ens["per_curve_aggregate_headline"].values():
            percurve_max = max(percurve_max, abs(v["per_curve_max_abs_z"]))

    if pooled_max > 4 or percurve_max > 5:
        verdict = "SURPRISING -- non-null signal beyond public identities"
    elif pooled_max > 3:
        verdict = "WEAK-SIGNAL -- marginal, needs replication"
    else:
        verdict = ("DEAD -- pooled MI beyond public identities sits at the "
                   "shuffle null floor")

    report["verdict"] = {
        "pooled_max_abs_z_beyond_H1": round(pooled_max, 3),
        "per_curve_max_abs_z_beyond_H1": round(percurve_max, 3),
        "thresholds": {"dead": "<3", "weak": "3-4 pooled",
                       "hit": ">4 pooled or >5 per-curve"},
        "verdict": verdict,
    }
    print("\n=== VERDICT ===")
    print(json.dumps(report["verdict"], indent=2))
    (out_dir / "phase44_aggregation.json").write_text(
        json.dumps(report, indent=2, default=str))
    print("\nWrote results/phase44_aggregation.json")


if __name__ == "__main__":
    main()
