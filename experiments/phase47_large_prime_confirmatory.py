"""
Phase 47 - Large-prime-order confirmatory curve (referee point S5).

The clean Phase-43 set tops out at n = 829. This phase runs the identical
pre-registered protocol (e = 4, digits j in {1,2,3}, 300 draws per null,
seed 20260609, same exclusion rule) on ONE ordinary j=0 curve with
n ~ 10^4 prime: the smallest p >= 10^4 with p = 1 mod 3 and
N = #E(F_p) prime (hence ordinary, non-anomalous, gcd(n,p)=1).

Confirms, at a length scale >10x the previous maximum:
  * Theorem 1(b),(d): obstruction C_n = z([n]tau(G)) nonzero, valuation 1,
    additive order p^(e-1);
  * identity (3): z(delta(k)) + z(delta(n-k)) single-valued = C_n;
  * the U^2 control pattern of Table 1: real/inc elevated and growing
    with n, inc/rand and half/inc and real/syn at noise, syn/rand
    tracking real/inc.

This is a CONFIRMATORY run, not part of the pre-registered Phase-43 set;
it is reported separately in the paper for exactly that reason.

Output: results/phase47_large_prime.json
"""
from __future__ import annotations
import json
import math
import pathlib
import random
import time

import numpy as np

T0 = time.time()
SEED = 20260609
random.seed(SEED)
RNG = np.random.default_rng(SEED)
NREP = 300
E_PREC = 4          # paper precision e
WORK_PREC = 10      # internal working precision


# ---------------------------------------------------------------- utilities
def legendre(a: int, p: int) -> int:
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def is_prime(m: int) -> bool:
    if m < 2:
        return False
    for q in range(2, int(m ** 0.5) + 1):
        if m % q == 0:
            return False
    return True


def count_points(p: int, a: int, b: int) -> int:
    return 1 + sum(1 + legendre(x * x * x + a * x + b, p) for x in range(p))


def fp_add(P, Q, p, a):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)


class Lift:
    """x-Teichmuller / Hensel section and exact E(Z/p^prec) arithmetic."""

    def __init__(self, p, a, b, prec):
        self.p, self.a, self.b, self.prec = p, a, b, prec
        self.mod = p ** prec
        self._inv2 = pow(2, -1, self.mod)

    def teich_x(self, x0):
        t = x0 % self.mod
        for _ in range(self.prec + 3):
            t = pow(t, self.p, self.mod)
        return t

    def hensel_y(self, x, y0):
        c = (x * x * x + self.a * x + self.b) % self.mod
        y = y0 % self.p
        for _ in range(self.prec + 3):
            y = (y + c * pow(y, -1, self.mod)) * self._inv2 % self.mod
        assert y * y % self.mod == c and y % self.p == y0 % self.p
        return y

    def tau(self, P):
        x = self.teich_x(P[0])
        return (x, self.hensel_y(x, P[1]))

    def add(self, P, Q):
        x1, y1 = P
        x2, y2 = Q
        lam = (y2 - y1) * pow(x2 - x1, -1, self.mod) % self.mod
        x3 = (lam * lam - x1 - x2) % self.mod
        return (x3, (lam * (x1 - x3) - y1) % self.mod)

    def dbl(self, P):
        x1, y1 = P
        lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, self.mod) % self.mod
        x3 = (lam * lam - 2 * x1) % self.mod
        return (x3, (lam * (x1 - x3) - y1) % self.mod)

    def zdiff(self, A, B):
        """z(A - B) for A = B mod p: projective A + (-B), z = -X3/Y3."""
        X1, Y1 = A
        X2, Y2 = B[0], (-B[1]) % self.mod
        u = (Y2 - Y1) % self.mod
        v = (X2 - X1) % self.mod
        w = (u * u - v ** 3 - 2 * v * v * X1) % self.mod
        X3 = v * w % self.mod
        Y3 = (u * (v * v * X1 - w) - v ** 3 * Y1) % self.mod
        return (-X3) * pow(Y3, -1, self.mod) % self.mod


def val(x, p, prec):
    if x == 0:
        return prec
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


# ------------------------------------------------------------ curve search
def find_curve(start=10_000):
    p = start
    while True:
        p += 1
        if p % 3 == 1 and is_prime(p):
            N = count_points(p, 0, 7)
            if is_prime(N) and N not in (p, p + 1):
                return p, N


# ------------------------------------------------------------ U2 machinery
def U2(fv):
    F = np.fft.fft(fv) / len(fv)
    return float((np.abs(F) ** 4).sum() ** 0.25)


def phases(s_int64, j, p):
    return np.exp(2j * np.pi * ((s_int64 // p ** j) % p) / p)


def main():
    p, n = find_curve()
    print(f"[{time.time()-T0:6.1f}s] curve: p={p}, n=N={n} (prime), "
          f"p%3={p % 3}, gcd(n,p)={math.gcd(n, p)}")
    a, b = 0, 7
    pe = p ** E_PREC
    L = Lift(p, a, b, WORK_PREC)

    # generator
    G = None
    for x in range(p):
        c = (x ** 3 + b) % p
        if legendre(c, p) == 1:
            for yy in range(1, p):
                if yy * yy % p == c:
                    G = (x, yy)
                    break
            break

    # [m] tau(G) chain and delta(k)
    T = L.tau(G)
    mult_prev = T
    zd = [0] * (n - 1)          # zd[k-1] = z(delta(k)) mod p^e
    Pk = G
    zd[0] = 0                   # delta(1) = O
    for k in range(2, n):
        mult_prev = L.dbl(T) if k == 2 else L.add(mult_prev, T)
        Bk = L.tau(Pk := fp_add(Pk, G, p, a))
        zd[k - 1] = L.zdiff(mult_prev, Bk) % pe
    # obstruction C_n = z([n]tau(G)) = z(mult[n-1] - (-T))
    Cn = L.zdiff(mult_prev, (T[0], (-T[1]) % L.mod)) % pe
    print(f"[{time.time()-T0:6.1f}s] lift errors done")

    # exclusion rule (pre-registered): drop degenerate values
    drops = sum(1 for z in zd if z % p != 0)
    assert drops == 0, f"{drops} degenerate values"

    # Theorem 1 checks
    sums = {(zd[k - 1] + zd[n - k - 1]) % pe for k in range(1, n)}
    checks = {
        "antisymmetry_single_value": len(sums) == 1,
        "antisymmetry_equals_Cn": sums == {Cn},
        "Cn_valuation": val(Cn, p, E_PREC),
        "Cn_additive_order_is_p^(e-1)":
            (Cn * p ** (E_PREC - 1)) % pe == 0
            and (Cn * p ** (E_PREC - 2)) % pe != 0,
        "L": n - 1, "drops": drops,
    }
    print(f"[{time.time()-T0:6.1f}s] theorem checks: {checks}")

    # ---- U2 controls, identical protocol to Phase 43 ----
    s_np = np.array(zd, dtype=np.int64)
    Lseq = n - 1
    h = (n - 1) // 2

    rand_u2 = [U2(np.exp(2j * np.pi * RNG.random(Lseq))) for _ in range(NREP)]
    mu_r, sd_r = float(np.mean(rand_u2)), float(np.std(rand_u2))

    d = [(zd[i + 1] - zd[i]) % pe for i in range(Lseq - 1)]
    dh = [(zd[i + 1] - zd[i]) % pe for i in range(h - 1)]

    def resum(base, incs):
        out = [base]
        acc = base
        for x in incs:
            acc = (acc + x) % pe
            out.append(acc)
        return np.array(out, dtype=np.int64)

    inc_seqs, inc_half, syn_seqs = [], [], []
    for _ in range(NREP):
        dd = d[:]
        random.shuffle(dd)
        inc_seqs.append(resum(zd[0], dd))
        dd = dh[:]
        random.shuffle(dd)
        inc_half.append(resum(zd[0], dd))
        first = (RNG.integers(0, pe // p, h) * p).astype(np.int64)
        full = np.zeros(Lseq, dtype=np.int64)
        full[0:h] = first                       # k = 1..h
        full[n - 1 - h: n - 1] = (Cn - first[::-1]) % pe  # k = n-h..n-1
        syn_seqs.append(full)
    print(f"[{time.time()-T0:6.1f}s] null ensembles built")

    half_np = s_np[:h]
    table = {}
    for j in (1, 2, 3):
        u_real = U2(phases(s_np, j, p))
        u_inc = [U2(phases(ss, j, p)) for ss in inc_seqs]
        u_syn = [U2(phases(ss, j, p)) for ss in syn_seqs]
        u_hr = U2(phases(half_np, j, p))
        u_hi = [U2(phases(ss, j, p)) for ss in inc_half]
        table[j] = {
            "real_vs_inc": (u_real - np.mean(u_inc)) / np.std(u_inc),
            "inc_vs_rand": (np.mean(u_inc) - mu_r) / sd_r,
            "half_vs_inc": (u_hr - np.mean(u_hi)) / np.std(u_hi),
            "syn_vs_rand": (np.mean(u_syn) - mu_r) / sd_r,
            "real_vs_syn": (u_real - np.mean(u_syn)) / np.std(u_syn),
        }
        print(f"[{time.time()-T0:6.1f}s] digit {j}: "
              + "  ".join(f"{k}={v:+.2f}" for k, v in table[j].items()))

    out = {
        "phase": 47, "seed": SEED, "e": E_PREC, "nrep": NREP,
        "curve": {"p": p, "a": a, "b": b, "n": n, "G": list(G)},
        "Cn_mod_pe": Cn, "theorem_checks": checks,
        "u2_z_scores_by_digit": {str(j): {k: round(float(v), 2)
                                          for k, v in table[j].items()}
                                 for j in table},
        "elapsed_s": round(time.time() - T0, 1),
    }
    res = pathlib.Path(__file__).resolve().parent.parent / "results"
    res.mkdir(exist_ok=True)
    (res / "phase47_large_prime.json").write_text(json.dumps(out, indent=2))
    print(f"[{time.time()-T0:6.1f}s] wrote results/phase47_large_prime.json")


if __name__ == "__main__":
    main()
