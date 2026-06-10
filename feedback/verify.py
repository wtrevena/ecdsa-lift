"""Independent verification of key claims in Trevena, 'When a Control Does Not Control'.
Pure Python + numpy. Reproduces: curve orders (Table 1), Teichmuller lift error delta(k),
antisymmetry identity (3), obstruction [n]tau(G) (Thm 1d/b), Gowers U2 controls (Table 1 row p=67),
and Smart's attack boundary (Table 2).
"""
import random, math
import numpy as np

random.seed(12345)

def legendre(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def count_points(p, a, b):
    N = 1
    for x in range(p):
        N += 1 + legendre(x * x * x + a * x + b, p)
    return N

def is_prime(m):
    if m < 2: return False
    for q in range(2, int(m ** 0.5) + 1):
        if m % q == 0: return False
    return True

# ---------- Fp affine arithmetic ----------
def fp_add(P, Q, p, a):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0: return None
    if P == Q:
        lam = (3 * x1 * x1 + a) * pow(2 * y1, -1, p) % p
    else:
        lam = (y2 - y1) * pow(x2 - x1, -1, p) % p
    x3 = (lam * lam - x1 - x2) % p
    return (x3, (lam * (x1 - x3) - y1) % p)

def fp_mul(k, P, p, a):
    R = None; Q = P
    while k:
        if k & 1: R = fp_add(R, Q, p, a)
        Q = fp_add(Q, Q, p, a); k >>= 1
    return R

# ---------- Z/p^prec lift arithmetic ----------
class Lift:
    def __init__(self, p, a, b, prec):
        self.p, self.a, self.b, self.prec = p, a, b, prec
        self.mod = p ** prec

    def teich_x(self, x0):
        t = x0 % self.mod
        for _ in range(self.prec + 3):
            t = pow(t, self.p, self.mod)
        return t

    def hensel_y(self, x, y0):
        c = (x * x * x + self.a * x + self.b) % self.mod
        y = y0 % self.p
        inv2 = pow(2, -1, self.mod)
        for _ in range(self.prec + 3):
            y = (y + c * pow(y, -1, self.mod)) * inv2 % self.mod
        assert y * y % self.mod == c, "hensel fail"
        assert y % self.p == y0 % self.p
        return y

    def tau(self, P):
        x = self.teich_x(P[0])
        return (x, self.hensel_y(x, P[1]))

    def add(self, P, Q):
        x1, y1 = P; x2, y2 = Q
        lam = (y2 - y1) * pow(x2 - x1, -1, self.mod) % self.mod
        x3 = (lam * lam - x1 - x2) % self.mod
        return (x3, (lam * (x1 - x3) - y1) % self.mod)

    def dbl(self, P):
        x1, y1 = P
        lam = (3 * x1 * x1 + self.a) * pow(2 * y1, -1, self.mod) % self.mod
        x3 = (lam * lam - 2 * x1) % self.mod
        return (x3, (lam * (x1 - x3) - y1) % self.mod)

    def zdiff(self, A, B):
        """z(A - B) for A,B with the same reduction mod p (difference in kernel).
        Projective addition A + (-B), no division until final unit inverse; z = -X3/Y3."""
        X1, Y1 = A
        X2, Y2 = B[0], (-B[1]) % self.mod
        u = (Y2 - Y1) % self.mod
        v = (X2 - X1) % self.mod
        w = (u * u - v ** 3 - 2 * v * v * X1) % self.mod
        X3 = v * w % self.mod
        Y3 = (u * (v * v * X1 - w) - v ** 3 * Y1) % self.mod
        return (-X3) * pow(Y3, -1, self.mod) % self.mod

def val(x, p, prec):
    if x == 0: return prec
    v = 0
    while x % p == 0: x //= p; v += 1
    return v

# ============ PART 1: Table 1 curve orders ============
print("=== Part 1: Table 1 curve check (y^2=x^3+7) ===")
for p in [67, 163, 211, 349, 433, 577, 823]:
    N = count_points(p, 0, 7)
    print(f"p={p:4d}  p%3={p%3}  N={N:4d}  N prime={is_prime(N)}  gcd(N,p)=1={math.gcd(N,p)==1}")
print("paper claims n = 79,139,199,313,397,613,829")

print("\n=== Part 1b: sec.3 example p=97, and anomalous curve order ===")
print("p=97: N =", count_points(97, 0, 7))
print("F53, y^2=x^3+4x+7: N =", count_points(53, 4, 7), "(anomalous iff N=53)")

print("\n=== Part 1c: Phase 41 orders, y^2=x^3+7, p=31..163 (p=1 mod 3) ===")
for p in [x for x in range(31, 164) if is_prime(x) and x % 3 == 1]:
    print(f"p={p:4d} N={count_points(p,0,7)}", end='   ')
print()

# ============ PART 2: lift error, antisymmetry, obstruction (p=67) ============
print("\n=== Part 2: delta(k), antisymmetry (3), obstruction (Thm 1) on p=67 ===")
p, a, b, e = 67, 0, 7, 4
prec = 10
pe = p ** e
n = count_points(p, a, b)  # 79, prime
# find a generator point
G = None
for x in range(p):
    if legendre(x**3 + a*x + b, p) == 1:
        y = pow(x**3 + a*x + b, (p + 1) // 4, p) if p % 4 == 3 else None
        if y is None:
            for yy in range(1, p):
                if yy * yy % p == (x**3 + a*x + b) % p: y = yy; break
        G = (x, y); break
assert fp_mul(n, G, p, a) is None, "G order check"
L = Lift(p, a, b, prec)
T = L.tau(G)
# mult[m] = [m] tau(G), m=1..n-1, by repeated affine addition
mult = {1: T, 2: L.dbl(T)}
for m in range(3, n):
    mult[m] = L.add(mult[m - 1], T)
# obstruction C_n = z([n] tau(G)) = z(mult[n-1] - (-T))
negT = (T[0], (-T[1]) % L.mod)
Cn = L.zdiff(mult[n - 1], negT)
print(f"C_n = z([n]tau(G)) mod p^4 = {Cn % pe},  valuation = {val(Cn % pe, p, e)} (claim: 1)")
print(f"additive order p^(e-1)? p^3*Cn=0: {Cn*p**3 % pe == 0}, p^2*Cn!=0: {Cn*p**2 % pe != 0}")
# delta(k) for all k
zd = {}
for k in range(1, n):
    Pk = fp_mul(k, G, p, a)
    B = L.tau(Pk)
    zd[k] = L.zdiff(mult[k], B) % pe
ok_kernel = all(zd[k] % p == 0 for k in range(1, n))
sums = set((zd[k] + zd[n - k]) % pe for k in range(1, n))
print(f"all z(delta(k)) = 0 mod p (kernel): {ok_kernel}")
print(f"antisymmetry: #distinct values of z(d(k))+z(d(n-k)) = {len(sums)} (claim: 1); equals C_n: {sums == { Cn % pe }}")
print(f"delta(1) = O: {zd[1] == 0}")

# ============ PART 3: Gowers U2 controls (Table 1, row p=67) ============
print("\n=== Part 3: Gowers U2 controls, p=67 row of Table 1 ===")
Lseq = n - 1
s_real = [zd[k] for k in range(1, n)]  # s(k), k=1..n-1

def U2(fv):
    F = np.fft.fft(fv) / len(fv)
    return float((np.abs(F) ** 4).sum() ** 0.25)

def phases(s, j):
    return np.exp(2j * np.pi * np.array([(x // p ** j) % p for x in s]) / p)

NREP = 300
rng = np.random.default_rng(20260609)
rand_u2 = [U2(np.exp(2j * np.pi * rng.random(Lseq))) for _ in range(NREP)]
mu_r, sd_r = np.mean(rand_u2), np.std(rand_u2)

def incshuf_seqs(s):
    d = [(s[i + 1] - s[i]) % pe for i in range(len(s) - 1)]
    out = []
    for _ in range(NREP):
        dd = d[:]; random.shuffle(dd)
        ss = [s[0]]
        for x in dd: ss.append((ss[-1] + x) % pe)
        out.append(ss)
    return out

inc_seqs = incshuf_seqs(s_real)
half = s_real[:(n - 1) // 2]
inc_half_seqs = incshuf_seqs(half)

def syn_seq():
    h = (n - 1) // 2
    first = [random.randrange(0, pe, p) for _ in range(h)]
    full = [0] * Lseq
    for i, k in enumerate(range(1, h + 1)):
        full[k - 1] = first[i]
        full[(n - k) - 1] = (Cn - first[i]) % pe
    return full

syn_seqs = [syn_seq() for _ in range(NREP)]

print("digit |  real/inc  inc/rand  half/inc  syn/rand  real/syn   (paper p=67: +3.5..6.7, -0.3..0.2, -0.6..1.5, +3.8..4.1, -0.4..1.2)")
for j in (1, 2, 3):
    u_real = U2(phases(s_real, j))
    u_inc = [U2(phases(ss, j)) for ss in inc_seqs]
    u_syn = [U2(phases(ss, j)) for ss in syn_seqs]
    u_half_real = U2(phases(half, j))
    u_half_inc = [U2(phases(ss, j)) for ss in inc_half_seqs]
    z_real_inc = (u_real - np.mean(u_inc)) / np.std(u_inc)
    z_inc_rand = (np.mean(u_inc) - mu_r) / sd_r
    z_half = (u_half_real - np.mean(u_half_inc)) / np.std(u_half_inc)
    z_syn_rand = (np.mean(u_syn) - mu_r) / sd_r
    z_real_syn = (u_real - np.mean(u_syn)) / np.std(u_syn)
    print(f"  {j}   |  {z_real_inc:+7.2f}  {z_inc_rand:+7.2f}  {z_half:+7.2f}  {z_syn_rand:+7.2f}  {z_real_syn:+7.2f}")

# ============ PART 4: Smart's attack boundary (Table 2) ============
print("\n=== Part 4: anomalous boundary, y^2=x^3+4x+7 / F_53 ===")
p2, a2, b2 = 53, 4, 7
N2 = count_points(p2, a2, b2)
print(f"N = {N2} = p: {N2 == p2}")
L2 = Lift(p2, a2, b2, 8)

def smart_phi(P):
    Tp = L2.tau(P)
    m = {1: Tp, 2: L2.dbl(Tp)}
    for mm in range(3, p2):
        m[mm] = L2.add(m[mm - 1], Tp)
    negTp = (Tp[0], (-Tp[1]) % L2.mod)
    zP = L2.zdiff(m[p2 - 1], negTp)  # z([p] tau(P))
    assert zP % p2 == 0
    return (zP // p2) % p2

# pick G2
G2 = None
for x in range(p2):
    if legendre(x**3 + a2*x + b2, p2) == 1:
        for yy in range(1, p2):
            if yy*yy % p2 == (x**3 + a2*x + b2) % p2: G2 = (x, yy); break
        break
phiG = smart_phi(G2)
succ = 0
for trial in range(12):
    k_secret = random.randrange(2, p2)
    Q = fp_mul(k_secret, G2, p2, a2)
    k_rec = smart_phi(Q) * pow(phiG, -1, p2) % p2
    succ += (k_rec == k_secret)
print(f"Smart recovery: {succ}/12 (paper: 12/12); phi(G) nonzero: {phiG != 0}")

# obstruction vanishes on anomalous curve:
Tp = L2.tau(G2)
m = {1: Tp, 2: L2.dbl(Tp)}
for mm in range(3, p2):
    m[mm] = L2.add(m[mm - 1], Tp)
zob = L2.zdiff(m[p2 - 1], (Tp[0], (-Tp[1]) % L2.mod))
print(f"[n]tau(G) on anomalous curve: z = {zob % p2**2} mod p^2 -> vanishes mod p^2 at e=2: {zob % p2**2 == 0}")
print("\nDONE")
