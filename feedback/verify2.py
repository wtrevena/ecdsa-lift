"""Part 2 of verification: (A) pin down the Sec.8/Table 2 'obstruction vanishes' question;
(B) reproduce full Table 1 (all 7 curves) and the beta~0.4 Fourier decay claim."""
import random, math, time
import numpy as np
from verify import Lift, fp_mul, fp_add, legendre, count_points, val

random.seed(777)
t0 = time.time()

# ---------- (A) anomalous curve: what exactly vanishes? ----------
print("=== (A) Anomalous y^2=x^3+4x+7/F53: section-(in)dependence of [p]tau(G) ===")
p2, a2, b2 = 53, 4, 7
L2 = Lift(p2, a2, b2, 8)
G2 = None
for x in range(p2):
    if legendre(x**3 + a2*x + b2, p2) == 1:
        for yy in range(1, p2):
            if yy*yy % p2 == (x**3 + a2*x + b2) % p2: G2 = (x, yy); break
        break

def p_mult_z(L, P0, lift_pt, n_):
    """z([n_] lift_pt) via affine chain + projective last step."""
    T = lift_pt
    m = {1: T, 2: L.dbl(T)}
    for mm in range(3, n_):
        m[mm] = L.add(m[mm - 1], T)
    return L.zdiff(m[n_ - 1], (T[0], (-T[1]) % L.mod))

# section 1: Teichmuller-x
tau1 = L2.tau(G2)
z1 = p_mult_z(L2, G2, tau1, p2)
# section 2: naive-x lift (x = x0 as integer), Hensel y
x_naive = G2[0] % L2.mod
tau2 = (x_naive, L2.hensel_y(x_naive, G2[1]))
z2 = p_mult_z(L2, G2, tau2, p2)
# section 3: x0 + p*17 lift
x3 = (G2[0] + p2 * 17) % L2.mod
tau3 = (x3, L2.hensel_y(x3, G2[1]))
z3 = p_mult_z(L2, G2, tau3, p2)
print(f"z([p]tau(G)) mod p^2, three different sections: {z1%p2**2}, {z2%p2**2}, {z3%p2**2}")
print(f"  equal mod p^2 (section-independent): {z1%p2**2 == z2%p2**2 == z3%p2**2}")
print(f"  nonzero mod p^2 (i.e. [p]tau(G) != O at e=2): {z1%p2**2 != 0}; val = {val(z1%p2**2, p2, 2)}")
print(f"  mod p^4: {z1%p2**4}, val = {val(z1%p2**4, p2, 4)}")
# and delta is NOT well-defined: delta(k+p)-delta(k) = [p]tau(G) != 0
# also: [p]*(kernel element) = O at e=2:
# kernel elt: tau1 - tau2 reduced... compute z(tau1 - tau2), then [p] of it: log mult by p -> val+1 >= 2
zk = L2.zdiff(tau1, tau2)
print(f"kernel elt c=tau1(G)-tau2(G): z={zk % p2**2} (val {val(zk%p2**2,p2,2)}); p*z(c) mod p^2 = {p2*zk % p2**2} (=> [p]c=O at e=2)")

# ---------- (B) full Table 1 reproduction + beta fit ----------
print("\n=== (B) Full Table 1 (real/inc, inc/rand, half/inc, syn/rand, real/syn) ===")
NREP = 300

def U2(fv):
    F = np.fft.fft(fv) / len(fv)
    return float((np.abs(F) ** 4).sum() ** 0.25)

def curve_run(p, nexp):
    a, b, e, prec = 0, 7, 4, 10
    pe = p ** e
    n = count_points(p, a, b)
    assert n == nexp
    G = None
    for x in range(p):
        c = (x**3 + b) % p
        if legendre(c, p) == 1:
            for yy in range(1, p):
                if yy*yy % p == c: G = (x, yy); break
            break
    L = Lift(p, a, b, prec)
    T = L.tau(G)
    mult = {1: T, 2: L.dbl(T)}
    for m in range(3, n):
        mult[m] = L.add(mult[m - 1], T)
    Cn = L.zdiff(mult[n - 1], (T[0], (-T[1]) % L.mod)) % pe
    # delta via tau(kG): walk kG over Fp incrementally
    zd = []
    Pk = G
    for k in range(1, n):
        B = L.tau(Pk)
        zd.append(L.zdiff(mult[k], B) % pe)
        Pk = fp_add(Pk, G, p, a)
    assert all(z % p == 0 for z in zd)
    assert len(set((zd[k - 1] + zd[n - k - 1]) % pe for k in range(1, n))) == 1
    return n, pe, Cn, zd

def phases(s, j, p, pe):
    return np.exp(2j * np.pi * np.array([(x // p ** j) % p for x in s]) / p)

def stats_for(p, nexp):
    n, pe, Cn, s_real = curve_run(p, nexp)
    Lseq = n - 1
    rng = np.random.default_rng(p)
    rand_u2 = [U2(np.exp(2j * np.pi * rng.random(Lseq))) for _ in range(NREP)]
    mu_r, sd_r = np.mean(rand_u2), np.std(rand_u2)
    d = [(s_real[i + 1] - s_real[i]) % pe for i in range(Lseq - 1)]
    inc_seqs = []
    for _ in range(NREP):
        dd = d[:]; random.shuffle(dd)
        ss = [s_real[0]]
        for x in dd: ss.append((ss[-1] + x) % pe)
        inc_seqs.append(ss)
    h = (n - 1) // 2
    halfs = s_real[:h]
    dh = [(halfs[i + 1] - halfs[i]) % pe for i in range(h - 1)]
    inc_half = []
    for _ in range(NREP):
        dd = dh[:]; random.shuffle(dd)
        ss = [halfs[0]]
        for x in dd: ss.append((ss[-1] + x) % pe)
        inc_half.append(ss)
    syn_seqs = []
    for _ in range(NREP):
        first = [random.randrange(0, pe, p) for _ in range(h)]
        full = [0] * Lseq
        for i, k in enumerate(range(1, h + 1)):
            full[k - 1] = first[i]
            full[(n - k) - 1] = (Cn - first[i]) % pe
        syn_seqs.append(full)
    rows = {}
    fhat_max = {}
    for j in (1, 2, 3):
        u_real = U2(phases(s_real, j, p, pe))
        u_inc = [U2(phases(ss, j, p, pe)) for ss in inc_seqs]
        u_syn = [U2(phases(ss, j, p, pe)) for ss in syn_seqs]
        u_hr = U2(phases(halfs, j, p, pe))
        u_hi = [U2(phases(ss, j, p, pe)) for ss in inc_half]
        rows[j] = ((u_real - np.mean(u_inc)) / np.std(u_inc),
                   (np.mean(u_inc) - mu_r) / sd_r,
                   (u_hr - np.mean(u_hi)) / np.std(u_hi),
                   (np.mean(u_syn) - mu_r) / sd_r,
                   (u_real - np.mean(u_syn)) / np.std(u_syn))
        F = np.fft.fft(phases(s_real, j, p, pe)) / Lseq
        fhat_max[j] = float(np.abs(F)[1:].max())  # leading nontrivial coefficient
        Fs = np.fft.fft(phases(syn_seqs[0], j, p, pe)) / Lseq
        fhat_max[(j, 'syn')] = float(np.mean([np.abs(np.fft.fft(phases(ss, j, p, pe)) / Lseq)[1:].max() for ss in syn_seqs[:30]]))
    return n, rows, fhat_max

paper = {67: "+3.5..6.7", 163: "+3.7..7.4", 211: "+5.1..7.8", 349: "+5.4..7.1",
         433: "+8.2..9.6", 577: "+7.8..11.5", 823: "+10.5..12.5"}
results = {}
print(" p    n   j |  real/inc inc/rand half/inc syn/rand real/syn | paper real/inc")
for p, nexp in [(67, 79), (163, 139), (211, 199), (349, 313), (433, 397), (577, 613), (823, 829)]:
    n, rows, fh = stats_for(p, nexp)
    results[p] = (n, rows, fh)
    for j in (1, 2, 3):
        r = rows[j]
        print(f"{p:4d} {n:4d}  {j} |  {r[0]:+7.2f} {r[1]:+7.2f} {r[2]:+7.2f} {r[3]:+7.2f} {r[4]:+7.2f} | {paper[p]}")

# beta fit: |fhat|_max ~ C n^-beta across curves, per digit
print("\nbeta fit (|f^|_lead ~ C n^-beta), real vs synthetic (paper: real ~{0.34,0.45,0.36}, syn ~0.41 flat):")
ns = np.array([results[p][0] for p in results])
for j in (1, 2, 3):
    yr = np.array([results[p][2][j] for p in results])
    ys = np.array([results[p][2][(j, 'syn')] for p in results])
    br = -np.polyfit(np.log(ns), np.log(yr), 1)[0]
    bs = -np.polyfit(np.log(ns), np.log(ys), 1)[0]
    print(f"  digit {j}: beta_real = {br:.2f}   beta_syn = {bs:.2f}")

print(f"\nelapsed {time.time()-t0:.1f}s\nDONE2")
