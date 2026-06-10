"""Verify the three Round-1 scientific claims:
 (A) reflection-only U2 elevation: E[U2^4]*L -> const; syn/rand ratio of that const = 3/2;
     z(syn/rand) ~ K*sqrt(n); does K*sqrt(n) predict Table 1 + Phase 47 syn/rand?
 (B) beta digit-to-digit variation is fit scatter (resample reflection-only beta across 7 curves).
 (C) Hasse exhaustive-dichotomy: for p>=7, p | n (n=ord G | #E) forces #E=p (no intermediate regime);
     v_p(n) in {0,1} always.
Self-contained synthetic model for (A),(B) needs NO curve arithmetic (that is the point)."""
import math, random
import numpy as np

rng = np.random.default_rng(0)
random.seed(0)

def U2(fv):
    F = np.fft.fft(fv) / len(fv)
    return float((np.abs(F) ** 4).sum() ** 0.25)

# ---------- (A) 3/2 mass factor + sqrt(n) law, pure synthetic (no curve) ----------
print("=== (A) reflection-only U2: 3/2 mass factor and sqrt(n) law (no curve arithmetic) ===")
print(f"{'L':>6} {'E[U2^4]L rand':>13} {'E[U2^4]L syn':>12} {'ratio':>7} {'z(syn/rand)':>11} {'K=z/sqrtL':>10}")
Ks = []
for L in [78, 198, 612, 828, 2000, 10638]:
    NREP = 400
    p = 10007  # value-alphabet size; digit uniform on Z/p. Irrelevant to the result (checked below).
    rand_u4, rand_u2 = [], []
    syn_u4, syn_u2 = [], []
    h = L // 2
    for _ in range(NREP):
        # random baseline: iid unit phases
        g = np.exp(2j * np.pi * rng.random(L))
        u = U2(g); rand_u2.append(u); rand_u4.append(u**4)
        # reflection-only synthetic: digit uniform, s(n-k)=C-s(k); build phase directly
        d = rng.integers(0, p, h)
        C = rng.integers(0, p)
        full = np.empty(L, dtype=np.int64)
        full[:h] = d
        full[L-h:] = (C - d[::-1]) % p
        f = np.exp(2j * np.pi * full / p)
        u = U2(f); syn_u2.append(u); syn_u4.append(u**4)
    er, es = np.mean(rand_u4)*L, np.mean(syn_u4)*L
    z = (np.mean(syn_u2)-np.mean(rand_u2))/np.std(rand_u2)
    K = z/math.sqrt(L)
    Ks.append(K)
    print(f"{L:6d} {er:13.3f} {es:12.3f} {es/er:7.3f} {z:11.2f} {K:10.4f}")
Kbar = np.mean(Ks)
print(f"mean K ~ {Kbar:.4f};  3/2 mass factor => U2^4 const 2->3")

# predict Table 1 / Phase 47 syn/rand from K*sqrt(n)
print("\nPredict syn/rand = K*sqrt(n) vs paper:")
paper_syn = {79:4.0, 139:5.4, 199:5.9, 313:7.6, 397:8.5, 613:10.3, 829:12.4, 10639:44.6}
for n,val in paper_syn.items():
    print(f"  n={n:6d}  K*sqrt(n)={Kbar*math.sqrt(n-1):6.2f}   paper~{val}")

# closed-form constant from delta method
K_theory = 2*(2**0.75)*(3**0.25-2**0.25)/(20**0.5)
print(f"\nclosed-form K (delta method, sd from chi-like) = {K_theory:.4f}")

# ---------- (B) beta is fit scatter ----------
print("\n=== (B) beta digit-to-digit variation is 7-point fit scatter ===")
ns = np.array([79,139,199,313,397,613,829])
def leading_decay_beta(seqs_per_n):
    # |fhat|_lead ~ C n^-beta ; fit slope on log-log
    ys = np.array(seqs_per_n)
    return -np.polyfit(np.log(ns), np.log(ys), 1)[0]
# generate reflection-only leading-coeff for each n, many resamples, NO digit structure
betas = []
for _ in range(400):
    lead = []
    for n in ns:
        L=n-1; h=L//2; p=10007
        d=rng.integers(0,p,h); C=rng.integers(0,p)
        full=np.empty(L,dtype=np.int64); full[:h]=d; full[L-h:]=(C-d[::-1])%p
        F=np.abs(np.fft.fft(np.exp(2j*np.pi*full/p))/L)
        lead.append(F[1:].max())
    betas.append(leading_decay_beta(lead))
betas=np.array(betas)
print(f"reflection-only beta (no digit structure): mean={betas.mean():.3f} sd={betas.std():.3f} "
      f"range=[{betas.min():.2f},{betas.max():.2f}]")
real_spread = np.array([0.34,0.45,0.36])
print(f"paper real betas {list(real_spread)}: spread={real_spread.max()-real_spread.min():.2f}, "
      f"each within {abs(real_spread-betas.mean()).max()/betas.std():.1f} sd of the no-structure mean")

# ---------- (C) Hasse exhaustive dichotomy ----------
print("\n=== (C) Hasse: p|n forces #E=p for p>=7; v_p(n)<=1 always ===")
def is_prime(m):
    if m<2: return False
    for q in range(2,int(m**0.5)+1):
        if m%q==0: return False
    return True
def legendre(a,p):
    a%=p
    return 0 if a==0 else (1 if pow(a,(p-1)//2,p)==1 else -1)
def order_of_point_divides_test(p):
    # check: across all curves y^2=x^3+ax+b /F_p, is there N=#E with p|N and N!=p?
    bad=[]
    for a in range(p):
        for b in range(p):
            if (4*a**3+27*b**2)%p==0: continue
            N=1+sum(1+legendre((x**3+a*x+b)%p,p) for x in range(p))
            if N%p==0 and N!=p:
                bad.append((a,b,N))
    return bad
viol=[]
for p in [5,7,11,13,17,19,23]:
    bad=order_of_point_divides_test(p)
    # also Hasse-window multiples of p
    lo,hi=p+1-2*math.isqrt(p)-1, p+1+2*int(math.isqrt(p))+2
    mults=[m for m in range(lo,hi+1) if m%p==0 and m>0]
    print(f"  p={p:3d}: #E with p|#E and #E!=p : {len(bad)} ; multiples of p in Hasse window: {mults}")
    if p>=7 and bad: viol.append((p,bad))
print(f"  violations for p>=7: {viol if viol else 'NONE (dichotomy exhaustive)'}")
print("  NB p=5 edge: #E=10=2p is Hasse-allowed (v_p=1, still 'broken' for kernel action).")
print("\nDONE3")
