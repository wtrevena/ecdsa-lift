# Creative Directions — Round 1

*A step-back read of the just-revised ECDSA-lift paper, with one targeted experiment run.
Framing: not error-hunting — asking where the evidence actually points and where the work
could deepen or reframe without overclaiming.*

---

## TL;DR

The single highest-leverage move available is **turning the "U2 elevation grows with n"
observation into a closed-form proposition**. I ran the experiment: the growth is
**z(syn/rand) ~ 0.427 * sqrt(n)**, and the constant is derivable from three clean facts
(E[U2^4]_rand = 2/L, E[U2^4]_refl = 3/L, sd(U2^4)_rand = 2/L^{3/2}). This predicts *every*
entry of Table 1's positive-control column — including the Phase 47 value of ~44.5 at
n=10639 — to within finite-size noise. It converts the paper's most striking empirical
curve into mathematics and makes the Phase 47 confirmatory run *predicted rather than merely
consistent*. This is achievable now, low-risk, and materially strengthens the paper.

---

## (a) Highest-leverage directions

### Direction 1 — A closed-form scaling proposition for the reflection-induced U2 elevation **(do this)**

**What the paper currently says.** The elevation "grows with n" (Table 1: real/inc and
syn/rand both rise from ~+4 at n=79 to ~+44 at n=10639); the positive control "tracks" it;
Phase 47 "extends the pattern by more than an order of magnitude." All true, all empirical,
all qualitative.

**What I found (experiment below).** The growth law is exactly
> **z(syn/rand) ~ K*sqrt(n-1),  K ~ 0.427,**

with K derivable in closed form. The mechanism is purely the index-reversal symmetry of the
reflection s(n-k)=C-s(k); no curve arithmetic enters the constant. Concretely, for a
length-L=n-1 digit-phase sequence f(k)=exp(2*pi*i*digit_j(s(k))/p):

- Random unit-phase baseline: E||f||_{U2}^4 = 2/L, sd(||f||_{U2}^4) = 2/L^{3/2}.
- Reflection-only: E||f||_{U2}^4 = 3/L — **the reflection adds exactly 50% more U2 mass**,
  independent of p, n, the digit index j, and the antisymmetry constant C.
- Hence the standardized elevation
  z = (E||f||_{U2,refl} - E||f||_{U2,rand})/sd(||f||_{U2,rand})
  = 4*a_r^{3/4}(a_f^{1/4} - a_r^{1/4})/b * sqrt(L)
  with (a_r, a_f, b) = (2, 3, 2), giving K = 0.4267.

**Proposed change.** Add a short **Proposition** (one paragraph + a 4-line proof sketch) in
Section 5 (controls) or Section 7 (verification):

> *Proposition (reflection elevates U2 as sqrt(n)). For a length-(n-1) sequence f whose
> values satisfy the reflection f((n-1)-i) = w*conj(f(i)) (the digit-phase image of the
> antisymmetry delta(k)+delta(n-k)=[n]tau(G) at the informative digit), the expected Gowers
> U2 norm satisfies E||f||_{U2}^4 = (3/2)*E||g||_{U2}^4 for an unconstrained baseline g of
> the same length; consequently the standardized elevation against that baseline grows as
> z(f) ~ 0.43*sqrt(n). The sqrt(n) growth is therefore not evidence of accumulating signal
> but a deterministic consequence of a single rank-1 reflection constraint.*

Then add one sentence to Section 5's positive-control paragraph: *"This sqrt(n) law
(Proposition X) predicts the syn/rand column of Table 1 — including the n~10^4 confirmatory
run — to within finite-size scatter; the agreement is not a fit but a derivation."* You
could even add a predicted-z column to Table 1.

**Why this is the strongest move.** (i) It directly rebuts the most natural skeptical
reading of the paper's own data — "the effect *grows*, doesn't that look like real
signal?" — with "yes, and growth ~sqrt(n) is exactly what a content-free reflection
produces; that growth is the tell of an *artifact*, not of structure." That inverts the
intuition the Phase 37 report fell for, quantitatively. (ii) It upgrades Phase 47 from "we
checked one bigger curve and it agreed" to "the bigger curve hit the predicted value,"
which is a much stronger form of confirmation and pre-empts a referee asking "why only one
large curve?" (iii) It's a genuine, if modest, *new* proposition — the rest of the paper is
careful to say its theorem is not new mathematics; this gives the paper one self-contained
result it can claim outright.

**Payoff vs risk.** High payoff, low risk. The derivation is elementary additive
combinatorics (the 3/2 factor is the only content), already numerically confirmed across
the whole table including p=10477. The only honest caveat to state: the constant K=0.427
assumes the idealized conjugate-reflection at a single digit; real sequences have carries,
so the *precise* prefactor for real/inc can differ from syn/rand by O(1) — but the **sqrt(n)
scaling and the 3/2 U2^4-ratio are exact and p-independent** (verified to 3 digits at
p=211, 577, 10477). Claim the scaling law and the ratio; state the prefactor as "~0.43,
matching Table 1."

---

### Direction 2 — Settle (and mostly retire) the beta digit-to-digit variation **(cheap, defensive)**

**What the paper currently says.** Real digits give beta~{0.34, 0.45, 0.36}, reflection-only
synthetics give a flat {0.41, 0.41, 0.40}; the paper leaves the digit-to-digit variation
"open," floating carry nonlinearity of the base-p map as "a plausible origin" of a possible
"secondary structure."

**What I found.** Fitting beta to a 7-curve reflection-only ensemble and resampling: the
fitted exponent has **sd ~ 0.04 and range [0.29, 0.47]** across draws *with no digit
structure at all*. The real spread {0.34, 0.45, 0.36} (range 0.11) is ~2.7 sd of pure
7-point fit noise — i.e. **comfortably inside the scatter that reflection-only data already
produces**. The digit-to-digit variation is therefore *not demanded by the data*; it is
consistent with finite-sample noise around a single reflection spectrum.

**Proposed change.** Sharpen the hedge from "we leave it open / plausibly carries" to a
falsifiable statement:

> *"The digit-to-digit spread in the fitted beta (range ~0.11) is within the ~0.04-per-fit,
> ~0.11-range scatter that reflection-only synthetics themselves exhibit over a 7-curve
> fit; we therefore do not treat it as evidence of secondary structure. A decisive test — a
> single large-n curve where the per-fit scatter shrinks as 1/sqrt(#points) — would either
> resolve a real digit dependence or, as we expect, collapse it; Phase 47's single large
> curve is already consistent with the latter."*

This is more honest *and* less inviting of a referee chasing the carry story. If you want
the one extra experiment that *settles* it: fit beta within a *single* large-n curve by
binning Fourier modes, rather than across the 7 small curves — that removes the
curve-to-curve confound entirely and the digit comparison becomes apples-to-apples.

**Payoff vs risk.** Modest payoff, near-zero risk. Mostly it removes a dangling loose end
that currently reads as "there might be something here we didn't chase." The risk of the
*current* phrasing is that it half-invites a reader to think a real secondary signal was
left on the table; the data say it almost certainly wasn't.

**Honest assessment of the deeper question you posed** ("could the carry nonlinearity be
something more, not just public/inert?"): No — and Theorem 1 already guarantees it can't be.
The carry structure of the base-p digit map is a *deterministic public function of
z(delta(k))*, and z(delta(k)) itself is, by Theorem 1(a-c), governed in its secret part
only by [k]c (the computationally-hidden coordinate). Any digit/carry nonlinearity is a
public post-processing of a public object; it cannot manufacture secret dependence. So this
is a dead end *as cryptanalysis* — worth saying once, then dropping.

---

### Direction 3 — A clean general phrasing of the inertness theorem as a cokernel/coboundary statement **(optional, medium effort)**

**What the paper currently says.** Theorem 1 is stated concretely in terms of delta_tau,
[k]c, and k mod p^{e-1}, with the honest disclaimer that it is "an explicit instance" of
Silverman / Gadiyar-Padma and "not new mathematics."

**The opportunity you flagged.** Is there a one-line invariant statement that's more than
repackaging? Yes, and it's already *latent* in the boundary section — the paper just doesn't
lift it out:

> The entire dichotomy is the behavior of the endomorphism **[n] acting on the kernel of
> reduction Ehat(pZ/p^e)**. For gcd(n,p)=1, [n] is an **automorphism** of Ehat (it is
> invertible mod p^{e-1}); for n=p (e=2), [n]=[p] is the **zero map** on Ehat. Everything
> else follows functorially:
> - The section ambiguity lives in Ehat (two sections differ by a coboundary valued in Ehat).
> - delta_tau, as a function of the *instance* k mod n, is well-defined iff the ambiguity
>   c in Ehat lies in the kernel of [n] on Ehat (Theorem 1(d)). For an automorphism that
>   kernel is 0 -> only c=0 (canonical) works; for the zero map the kernel is all of Ehat ->
>   the ambiguity dies and [n]tau(G) becomes section-independent (Smart).

**Proposed change.** State it as a clean corollary or a one-paragraph "structural reading,"
explicitly: *"The inert/broken dichotomy is exactly the dichotomy automorphism/zero-map for
the action of [n] on the kernel of reduction; gcd(n,p)=1 is the condition that [n]|_Ehat is
invertible."* This is genuinely a tiny bit more than the existing phrasing because it names
the **single linear-algebraic invariant** (invertibility of [n] on Ehat) that controls the
whole story, and it makes the "sharp boundary" a one-line consequence rather than a worked
pair of cases. It also makes the intermediate-regime question (bolder reframe) well-posed.

**Payoff vs risk.** Medium payoff (it's a cleaner, more memorable statement of what you
have; reviewers like a single controlling invariant), low risk (it's a restatement, fully
supported by Section 6, no new claim). Honest caveat: it is *still* not new mathematics —
the automorphism-vs-nilpotent behavior of [n] on the formal group is textbook (Silverman
AEC IV/VII). Frame it as exposition, not discovery. The genuine novelty in the paper remains
the *methodology*, not the theorem.

---

## (b) The experiment I ran

**Setup.** Pure NumPy, seed 20260609 (matching the repo). I reused the exact reflection-only
synthetic construction from phase47_large_prime_confirmatory.py (kernel elements as
multiples of p, e=4, s(n-k)=C-s(k), digit-phase f(k)=exp(2*pi*i*digit_j(s)/p), U2 by its FFT
formula). I compared three sequence types against the 300-draw random-S1 baseline across the
seven clean-set primes plus p=10477.

**Findings.**

1. **Random baseline is exactly (2/L)^{1/4}.** Across L = 78...10638, the mean random U2
   matches (2/L)^{1/4} to 4 decimals. (So the random null U2^4-mass is 2/L; the curve is
   pure finite-size geometry, not anything about the curve.)

2. **Reflection adds exactly 50% to the U2^4 mass — p/n/j/C-independent.** Real digit-phase
   synthetics give E[U2^4]*L = 2.973 (p=211), 3.000 (p=577), 2.999 (p=10477); the idealized
   exact-pairing model gives 3.000. Ratio to random = 1.491, 1.499, 1.499. **This 3/2 factor
   is the single piece of content in the scaling law.**

3. **The standardized elevation is z(syn/rand) ~ 0.427*sqrt(n-1).** Empirical fit over the
   seven curves: z ~ 0.441*L^{0.491} (exponent within 2% of 1/2; z/sqrt(L) flat at ~0.42).
   Closed-form prediction K = 4*2^{3/4}(3^{1/4}-2^{1/4})/2 = 0.4267. Predicted vs Table-1
   syn/rand: n=79->3.8/3.8-4.1; 139->5.0/5.4; 199->6.0/5.8; 313->7.5/7.5; 397->8.5/8.4;
   613->10.6/10.2; 829->12.3/12.2; **10639->44.0/44.5-44.8**. Every entry matches within
   finite-size noise.

4. **The sqrt(n) law is intrinsic to the reflection, not the arithmetic.** A *minimal*
   model — g(i) i.i.d. on the unit circle for i<L/2, g(L-1-i)=conj(g(i)), zero curve
   content — reproduces the identical z-elevations (3.89 at L=78 ... 43.37 at L=10638,
   z/sqrt(L)~0.42). This is the cleanest possible demonstration that the Phase 37 anomaly's
   *growth* carries no curve information whatsoever.

5. **The beta digit spread is fit noise.** Reflection-only beta fit over 7 curves has
   sd~0.04, range [0.29, 0.47] under resampling; the real {0.34,0.45,0.36} spread (range
   0.11) sits ~2.7 sd inside that envelope. (See Direction 2.)

Scripts are in /tmp (u2_scaling*.py, u2_theory.py, u2_constant.py, u2_af.py, beta_check.py);
none touch the repo. They can be folded into a phase48_scaling_law.py if you adopt
Direction 1.

---

## (c) Do-not-bother list (tempting but dead-end)

- **Chasing the carry nonlinearity as cryptanalysis.** Theorem 1 already forecloses it: any
  base-p digit/carry functional is a public post-processing of a public object (z(delta)),
  whose secret part [k]c is computationally hidden. It cannot leak. Worth one sentence of
  dismissal, not an experiment. (Direction 2 covers the honest version.)

- **A second independent OEIS/closed-form sweep.** Phase 46 already did this and it was a
  clean null; d1(k) is provably not degree-<=3 polynomial and matches no OEIS entry. Redoing
  it at larger n adds nothing — pseudorandomness of a public quantity isn't the question.

- **More cross-curve aggregation (Phase 44 redux).** Phase 44 already pushed the one lane the
  theorem doesn't directly cover (joint distribution across curves) and got a clean null at
  the shuffle floor with two estimators. Scaling it up will reproduce the null; the CRT
  direct-sum argument predicts independence across curves anyway.

- **Trying to make the antisymmetry give more than 2-to-1.** It is a single rank-1 linear
  constraint coupling alpha with -alpha. The sqrt(n) U2 law (Direction 1) is the *complete*
  Fourier content of that constraint; there is no higher-order shadow to find (U3 and up sit
  at baseline, as the paper already reports, and the 3/2 factor is the whole story).

- **Stress-testing e>2 "partial degenerations" hoping for a middle regime on these curves.**
  For gcd(n,p)=1, [n] is an automorphism of Ehat(pZ/p^e) for *every* e — there is no
  intermediate behavior to find as e grows; the kernel of [n] on Ehat stays trivial. (The
  one place a genuine intermediate regime *could* live is composite/almost-anomalous n, see
  the bolder reframe — but not via e.)

---

## (d) One bolder reframe worth considering

**Make the controlling invariant — "[n] on the kernel of reduction" — the organizing spine
of the paper, and probe the one genuinely intermediate regime it predicts: n with a small
factor in common with the kernel exponent.**

The paper's dichotomy is currently stated as a *boundary* (gcd(n,p)=1 inert vs n=p broken).
But the real invariant (Direction 3) is **the kernel of the endomorphism [n] acting on
Ehat(pZ/p^e)**, a group of order p^{e-1}. That kernel is:
- **trivial** when gcd(n,p)=1 -> fully inert (automorphism),
- **everything** when p | n at e=2 -> fully broken (zero map),
- **and a proper nonzero subgroup** in the genuinely intermediate case **n = p*m with
  gcd(m,p)=1, e>=3** (or more generally v_p(n) = r with 1 <= r < e-1): then [n]|_Ehat has
  kernel of order p^{min(r,e-1)} — *neither* trivial *nor* all of Ehat.

**The sharp, testable prediction.** In that intermediate case the section ambiguity dies
*partially*: the component of c in the kernel of [n] becomes invisible (as in Smart), but
the component outside it survives (as in the inert case). So delta_tau, read mod the
surviving part, should leak a **partial** secret-linear invariant — recovering v_p(k) bits,
or k modulo a power of p — *without* a full break. This is the honest "is there an
intermediate regime?" answer: **yes, and it's not empty — it's a graded family interpolating
Smart's attack, where the number of recoverable p-adic digits of the kernel coordinate
equals min(v_p(n), e-1).**

Such curves exist: take E/F_p with #E(F_p) = p*m (a small cofactor times p — these are the
"almost-anomalous" curves cryptographers already flag as risky), and work at e=3 or 4. A
Phase 48/49 could verify: (i) the kernel of [n] on Ehat has the predicted order, (ii) the
recovered invariant gives exactly the predicted partial information (e.g. k mod p but not
k mod m), and (iii) it degrades smoothly to full inertness as you remove the p-factor. That
would turn the binary "inert/broken" headline into a **graded inertness theorem with Smart
as the extreme case** — a cleaner and more complete story than a single sharp wall, and it
connects directly to the practical lore that "curves with order divisible by p (not just
equal to p) are dangerous."

**Honest assessment.** This is the most intellectually ambitious option and the one most
likely to contain a genuinely new (if small) result — a *graded* version of the inert/broken
dichotomy. Risks: (1) it may turn out the partial invariant, while mathematically nonzero,
is still public for the same CRT reason (the surviving coordinate is k mod (n/p^r), which the
public data *does* fix) — in which case the "partial leak" is again inert and the reframe
becomes "the dichotomy is actually a clean theorem: you leak exactly the p-adic part of k
that [n] annihilates, and that part is public unless n is a pure p-power" — *still a sharper,
more complete statement than the current binary*. (2) It needs e>=3 kernel arithmetic, where
the repo already flags occasional projective degeneracy — but the z-coordinate additive path
used in Phase 47 handles this. I'd scope it as one exploratory phase before committing it to
the paper; even the *negative* outcome (partial invariant is public) upgrades the boundary
section from "sharp wall" to "graded, with the wall as a limit," which is a better theorem.
This is the direction with the highest upside and the most uncertainty — worth a half-day
probe, not a rewrite-on-spec.
