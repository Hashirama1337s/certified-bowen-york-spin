# Certified spin, mass and horizons of the Bowen-York spinning puncture family

**Moki & Julio**

Version 1.0.1 (2026-11-17). Every number below is traced to a result file in `CLAIMS.md`. Certified
lower ends are printed rounded down and upper ends rounded up, from the stored exact rationals.

---

## Abstract

We study the single, non-boosted, conformally flat, maximal Bowen-York spinning puncture (puncture boundary
conditions, J = J_ADM) with bare mass μ > 0, and its extreme limit μ = 0. With exact rational arithmetic and outward
ball arithmetic we certify sub/super-solution pairs of the Hamiltonian constraint and two-sided ADM mass enclosures for
every spin parameter j = 4J/μ² in [0, 11994.8], and a one-sided bound for every j > 0. The supremum of J/M_ADM² over
the family lies in [0.9282004538, 0.9282004687] (computer-assisted); it rounds to 0.928200, and this rounding
decision was replicated by a blind, independently formulated certificate. At fixed J ≠ 0 the mass M(μ) is strictly
increasing with difference quotients in (0, 1), strictly convex, and C^1, with M'(μ) equal to the exact tangent slope.
Monotonicity, slopes in [0, 1] and continuity at μ = 0 are corollaries of Dain and Gabach Clement's Lemmas 3.1-3.2,
which we record. It follows that the supremum is not attained and equals the spin of the extreme data, which we
certify directly on the throat: χ(0) = J/M_0² lies in [0.92820045875619, 0.92820045878266] (blind-replicated). For the
extreme data we prove a mass formula from the throat profile V alone, M_0 = (1/(6π)) ∫ (4|∇V|² + V²) dΩ =
(2/3) Σ (2l+1) a_l² (at J = 1), a Lichnerowicz analogue of the radial Pohozaev invariant of Korevaar, Mazzeo, Pacard
and Schoen. For μ > 0 the throat is the only closed connected MOTS and minimal surface and the strict area minimiser among
coordinate spheres. The extreme data contain no closed minimal surface, MOTS or apparent horizon, which answers a
question Dain and Gabach Clement left open ("we were unable to show this"); we give two independent proofs.
Limitations are listed in Section 11.

---

> **In plain words**
>
> A black hole that spins carries angular momentum J and mass M. In Einstein's theory the ratio J/M² can never exceed
> 1 for a black hole that has settled down. A popular way to set up a spinning black hole on a computer, the
> Bowen-York construction, has a lower ceiling. Numerical work since 1990 pointed to a ceiling near 0.928 (printed as 0.928 in 2002, and as
> 0.928200 in 2008).
>
> We proved where the ceiling is, with a computer that carries a guaranteed error bar on every number. No member of
> the family reaches it. The family only gets closer and closer as one parameter, the bare mass, goes to zero. In
> that limit the black hole turns into an infinitely long throat. We computed the spin of that limiting object
> directly: J/M² lies between 0.92820045875619 and 0.92820045878266.
>
> We also proved four facts about the family by hand. The mass changes without jumps or kinks (it is continuously differentiable) and is a convex curve in the bare
> mass. The mass of the limiting throat follows from one curve on a sphere. For a finite bare mass the black hole's
> surface sits exactly at the throat, with no other candidate anywhere. For the limiting throat there is no horizon at
> all.
>
> The idea in one picture: a ladder whose rungs get closer and closer to a ceiling. We weigh every rung with a
> certified scale, prove that the rungs only go up, and then measure the ceiling itself. Where the picture breaks:
> the ceiling is not a rung. It is a different object, the infinite throat, and we measure it with a separate
> certificate.

---

## 1. Results

Notation (Section 2): Φ = 1 + μ/(2r) + u, Δu = −S Φ^-7 with S = (9/4) J² sin²θ / r^6, M_ADM = μ + 2 lim r u,
j = 4J/μ², δ(j) = M_ADM/μ − 1 and χ = J/M_ADM² = j/(4(1 + δ)²). All statements concern the family qualified in the
abstract. "Computer-assisted" means exact rational residual polynomials and outward ball (arb) arithmetic.

**Theorem 1 (certified data, computer-assisted).** For every j in [0, 11994.815017806193245] an explicit pair
u_- ≤ u ≤ u_+ is certified, with certified enclosures of δ(j) and χ(j); for every j > 0, χ(j) is bounded above by a
certified piece (Section 3). Up to j = 1030 the point certificates have relative width 1e-7 or better (Tables 1-3).
Four independent blind certifiers reproduced the masses inside the certified intervals at eight spins, to 16-26
significant digits.

**Theorem 2 (supremum spin, computer-assisted).**

  sup over μ > 0 of J/M_ADM² lies in [0.9282004538, 0.9282004687]

(ends rounded outward to 20 decimals: 0.92820045389621230212 and 0.92820046861667179726). Since the upper end is below 0.9282005, the supremum
rounds to 0.928200 at six decimals. A blind certificate with a different formulation gives the tail bound
χ ≤ 0.9282004721 < 0.9282005 (Section 4). The supremum is not attained (Corollary 5).

**Proposition 3 (we record: corollaries of Dain-Gabach Clement Lemmas 3.1-3.2).** At fixed J, M(μ) is nondecreasing on
(0, ∞), its difference quotients lie in [0, 1], and M(μ) − μ ≤ M_0 ≤ M(μ), so M(μ) → M_0 as μ → 0 (Section 5). These
lines are not printed in the sources we opened, but they follow from the published lemmas, and the credit is theirs.

**Theorem 4 (mass function at fixed J ≠ 0).**
- (a) Strict monotonicity: 0 < (M(μ1) − M(μ2))/(μ1 − μ2) < 1 for μ1 > μ2 > 0.
- (b) Strict convexity: M((1−t)μ_a + tμ_b) < (1−t)M(μ_a) + tM(μ_b) for 0 < μ_a < μ_b and t in (0, 1). To our
  knowledge this has not been stated before (search to 2026-11-14); the mechanism is classical in spirit (Laetsch 1975,
  multiplier form).
- (c) M is C^1 on (0, ∞), and M'(μ) = 1 + 2 lim r w, where w solves the tangent problem Δw = 7SΦ^-8 (1/(2r) + w),
  w → 0. This is the exact tangent slope; 0 < M'(μ) < 1.

**Corollary 5.** χ is strictly increasing in j; the supremum over μ > 0 is not attained and equals χ(0) = J/M_0², the
spin of the extreme (μ = 0) data u_0 = lim u_μ (J_ADM = J at μ = 0; M_0 is the ADM mass of u_0, Appendix A).

**Theorem 6 (extreme spin, computer-assisted).** χ(0) lies in [0.92820045875619, 0.92820045878266] and
M_0/√J lies in [1.0379563989038, 1.0379563989187]. A blind certifier with a different formulation gives
χ(0) in [0.928200458542311, 0.928200458988139], which contains this window.

**Theorem 7 (throat flux mass formula).** For the extreme data at J = 1, with V the throat profile (the positive
solution of L_S V − V/4 + (9/4) sin²θ V^-7 = 0 on S²) and V = Σ a_l P_l(cos θ),

  M_0 = (1/(6π)) ∫_{S²} (4|∇V|² + V²) dΩ = (2/3) Σ_l (2l+1) a_l².

**Lemma 8 (μ > 0).** The throat r = μ/2 is the only closed connected MOTS and the only closed connected minimal
surface of the data (immersed surfaces cover it),
every closed embedded weakly outer-trapped surface lies in {r ≤ μ/2}, and the area of the coordinate spheres has its
strict global minimum at r = μ/2 (Dain's eq. (88)).

**Theorem 9 (no horizon at μ = 0).** For J ≠ 0 the extreme data contain no closed immersed minimal surface, no closed
MOTS for either normal, no closed embedded weakly outer-trapped surface, and hence no apparent horizon. Every
coordinate sphere is strictly mean-convex. Two independent proofs are given (Section 8).

**Probe only (not certified).** The extreme throat area gives A_0/(8πJ) ≈ 1.1991. This is a probe value. The
inequality A ≥ 8π|J| for the μ > 0 throats is due to Dain and Reiris (2011); it does not directly cover the μ = 0
data (cylindrical end).

| Item | Status | Section |
|---|---|---|
| Theorem 1 | certified (computer-assisted), refereed, blind-replicated | 3 |
| Theorem 2 | certified, refereed, decision blind-replicated | 4 |
| Proposition 3 | proved; credit Dain-Gabach Clement | 5 |
| Theorem 4 (a) strictness, (b), (c) | proved, refereed | 5 |
| Corollary 5 | proved, refereed | 5, 6 |
| Theorem 6 | certified, refereed, blind-replicated | 6 |
| Theorem 7 | proved, refereed; sealed float test held | 7 |
| Lemma 8 | proved, refereed | 8 |
| Theorem 9 | proved twice (one blind), refereed | 8 |
| A_0/(8πJ) ≈ 1.1991 | probe | 8.4 |

## 2. The data, conventions and the comparison lemma

**Data.** On R³ \ {0}, h = Φ^4 δ and K_ij = Φ^-2 Â_ij with the Bowen-York spin term
Â_ij = (3/r³)(ε_kil S^l n^k n_j + ε_kjl S^l n^k n_i), equivalently Â_ij = (6/r³) n_(i ε_j)kl J^k n^l as in
[DGC] eqs. (5)-(6). Then Â_ij Â^ij = 18 J² sin²θ / r^6 (checked symbolically; [DGC] eq. (8)), so the Hamiltonian constraint reads Δu = −S Φ^-7 with
S = (9/4) J² sin²θ r^-6. The data are maximal (tr K = 0) and Â(n, n) = 0 for the radial direction n.

**Scaling.** u'(x') = u(x'/λ) solves the problem with (λμ, λ²J), so δ = M_ADM/μ − 1 and χ depend on j = 4J/μ² only. We
compute at μ = 1, J = j/4. In Cook and York's units, a = μ/2 and J/a² = j (Section 10).

**Facts about the true solution (E1-E3).**
- E1. For every μ > 0 there is a solution u_μ with u → 0 at infinity, u in C^{2,α}(R³) (Dain-Gabach Clement p. 4-5,
  crediting Brandt-Brügmann and Cantor; Brandt-Brügmann p. 2). Uniqueness is not used by the certificates.
- E2. u_μ ≥ 0: Δu ≤ 0, u is continuous and tends to 0, so the minimum principle applies.
- E3. lim r u exists uniformly in angle, lim r u = (1/(4π)) ∫ S Φ^-7 dx, M_ADM = μ + 2 lim r u, and J_ADM = J
  (Dain-Gabach Clement p. 5). Equivalently M_ADM = μ + (1/(2π)) ∫ S Φ^-7 dx. (Dain-Gabach Clement eq. (13) prints the
  prefactor 1/(4π); with M = μ + 2 lim r u it must be 1/(2π), as in Brandt-Brügmann eq. (16) at first order; Section 10.)

**Lemma G (Green representation).** For μ > 0, u_μ = N[g_μ] with g_μ = S Φ_μ^-7 and N[g](x) = (1/(4π)) ∫ g(y)/|x − y| dy;
lim r u_μ = (1/(4π)) ∫ g_μ. Proof: g_μ ≤ S (2r/μ)^7 = O(r) near 0 (Lipschitz) and g_μ ≤ S = O(r^-6) at infinity, so
N[g_μ] is C², decays and solves ΔN = −g_μ; u − N[g_μ] is harmonic on R³ and tends to 0, hence vanishes. For the tail,
split at |y| = |x|/2: dominated convergence on the inner part and g ≤ C|y|^-6 on the outer part (O(|x|^-3)). At leading
order in J this gives δ = j²/40.

**Comparison class and Lemma C.** A function w is in class K if it is continuous on R³, C² on R³ \ {0}, tends to 0
uniformly at infinity, and satisfies |∇w| ≤ C r^-1/2, |Δw| ≤ C r^-3/2 near 0. For w in K the distributional Laplacian
has no point mass at 0 (the flux through small spheres is O(ε^{3/2})), and Kato's inequality Δw^+ ≥ 1_{w>0} Δw holds;
a continuous subharmonic function that tends to 0 at infinity is ≤ 0. Hence (Lemma C-REP4): if u_± are in K,
Φ_μ[u_±] > 0, Δu_+ ≤ −SΦ_μ[u_+]^-7 and Δu_- ≥ −SΦ_μ[u_-]^-7 on R³ \ {0}, then u_- ≤ u ≤ u_+, and the masses are bracketed
by 2 lim r u_±. The sub-side variant needs the inequality only on {u_- > u} (Lemma C-REP4').

**Lemma TL (tangent lemma).** Fix J and μ1. Let u_- be a certified sub-solution at μ1 with Φ_1 = Φ_{μ1}[u_-] > 0, let
V1 = 7SΦ_1^-8, and let ψ_+ in K satisfy Δψ_+ ≤ V1 (1/(2r) + ψ_+). Then for every μ' in (0, μ1],
u_- − (μ1 − μ')ψ_+ ≤ u_{μ'} and

  M_ADM(μ') ≥ M_-(μ1) − (μ1 − μ')(1 + 2 lim r ψ_+),   M_-(μ1) = μ1 + 2 lim r u_-.

Proof: the tangent inequality of the convex map x ↦ x^-7 gives −SΦ'^-7 ≤ −SΦ_1^-7 − (μ1 − μ')V1 (1/(2r) + ψ_+) for
Φ' = Φ_1 − (μ1 − μ')(1/(2r) + ψ_+), the bare-mass-μ' conformal factor of u~ = u_- − (μ1 − μ')ψ_+; on {u~ > u_{μ'}}
we have Φ' > 0 by E2, so Δu~ ≥ −SΦ'^-7 there and Lemma C-REP4' applies. The certified
coefficient 2L_hi = 1 + 2 lim r ψ_+ bounds M from below on (0, μ1]; it is built on the sub-solution, so it is not the
slope M'(μ1) and does not bound it (Section 5). With ψ_+ = 0 the lemma reproduces the plain monotone cell without
Dain-Gabach Clement Lemma 3.1. The inequality is elementary; what we claim is its use, a certified one-sided global
mass bound from one node.

## 3. Certified data for every spin (Theorem 1)

**P1: j in [0, 1] (sealed 565736c4c4).** Nineteen slabs in j; on each, an explicit sub/super-solution pair is certified
on the whole compactified square. The "for every J" statement holds at the slope-route accuracy (about 2e-11
relative); the balls below are at the evaluated spins.

Table 1. P1, μ = 1 (arb balls from `result_p1.json`).

| J | δ = M_ADM − 1 | χ |
|---|---|---|
| 0.0025 | [2.49998740540817e-6 ± 5.43e-21] | [0.002499987500109847 ± 3.31e-19] |
| 0.025 | [0.00024987415806877 ± 5.04e-18] | [0.02498751097331902 ± 2.15e-18] |
| 0.0625 | [0.0015576057054402 ± 5.39e-17] | [0.06230575324432787 ± 8.96e-18] |
| 0.125 | [0.006172885617036 ± 8.62e-16] | [0.123470951082942 ± 2.50e-16] |
| 0.25 | [0.0238361051533 ± 3.31e-14] | [0.23849491691357 ± 6.40e-15] |

Small spin, certified on the first slab: δ/j² → [0.025000000000 ± 2.13e-13], containing the published 1/40, and
(δ − j²/40)/j⁴ → [−0.00125946969697 ± 8.30e-15], containing −133/105600. The coefficient −133/105600 was derived
symbolically and re-derived by a referee seat by a different route and by a blind seat; in bare-mass units
M_ADM = μ + 2J²/(5μ³) − (266/825) J⁴/μ⁷ + O(J⁶). We did not find the J⁴ term in print (search to 2026-11-14); it follows
from the published first-order solution of Gleiser, Nicasio, Price and Pullin by one quadrature. Existence of the
limit defining it is not proved by the certificate.

**P2a: j in [1, 10.144] (sealed d62bfc0e0b).** 123 point certificates (grid ratio 51/50) and 122 monotone cells.
Cells use δ(j)/√j nondecreasing, a two-line corollary of Dain-Gabach Clement Lemma 3.1 plus scaling (credit theirs);
the cell inequality is checked in Lean (commit 036ea0f37b, three theorems; Lean covers only that inequality, not
Lemma 3.1 or the scaling).

Table 2. P2a and P2b point certificates, μ = 1 (decimals regenerated from the exact rationals; upper ends rounded up;
identical to `results/DISPLAY-TABLES.md`).

| j | certified δ = M_ADM − 1 | rel. width | certified χ |
|---|---|---|---|
| 1 | [0.023836105153325914, 0.023836105153325915] | 2.68e-22 | [0.23849491691356, 0.23849491691357] |
| 2 | [0.084685104039633339, 0.084685104039633341] | 5.36e-18 | [0.42497428536413, 0.42497428536414] |
| 4 | [0.250423126410051656, 0.250423126410055941] | 1.71e-14 | [0.63956693845623, 0.63956693845625] |
| 6 | [0.423044107116482903, 0.423044107116483036] | 3.12e-16 | [0.74072079402632, 0.74072079402633] |
| 8 | [0.586331612310199175, 0.586331612310201659] | 4.23e-15 | [0.79477104734360, 0.79477104734361] |
| 10 (P2a) | [0.738555349918561234, 0.738555349918583462] | 3.01e-14 | [0.82710941638178, 0.82710941638181] |
| 10.144253562168539132 | [0.749129523043965742, 0.749129523043991111] | 3.39e-14 | [0.82892676839664, 0.82892676839667] |
| 30 | [1.879329556666385383, 1.879329557821432295] | 6.15e-10 | [0.90464567733076, 0.90464567805657] |
| 100 | [4.202432956331064025, 4.202432958278665086] | 4.63e-10 | [0.92369166339863, 0.92369166409023] |
| 300 | [7.993758667606936304, 7.993758727426259015] | 7.48e-9 | [0.92721147547924, 0.92721148781339] |
| 1000 | [15.413187550627116440, 15.413188393220564507] | 5.47e-8 | [0.92801308498934, 0.92801318027100] |
| 1030.3467644609451531 | [15.660301597447063058, 15.660302819734156599] | 7.81e-8 | [0.92802064647561, 0.92802078264455] |

**P2b: j in [10, 1030.35] (sealed eba765ebe0).** 100 point certificates (ratio 21/20), 99 cells, degrees from (36, 10)
to (112, 14). It contains P2a's δ(10). The largest certified cell upper bound on χ over [10, 1030.35] is
0.9326738205130032 (cell [153.67, 161.36]).

**P2c: a uniform ceiling (sealed 3031217ac8).** Fourteen nodes (ratio 5/4) to j1 = 11994.815017806193245 and a tail
χ ≤ j1/(4 δ_lo(j1)²) for j ≥ j1 (from δ(j)/√j nondecreasing) give J/M_ADM² ≤ c = 0.96223488883111534543 for every
μ > 0. Per piece: P1 0.25, P2a cells 0.83394087233933423687, P2b cells 0.93267382051300310048, P2c cells
0.93989487162744096378, tail 0.96223488883111534543 (binding). A uniform ceiling below 1 for the whole family was not
found in print; per-member χ < 1 is published (Dain-Lousto-Zlochower 2008; Dain-Lousto-Takahashi 2002), and the
analytic barriers of Dain and Gabach Clement give only χ ≤ 675/256.

Table 3. P2c point certificates, μ = 1.

| j | certified δ | rel. width | certified χ |
|---|---|---|---|
| 3000 | [27.425202800381563409, 27.427267383083642380] | 7.53e-5 | [0.92809201147916, 0.92822683495341] |
| 10000 | [50.896029437271015641, 50.900045039378384190] | 7.89e-5 | [0.92812087160051, 0.92826450911423] |
| 11994.815017806193245 | [55.824678669721996208, 55.853551146282005455] | 5.17e-4 | [0.92772309989694, 0.92866608695635] |

**Blind replication.** Four fenced blind seats wrote independent certifiers from the equations (different radial maps,
bases, comparison functions and positivity methods). All sealed bars held:
- REP1 (e44dd02296): M_ADM(J = 0.25) in [1.02383610515332591435335055921, 1.02383610515332591435338447325], inside P1's
  ball; δ(0.0625) to 26 digits.
- REP2 (227cbcafc2): χ(10) in [0.827109416381791178798080260634, 0.827109416381791178800503789071]; δ(6), δ(10) inside
  P2a.
- REP3 (a669f700df): χ(1000) in [0.928013132630170991794, 0.928013132630171290517]; δ(300), δ(1000) inside P2b.
- REP4 (ad62cca292): δ(3000) in [27.42623509173264567580156405725267, 27.42623509173264567580156540633947] and
  δ(10000) in [50.89803723832470319819716946884112, 50.89803723832470319820952343563975], inside P2c.

Every blind seat also passed a small-spin control (δ/j² contains 1/40) and a null control (source 9/4 × (1 + 1e-3)
moves the enclosure off). REP4's comparison lemma, derivations and rounding were checked by a referee seat, and its
method became the second engine module (LEMMA-REP4).

**External checks.** Cook and York (1990), Table I, list Energy/a = 2(1 + δ) at J/a² = j. Our certified values round
to their printed entries at j = 1, 10, 30, 300, 1000 and 10000 (2.048, 3.477, 5.759, 17.99, 32.83, 103.8). At j = 100
our certified 10.404865-10.404866 is one unit below their last printed digit (10.41); at j = 3 an unsealed probe
enclosure (2.32846) is likewise one unit off (2.329). Both are far inside the 2% level at which Cook and York report agreement with earlier codes. Their
M_AH/a column agrees with our probe enclosures at all eight spins (Section 8.4). The Bowen-York row of Lousto, Nakano,
Zlochower, Mundim and Campanelli (2012), Table I, prints S/M_ADM² = 0.9282; REP4's δ(10⁴) gives
χ(10⁴) in [0.9281926861896467, 0.9281926861896468], which agrees at the printed digits (the spin of their row is our
inference, not stated in the text we read).

**Novelty (gate, to 2026-11-14).** To our knowledge this is the first computer-assisted, interval-arithmetic enclosure
of black-hole initial data solving the Einstein constraint equations of four-dimensional general relativity.
Computer-assisted proofs in gravity so far concern self-similar collapse (Reiterer-Trubowitz; a Zenodo mode count, 10.5281/zenodo.23270506),
a static black hole in Einstein-Weyl gravity (Goldstein-Jejjala 2026) and Riemannian Einstein and Ricci-flat metrics.
Analytic sub/super-solution bounds for this family, uniform but crude, go back to Dain and Gabach Clement. We do not
claim a "first certified black hole".

## 4. The supremum spin (Theorem 2)

**Construction.** Nodes j = 2^k, k = 0..21, each with a point certificate (Lemma C-REP4, sqrt-map engine) and a
tangent certificate (Lemma TL). Cover: [0, 1] by P1 (χ ≤ j/4 ≤ 1/4); each cell [2^k, 2^(k+1)] by Lemma TL from node
2^k, with μ' = √(j_a/j) in [1/√2, 1] and a rational lower bound of √(1/2) (the safe direction); and [2^21, ∞) by the
tail from node 2^21 (binding). The certified upper bound is c, and the certified lower end of the supremum is
X3 = χ_lo(2^21).

| piece | certified χ bound (rounded up) |
|---|---|
| tail [2097152, ∞) (binding) | 0.928200468617 |
| cell [1048576, 2097152] | 0.928200456416 |
| cell [524288, 1048576] | 0.928200454324 |

At the top node, δ(2^21) lies in [750.5607143681, 750.5607144170] and 2L_hi ≤ 5.911e-6.

**History of the window.** P2c [0.928146, 0.962235] (width 0.0340884); P2d (sealed 9677804a51, nodes to 2^18)
[0.9282003665, 0.9282009626], width 5.95979e-7; P2e (sealed 560847f6a4) [0.9282004538, 0.9282004687], width
1.47205e-8, about 2.3 million times narrower than P2c. The bar 0.9282005 lies 3.1383e-8 above c and 4.6103e-8 above
X3; to flip the decision, δ_lo(2^21) would have to be about 260 times its certified width too high, or the slope about
3.1 times its certified bound.

**Blind replication.** REP5 (9e0df55e7a) certified the tangent slopes independently at j = 64 and 4096:
s(64) in [0.0520407975655210747020518703815, 0.0520407975655210747020518836049], within 7e-10 of P2d's coefficient,
and s(4096) within 3.5e-10. It also caught a direction error in its own brief: Lemma TL bounds M directly, and its
coefficient is not a true-slope bound; the seat built correct one-sided bounds and a secant route that uses convexity
(Theorem 4 (b)). REP6 (0405395ed8) FAILED as sealed: its blind δ(2^21) in [750.5607142022, 750.5607145913] contains
ours, but its blind tangent coefficient (0.0140665) was 2380 times looser and its tail (0.928235200531265) did not
decide the digit. REP6b (3e929bb723) HELD: the same seat, tightening its own tangent, certified
δ(2^21) in [750.5607143858, 750.5607143991] (nested inside ours), 2L_hi ≤ 7.3143e-6 and χ_tail ≤ 0.9282004721 < 0.9282005
with a quarter-power/log radial map and polyellipse-Cauchy positivity. The margin is thin: the cover absorbs 2L_hi up
to about 1.8e-5.

**A second, tangent-free route.** With convexity (Theorem 4 (b)), M'(μ1) is at most the right chord, so two point
certificates at 2^20 and 2^21 bound the tail. On the certified P2e data this gives χ ≤ 0.92820047282 < 0.9282005
(sealed test B2V, a580a8ad0f).

**Comparison with print.** The value 0.928200 appears in Lovelace, Owen, Pfeiffer and Chu (2008), eq. (5), as
"ε_J ≤ 0.928200" from extrapolation; among the sources we opened it is the earliest six-digit print (Dain, Lousto and
Takahashi 2002, eq. (5), print 0.928). Our certified supremum rounds to the same six decimals. Read literally as a
strict bound, the printed inequality is exceeded by certified members: χ(2^21) ≥ 0.9282004538, about 4.5e-7 above it.
This is a remark on wording, not a defect.

## 5. The mass function at fixed J (Proposition 3, Theorem 4, Corollary 5)

All maximum principles used here are proved in Appendix C. Fix J ≠ 0 (at J = 0, M = μ and the strict statements fail).

**A1 (Dain-Gabach Clement Lemma 3.2, with a complete proof).** For μ1 > μ2 > 0, Φ_μ1 ≥ Φ_μ2 on R³ \ {0}. Let
D = Φ_μ1 − Φ_μ2 = (μ1 − μ2)/(2r) + (u_μ1 − u_μ2). D → +∞ at 0 and D → 0 at infinity. On {D < 0}, whose closure misses 0,
Φ_μ1^-7 > Φ_μ2^-7, so ΔD = −S(Φ_μ1^-7 − Φ_μ2^-7) ≤ 0; D = 0 on the boundary and tends to 0 at infinity, so the weak
minimum principle (Lemma C.1) gives D ≥ 0 there, a contradiction unless the set is empty.

**Proposition 3 (we record).** By Lemma G, r(Φ_μ − 1) → M(μ)/2, so lim r D = (M(μ1) − M(μ2))/2 ≥ 0. Dain-Gabach Clement
Lemma 3.1 (u_μ1 ≤ u_μ2) gives M(μ1) − μ1 ≤ M(μ2) − μ2. Hence 0 ≤ (M(μ1) − M(μ2))/(μ1 − μ2) ≤ 1. Letting μ' → 0 in
u_μ ≤ u_μ' and Φ_μ' ≤ Φ_μ gives u_μ ≤ u_0 and Φ_0 := 1 + u_0 ≤ Φ_μ, so (M(μ) − μ)/2 ≤ liminf r u_0 ≤ limsup r u_0 ≤
M(μ)/2 uniformly in angle; since M is monotone it has a limit M_0 = 2 lim r u_0 and M(μ) − μ ≤ M_0 ≤ M(μ). Dain and
Gabach Clement report that numerical calculations show the mass decreasing as μ → 0 and add "we did not prove this
analytically"; Gabach Clement's 2012 revision still says "we expect". The lines above use only their lemmas, the existence of lim r u
(Lemma G) and the mass read-off; we record them and credit them to Dain and Gabach Clement.

**Theorem 4 (a), strictness.** D ≥ 0 solves ΔD = cD on R³ \ {0}, with c = 7Sξ^-8 for ξ between Φ_μ2 and Φ_μ1, so
0 ≤ c ≤ 7SΦ_μ2^-8 ≤ C min(r², r^-6) and c ≤ C_V r^-6 for large r with C_V = 63J²/4. D blows up at 0, so it is not
identically zero, and the strong maximum principle (Lemma C.3) gives D > 0; let m_R = min over |x| = R of D. Choose
R^4 ≥ C_V and ψ = A(1/r + R/r²) on r ≥ R with A = m_R R/2. Then Δψ = 2AR r^-4 ≥ cψ (since cψ ≤ 2A C_V r^-7), ψ ≤ D on
|x| = R, and (Δ − c)(ψ − D) ≥ 0; on {ψ > D}, ψ − D is subharmonic, nonpositive on the boundary and at infinity, so
ψ ≤ D and lim r D ≥ A > 0. For the upper bound, M(μ1) − M(μ2) = (μ1 − μ2) + (1/(2π)) ∫ S(Φ_μ1^-7 − Φ_μ2^-7) dx, and the
integrand is negative wherever S > 0 because D > 0.

**Theorem 4 (b), strict convexity.** Let Φ_t = (1−t)Φ_a + tΦ_b = 1 + μ_t/(2r) + u_t with u_t = (1−t)u_a + tu_b. Since
x^-7 is convex, (1−t)Φ_a^-7 + tΦ_b^-7 ≥ Φ_t^-7, so Δu_t ≤ −SΦ_t^-7: u_t is a super-solution at μ_t. Put
d = u_t − u_{μ_t}, which is C² on R³ (the 1/(2r) terms cancel) and tends to 0. On {d < 0}, Φ_t < Φ_{μ_t} and
Δd ≤ S(Φ_{μ_t}^-7 − Φ_t^-7) ≤ 0, so d ≥ 0 by Lemma C.1. Then lim r d ≥ 0 reads (1−t)M(μ_a) + tM(μ_b) ≥ M(μ_t) (the μ terms
cancel). For strictness, Δd = cd − SB with 0 ≤ c ≤ 7S and B = (1−t)Φ_a^-7 + tΦ_b^-7 − Φ_t^-7 ≥ 0, and SB > 0 on an open
set (Φ_a ≠ Φ_b near 0, and S > 0 off the axis when J ≠ 0). So d is not identically zero, d > 0 by Lemma C.3, and the
ψ-comparison above (now (Δ − c)(ψ − d) = Δψ − cψ + SB ≥ 0) gives lim r d > 0. Convexity does not follow from the mass
integral alone: pointwise comparison of the integrands in E3 does not decide the sign; the 1/r tail is needed.

**Theorem 4 (c), differentiability.** For μ > 0 let V = 7SΦ_μ^-8 ≤ C min(r², r^-6).
- The tangent problem Δw = V(1/(2r) + w), w → 0, has a unique solution: existence by Lax-Milgram on D^{1,2}(R³) with
  the form ∫∇v·∇φ + Vvφ, local regularity because V and V/r are Lipschitz near 0, and w = −N[Vz] with z = 1/(2r) + w,
  so w decays; uniqueness among solutions that tend to 0 by Lemma C.1 applied to the difference.
- For 0 < |h| < μ/2 put q_h = (u_{μ+h} − u_μ)/h and z_h = 1/(2r) + q_h = (Φ_{μ+h} − Φ_μ)/h. A1 gives z_h ≥ 0 and
  Dain-Gabach Clement Lemma 3.1 gives q_h ≤ 0, so 0 ≤ z_h ≤ 1/(2r). z_h solves Δz_h = V_h z_h with
  V_h ≤ 7S(1 + μ/(4r))^-8 and V_h → V pointwise. The error e = q_h − w solves (Δ − V)e = f_h := (V_h − V)z_h with
  |f_h| ≤ C min(r, r^-7) uniformly in h and f_h → 0 pointwise.
- With P = N[|f_h|] (C² because |f_h| is locally Hölder), (Δ − V)(±e − P) = ±f_h + |f_h| + VP ≥ 0, so |e| ≤ P by
  Lemma C.1 (a) on balls B_R with R → ∞: V is bounded, e and P are continuous at 0 and tend to 0 at infinity, so the
  boundary term vanishes in the limit. By Lemma G, |lim r e| ≤ (1/(4π)) ∫ |f_h| → 0 (dominated convergence), so (M(μ+h) − M(μ))/h →
  1 + 2 lim r w. Also q_h → w pointwise.
- Signs: w ≤ 0 and z ≥ 0 as limits; z > 0 by Lemma C.3; M'(μ) = 2 lim r z > 0 by the ψ-comparison with c = V; and
  M'(μ) = 1 − (1/(2π)) ∫ Vz < 1. A convex differentiable function on an open interval is C^1.

The exact tangent slope 1 + 2 lim r w equals M'(μ1). The certified coefficient 2L_hi of Lemma TL is a different
object (built on the sub-solution's potential V1 ≥ V) and does not even imply 2L_hi ≥ M'(μ1), because M_-(μ1) < M(μ1).
At small spin the series identity 2L = 1 − (3/40) j² is dM/dμ at that order.

**Corollary 5.** χ(μ) = J/M(μ)² strictly decreases in μ (Theorem 4 (a)), so χ(j) = j/(4(1 + δ(j))²) strictly increases
in j, no member attains the supremum, and sup χ = lim_{μ→0} J/M(μ)² = J/M_0². It is the spin of the extreme data
because (i) u_0 = lim u_μ exists and is the extreme solution of Dain and Gabach Clement (their Thm 2.1); (ii)
J_ADM = J at μ = 0, since the angular-momentum integrand is conformally invariant for these data (Φ^-2 from K_ij,
Φ^-2 from the unit normal, Φ^4 from the area element); (iii) M_0 is the ADM mass of u_0 including derivative decay
(Appendix A). Uniqueness of u_0 is proved by Gabach Clement (Lemma 4.3, in H²) and by Waxenegger, Beig and Ó Murchadha
(Thm 3.1); see Section 10 for a remark on the former. Our statements are about u_0 = lim u_μ and use uniqueness only
to call it "the" extreme data.

**Sealed check on certified data (B2V, a580a8ad0f, HOLD).** On the 22 certified P2e nodes: all 20 three-node chord
inequalities decided in the convex direction, all 21 chords in (0, 1), χ increasing on all 21 pairs, and the
tangent-free tail above.

**Conjecture E (heuristic, not proved).** Matched asymptotics at the throat suggest M(μ) − M_0 ~ μ^{2ν} with
ν = √(λ_1 + 1/4), λ_1 the ground state of −L_S + 7(9/4) sin²θ V^-8 on the sphere. Probe: λ_1 = 1.65997689251753810055,
ν = 1.38201913608948921114 (arb midpoint, not an enclosure), so the local j-exponent of the chord slopes should approach
ν − 1/2 = 0.882019136. B2V's sealed test: at the highest pair whose interval exponent has width below 2e-4, the local
exponent is 0.882059 (width 5.1e-5), within 2e-3 of the prediction and excluding Lovelace et al.'s spherical-average
value √2 − 1/2 = 0.914213562. Lovelace et al. fitted a power law with exponent "−1.4" in S/m_p² (proportional to j); in j
this corresponds to ν = 1.382 (equivalently μ^2.8 against μ^{2ν} = μ^2.764); whether it is a two-digit shadow is not
established.

## 6. The extreme data (Theorem 6)

**Comparison at μ = 0 (B3-5).** Let w be C² on R³ \ {0} with 1 + w > 0, w → 0 uniformly at infinity, w → +∞ at 0 and
Δw ≤ −S(1 + w)^-7. Then u_μ ≤ w for every μ > 0 (on {u_μ > w}, Φ_μ[u_μ] > Φ_μ[w] > 1 + w > 0, so u_μ − w is subharmonic
there, and the set stays away from 0), hence u_0 ≤ w, and

  (1/(2π)) ∫ S(1 + w)^-7 dx ≤ M_0 ≤ 2 limsup r w.

The upper bound follows from M(μ) ≤ μ + 2 limsup r w and Proposition 3; the lower bound from E3 in integral form,
monotone convergence as μ → 0 and Proposition 3. Integrability is automatic (the integral is at most 2πM_0). One
super-solution gives both ends; the gap between them is its integrated residual.

**B3T (sealed 10a65e940a, HOLD).** On the throat, σ = log r, x = tanh(σ/8), so r^{1/4} = (1+x)/(1−x) and
√r = ((1+x)/(1−x))². The ansatz w = r^-1/2 G/(1 + √r), G = g0 + (1−x)² H(x, cos θ) with Chebyshev times even Legendre
modes (N = 112, K = 12), is solved in floating point and then certified on an exact cover of [−1, 1] × [0, 1] by 17
patches: tensor Chebyshev interpolants at 65 × 65 Lobatto nodes in arb (256 bits), with the interpolation error bounded
through a Cauchy estimate on a Bernstein polyellipse and analyticity checked by pluriharmonic real parts. The shift
G → G + η lowers the normalised residual by exactly η. Certified: sup Q_N ≤ 2.16576e-12 < η = 2.51977e-12,
inf G_N ≥ 0.5183644946. Then M_up = 2(g0 + η) exactly, and M_lo is a rigorous lower quadrature (loss ≤ 8.9e-16):

  M_0/√J in [1.0379563989038, 1.0379563989187],   χ(0) in [0.92820045875619, 0.92820045878266],

width 2.647e-11, equal to the integrated residual. Controls: extreme Reissner-Nordström (exact M = 1) inside its window;
source × 1.001 moves the window off P2e's; a forced shift below the certified residual fails. The window lies inside
P2e's, as Corollary 5 requires. A referee seat re-derived the residual, the shift direction, the mass integrand and
the polyellipse bound, and reran the certificate byte-identically.

**REP7 (sealed 133b99d241, HOLD; blind).** A different formulation (ψ = √(1 + F/r), a far element in 1/r, Lobatto
elements in t = −log r, a Newton-potential comparison) gives M_0 in [1.037956398789003, 1.037956399038277] and
χ(0) in [0.928200458542311, 0.928200458988139]; a second discretisation gives [0.928200458680027, 0.928200458850410].
Both contain B3T's window. Its RN control contains 1.

**Wording (referee seat).** Certified: sup_{μ>0} J/M² = χ(0) of the extreme (μ = 0) Bowen-York data (J_ADM = J at
μ = 0; u_0 unique by Gabach Clement Lemma 4.3, with u_0 = lim u_μ in that class by Dain-Gabach Clement Thm 2.1) lies in
[0.92820045875619, 0.92820045878266]. In print: 0.928 at μ = 0 (Hannam, Husa and Ó Murchadha 2009, numerical) and
0.928200 extrapolated (Lovelace et al. 2008, numerical). "Directly at μ = 0" describes the computation only; the proof
of the bounds runs through the μ > 0 family (B3-5). A certified value at μ = 0 was not found in print (gate to
2026-11-14).

## 7. The throat flux mass formula (Theorem 7)

**The autonomous throat equation.** Put σ = log r and Ω = r^{1/2} Φ. Since Δ(r^-1/2 X) = r^-5/2 (X_σσ − X/4 + L_S X),
with L_S the Laplacian of the round sphere, the constraint becomes, for μ = 0 and μ > 0 alike,

  (E)  Ω_σσ + L_S Ω − Ω/4 + (9/4) J² sin²θ Ω^-7 = 0   on R × S².

The equation has no explicit σ. For the extreme data Ω = e^{σ/2} + r^{1/2} u_0.

**Flux.** At J = 1 let E(σ) = ∫_{S²} [Ω_σ²/2 − |∇_S Ω|²/2 − Ω²/8 − (3/8) sin²θ Ω^-6] dΩ. Then
dE/dσ = ∫ Ω_σ [Ω_σσ + L_S Ω − Ω/4 + (9/4) sin²θ Ω^-7] dΩ = 0 after one integration by parts on the closed sphere (Ω is
smooth on S² because u_0 is C^∞ on R³ \ {0}, Gabach Clement Thm 2.1).

**Ends.** At σ → +∞, Ω = r^{1/2}(1 + M_0/(2r) + ρ) with ρ = o(r^-1), ∂_r ρ = o(r^-2), ∂_θ ρ = o(r^-1/2) (Appendix A); the
ρ = 0 integrand is exactly −M_0/4, so E → −πM_0. At σ → −∞, Ω → V in C^1 uniformly (Appendix B), so
E → −∫ [|∇V|²/2 + V²/8 + (3/8) sin²θ V^-6] dΩ. Multiplying the sphere equation by V and integrating gives
(9/4) ∫ sin²θ V^-6 = ∫ |∇V|² + ∫ V²/4. Substituting,

  M_0 = (1/(6π)) ∫ (4|∇V|² + V²) dΩ = (2/3) Σ (2l+1) a_l²,

using ∫ P_l² dΩ = 4π/(2l+1), ∫ |∇P_l|² dΩ = 4π l(l+1)/(2l+1) and 4l(l+1) + 1 = (2l+1)². For a sphere equation
L_S V − V/4 + C V^-p = 0 the same steps give M = (1/π) ∫ [|∇V|²/2 + V²/8 + C V^{1−p}/(p−1)]; extreme Reissner-Nordström
(p = 3, C = q²/4, V = √q) returns M = q.

**Status.** A referee seat re-derived every step symbolically and required the two limit arguments, now written in
Appendices A and B. A sealed slice test (FLUXT, 0ce319a05a, HOLD) evaluated E(σ) at ten slices σ in [−6, 6] on the
independent two-dimensional float solution of B3T, which was computed with no knowledge of the identity: E was
constant to a spread of 5.773e-15 against a float tolerance of 2.562e-12 stated before evaluation, and |E/(−π) − M_0|
≤ 8.882e-16; the null (gradient coefficient 1/4 instead of 1/2) drifted by 7.098e-3 (2.77e9 tolerances). The tolerance
is a float quantity, not a rigorous one.

**Probe digits (not certified).** Solving the one-dimensional sphere problem by Legendre collocation and Newton at
high precision gives M_0 = 1.03795639891365137036534431609 and χ(0) = 0.9282004587652052045098878 (truncated). A
referee seat's independent shooting solver reproduced M_0 to 25 digits. Two truncations of one code agree to about 180
digits; these are probe digits. Both values lie inside the certified windows of Theorem 6 (a post hoc check: the
derivation came after those windows were known). A certified one-dimensional enclosure is designed but not built.

**Credit.** The mechanism is not new. Korevaar, Mazzeo, Pacard and Schoen (1999) built the same conserved quantity,
the "radial Pohozaev invariant", for the Yamabe equation on the cylinder t = −log|x| with v = |x|^{(n−2)/2} u (for
n = 3 exactly our Ω); for n = 3 their kinetic, gradient and −v²/8 terms are ours, and only the potential differs.
They also evaluate it at the singular end from limits of translates. The Pohozaev integral at an asymptotically flat
end equals a constant times the mass (Schoen; see Marques 2004, Corvino-Pollack 2011). Chruściel and Mazzeo (2012)
wrote the conserved Hamiltonian for x-only Lichnerowicz solutions on a cylinder, without a mass link. The
one-dimensional throat equation and the area A_0 = 2π ∫ D⁴ sin θ dθ are due to Hannam, Husa and Ó Murchadha (2009,
eqs. (22), (25)-(26)), who note it can be computed "to arbitrary accuracy", and to Gabach Clement (eq. (56)). What is
ours: the Bowen-York potential and the observation that the r^-6 sin²θ source keeps the equation autonomous, the
angle-dependent throat with its gradient term (Korevaar et al. Prop. 1 forces radial profiles in the Yamabe case), the
virial reduction and the Legendre sum of squares. To our knowledge this mass formula, and hence a one-dimensional
route to J/M², has not been written before (search to 2026-11-14). Dain (2010) wrote that there is "no obvious
relation" between his (different) throat mass functional and the total mass.

## 8. Horizons (Lemma 8, Theorem 9)

Throughout, H^N = div_Σ N (round spheres have H > 0 for the outward normal), θ^±_N = ±H^N + tr_Σ K, and
tr_Σ K = tr K − K(N, N). For these data tr K = 0 and K(ν, ν) = 0 for the radial unit normal ν, for either sign
convention of K. For the coordinate spheres S_r,

  H(S_r) = Φ^-2 (2/r + 4Φ_r/Φ) = 4 r^-3/2 Φ^-3 Y = 4w/(rΦ³),   Y = ∂_σ Ω,   w = Φ/2 + r∂_r Φ = r^-1/2 Y,

and the area is A(r) = ∫ Ω⁴ dΩ (since r²Φ⁴ = Ω⁴), with dA/dσ = 4 ∫ Ω³ Y dΩ. Differentiating (E) gives
(∂_σσ + L_S − 1/4 − V)Y = 0 with V = (63/4) J² sin²θ Ω^-8 ≥ 0.

**Touching lemma.** Let Σ be a C² immersed closed surface and p a point where r|Σ attains its maximum r_0. Then
H^ν_Σ(p) ≥ H(S_{r_0})(p). Proof: Δ_Σ(F|Σ) = tr_Σ Hess F − H^N dF(N) for F = r and any unit normal N; at p the tangent
planes agree, N = ν, and on S_{r_0} the same identity gives tr Hess F = H(S_{r_0}) |∇F|; since Δ_Σ(F|Σ)(p) ≤ 0, the claim
follows. The same holds at a minimum of r with F = −r and the normal −ν. This is the standard barrier argument used by
Benko and Stavrov Allen (2021, Prop. 3.4) and Reiris (arXiv:1209.4530, proof of Cor. 3) in other data.

### 8.1 Lemma 8 (μ > 0)

Let a = μ/2, σ_a = log a and Ω = e^{σ/2} + a e^{−σ/2} + r^{1/2} u_μ. **Claim:** Y > 0 for σ > σ_a, Y < 0 for σ < σ_a,
Y = 0 on σ = σ_a.

Proof (reflection, needing neither inversion symmetry nor derivative bounds). For λ > σ_a set
w(σ) = Ω(2λ − σ) − Ω(σ) on {σ < λ} × S². The reflection preserves (E), so
w_σσ + L_S w − w/4 = −(9/4) J² sin²θ [Ω(2λ − σ)^-7 − Ω(σ)^-7]; w = 0 on σ = λ; and as σ → −∞,
w ≥ (e^λ − a)e^{−σ/2} − O(e^{σ/2}) → +∞ (using u_μ ≥ 0 and u_μ bounded). At a negative interior minimum,
w_σσ ≥ 0, L_S w ≥ 0 and −w/4 > 0 make the left side positive, while Ω(2λ − σ) < Ω(σ) makes the right side ≤ 0.
So w ≥ 0, and ∂_σ w(λ) = −2Y(λ) ≤ 0. The mirror argument for λ < σ_a uses the upper barrier
u_μ ≤ √(1 + Q/r) − 1 (Gabach Clement eqs. (29)-(31)) at σ → +∞. Strictness on each side follows from Lemma C.3
(transferred to the cylinder, Appendix C), since Y ≡ 0 on either side is impossible.

Consequences: H(S_r) > 0 toward larger r for r > a and toward smaller r for r < a. At the max-r point of any closed C²
surface reaching r > a, and at the min-r point of one reaching r < a, the touching lemma gives θ^N ≠ 0 and H ≠ 0. A
closed MOTS or minimal surface therefore lies in S_a; a closed surface immersed in S² covers it, so it is S_a. S_a is a
MOTS (Y = 0 there, Cook-York sec. V). The area has the sign of σ − σ_a in its derivative, so r = a is the strict global
minimum among coordinate spheres; Dain (2010, sec. VI, eq. (88)) took this "by definition of minimal surface".

Credit and gate: that r = a is a MOTS is known (Cook-York 1990; "well-known" in Lovelace et al., crediting York and
Piran 1982 [not opened; as cited by LOPC]). Cook and York wrote that this "does not preclude the existence of other
apparent horizons". Uniqueness of the closed connected MOTS and minimal surface is new to our knowledge; least area among coordinate spheres was previously
asserted (Dain 2010) and seen numerically (Lovelace et al. 2008, Fig. 3), and is now proved.

### 8.2 Theorem 9, first proof (sliding)

Inputs (cited, opened): (A0) Gabach Clement Thm 2.1: u_0 = V/(√r (1 + b√r)) + U with V in C^∞(S²) positive, U in
H'²_{−1/2}; (A1) Dain-Gabach Clement 2010, Lemma A.1, eqs. (78)-(79) and the annulus argument (82)-(89): for f in
H'²_δ, r^{−δ}|f| → 0 at 0 and at infinity uniformly in angle; (A3) 0 ≤ u_0 ≤ √(1 + Q/r) − 1 (Gabach Clement
(29)-(31)). The proof of (A0) in print contains a step that needs the repair of Appendix D.

- A2: |Ω_0 − V| ≤ √r + b c_+ √r + sup r^{1/2}|U| → 0 uniformly as σ → −∞; and e^{σ/2} ≤ Ω_0 ≤ √(e^σ + Q).
- Sliding: fix t > 0 and w(σ) = Ω_0(σ + t) − Ω_0(σ). w → 0 uniformly at σ → −∞ (A2) and w → +∞ at σ → +∞ (A3). If
  inf w < 0, it is attained, and at that point the left side of
  w_σσ + L_S w − w/4 = −(9/4) J² sin²θ [Ω_0(σ+t)^-7 − Ω_0(σ)^-7] is positive while the right side is ≤ 0. So Ω_0 is
  nondecreasing in σ and Y_0 ≥ 0.
- Strictness: Y_0 solves the linear equation with zeroth-order coefficient −(1/4 + V) ≤ 0; if Y_0(p) = 0, Lemma C.3
  forces Y_0 ≡ 0, which contradicts Ω_0 → V at one end and Ω_0 ≥ e^{σ/2} → ∞ at the other.

So H(S_r) > 0 at every point of every coordinate sphere. Note that H → 0 as r → 0 is not attained; the argument
needs only pointwise positivity, because a closed surface stays at r ≥ r_min > 0.

The sliding method is classical (Berestycki and Nirenberg 1991, as described in Rivière's survey of Nirenberg's work,
Sec. VII.2: "comparing the solution with its translations along the axis of the cylinder"). To our knowledge this is
its first use for the Lichnerowicz equation on asymptotically cylindrical initial data.

### 8.3 Theorem 9, second proof (dilation barrier; blind)

A blind seat, given only the statement and the inputs, proved positivity of the same function w = r^-1/2 Y by a
different route.
- Dilations x ↦ λx map solutions of the extreme equation to solutions (Φ_λ = λ^{1/2}Φ(λx)), and w = dΦ_λ/dλ at λ = 1
  solves (Δ − Q)w = 0 on R³ \ {0} with Q = 7 · (9/4) J² sin²θ r^-6 Φ^-8 ≥ 0 (symbolic identity, with negative controls).
- Two-sided throat bound: (c_-/2) r^-1/2 ≤ Φ ≤ 2c_+ r^-1/2 near 0, from (A0) and (A1), or independently from Dain-Gabach
  Clement's barriers (eqs. (16), (96)). An interior gradient estimate then gives w ≥ −C_1 r^-1/2 near 0, and w → 1/2 at
  infinity.
- Barrier: v_ε = w + ε r^{-3/4} satisfies (Δ − Q)v_ε = ε[(3/4)(−1/4) r^{-11/4} − Q r^{-3/4}] < 0. On an annulus
  {ρ < r < R} with ρ^{1/4} < ε/C_1 and R large, v_ε > 0 on the boundary, and an interior negative minimum would give
  (Δ − Q)v_ε ≥ 0. So v_ε ≥ 0, hence w ≥ 0, and w > 0 by Lemma C.3.

Mechanism: Schwarzschild and the μ > 0 data have w ~ −r^-1 near 0 (the borderline exponent, where minimal surfaces
exist), while the extreme throat has w ≥ −C r^-1/2, strictly below it. This proof avoids Gabach Clement's decomposition
if the barrier bounds are used, and it does not use uniqueness.

### 8.4 The theorem, consequences, and what is probe

**Proof of Theorem 9.** Let Σ be closed (compact in M = R³ \ {0}) and p its max-r point.
- Minimal: 0 = H^ν_Σ(p) ≥ H(S_{r_0})(p) > 0 is impossible.
- MOTS for any unit normal N: N(p) = ±ν and tr_Σ K(p) = 0, so θ^+_N(p) = ±H^ν_Σ(p) ≠ 0; θ^- = 0 is excluded since
  θ^-_N = −θ^+_{−N}.
- Weakly outer-trapped (embedded, θ^+ ≤ 0, outer toward the asymptotically flat end): Σ separates R³; the ray beyond p
  lies in {r > r_0}, misses Σ and is unbounded, so the outer normal at p is ν and θ^+(p) = H^ν_Σ(p) > 0.
- An apparent horizon would be a closed MOTS.

Surfaces that wrap the cylinder are covered, because a closed surface cannot reach the end. Scope: P = 0 and no A-term
(otherwise K(ν, ν) ≠ 0 at p), data only (nothing about an event horizon).

Dain and Gabach Clement (Sec. 2) state that they believe the extreme solution has no minimal surface or apparent
horizon, by analogy with extreme Kerr-Newman and as indicated by numerics, and add: "we were unable to show this". No
later source we found closes it (gate to 2026-11-14). The physical expectation is already
in print: Dain, Lousto and Zlochower (2008) saw numerically that "Initially the horizon has a coordinate radius of
zero", Hannam, Husa and Ó Murchadha (2009) write that "the horizon is located at the puncture", and the exact extreme
Kerr-Newman and Reissner-Nordström analogues are classical. Theorem 9 is a proof of that expectation, not a surprise:
the horizon of these data exists only as the asymptotic cylindrical end, at infinite distance.

**Area corollaries (proved).** dA_0/dσ = 4 ∫ Ω_0³ Y_0 > 0, so the area of the coordinate spheres of the extreme data
decreases strictly to A_0 = ∫ V⁴ dΩ as r → 0 (by A2). The limit value is Hannam et al.'s eq. (26); the strict monotone
approach is new to our knowledge. With Φ_0 < Φ_μ strictly, every μ > 0 throat has area A_μ(μ/2) > A_0(μ/2) ≥ A_0.

**Probe values (not certified).** The one-dimensional probe gives A_0/J = 30.13745891933985963534 and hence
c_BY = A_0/(8πJ) ≈ 1.1991 (only the lower bound A_μ(μ/2) > A_0 is proved here; that the constant is sharp,
A_μ(μ/2) → A_0, is not proved here). It is not certified; a certified A_0 window is designed but not built. The
inequality A ≥ 8π|J| for the μ > 0 throats is Dain and Reiris's theorem (2011), which does not directly cover the μ = 0
data (cylindrical end); we claim only the probe value of c_BY. That Bowen-York horizons stay about 20% away from the extreme-Kerr equality is known numerically (Lovelace et
al. 2008, χ_H ≤ 0.9837; Hannam et al. 2009, 0.9837). Probe enclosures of M_AH/a at Cook and York's eight spins
(engine sandwich on r = a, exact polynomial area integral, unsealed) match all eight printed values: 2.033, 2.227,
3.034, 4.673, 8.071, 13.67, 24.69, 77.56.

## 9. Sealed tests and verification

Every test below was sealed (prediction and bar hashed) before its run and scored as written.

| test | seal | what | verdict |
|---|---|---|---|
| P1 | 565736c4c4 | every J in [0, 1/4]; small-spin 1/40 | HOLD (process slip disclosed) |
| P2a | d62bfc0e0b | every j in [1, 10.144] | HOLD |
| P2b | eba765ebe0 | every j in [10, 1030.35] | HOLD |
| P2c | 3031217ac8 | uniform ceiling c < 1 | HOLD |
| P2d | 9677804a51 | tangent-lemma window | HOLD |
| P2e | 560847f6a4 | sixth digit decided | HOLD |
| REP1-REP4 | e44dd02296, 227cbcafc2, a669f700df, ad62cca292 | blind point certificates | HOLD (all four) |
| REP5 | 9e0df55e7a | blind tangent slopes | HOLD |
| REP6 | 0405395ed8 | blind top node incl. tail | FAILED (tangent too loose) |
| REP6b | 3e929bb723 | blind top node, tightened | HOLD |
| B2V | a580a8ad0f | convexity, chords, tail, exponent | HOLD |
| B3T | 10a65e940a | χ(0) on the throat | HOLD |
| REP7 | 133b99d241 | blind χ(0) | HOLD |
| FLUXT | 0ce319a05a | flux constant on 10 slices | HOLD |

Disclosures. P1: the controls check returned NOT GO (untested resume path) and the launch went ahead; the run used one
invocation, so the untested path was never exercised. P2c: the build seat had read a logged estimate of c (nothing
tuned). P2d/P2e: degree schedules and the N-TAN control table were chosen after rigorous probes of the same nodes (not
blind; validity unaffected); the sealed CEIL control of P2d ran on a mock grid, and the real-data recomputations are
post hoc (P2e runs them inside the instrument). REP1-REP4 seats saw partial context (disclosed per verdict); REP4's
instrument was pinned after its run; REP5 saw earlier summaries through a stray scratch file (no slope values). B3T:
N, K and the η factor were tuned on float probes of the same object. An extreme-data monotonicity probe (pre-stated,
unsealed) failed because its bar sat below the solve's noise floor; it tests nothing about the theorem.

Referees. "Seats" are separate internal working sessions of the lab, each with its own files; referee seats and blind
seats are internal checks, not external peer review. A blind seat did not see our code or the values it was asked to reproduce; partial context is disclosed per test (Section 9). Fifteen internal
hostile referee reports cover the equation and normalisation (P1), the monotone-cell lemmas, each
certificate run (P2a-P2e), the sqrt-map certifier and LEMMA-REP4, the tangent lemma, the Wing-1 theorems (Sections 5-6),
B3T, the horizon lemma, the flux identity and the no-horizon proof. No report found an error in a certified value or in
a proof used here. Items judged WRONG or CONDITIONAL (a prior-art framing, the then-unproved extreme variant of Lemma 8,
a proposed re-bar of a probe, and a K = 0 engine branch that stamped a false certificate and was never used in a run)
were fixed or dropped; the wording fixes are listed in Appendix E. A sixteenth internal report (REFEREE-RELEASE) checked
this text. Lean (an interactive proof assistant whose kernel checks each step) covers the monotone-cell inequality only
(036ea0f37b).

## 10. Fairness notes and corrections

1. **Gabach Clement 2012, Lemma 4.2.** U is defined by Δ U = Δ(u − V/(√r (1 + b√r))) (eq. (65)), while Lemma 4.2
   produces the unique solution of (Δ − h/(4r²))U = f in H'²_{−1/2} (eq. (66)). The identification (71)-(72),
   "Δ[u − (V/... + U)] = 0", does not follow as printed: the difference D solves (Δ − h/(4r²))D = 0, not ΔD = 0. A
   three-line barrier argument repairs it (Appendix D). The boundedness of h needs a lower bound √r u ≥ c > 0; for u_0
   it holds by Gabach Clement's barrier (eqs. (38), (43)-(44)) and by Dain-Gabach Clement (eqs. (16), (96)). For an
   arbitrary H² solution, which the uniqueness proof of Lemma 4.3 needs, the text as read does not establish it.
   Uniqueness is independently proved by Waxenegger, Beig and Ó Murchadha (Thm 3.1). None of our results uses
   uniqueness except to name the extreme data. Both the referee seat and the blind seat flagged this step
   independently.
2. **Lovelace et al. 2008, "ε_J ≤ 0.928200".** Read literally, certified members exceed it (Section 4); as a six-decimal
   rounding of the limit it is correct, and our certified supremum confirms the rounding.
3. **Cook and York 1990.** Their Table I agrees with our certified masses at every printed digit except j = 100 (and
   j = 3, probe), where the last digit is one unit high, within their stated accuracy. An earlier lab log line claimed
   agreement at every digit; that was wrong.
4. **Dain and Gabach Clement 2009, eq. (13).** The printed prefactor 1/(4π) gives j²/80 at leading order; with
   M = μ + 2 lim r u the prefactor is 1/(2π) (a factor-2 discrepancy consistent with a misprint).
5. **Gabach Clement eq. (23)** prints V/r²; eq. (24) and Waxenegger-Beig-Ó Murchadha show it should read V⁴/r²
   (not load-bearing). Her condition (28) with 9J is sufficient; for A = P = 0 the exact super-solution condition is
   Q² ≥ 3|J|.
6. Dain-Gabach Clement's sentence "In the extreme case they do not contain any apparent horizons or minimal surfaces"
   (their Sec. 1) is about Kerr-Newman slices, not Bowen-York data; it must not be cited for Theorem 9.
7. **Our own corrections, applied here.** Upper ends in the P2a, P2b and P2c verdict tables were copied from an engine log
   line that truncates; Tables 2-3 are regenerated from the exact rationals. P2e display slips (a lower end
   0.92820045390 that should read 0.92820045389; a ledger row with both ends rounded the wrong way) are corrected to
   [0.9282004538, 0.9282004687]. A design file printed √r = ((1+ξ)/(1−ξ))⁴; the correct form, used by every instrument,
   is ((1+ξ)/(1−ξ))². A build note described the B3T width as "about 2η"; it equals the integrated residual.

## 11. Limitations (what we do not claim)

- Floating-point and library layer: python-flint (arb) ball arithmetic with outward rounding is trusted.
- E1 (existence for μ > 0) is cited; Cantor's 1979 paper, on which the cited proof rests, was not opened. The
  Riemannian Penrose inequality (Huisken-Ilmanen, Bray) was not opened, so we do not claim M_ADM ≥ M_AH.
- Standard analysis used without opening a textbook: the interior W^{2,p} (Calderón-Zygmund) estimate and the Morrey and
  Sobolev embeddings (Appendix B; Gabach Clement's weighted embedding via Dain-Gabach Clement Lemma A.1, which was opened),
  Newtonian-potential regularity, Jordan-Brouwer separation. The maximum principles are proved here (Appendix C).
- Conjecture E (the throat exponent) is unproved; c_BY ≈ 1.1991 and all flux digits beyond Theorem 6 are probe values.
- "For every J" in P1 holds at the slope-route accuracy (about 2e-11 relative); every other "for every j" rests on
  monotone or tangent-lemma cells and is a statement about bounds, not about point values.
- Theorem 2's rounding decision rests on a single top node; it is blind-replicated (REP6b) with a thin margin. The PDE
  constants (the 9/4 source and the exponent 7) are a common-mode input shared by all stages, guarded by the 1/40
  control, the symbolic check of the source term and the agreement with independent numerical literature.
- No claim that every χ value below the supremum is attained by our certificates (continuity of j ↦ δ is proved in
  Section 5 but was not used by the cell certificates).
- The final text had one internal referee round (REFEREE-RELEASE); novelty gates must be re-run within a week of release (arXiv,
  INSPIRE full text, Zenodo). Gate blind spots: pure-mathematics full texts are not searchable through INSPIRE, and
  numeric full-text searches miss known hits.
- Scope: one spinning puncture, conformally flat, maximal, no boost. Kerr is not conformally flat; boosted data break
  K(ν, ν) = 0.

## 12. Related work

- Data and puncture method: Bowen and York (1980, not opened); Brandt and Brügmann (1997); Cook and York (1990).
- Extreme limit: Dain and Gabach Clement (2009), Gabach Clement (2009/2012), Waxenegger, Beig and Ó Murchadha (2011),
  Hannam, Husa and Ó Murchadha (2009), Dain, Lousto and Zlochower (2008), Aceña and Gabach Clement (arXiv 2015).
- Numerical spin ceilings: Dain, Lousto and Takahashi (2002), Lovelace, Owen, Pfeiffer and Chu (2008), Lousto, Nakano,
  Zlochower, Mundim and Campanelli (2012).
- Perturbative masses: Gleiser, Nicasio, Price and Pullin (1997), eq. (19); Cabero and Krishnan (2014); Altas and Tekin
  (2021).
- Area-spin: Dain (2010); Dain and Reiris (2011); Jaramillo, Reiris and Dain (2011); Aceña, Dain and Gabach Clement (2011).
- Methods: Korevaar, Mazzeo, Pacard and Schoen (1999); Marques (2004); Corvino and Pollack (2011); Chruściel and Mazzeo
  (2012); Laetsch (1975); Berestycki and Nirenberg (1991); Andersson, Mars and Simon (2007).
- Nearest numerical neighbour: De Ceuster and Li (2026) solve the bare-mass tangent problem in floating point with a
  1e-10 residual gate; they use "certified" for that gate, not for a proof.

## 13. Open doors

What this work leaves open, and what we noticed and did not follow up. "Open according to" quotes the literature; "open as
far as our dated search found" reports our own search; "UNTESTED" marks an observation of ours that we did not test.

- **Horizons in perturbed extreme data.** Theorem 9 proves the first half of the expectation in [DGC, Sec. 6]; the second
  half, that "a small perturbation of an extreme solution (in particular, the extreme Bowen-York data) will always have
  an horizon", is left there as a conjecture to "prove or disprove". Along the family it holds (Lemma 8: for μ > 0 the
  throat is a closed minimal surface). [DGC10, Sec. 4] asks whether their compactly supported deformations of extreme
  Kerr data contain trapped surfaces. UNTESTED: the sliding argument of Section 8.2 does not use an explicit spacetime.
  *Status:* open according to Dain and Gabach Clement 2009 [DGC, Sec. 6] (no later answer found in our search of
  2026-11-17).
- **Conjecture E.** The throat exponent ν = √(λ_1 + 1/4) of Section 5 is a heuristic; the sealed test B2V agrees
  numerically (0.882059 vs 0.882019) and excludes the spherical-average value √2 − 1/2. A proof, and the relation to the
  fitted exponent of [LOPC, eq. (50)], are open in our work. *Status:* UNTESTED observation from this work (a heuristic
  with one sealed numerical check).
- **The area-spin constant c_BY.** Section 8.4 proves A_μ(μ/2) > A_0 and gives the probe c_BY = A_0/(8πJ) ≈ 1.1991;
  sharpness (A_μ(μ/2) → A_0) and a certified A_0 are not proved. A ≥ 8π|J| for the μ > 0 throats is [DR]; that
  Bowen-York horizons stay about 20% from the extreme-Kerr equality is numerical in [LOPC] and [HHO, Sec. VI].
  *Status:* open as far as our dated search found (2026-11-17; arXiv, INSPIRE, Zenodo).
- **Flux masses for other cylindrical ends.** The identity of Section 7 appears to extend to the maximally sliced spinning
  trumpets of [HHO, eqs. (20), (22)], giving M(C, S) = (1/(6π)) ∫ (4|∇D|² + D²) dΩ with the Schwarzschild trumpet as an
  exact control (M = m); boosted and dipole data are not autonomous along the throat. The two limit steps (Appendices A-B)
  have not been redone for these cases. *Status:* UNTESTED observation from this work.

## References

Only sources opened by the lab are listed with locators. Items marked [abstract] or [not opened] are cited for
attribution only.

- [AMS] L. Andersson, M. Mars, W. Simon, arXiv:0704.2889 (ATMP 12, 2008), Thm 7.1.
- [AT] E. Altas, B. Tekin, arXiv:2012.04450, eqs. (69)-(70).
- [AGC] A. Aceña, M. E. Gabach Clement, arXiv:1510.02079, Sec. 1, eq. (30).
- [ADG] A. Aceña, S. Dain, M. E. Gabach Clement, arXiv:1012.2413 (CQG 28, 105014), Sec. 2, Cor. 2.2.
- [BB] S. Brandt, B. Brügmann, Phys. Rev. Lett. 78, 3606 (1997), arXiv:gr-qc/9703066 v2, p. 2 right column (existence
  and uniqueness), p. 3 eq. (16).
- [BN] H. Berestycki, L. Nirenberg, Bol. Soc. Brasil. Mat. 22 (1991) 1-37, doi:10.1007/BF01244896 [not opened; as
  described in Riv].
- [Bray] H. L. Bray, arXiv:math/9911173 (Riemannian Penrose inequality) [abstract seen only as a paraphrase; cited only in Section 11].
- [Can] M. Cantor, J. Math. Phys. 20, 1741 (1979) [not opened; cited by DGC and BB as the basis of E1].
- [GT] D. Gilbarg, N. S. Trudinger, Elliptic Partial Differential Equations of Second Order, Classics in Mathematics
  (2001), doi:10.1007/978-3-642-61798-0, Ch. 3 [not opened; "see also" only, Appendix C].
- [Gon] Gonçalves, Zenodo record 10.5281/zenodo.23270506 (2026), computer-assisted Choptuik mode count [abstract].
- [HI] G. Huisken, T. Ilmanen (2001) [not opened; cited only in Section 11].
- [BSA] T. Benko, I. Stavrov Allen, arXiv:2103.11036, Sec. 1.2, Thm 1.11, Prop. 3.4.
- [BY] J. M. Bowen, J. W. York, Phys. Rev. D 21, 2047 (1980), doi:10.1103/PhysRevD.21.2047 [not opened].
- [CK] M. Cabero, B. Krishnan, arXiv:1407.7656, eqs. (36), (47).
- [CM] P. T. Chruściel, R. Mazzeo, arXiv:1201.4937, Sec. 3 after eq. (3.1).
- [CP] J. Corvino, D. Pollack, arXiv:1102.5050, Prop. 3.9 and the mass formula following it.
- [CY] G. B. Cook, J. W. York, Phys. Rev. D 41, 1077 (1990), doi:10.1103/PhysRevD.41.1077, p. 1078 col. 2, Sec. V
  (p. 1081), Table I (p. 1082).
- [Dain] S. Dain, arXiv:1008.0019 v2 (PRD 82, 104010), Sec. IV, Sec. VI eqs. (87)-(92).
- [DCL] F. De Ceuster, T. G.-F. Li, arXiv:2609.10173, Sec. IV C.
- [DGC] S. Dain, M. E. Gabach Clement, "Extreme Bowen-York initial data", CQG 26, 035020 (2009), arXiv:0806.2180 v1:
  pp. 4-5 (eqs. (10)-(13), positivity, J_ADM), Sec. 2 (the mass remark and the horizon question), Lemmas 3.1-3.2
  (p. 8), Thm 2.1 and eqs. (16), (19)-(20), (96), Sec. 6 (Final comments).
- [DGC10] S. Dain, M. E. Gabach Clement, arXiv:1001.0178 (CQG 27, 075003), Lemma A.1, eqs. (78)-(79), (82)-(89), Sec. 4
  (Final Comments).
- [DLT] S. Dain, C. O. Lousto, R. Takahashi, arXiv:gr-qc/0201062, eq. (5).
- [DLZ] S. Dain, C. O. Lousto, Y. Zlochower, arXiv:0803.0351 (PRD 78, 024039), Sec. II B, Sec. III.
- [DR] S. Dain, M. Reiris, arXiv:1102.5215 (PRL 107, 051101), abstract and eq. (10).
- [GC] M. E. Gabach Clement, arXiv:0911.0258 v2: Thm 2.1 (eq. (21)), Lemma 3.1 (eqs. (28)-(31)), eqs. (38), (43)-(44),
  Lemma 4.1 (eqs. (56), (59), (64)), Lemma 4.2 (eqs. (65)-(72)), Lemma 4.3, Sec. 5.
- [GJ] K. Goldstein, V. Jejjala, arXiv:2609.21797 v2, abstract and p. 1-2.
- [GNPP] R. J. Gleiser, C. O. Nicasio, R. H. Price, J. Pullin, arXiv:gr-qc/9710096, eqs. (16)-(17), (19).
- [HHO] M. Hannam, S. Husa, N. Ó Murchadha, arXiv:0908.1063 v4 (PRD 81, 064020), eqs. (20), (22), Sec. VI, eqs. (25)-(26),
  Fig. 8.
- [JRD] J. L. Jaramillo, M. Reiris, S. Dain, arXiv:1106.3743, abstract.
- [KMPS] N. Korevaar, R. Mazzeo, F. Pacard, R. Schoen, Invent. Math. 135 (1999) 233-272, doi:10.1007/s002220050285,
  arXiv:math/9807038: Prop. 1, the "radial Pohozaev invariant", Secs. 3-4, Prop. 5, Sec. 6.2.
- [Lae] T. Laetsch, J. Math. Anal. Appl. 51 (1975) 653-669, doi:10.1016/0022-247X(75)90116-X [abstract].
- [LNZMC] C. O. Lousto, H. Nakano, Y. Zlochower, B. C. Mundim, M. Campanelli, arXiv:1203.3223 v2, Table I.
- [LOPC] G. Lovelace, R. Owen, H. P. Pfeiffer, T. Chu, arXiv:0805.4192, eq. (5), Sec. II and Fig. 3, eqs. (49)-(50),
  (53)-(54).
- [Mar] F. C. Marques, arXiv:math/0408063, eq. (6.9) and Sec. 7.
- [Rei] M. Reiris, arXiv:1209.4530 (CQG 31, 025001), proof of Cor. 3.
- [RT] M. Reiterer, E. Trubowitz, arXiv:1203.3766, CMP 368, 143 (2019) [abstract].
- [Riv] T. Rivière, "Exploring the unknown: the work of Louis Nirenberg on partial differential equations",
  arXiv:1505.04930, Sec. VII.2.
- [YP] J. W. York, T. Piran, in Spacetime and Geometry (eds. Matzner, Shepley; Univ. of Texas, Austin, 1982),
  pp. 147-176 [not opened; as cited by LOPC ref. 48].
- [WBO] G. Waxenegger, R. Beig, N. Ó Murchadha, arXiv:1107.3083 v2, Thm 3.1.

---

## Appendix A. E3_0: the ADM mass of the extreme data, with derivative decay

**Statement.** Let u_0 = lim u_μ and g_0 = SΦ_0^-7. Then (i) g_0 is in L¹ and (1/(2π)) ∫ g_0 = M_0; (ii) u_0 = N[g_0];
(iii) r²∂_r u_0 → −M_0/2 and |∇_T u_0| = O(r^-3) uniformly in angle. Hence Φ_0 = 1 + M_0/(2r) + ρ with ρ = o(r^-1),
∂_r ρ = o(r^-2), ∂_θ ρ = o(r^-1/2), the flux ∫ r² ∂_r ρ dΩ → 0, and M_0 is the ADM mass of u_0.

**Proof.** (i) By A1 (Section 5), Φ_μ decreases to Φ_0 as μ decreases, so g_μ = SΦ_μ^-7 increases to g_0. By Lemma G,
(1/(2π)) ∫ g_μ = M(μ) − μ ≤ M(1) for μ ≤ 1, and M(μ) − μ → M_0 (Proposition 3). Monotone convergence gives the claim.
No lower bound r^-1/2 on Φ_0 is needed.

(ii) g_0 is smooth on R³ \ {0} (u_0 is C^∞ there), integrable, and g_0 ≤ S = O(r^-6) at infinity. Near 0, Φ_0 ≥ c r^-1/2
(the lower barrier of Dain-Gabach Clement, eqs. (16), (96)) gives g_0 ≤ C r^-5/2. So N[g_0] is C² on R³ \ {0}, solves
ΔN = −g_0 there, is O(r^-1/2) near 0 and tends to 0 at infinity. h = u_0 − N[g_0] is harmonic on R³ \ {0},
h = O(r^-1/2) = o(r^-1) at 0 (u_0 = O(r^-1/2) by the upper barrier), and h → 0 at infinity. By Bôcher's theorem the
singularity at 0 is removable (a harmonic function in a punctured ball that is o(r^-1) extends harmonically), and by
Liouville h ≡ 0. (Alternatively, u_μ = N[g_μ] and monotone convergence give u_0 = N[g_0] directly.)

(iii) ∇N[g](x) = −(1/(4π)) ∫ g(y)(x − y)/|x − y|³ dy. For |y| < |x|/2, |x|² x̂·(x − y)/|x − y|³ → 1 boundedly, so dominated
convergence gives |x|² x̂·∇N → −(1/(4π)) ∫ g_0 = −M_0/2. For |y| ≥ |x|/2, g_0(y) ≤ C|y|^-6 gives a contribution
O(|x|^-5) to |∇N|, hence O(|x|^-3) after multiplying by |x|². The tangential part is bounded by
(C/|x|³) ∫ g_0(y)|y| dy + O(|x|^-5), and ∫ g_0|y| < ∞. Finally r u_0 → M_0/2 uniformly (Proposition 3), which gives the
three statements on ρ. ∎

## Appendix B. C^1 convergence at the throat (W^{2,p} and Morrey)

**Statement.** Ω_0(σ, ·) → V in C^1(S²) uniformly as σ → −∞, together with ∂_σ Ω_0 → 0.

**Proof.** C^0 convergence is A2 (Section 8.2). Let W = Ω_0 − V. Subtracting the sphere equation from (E),
(∂_σσ + L_S − 1/4)W = F := −(9/4) sin²θ (Ω_0^-7 − V^-7). Choose σ_1 with Ω_0 ≥ c_-/2 for σ ≤ σ_1 (possible by A2 and
V ≥ c_-). By the mean value theorem, |F| ≤ 7 · (9/4) (c_-/2)^-8 |W| there. Let Q_k = [σ − k, σ + k] × S² with σ + 2 ≤ σ_1.
The operator is uniformly elliptic with smooth coefficients and invariant under σ-translation, so the interior
W^{2,p} estimate (Calderón-Zygmund, p > 3) holds with a constant C independent of σ:
||W||_{W^{2,p}(Q_1)} ≤ C (||W||_{L^p(Q_2)} + ||F||_{L^p(Q_2)}) ≤ C' sup_{Q_2} |W|. Morrey's embedding (p > 3 in three
dimensions) gives ||W||_{C^{1,α}(Q_1)} ≤ C'' sup_{Q_2} |W| → 0. Hence ∂_σ Ω_0 → 0 and ∇_S Ω_0 → ∇V uniformly, and the
integrand of E(σ) converges uniformly. A Schauder estimate alone would need a Hölder norm of F tending to 0, which is
not available a priori; the L^p route avoids that. ∎

## Appendix C. Maximum principles, proved here

Let Ω ⊂ R³ be open and connected, c ≥ 0 a locally bounded function on Ω, and L = Δ − c. All functions are C² on Ω.
No textbook page is load-bearing; see also Gilbarg and Trudinger, Ch. 3 [not opened].

**Lemma C.1 (weak maximum principle).** (a) If Ω is bounded, c is bounded, Lw ≥ 0 in Ω and w is continuous on the
closure, then sup_Ω w ≤ max(0, sup_{∂Ω} w). (b) In particular, if v satisfies Δv ≤ 0 on an open set U (possibly
unbounded) and liminf v ≥ 0 at every boundary point of U and as |x| → ∞, then v ≥ 0 on U.

*Proof.* (a) First let Lw > 0. If w had a positive maximum at an interior point p, then Δw(p) ≤ 0 and −c w(p) ≤ 0, so
Lw(p) ≤ 0, a contradiction. In general put w_η = w + η e^{αx_1} with α² > sup c; then
L e^{αx_1} = (α² − c) e^{αx_1} > 0, so Lw_η > 0 and sup w_η ≤ max(0, sup_{∂Ω} w_η); let η → 0. (b) If v(x_0) < −ε < 0,
the set {v < −ε/2} containing x_0 has compact closure inside U by the liminf hypotheses; on it Δv ≤ 0 and v = −ε/2 on
its boundary, so (a) with c = 0 applied to −v gives v ≥ −ε/2, a contradiction. ∎

**Lemma C.2 (Hopf boundary lemma, the case used).** Let B = B_R(y) with closure in Ω, v ≥ 0 in Ω, Lv ≤ 0, v > 0 in B,
and v(x_0) = 0 for some x_0 on ∂B. Then ∂v/∂n(x_0) < 0, with n the outward unit normal of B at x_0.

*Proof.* Let ρ = |x − y| and h = e^{−αρ²} − e^{−αR²}. In three dimensions Δe^{−αρ²} = (4α²ρ² − 6α) e^{−αρ²}, and
−ch ≥ −c e^{−αρ²} since c ≥ 0. On the annulus A = {R/2 < ρ < R}, Lh ≥ e^{−αρ²}(α² R² − 6α − sup_A c) > 0 for α large. On
the inner sphere ρ = R/2, v ≥ m > 0 (a compact subset of B); on the outer sphere h = 0. With ε = m/max h, the function
z = v − εh satisfies z ≥ 0 on ∂A and Lz = Lv − εLh < 0 in A. If z had a negative minimum at an interior point p, then
Δz(p) ≥ 0 and −c z(p) ≥ 0, so Lz(p) ≥ 0, a contradiction. So z ≥ 0 in A, and z(x_0) = 0. Hence ∂z/∂n(x_0) ≤ 0, that is
∂v/∂n(x_0) ≤ ε ∂h/∂n(x_0) = −2εαR e^{−αR²} < 0. ∎

**Lemma C.3 (strong maximum principle).** If v ≥ 0 and Lv ≤ 0 in Ω and v(p) = 0 at some p in Ω, then v ≡ 0.

*Proof.* Let Z = {v = 0}, closed in Ω and nonempty. Suppose Z ≠ Ω. Since Ω is connected, Ω \ Z is open and not closed
in Ω, so there is a point y in Ω \ Z with dist(y, Z) < dist(y, ∂Ω). Take the largest ball B_R(y) contained in Ω \ Z;
its boundary meets Z at a point x_0, and its closure lies in Ω. Replacing it by a smaller ball inside it that touches
∂B_R(y) only at x_0, we may assume v > 0 in B. At x_0, v attains its minimum 0 over Ω, so ∇v(x_0) = 0. Lemma C.2 gives
∂v/∂n(x_0) < 0, a contradiction. ∎

**Transfer to the cylinder.** If Y solves (∂_σσ + L_S − 1/4 − V)Y = 0 on R × S² with V ≥ 0 locally bounded, then
y = r^-1/2 Y solves (Δ − V r^-2) y = 0 on R³ \ {0} (by Δ(r^-1/2 X) = r^-5/2 (X_σσ − X/4 + L_S X)), with c = V r^-2 ≥ 0
locally bounded, and y has the same sign and zeros as Y. So Lemma C.3 applies to the cylinder uses in Section 8.

**Uses.** Lemma C.1 in A1, Theorem 4 (a)-(c), B3-5 and Appendix D; Lemma C.3 in Theorem 4 (strictness, z > 0), Lemma 8
and both proofs of Theorem 9. In each use c is locally bounded: c ≤ 7SΦ^-8 ≤ C min(r², r^-6) for μ > 0, c = V r^-2 on the
cylinder, Q = 7 (9/4) J² sin²θ r^-6 Φ^-8 with Φ ≥ 1 in the dilation proof, and c = h/(4r²) with h bounded in Appendix D.

## Appendix D. The barrier repair of Gabach Clement (71)-(72)

**Setting.** D := u − V/(√r (1 + b√r)) − U solves (Δ − c)D = 0 on R³ \ {0} with c = h/(4r²) ≥ 0, D is smooth off 0,
D = o(r^δ) at 0 for some δ > −1 (from the weighted-space memberships via Dain-Gabach Clement Lemma A.1), and D → 0 at
infinity. h is bounded, which needs √r u ≥ c_0 > 0 near 0; for u_0 this follows from u_0 ≥ u_0^- with
√r u_0^- → V^- > 0 (Gabach Clement eqs. (38), (43)-(44)).

**Claim.** D ≡ 0.

**Proof.** Fix ε > 0. Since Δ(1/r) = 0 off 0, L(D − ε/r) = cε/r ≥ 0 with L = Δ − c. Near 0, D − ε/r → −∞ because
D = o(r^δ) and r^δ = o(r^-1). On the annulus {ρ < r < R}, Lemma C.1 (a) gives
sup (D − ε/r) ≤ max(0, sup over |x| = ρ, sup over |x| = R of (D − ε/r)). For ρ small the inner value is negative; the
outer value is at most sup over |x| = R of |D|, which tends to 0 as R → ∞. Hence D ≤ ε/r on R³ \ {0}, and letting
ε → 0, D ≤ 0. The same argument for −D gives D ≥ 0. ∎

## Appendix E. Referee wording fixes applied in this text

- J ≠ 0 is a hypothesis of Theorem 4 (at J = 0, M = μ).
- "The tangent-lemma slope 2L is exactly M'" is false; the text says: the exact tangent slope 1 + 2 lim r w equals M',
  and the series identity 2L = 1 − (3/40) j² is dM/dμ at that order.
- Monotonicity, slopes in [0, 1] and continuity at μ = 0 are recorded as corollaries of Dain-Gabach Clement Lemmas
  3.1-3.2, with credit; never "an open point closed".
- The uniqueness of u_0 is cited (Gabach Clement Lemma 4.3; Waxenegger-Beig-Ó Murchadha Thm 3.1), not called open.
- B3-5 has no r^-1/2 hypothesis; integrability is automatic.
- Conjecture E: the word "exact" is removed; the throat-profile input is cited, not assumed.
- Supremum, not maximum or attainment, everywhere; "computer-assisted" on every certified statement; the family
  qualifiers (single, non-boosted, conformally flat, maximal, puncture boundary conditions, J = J_ADM, μ > 0) on every
  family statement; μ = 0 excluded from the P2c-P2e statements.
- "Certified spins from 0 to 0.8289" (P2a) replaced by enclosures for every j; "the largest certified χ" replaced by
  the certified lower bound; window widths reported from the runs, never "within 1e-6 of the true supremum".
- "P2e decides the sixth digit" names the family and says that the bar lies 3.13e-8 above c and 4.61e-8 above X3.
- The small-spin anchor 3/40 is three times 1/40 for any exponent (a scaling identity), so it is a consistency check of
  the tangent problem, not an independent guard of the constant 7.
- Theorem 7: "new identity" replaced by the Korevaar-Mazzeo-Pacard-Schoen credit; the spherical-model control is an
  algebra check that cannot fail; the 180 digits are probe digits.
- Lemma 8: "outermost MOTS" strengthened to "only closed connected MOTS"; the mean-curvature formula has the power r^-3/2.
- Theorem 9: "new trick" replaced by the Berestycki-Nirenberg sliding method; the touching step credited to Benko-Stavrov
  Allen and Reiris; "horizon infinitely far away" marked as a heuristic picture; the extreme-data probe that failed
  is described as pre-stated and unsealed.
- c_BY is a probe value; A ≥ 8π|J| is credited to Dain and Reiris.
- Lean covers only the monotone-cell inequality.
- The convexity proof says that pointwise comparison of the integrands does not decide the sign.

---

*Prepared with AI assistance; what was checked, and by what, is stated above.*
