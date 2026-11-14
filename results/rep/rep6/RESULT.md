# REP6 (BHSLICE blind replication, seal 0405395ed8): result

## Disclosures (brief hygiene, read first)
- The brief gave one target number: the K1 leading coefficient 1/40. I derived it myself before using it (see Controls).
- The brief also fixed the normalisations (M_ADM = mu + 2 lim r u, j = 4J/mu^2). These define the object. I checked them myself (derivation below).
- From general literature memory, not from any file, I knew that Bowen-York J/M^2 tends to about 0.93 at large spin. I did not use this anywhere in the computation.
- Forbidden paths: none opened. I read DOCTRINE.md and LEMMA-REP4.md. I did not open lichnerowicz_axisym.py; it was allowed but not needed.
- Map family: my radial variable is xi = tanh(log(2r)/8) = ((2r)^(1/4) - 1)/((2r)^(1/4) + 1). It is a log/tanh compactification, and it is also a quarter-power map, so it belongs to the same power-map family as the sqrt map (exponent 1/4 instead of 1/2). The rest of the formulation is different: the sigma = log(2r) equation, the use of inversion symmetry, the 1/(2r+1) prefactor, the Legendre angular basis, and the positivity method.

## Method
1. Formulation, derived and checked in ctrl_formula.py, which agrees with brute-force mpmath Laplacians to 1e-40. Put mu = 1, sigma = log(2r), u = e^{-sigma/2} W, Omega = 2cosh(sigma/2) + W. Then the equation becomes
   W_ss - W/4 + Lap_S W = -(9/4) j^2 sin^2(theta) Omega^-7.
   This has constant coefficients in sigma and is invariant under sigma -> -sigma (inversion symmetry), so W is even.
   - Ansatz: u = 2G/(2r+1) with G = g0 + w^4 H(xi, c), where w = 1 - xi^2, H uses even Chebyshev polynomials T_2n(xi) times even Legendre polynomials P_2k(c), and c = cos(theta).
   - Read-off: delta = M_ADM - 1 = 2 g0.
   - Normalised residual: Rhat = r^2 (Delta u + S Phi^-7) / (e^{-sigma/2} sech^3(sigma/2)). It is polynomial in (xi, c) plus CS (1-c^2) A^10 w^8 / E^7, with E = A^2 + 2 w^4 G, A = 2(1 + 6xi^2 + xi^4) and CS = (9/32) j^2.
2. Sandwich. Shifting G -> G + eta adds -eta/2 to the linear part of Rhat. The source term is monotone in G. Hence:
   - u_+ = 2(G_N + eta_p)/(2r+1) is a strict super-solution if eta_p > 2 sup Rhat_N;
   - u_- = 2(G_N - eta_m)/(2r+1) is a strict sub-solution if eta_m > -2 inf Rhat_N and E_- > 0.
   - Both lie in class K (LEMMA-REP4 sec. 2), structurally: w^4 = 512 r/(1 + (2r)^{1/4})^8, so u_+- - u_+-(0) = r * (a smooth function of ((2r)^{1/4}, c)). That gives continuity at 0 (u(0) = 2(g0 +- eta)), |grad| = O(1), |Lap| = O(1/r) and 1/r decay.
   - Lemma C-REP4 then gives delta in [2(g0 - eta_m), 2(g0 + eta_p)].
3. Tangent. psi_+ = 2Z/(2r+1) with Z = Z_D + zeta0 + kappa + w^4 H_e. Here 2 Z_D/(2r+1) = -r d_r u_N is the dilation part, and H_e uses the full Chebyshev basis (Ne = 120, Ke = 11).
   - Because the sigma operator has constant coefficients, the residual of Z_D is (1+3t)/2 Lin_N - D Lin_N. It is computed with dual numbers, with no 0/0 at the ends (ctrl_tangent.py: agrees with brute force to 1e-40).
   - TL-H is equivalent to That <= 0, where That is normalised as above and V1 is built from MY certified u_-.
   - Adding kappa lowers That by at least kappa/2, so kappa > 2 sup That_N.
   - 2L_hi = 1 + 2(g0 + zeta0 + kappa).
   - Class K holds structurally: 2r/(2r+1) is smooth in r, and every other non-constant piece has the form w^4 times a polynomial.
4. Positivity certificate (my own; not a box-polynomial bound). The domain is covered by dyadic patches; for each patch:
   - Sample f = Rhat (or That, or E) exactly in arb at 65 x 65 Chebyshev-Lobatto points and form the tensor interpolant p. Then sup p <= c00 + sum'|c_jk|.
   - Error bound: on the Bernstein polyellipse E_rho1 x E_rho2 with |f| <= M, the Chebyshev coefficients satisfy |a_jk| <= 4 M rho1^-j rho2^-k (Laurent/Cauchy bound on the Joukowski torus; derived, not looked up). Lobatto aliasing maps T_j to T_j' exactly, so ||f - p|| <= 8M(s1 s2 - t1 t2), with s = rho/(rho-1) and t = s(1 - rho^-(n+1)).
   - M comes from acb ball evaluation on balls covering the distinguished boundary, using the conj symmetry and the maximum modulus principle.
   - Analyticity: Re E > 0 and Re A > 0 on that boundary. Re E is pluriharmonic, so by the minimum principle it stays positive inside, hence E != 0 and A != 0 there.
   - An exact Fraction check confirms that the patches tile the domain.
   - Symmetry reduces the main case to [0,1]^2. The tangent case uses xi in [-1,1] and c in [0,1].

## Results (j = 2^21, mu1 = 1; copied from result_rep6.json)
- (a) delta = M_ADM - 1 in [750.5607142022, 750.5607145913] (exact dyadics in JSON; M_ADM in [751.5607142022, 751.5607145913]).
  - Ingredients: sup Rhat_N <= 4.96892e-8, inf >= -4.75791e-8, inf E_- >= 90.826342, 16 patches, candidate N = 84, K = 11.
- (b) 2L_hi = 0.0140665057657, rounded up (exact 129740615934886563/9223372036854775808).
  - Ingredients: sup That_N <= 0.00351532, kappa = 129692229038850723/18446744073709551616, 32 patches.
  - TL-H is certified against my u_- (Phi_1 > 0 via E_- >= 90.826342). Class K is structural (above).
- (c) chi_tail = (j/4)/(1 + delta_lo - max(2L_hi,0))^2 = 0.928235200531265, rounded up (exact fraction in JSON).

## Controls
- C0 and C0T: the residual formulas agree with brute-force mpmath Laplacians to 1e-40 (main and tangent).
- K1: linear theory gives lim r u = (3/2) J^2 * 128 int r^3/(2r+1)^7 dr = J^2/5, via a Beta integral 8 B(4,3) = 2/15. So delta = 2J^2/5 = j^2/40.
  - At j = 2^-20 the certified delta/j^2 is in [0.0249998807907073, 0.0250001192092865], which contains 1/40.
  - At j = 2^-6 the enclosure [0.0249996925180803, 0.0249996925190627] EXCLUDES 1/40. This is expected: it resolves the O(j^2) correction, about -1.26e-3 j^2, and shows the control could fail.
- KT: psi = 0 gives That = -forcing <= 0 analytically, so 2L = 1. I started an interpolation run of this case and stopped it at 597 patches (the bound is loose where That touches 0 at c = 1). It is not needed.
- KN (source 9/4 * 1.001): delta in [750.748516021, 750.748552358]. The shift is +0.1878, which matches M proportional to fac^{1/4} at large j (751.56 * 2.498e-4 = 0.1878).
- Toy (J = 0): G = 0 exactly; delta in [-2^-59, 2^-59]; 2L_hi = 1.
- Nulls on one patch: candidate with g0 + 1e-3 gives sup/inf 5.6e-2 / -0.33, and CS*1.001 gives 34.8 / -10.8. In both cases the certificate rejects at the production tolerance.

## Convergence
- Second truncation (N = 72, K = 11): delta in [750.5606962809, 750.5607325779]. It contains the N = 84 enclosure.
- 384-bit precision: identical delta enclosure to all printed digits.
- 48 samples instead of 64: [750.5607137287, 750.5607150557], which contains the main enclosure.
- Tangent, second candidate (mainB, Ne = 100): 2L_hi = 0.198916 (looser, consistent).
- Tangent at 384 bits: 2L_hi = 0.0140665057657 (it differs only near the 1e-15 level).
- Float residuals fall geometrically with N: max|Rhat| = 4.3, 0.14, 2.8e-3, 1.05e-4, 1.4e-6, 1.5e-8 for N = 30, 40, 50, 60, 72, 84.
- Digits I trust: all printed digits of the enclosures are rigorous. delta itself is pinned to about 9-10 significant digits (750.56071440 +- 2e-7).

## Times and memory
- Float solves: 37 s (N = 84), about 26 s (N = 72), about 27 s (KN).
- Main certificate: about 4 s (16 patches). Tangent certificate: about 20 s.
- KT partial run: about 10 min (stopped).
- Total compute about 25 min.
- Peak memory 6.35 GB, in one early float probe (a dense check grid). All production runs were at or below 0.82 GB.
- Resume path: exercised (the KT run resumed from its JSON checkpoint).

## Weaknesses
- The quarter-power map is in the same power family as the sqrt map (disclosed above).
- 2L_hi = 0.0141 is loose: the float tangent suggests 2L of about 5e-6. It is limited by the odd residual that comes from D(Rhat_N) and from using u_- rather than u_N in V1.
- E1-E3 and the lemmas are taken from LEMMA-REP4.md; I did not re-derive them. The polyellipse/Cauchy and aliasing bounds were derived by me and are not refereed.
- arb/python-flint ball arithmetic is trusted.
- The Lin/source formulas are verified numerically at sample points, not symbolically.
- The analytic KT control was not completed by the machinery.
