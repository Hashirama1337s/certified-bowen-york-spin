# REP7 blind replication: extreme Bowen-York spin, rigorous M_0 and chi(0)

Seat REP7, 2026-11-14. Lead seal 133b99d241 (its BAR was not read).

## Disclosures (read first)
- The brief gave no target number. It did name the first seat's formulation (w = r^-1/2 G/(1+sqrt r), xi = tanh(log r/8)); I used a different one.
- The RN control's exact value M = 1 is part of how the control is defined. I derived the source myself (below).
- From prior reading I remembered that the extreme Bowen-York spin ratio is "about 0.93". No digits of it were used. Every number below comes from free computation.
- I opened no forbidden path. No web search. No identification attempted (none was asked for).

## Object
Delta u = -S (1+u)^-7 with S = (9/4) sin^2(theta) / r^6. u -> 0 at infinity, u ~ r^-1/2 V(theta) at r = 0.
M_0 = 2 lim r u and chi(0) = 1/M_0^2.

## Formulation (REP7's own, derived here)
- psi = 1 + u = sqrt(1 + F/r). Equivalently g := r psi^2 = e^-t + F, with t = -log r and eta = cos(2 theta).
- Cylinder form. With phi = r^1/2 psi:  Delta psi = r^-5/2 (phi_tt + L phi - phi/4), and S psi^-p = K r^-5/2 phi^-p.
  - BY: p = 7, K = (9/4) sin^2 = (9/8)(1 - eta).
  - RN: p = 3, K = 1/4.
- Substituting phi = sqrt(g) and multiplying by 2 g^3/2 gives
  E(F) = e^-t (F_tt + F_t + LF) + F (F_tt + LF) - (F_t^2 + 4(1-eta^2) F_eta^2)/2 - F^2/2 + 2 K g^((3-p)/2) = 0,
  with L f = 4(1-eta^2) f_etaeta - 2(1+3 eta) f_eta.
- Residual: R := Delta u_a + S (1+u_a)^-p = r^-5/2 E / (2 g^3/2).
  - Far zone, x = 1/r: R = x^4 Q with Q = (xE) / (2 x (1+xF)^3/2).
  - F tends to V^2 at the throat and to M_0 at infinity, so the mass is read off directly as M = F(x=0).
- RN control derivation. psi = sqrt(1 + M/r) gives Delta psi = -(M^2/(4 r^4)) psi^-3. So the source 1/(4 r^4) has exact solution u = sqrt(1+1/r) - 1, with M = 1 (that is, F = 1).
- Algebra check (formula_check.py). The coded r^-5/2 Phi and x^4 Q were compared with a direct spherical Laplacian of psi (mpmath numerical differentiation, 40 digits), for an arbitrary test F at 4 points per case. Worst relative difference: 4.1e-40, PASS. The first attempt "failed" because of a sign slip in the test itself (t = log x); the code was right.

## Numerics
- One far element in x = 1/r on [0, x1]. It uses Gauss-Radau nodes, so x = 0 is never a node; the regular solution is picked out because the unknown is a polynomial in x.
- 12 Lobatto elements in t, up to T+ = 26 (BY2: 30). Chebyshev-Lobatto grid in eta.
- C^1 matching at element interfaces, Neumann condition at T+, Newton iteration in float64 (about 2 s).

## Certificate
1. u_a is defined as the exact interpolant of the float nodal values.
   - It is continuous across elements by construction, because the shared nodes sit at s = +-1 exactly.
   - The far polynomial is corrected by ((1-s)/2)^2 times its non-constant part at x = 0. This keeps the value and derivative at the interface, makes F_a(0, eta) exactly constant, and makes Ehat(0, .) = 0 exactly.
2. Cells: 4 x 4 per element (BY2: 5 x 3). On each cell, sup|Phi| <= sum|c_kl| + 8 M (t1 S2 + S1 t2), where:
   - c_kl are the first-kind Chebyshev interpolation coefficients at 24^2 nodes (BY2: 28^2), from ball evaluations.
   - The remainder uses the Bernstein polyellipse rho = 4 (BY2: 3.5). M is the max of |Phi| over ball covers of the distinguished boundary |w1| = rho1, |w2| = rho2 (half of the w1 circle, by conjugate symmetry).
   - Analyticity: Re g > 0 (far zone: Re(1+xF) > 0) is checked on every box. Re g is pluriharmonic, so this certifies it on the whole closed polyellipse.
   - Derivations, done by me in the notes and not taken from an opened source:
     - |a_kl| <= 4 M rho1^-k rho2^-l, from a Cauchy estimate on the Joukowski torus;
     - the first-kind aliasing T_{2m(n)+-k} -> (-1)^m T_k, which gives |f - p| <= 2 sum_tail |a_kl|.
3. Kinks at the interfaces (and at T+) are surface measures with density |[d_r psi]| = |[F_t]| / (2 r^3/2 sqrt g). They are bounded by the sum of |eta-coefficients|.
4. Throat tail, t > T+: F = F_T(eta), so |Phi| <= A + B e^-t. A comes from the 1-D version of step 2; B is a crude bound. This gives a contribution 2A sqrt(r+) + (2/3) B r+^3/2.
5. rho(r) >= sup_theta |R|, as a radial measure (piecewise r^-5/2 B in t-cells, r^-4 B in far cells, plus shells). h is its Newton potential, h(r) <= delta/r, and delta = int rho s^2 ds.
   - w = u_a + h is a super-solution (h >= 0, so S(1+w)^-7 <= S(1+u_a)^-7).
   - v = u_a - h is a sub-solution.
   - 1 + v > 0 is checked: three inequalities, all True.
6. Comparison, with hypotheses: u_0 is C^2 on R^3 \ {0}, u_0 -> 0 at infinity, and u_0 = O(r^-1/2) (o(1/r) suffices), as stated in the object.
   - Super side: d = u_0 - w. On the open set {d > 0}, d has nonnegative distributional Laplacian, so it is subharmonic. Also d+ <= u_0 + 1 = O(r^-1/2).
   - Then d - eps/r is subharmonic on {d > 0}, negative near 0, and has limsup <= 0 at infinity. So d <= eps/r for every eps, hence d <= 0.
   - Sub side: the same argument with d = v - u_0. Here u_0 > 0 by the minimum principle and v = O(r^-1/2).
   - Therefore M_a - 2 delta <= M_0 <= M_a + 2 delta.
   - This needs C^0 plus measure Laplacians, not C^2; the only interior singularities are the kink shells, and they are accounted for in step 3.

## Results (outward-rounded decimals, from exact binary endpoints; result_rep7.json)
| case | M window | width |
|---|---|---|
| BY (primary) | [1.037956398789003, 1.037956399038277] | 2.49e-10 |
| BY2 (independent grid, map and verification parameters) | [1.037956398866011, 1.037956398961277] | 9.5e-11 |
| RN control | [0.999999999936687, 1.000000000063303] (contains 1) | 1.27e-10 |
| KN control (9/4 -> 9/4 (1+1e-3)) | [1.038215790677670, 1.038215790845716] | 1.68e-10 |
| BYP (BY spoiled by +1e-7 at every node) | [1.037954934262893, 1.037958063564387] | 3.1e-6 |

- chi(0) from BY: [0.928200458542311, 0.928200458988139]. From BY2: [0.928200458680027, 0.928200458850410].
- Width target of 1e-8: met (BY 2.5e-10).
- Trusted digits:
  - M_0 = 1.0379563989 (BY2 lies inside BY, and its window is 1.0379563988660 to 1.0379563989613);
  - chi(0) = 0.92820045876 +- 9e-11.
- Convergence: two discretisations whose rigorous windows nest. The numerical M_a values agree to 2e-14 (N = 16, 20, 24 probes: 1.0379563989136602 / 1.0379563989136400 / 1.0379563989136529; BY vs BY2 M_a differ by 3.8e-15).

## Controls
- RN: started away from the exact solution (F0 = 1.3 + 0.2 eta). Its window contains 1.
- KN: the scaling r -> b r with b^4 = J^2 predicts M_KN = 1.001^(1/4) M_0. The BY window scaled this way, [1.038215790637013, 1.038215790886348], overlaps the KN window. KN is disjoint from BY, and KN chi = sqrt(1.001)/M^2 lies in [0.928200458614963, 0.928200458915442], as scaling predicts.
  - The KN grid is not shifted with the scaling, so this is a genuine discretisation-independent check.
- BYP: the certificate widens to delta = 7.8e-7 and still contains the BY window.

## Resources
- Each case takes 15-24 s wall (Newton about 2 s; 208 cells take about 10 s in python-flint at 128 bits).
- Peak memory 2.41 GB under the 8 GB memcap. Total compute, probes included, is under 4 minutes.
- Instrument: rep7.py, sha256 c473c57be1f92f963af66b4d27e236d1c3f86d9541d324c1b36af41eaa019c32. All five final cases ran on this hash.

## Weaknesses
- The polyellipse/aliasing bound and the "continuous with nonnegative distributional Laplacian => subharmonic" lemma are standard facts that I derived or recalled. I did not transcribe them from an opened source. No Lean check, no referee yet.
- Hypotheses on u_0 (existence, regularity, O(r^-1/2) throat growth) are taken from the object definition and are not proved here.
- Trust in the arb ball library (python-flint 0.9.0) is assumed.
- The arc-ball centres are floats; they are covered by an explicit radius margin.
- The comparison argument relies on M_0 being defined as a limit along rays; the bounds hold for both liminf and limsup.
