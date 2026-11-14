# REP5 - blind replication seat, BHSLICE slope (lead seal 9e0df55e7a)

## Disclosures first

- What the brief handed me: the equations of the object, including the tangent equation and the normalisation
  s = 1 + 2 lim r w. It gave no target numbers. I derived both again (section 2).
- Accidental exposure. My first `py` run started in the session <SCRATCH>. That folder is shared with other
  sessions and holds a stray `inspect.py`, which Python imported in place of the standard module. It printed
  summaries of BHSLICE result_p1.json, result_p2a.json, result_p2b.json, result_p2c.json and ceiling_p2c.json.
  - None of these files is on my forbidden list. None of them holds a slope value.
  - The printout did include a key named `delta4_contains_m133_105600`, which gives away the j^4 coefficient of
    delta. I saw it before I derived that coefficient. My derivation is exact sympy algebra with a check of the
    residual at each order, and it gives -133/105600.
  - The printout also showed P2a/P2c delta and chi enclosures near j = 1 and the chi ceiling 0.9622.
  - I used none of these numbers. After that I ran nothing from the <SCRATCH>.
- The gravity_engine `__init__` docstring says the first seat's tangent certificate lives in bowen_york_sqrt.py. It
  contains no numbers. I never opened that file, its tests (including test_bowen_york_tangent.py) or any other fenced
  path. I avoided a listing of <LAB>/BHSLICE; the one attempt was refused by the permission system before it ran.
- How independent this is. I reused the engine's coordinate s = a/(r+a), its ansatz and its exact Bernstein
  routine. The residual-bound code and the whole tangent certificate are my own. The secant/convexity check (4c)
  is a separate route that never touches the tangent problem.

## 1. Results

All values are rigorous and rounded outward. The exact rationals are in result_rep5.json.

| j | quantity | lower | upper |
|---|---|---|---|
| 64 | delta = M_ADM - 1 | 3.17060411682713773929617041116 | 3.17060411682713773929617041186 |
| 64 | slope s(j) | 0.0520407975655210747020518703815 | 0.0520407975655210747020518836049 |
| 4096 | delta | 32.2150822519650947360600776016 | 32.2150822519650947655917481627 |
| 4096 | slope s(j) | 0.00131983286607973853976366731638 | 0.00131983286607990373316029655626 |
| 64 | KN control (exponent 6 in the tangent operator) | - | 0.0712192164366958717014770166549 |
| 0 | toy, J = 0 | 1 (exact) | 1 (exact) |

The quantity the brief literally asked for (a super-solution with V = 7 S Phi_-^-8) is NOT an upper bound on s;
see 3. Its values are 0.0520407975655210747020518835376 at j = 64 and 0.00131983286607990368908184511072 at
j = 4096. Both lie below my true upper bounds, as the theory predicts.

## 2. Derivations

### Coordinate and divided operator

Take s = a/(r+a) with a = mu/2 = 1/2. Then r = (1-s)/(2s) and 1/(2r) = s/(1-s), so

    Phi = D/(1-s),   D = 1 + (1-s) u.

For a function f(r) times P_l, the radial chain rule (dr/ds = -a/s^2) gives

    Delta (f P_l) = (4 s^4/(1-s)) [ (1-s) f'' - 2 f' - l(l+1) f/(s^2 (1-s)) ] P_l.

Checks of this formula:
- Delta(a/(r+a)) = -2a^2/(r (r+a)^3) holds exactly at rational r.
- My Chebyshev-coefficient operator equals the engine's monomial Rop exactly (control C1).

### Source terms and the two equations

The source pieces are

    S = 144 J^2 (1-c^2) s^6/(1-s)^6
    S Phi^-7 = 144 J^2 (1-c^2) s^6 (1-s) D^-7
    S Phi^-8 = 144 J^2 (1-c^2) s^6 (1-s)^2 D^-8

With J = j/4 this gives

    point:   Lh u + tau src D^-7 = 0,                  tau src = (9/4) j^2 (1-c^2) s^2 (1-s)^2
    tangent: Lh w = 7 tau src (s + (1-s) w) D^-8

To get the tangent equation, differentiate Delta u = -S Phi^-7 in mu with J fixed. Since dPhi/dmu = 1/(2r) + w,
this gives Delta w = 7 S Phi^-8 (1/(2r) + w). The mass term is 2 lim r (s p(s)) = p(0), so delta = p_u(0) and
s = 1 + p_w(0).

### Scaling

M(mu, J) = mu G(4J/mu^2) with G = 1 + delta. Differentiating at fixed J:

    s(j) = 1 + delta(j) - 2 j delta'(j)

### Analytic bounds 0 <= s <= 1

Write V = 7 S Phi^-8 >= 0.
- (Delta - V) w = V/(2r) >= 0 with w -> 0, so w <= 0, which gives s <= 1.
- z = 1/(2r) + w solves Delta z = V z away from the origin, with z -> +inf at the origin and z -> 0 at infinity, so
  z >= 0. Since lim r z = s/2, this gives s >= 0.

### Small-j series (series.py, exact sympy)

Write u = sum_k j^(2k) U_k. Then

    Lh U_1 = -src
    Lh U_2 = 7 src (1-s) U_1
    Lh U_3 = src [7 (1-s) U_2 - 28 (1-s)^2 U_1^2]

Each Legendre mode has a unique polynomial regular solution (no log terms appear up to this order). The program
asserts the residual is exactly 0 at every order. The coefficients are

    delta = j^2/40 - (133/105600) j^4 + (8923/84864000) j^6 + O(j^8)
    s     = 1 - (3/40) j^2 + 7 (133/105600) j^4 - 11 (8923/84864000) j^6 + ...

## 3. Certificates, and the direction issue in the brief

### Residual bound

The brief asks for a super-solution of the tangent problem with V = 7 S Phi_-^-8. Before the direction issue, here
is how every residual is bounded. Each residual is an exact polynomial in (s, t = c^2):

    point:   Q = D^7 Lh u_N + tau src
    tangent: G = D^8 Lh psi - 7 tau src (s + (1-s) psi)

- The arb DCT gives an enclosure of every Chebyshev coefficient. It uses a Lobatto grid whose size is at least the
  degree bound worked out from the construction.
- The sum of their absolute values bounds sup|residual| on the closed square.
- The shift is e > S/2, using Lh(s) = -2.
- Exact Bernstein (the engine routine) certifies u_-/s > 0, which gives D_+-, D_N >= 1. It also certifies
  Z/s = 1 + (1-s) psi/s > 0, i.e. 1/(2r) + psi >= 0.

### Which way each inequality runs

For an upper bound on the TRUE slope I need Delta psi <= V (1/(2r) + psi) with the TRUE V.
- psi >= w then follows from the maximum principle for Delta - V, V >= 0: (Delta - V)(psi - w) <= 0 with decay.
- Because 1/(2r) + psi >= 0, the worst case over V in [V_lo, V_hi] is the SMALLEST V. That is
  V_lo = 7 S Phi_+^-8, built from the point SUPER-solution.
- So my upper bound fits psi against D_+ and certifies G_+[psi_N + e s] <= S - 2e < 0 and Z >= 0.

The brief's V_hi = 7 S Phi_-^-8 is the wrong direction for an upper bound.
- A super-solution with respect to V_hi bounds w_hi from above, where Delta w_hi = V_hi (1/(2r) + w_hi).
- But (Delta - V_hi)(w_hi - w) = (V_hi - V) z >= 0, so w_hi <= w.
- So the brief's object bounds a quantity that sits below s. It is not a bound on s.
- Used correctly, Phi_- gives the LOWER bound: a sub-solution phi_N - e s with respect to V_hi, plus Z >= 0, gives
  Delta phi >= V_hi z_phi >= V z_phi, hence phi <= w.

## 4. Controls

All controls passed. They are in controls_rep5.json and k1t_rep5.json.

**C0 - DCT routine.** It recovers known 2-D Chebyshev coefficients exactly; the bound came out 381.834... against
an exact value of 381.834.... An aliasing null, a polynomial of degree above the grid, is detected, so this control
can fail.

**C1 - operator algebra.** My Chebyshev-coefficient Rop, s^2(1-s) multiplication and d/ds agree exactly with the
engine's monomial forms. Lh(s) = -2 and the Laplacian identity both hold.

**C2 - positive control at j = 2.** I reproduced the engine's certified point enclosure at j = 2 with my own
residual code:
- mine: [0.0846851040396333398063188927148, 0.0846851040396333402603260100475]
- engine: [0.08468510403963333980631889, 0.08468510403963334026032602]

These are the same interval: my e is exactly tau times the engine's eta, so the two residual-bound algorithms
compute the same bound.

**T - toy, J = 0.** The point residual is S = 0 exactly. The tangent fit is all zeros, so the slope enclosure is
[1, 1] exactly, a non-strict but exact certificate. The point Bernstein check reports VIOLATION here. That is
expected and harmless: u is identically 0, so u_-/s = 0 is not strictly positive, but no comparison is needed when
u = 0 exactly.

**K1T - small-j check.** Certified slope enclosures at j = 1/4, 1/2 and 1 (N = 40, K = 8):

| j | slope lower | slope upper | (s - 2-term series)/j^6 | (s - 3-term)/j^8 | (delta - 3-term)/j^8 |
|---|---|---|---|---|---|
| 1/4 | 0.9953466586662603811810323 | 0.9953466586662603811810324 | -0.0011467 | 1.5812e-4 | -1.0561e-5 |
| 1/2 | 0.9817835480500189319327059 | 0.9817835480500189319327060 | -0.0011181 | 1.5406e-4 | -1.0346e-5 |
| 1 | 0.9327994339721661032993351 | 0.9327994339721661032995696 | -0.0010169 | 1.3974e-4 | -9.570e-6 |

How to read it:
- The column (s - 2-term)/j^6 tends to -11 c3 = -0.0011566 as j -> 0. Its gap from that value shrinks like j^2
  (9.9e-6 at j = 1/4, 3.85e-5 at j = 1/2).
- The slope residual is about -15 times the delta residual: 1.584e-4 against 1.581e-4. That is the s(j) = 1 + delta
  - 2j delta' relation at the 4th coefficient, checked between two independently certified objects.

**KN control.** With exponent 6 in the tangent operator only, the bound at j = 64 is <= 0.0712192164366958717014770166549.
That is 0.019 above the true enclosure, so the certificate does respond to the exponent.

**Secant / convexity check.** This is an independent rigorous route that uses only point certificates.
- Why it works: Phi_t = (1-t) Phi_0 + t Phi_1 is a super-solution for mu_t, because x^-7 is convex. So M(mu) at
  fixed J is convex.
- That gives [M(1) - M(1-h)]/h <= s <= [M(1+h) - M(1)]/h, with M(mu) = mu (1 + delta(j/mu^2)).

Results, both containing the tangent enclosure:
- j = 64, h = 2^-30: [0.05204079752240561813432952, 0.05204079760863653128994858]
- j = 4096, h = 2^-24: [0.001319832301227427371698210, 0.001319833430932115888786756]

## 5. Convergence

| j | main degrees | partner degrees | partner delta interval | partner slope interval |
|---|---|---|---|---|
| 64 | N = 160, K = 24 | N = 120, K = 20 | [3.17060411682713773929612960376, ...621121927] | [0.0520407975655210747013801603144, ...7027235936474] |
| 4096 | N = 340, K = 20 | N = 300, K = 18 | [32.2150822519650914925991767626, ...80090526494197] | [0.00131983286602438831788400072503, ...613525395503990435831] |

In both cases the main interval lies inside the partner interval, for delta and for the slope.

How many digits to trust. Every digit inside an interval is certified. Read as shared leading digits:
- delta(64) = 3.170604116827137739296170411..., about 28 significant digits.
- s(64) = 0.052040797565521074702051 8(7..8), about 26 significant digits.
- delta(4096) = 32.2150822519650947..., about 18 significant digits.
- s(4096) = 0.001319832866079(7..9), about 13 significant digits.

Earlier development runs agree, for example N = 210 at j = 4096 and N = 80 at j = 64.

## 6. Times and memory

All runs used the instrument with sha256 9b003feb... (the full hash is in the JSON).

| run | stages (s) | total | peak |
|---|---|---|---|
| j = 64, main | fit 10.4, point cert 6.0, tangent fits 2.0 + 2.1, tangent certs 8.3 + 7.8, KN 2.1 + 6.5 | about 45 s | 1.26 GB |
| j = 4096, main | fit 32.5, point cert 20.1, tangent fits 7.5 + 7.3, tangent certs 27.4 + 31.3 | about 126 s | 4.46 GB |
| j = 4096, partner | - | - | 3.37 GB |
| j = 4096, secant (3 point certificates) | - | about 1 min | 3.60 GB |
| controls | - | - | 0.24 GB |
| K1T | - | - | 0.18 GB |

## 7. Weaknesses

1. **Comparison lemmas are standard but not proved here.** These are maximum principles on R^3 minus the origin for
   bounded continuous functions: the origin is a removable point, functions decay at infinity, and V is bounded with
   V ~ r^2 near the origin. Polynomials in s have a linear-in-r kink at the origin, as in the engine.
2. **Assumed, not proved by me.** Existence of u, the fact that mu -> u is differentiable, and the identification
   dM/dmu = 1 + 2 lim r w. The secant route needs only existence of u and convexity, and it still brackets the
   slope (both one-sided derivatives).
3. **Floating-point arb itself is trusted.** That is python-flint 0.9.0 arb ball arithmetic and the engine's exact
   Bernstein routine. No Lean step was done.
4. **Wasted digits.** Multiplying by D^8 costs up to max D^8, about 1e4 to 1e5 at j = 4096. It costs digits, not
   validity.
5. **Not fully independent of the engine.** The coordinate, ansatz and abs-Chebyshev principle are shared with it;
   only the residual algorithm, the tangent certificate and the secant route are mine.
6. **Exposure incident.** See the disclosures at the top.

## Files

- <LAB>/BHSLICE/rep5/rep5.py
- <LAB>/BHSLICE/rep5/series.py
- <LAB>/BHSLICE/rep5/result_rep5.json
- <LAB>/BHSLICE/rep5/controls_rep5.json
- <LAB>/BHSLICE/rep5/k1t_rep5.json
- <LAB>/BHSLICE/rep5/ckpt/*.json
- <LAB>/BHSLICE/rep5/run_*.log
- <LAB>/BHSLICE/rep5/old/: superseded development runs, made with an earlier instrument hash.
