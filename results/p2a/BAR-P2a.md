# BAR-P2a: certified Bowen-York spinning-puncture data for every spin j = 4J/mu^2 in [1, 10]

Written 2026-11-14 by the lead before any P2 certifier code exists. Main quest, Level 2. Applies REFEREE-P2DESIGN fixes 3a, 3b, 3c, 6a, 6b and the S > 0 assertion.

instruments: p2a.py

## Object
DESIGN-P2.md sections 1-4.1 (engine: REP1's exact-rational residual moved to P1's s-coordinate and Legendre basis; R1 port, R3 Bernstein positivity of V_- = U_-/s, R4 margin eta = S (1 + 2^-20)/2 with an explicit S > 0 assertion, R5 slope route with exact rational endpoints, optional R6 integral route). Point certificates at every node of the geometric grid j_{k+1} = (51/50) j_k from j_0 = 1 up to the first node >= 10, plus the extra nodes j = 2, 4, 6, 8, 10. Cell enclosures K1 between neighbouring nodes using delta(j)/sqrt(j) nondecreasing (corollary of Dain-Gabach Clement Lemma 3.1 + scale invariance; Lean-checked cell inequality lean 036ea0f37b), evaluated in arb with outward rounding (fix 3b), with per-cell chi bounds chi_max = j_{k+1}/(4(1 + sqrt(j_{k+1}/j_k) delta_-(j_k))^2), chi_min = j_k/(4(1 + sqrt(j_k/j_{k+1}) delta_+(j_{k+1}))^2) (fix 3a). mu = 1, J = j/4.

## Predictions (expected values: TARGETS-P2.md item T1, sealed b18e905ff1)
- B1: every grid and extra node certifies (strict residual signs, V_- > 0 by Bernstein, S > 0).
- B2: at j = 2, 4, 6, 8, 10 the certified delta overlaps [ref (1 - 1e-9), ref (1 + 1e-9)] with ref from TARGETS T1.
- B3: relative width of the certified delta <= 1e-10 at j = 2, 4, 6, 8, 10.
- B4: every K1 cell on the 51/50 grid covering [1, 10] has relative delta width <= 3e-2, consecutive cells share their exact rational endpoint (no gap) (fix 3c: the grid is 51/50, not h = 0.1).
- B5: the j = 1 node overlaps P1's certified ball for delta(J = 1/4) and REP1's 21-digit enclosure.

## Controls (any failure = VOID; each run as its own command stopping on a non-zero exit)
- N1, N2: wrong-ansatz nulls (sub-solution slope +1e-3 relative; u+ - 1e-6 j^2 s) give certified VIOLATION.
- N4 (REFEREE-P1): eta' = 0.25 x the needed margin at j = 2 and j = 10 gives FAIL through the bound path (no node violation).
- N-pos: V_- shifted down past its certified Bernstein minimum gives FAIL; shifted up gives PASS.
- N-cell: swapped cell endpoints exclude a mid-cell certificate.
- T6a (fix 6a): exact identity at a nonzero rational j with U_N = v1: Q equals (9/4)(1-c^2) s^2 (1-s)^2 (1 - (1 + (1-s) j^2 v1)^7) exactly.
- T0b (fix 6b): j = 0 gives status REFUSED/TRIVIAL with S == 0 and p(0) == 1/40 exactly.
- T0c: every enclosure lies inside the stage-0 bracket and satisfies chi < 1.
- RESUME: an interrupted run resumed gives a byte-identical result file (endpoints stored as exact p/q strings; resume_tested true).
- NaN/inf guards; controls on the final file (instrument sha256 recorded); GO only through `lab.py go`.

## Verdict rule
P2a CERTIFIED if B1-B5 hold and every control behaves as stated. A certified enclosure disjoint from its reference window = investigate (instrument or reference), no retuning; a null that passes = VOID.

## Seal fields
- formula_ref: DESIGN-P2.md:1-120
- calibrated_on: 565736c4c4, e44dd02296
- depends_on: 565736c4c4, e44dd02296
- null_not_symmetry: N1, N2, N4, N-pos perturb the sandwich (they stop being sub/super-solutions); N-cell swaps an inequality; T0b/T6a are exact known answers.
- target_form: B1 boolean; B2 overlap with [ref(1-1e-9), ref(1+1e-9)]; B3 [0, 1e-10]; B4 [0, 3e-2] and boolean; B5 boolean overlap.
- ceiling: the Bowen-York spin limit chi_max ~ 0.928 (numerical; Lovelace et al. 0805.4192 eq (5): 0.928200) and Kerr's chi < 1; P2a reaches chi ~ 0.83 at j = 10 (TARGETS T1), distance to the BY ceiling ~ 0.10.
- toy_model: J = 0 (Schwarzschild puncture, u = 0: T0b) and the exact order-j^2 solution v1 through the full D^7 code path at a nonzero rational j (T6a).
- tolerances: Chebyshev/Legendre degrees from need for B3 with >= 10x headroom; Bernstein subdivision until positive or certified FAIL; queues must shrink.
- parameters: j grid 51/50 from 1 to >= 10 (about 117 nodes) plus j = 2, 4, 6, 8, 10
- files: BAR-P2a.md, p2a.py
- command: memcap 8 python p2a.py --main result_p2a.json
