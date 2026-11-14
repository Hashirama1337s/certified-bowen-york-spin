# BAR-P1: first certified spinning black-hole slice (Bowen-York puncture, mu = 1, every J in [0, 1/4] at once)

Written 2026-11-14 by the lead before any certifier code exists. Targets avoid rounded-constant containment (lesson from K3 WP2a): reference values are matched by OVERLAP with explicit windows.

instruments: p1.py

## Object
DESIGN.md (BHSLICE): Delta u = -9 J^2 sin^2(theta) / (4 r^6 (1 + mu/(2r) + u)^7), u -> 0 at infinity, u regular at the puncture; M_ADM = mu + 2 lim r u; mu = 1; j = 4J. Compactified s = a/(r + a), a = mu/2, c = cos theta, on the closed square. Certify explicit sub/super-solutions u- <= u <= u+ (comparison lemma of DESIGN section 2, weak maximum principle) on a perturbative slab J in [0, 0.1] and 18 slabs of width 0.05... covering [0, 1/4] in J (the seat may choose slabs; they must cover [0, 1/4]), with residual signs checked rigorously (arb) on the whole closed square for every J in each slab; then certified enclosures of delta(J) = M_ADM(J) - 1 and chi(J) = J / M_ADM(J)^2.

## Reference values (double precision, three resolutions agree to ~1e-15; DESIGN section 3)
delta(0.0025) = 2.4999874054e-6; delta(0.025) = 2.4987415807e-4; delta(0.0625) = 1.5576057054e-3; delta(0.125) = 6.1728856170e-3; delta(0.25) = 2.3836105153e-2.

## Predictions
- B1: every slab closes (strict residual signs for u+ and u-, u+- >= 0 where required, 1 + (1-s) u- >= 1/2).
- B2: at each of the 5 reference J, the certified delta enclosure OVERLAPS the window [ref (1 - 1e-8), ref (1 + 1e-8)].
- B3: relative width of the certified delta enclosure <= 1e-6 at each reference J (and on every slab at its endpoints).
- B4 (positive control): the first-slab certified enclosure of delta / j^2 as j -> 0 contains the exact published coefficient 1/40.

## Controls (any failure = VOID)
- N1 (null): raising the sub-solution's slope coefficient by 1e-3 relative must make the certifier report a sign violation.
- N2 (null): u+ - 1e-6 j^2 s must make the certifier report a sign violation.
- N3: the stage-0 closed-form bracket 0 <= u <= j^2 v1 (DESIGN stage 0) is reproduced by the certifier as a valid (looser) sandwich on [0, 1/4].
- No NaN/inf; controls run on the final sealed file (sha256 recorded; lab controls-check GO).

## Verdict rule
P1 CERTIFIED if B1-B4 hold and N1-N3 behave as stated. Any certified enclosure disjoint from its reference window = investigate (instrument or reference), no retuning; a null that passes = VOID.
Noise floor/ceiling: the reference values are double-precision numerics agreeing to ~1e-15 across resolutions, far inside the 1e-8 windows; the certified enclosure must contain the true value, so a correct certificate overlaps each window.

## Seal fields
- formula_ref: DESIGN.md:1-270
- calibrated_on: none - first BHSLICE test; reference values from the design seat's double-precision solve (probe0/solve.log)
- depends_on: none - new project
- null_not_symmetry: N1 and N2 perturb the sandwich functions (they stop being sub/super-solutions), N3 is a known looser sandwich.
- target_form: B1 boolean; B2 overlap with [ref(1-1e-8), ref(1+1e-8)]; B3 [0, 1e-6]; B4 interval containing exact 1/40.
- tolerances: Chebyshev degrees and slab widths chosen from need; refinement must leave >= 2x headroom on B3 (justified in BUILD-P1.md); queues must shrink.
- parameters: J in [0, 1/4] in >= 6 slabs; reference J = 0.0025, 0.025, 0.0625, 0.125, 0.25
- files: BAR-P1.md, p1.py
- command: memcap 8 python p1.py --main result_p1.json
