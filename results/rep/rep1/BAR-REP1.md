# BAR-REP1: blind replication of the P1 black-hole certificate

Written 2026-11-14 by the lead before the blind seat starts. The seat is fenced from p1.py, BUILD-P1.md, result_p1.json, controls_p1.*, run_p1.log, VERDICT-P1.md, BAR-P1.md, this file, probe0/ and LOG.md; it builds from DESIGN.md sections 1-2 only (the problem and the comparison lemma), with a different discretisation.

instruments: rep1/rep1.py

## Object
The Bowen-York spinning puncture (mu = 1) of DESIGN.md: an independent rigorous enclosure (sub/super-solutions certified by interval arithmetic, any method; or a validated Newton-Kantorovich enclosure) of delta(J) = M_ADM(J) - 1 at J = 0.0625 and J = 0.25.

## Lead's certified values (NOT shown to the seat)
P1: delta(0.0625) in [1.5576057054402e-3 +/- 5.39e-17]; delta(0.25) in [2.38361051533e-2 +/- 3.31e-14].

## Predictions
- Q1: the blind enclosure of delta(0.25) overlaps P1's ball [0.0238361051533 - 3.31e-14, 0.0238361051533 + 3.31e-14] and has relative width <= 1e-6.
- Q2: the blind enclosure of delta(0.0625) overlaps P1's ball [0.0015576057054402 - 5.39e-17, 0.0015576057054402 + 5.39e-17] and has relative width <= 1e-6.

## Controls (any failure = VOID)
- K1: the blind code reproduces the published small-spin coefficient: delta/j^2 -> 1/40 (j = 4J), from its own small-J run.
- K2: a perturbed source (coefficient 9/4 -> 9/4 (1 + 1e-3)) must move the blind delta(0.25) enclosure off P1's ball.

## Verdict rule
REPLICATED if Q1, Q2 hold and K1, K2 pass. Disjoint enclosures = a real discrepancy: stop and find the bug before any further claim.

## Seal fields
- formula_ref: DESIGN.md:1-270
- calibrated_on: 565736c4c4
- depends_on: 565736c4c4
- null_not_symmetry: K2 changes the equation's source coefficient (changes the answer).
- target_form: Q1, Q2 interval overlap with P1's certified balls and relative width in [0, 1e-6]; K1 interval containing exact 1/40.
- tolerances: from need (relative width target 1e-8 or better; 1e-6 is the bar).
- parameters: J = 0.0625, 0.25 (and a small J for K1)
- files: BAR-REP1.md, rep1/rep1.py
- command: memcap 8 python rep1/rep1.py
