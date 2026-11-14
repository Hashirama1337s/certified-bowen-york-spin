# BAR-REP5: blind replication of P2d's tangent (slope) certificates

Written 2026-11-14 by the lead while P2d main runs, before any P2d result exists and before the blind seat starts. Fenced from p2d.py, BUILD-P2d.md, result_p2d.json*, controls_p2d.*, run_p2d.log, VERDICT-P2d.md, BAR-P2d.md, AMEND-P2d-*.md, DESIGN-P2d.md, TARGETS-P2d.md, LEMMA-REP4.md, REFEREE-TL.md, REFEREE-LEMMA.md, GATE-P2d.md, seats/, LOG.md, every VERDICT-*.md, and the module gravity_engine/bowen_york_sqrt.py and its tests. The seat may use the s-coordinate engine gravity_engine/lichnerowicz_axisym.py (a different formulation from P2d's sqrt-map module) or its own code.

instruments: rep5/rep5.py

## Object
Bowen-York spinning puncture, J fixed, bare mass mu varied (j = 4J/mu^2). At mu1 = 1 and j in {64, 4096}: a rigorous upper bound on the slope s(j) := dM_ADM/dmu at mu = mu1 (fixed J), obtained by certifying a super-solution of the linearised (tangent) problem Delta w = 7 S Phi^-8 (1/(2r) + w) (with Phi from a certified sub-solution of the point problem), plus a rigorous enclosure of delta(j) at the same j.

## Predictions
- S1: the blind slope upper bound at j = 64 and j = 4096 lies within 1e-3 (absolute) of P2d's certified slope upper bound 2L_hi there (both bound the same true slope from above; agreement within 1e-3 tests that both certify the same object).
- S2: blind delta(64) and delta(4096) overlap P2d's point intervals.

## Controls (any failure = VOID)
- K1T: small-j slope consistent with 1 - (3/40) j^2 (consistency control).
- KN: replacing the exponent 7 by 6 in the tangent operator (only) moves the slope bound by more than 1e-3 at j = 64.

## Verdict rule
REPLICATED if S1, S2 hold and K1T, KN pass. A slope disagreement beyond 1e-3 = stop and find the bug before any P2d claim.

## Seal fields
- formula_ref: DESIGN-P2d.md:1-200
- calibrated_on: ad62cca292
- depends_on: 9677804a51, ad62cca292
- null_not_symmetry: KN changes the tangent operator (changes the answer).
- target_form: S1 slope bound difference in [0, 1e-3] in absolute value; S2 overlap booleans; K1T within the seat's stated tolerance.
- ceiling: P2c window [0.928147, 0.962235]; numerical 0.928200.
- toy_model: J = 0 gives slope exactly 1 (M = mu).
- tolerances: from need; memory cap 8 GB; one process.
- parameters: j = 64, 4096
- files: BAR-REP5.md, rep5/rep5.py
- command: memcap 8 python rep5/rep5.py
