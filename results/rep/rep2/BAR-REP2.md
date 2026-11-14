# BAR-REP2: blind replication of P2a at large spin

Written 2026-11-14 by the lead before the blind seat starts. Fenced from p2a.py, BUILD-P2a.md, result_p2a.json, controls_p2a.*, run_p2a.log, VERDICT-P2a.md, BAR-P2a.md, AMEND-P2a-1.md, TARGETS-P2.md, p2design/, seats/build-p2a/, this file, LOG.md. Built from DESIGN.md sections 1-2 and DESIGN-P2.md (which holds no target numbers).

instruments: rep2/rep2.py

## Object
The Bowen-York spinning puncture (mu = 1, j = 4J): an independent rigorous enclosure of delta(j) = M_ADM - 1 at j = 6 and j = 10, by a method other than P2a's (P2a: exact-rational residual in s = a/(r+a), Legendre in cos theta, Bernstein positivity).

## Lead's certified values (NOT shown to the seat)
P2a: delta(6) in [0.423044107116482903, 0.423044107116483035]; delta(10) in [0.738555349918561234, 0.738555349918583461].

## Predictions
- Q1: blind delta(10) overlaps P2a's interval and has relative width <= 1e-8.
- Q2: blind delta(6) overlaps P2a's interval and has relative width <= 1e-8.

## Controls (any failure = VOID)
- K1: the blind code at small j reproduces delta/j^2 -> 1/40.
- K2: source coefficient 9/4 -> 9/4 (1 + 1e-3) moves blind delta(10) off P2a's interval.

## Verdict rule
REPLICATED if Q1, Q2 hold and K1, K2 pass. Disjoint enclosures = a real discrepancy: stop and find the bug.

## Seal fields
- formula_ref: DESIGN.md:1-270
- calibrated_on: d62bfc0e0b
- depends_on: d62bfc0e0b, e44dd02296
- null_not_symmetry: K2 changes the equation's source coefficient (changes the answer).
- target_form: Q1, Q2 interval overlap with P2a's certified intervals and relative width in [0, 1e-8]; K1 interval containing exact 1/40.
- ceiling: Bowen-York chi_max ~ 0.928 (numerical); j = 10 gives chi ~ 0.827.
- toy_model: J = 0 (u = 0) through the blind code.
- tolerances: from need (relative width 1e-10 or better targeted).
- parameters: j = 6, 10 (and a small j for K1)
- files: BAR-REP2.md, rep2/rep2.py
- command: memcap 8 python rep2/rep2.py
