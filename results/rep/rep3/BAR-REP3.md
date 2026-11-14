# BAR-REP3: blind replication of P2b at high spin

Written 2026-11-14 by the lead before the blind seat starts. Fenced from p2b.py, BUILD-P2b.md, result_p2b.json*, controls_p2b.*, run_p2b.log, VERDICT-P2b.md, BAR-P2b.md, TARGETS-P2.md, p2bprobe/, seats/, P2C-SCOUT.md, this file, LOG.md, VERDICT-P2a.md, VERDICT-REP2.md, and the gravity_engine package (the seat must not reuse the certified engine).

instruments: rep3/rep3.py

## Object
Bowen-York spinning puncture (mu = 1, j = 4J): an independent rigorous enclosure of delta(j) = M_ADM - 1 at j = 300 and j = 1000 by a method other than the engine's (exact-rational residual in s = a/(r+a), Legendre in cos theta, Bernstein positivity); e.g. the REP2 method (y = 1/(1+r), cos^2 theta basis, Taylor-form positivity) extended to high degree, with radii controlled.

## Lead's certified values (NOT shown to the seat)
P2b: delta(300) in [7.993758667606936304, 7.993758727426259014]; delta(1000) in [15.41318755062711644, 15.41318839322056450].

## Predictions
- Q1: blind delta(1000) overlaps P2b's interval, relative width <= 1e-6.
- Q2: blind delta(300) overlaps P2b's interval, relative width <= 1e-6.

## Controls (any failure = VOID)
- K1: small-j delta/j^2 contains 1/40.
- K2: source 9/4 -> 9/4 (1 + 1e-3) moves blind delta(1000) off P2b's interval.

## Verdict rule
REPLICATED if Q1, Q2 hold and K1, K2 pass. Disjoint = real discrepancy: stop and find the bug.

## Seal fields
- formula_ref: DESIGN.md:1-270
- calibrated_on: eba765ebe0
- depends_on: eba765ebe0, 227cbcafc2
- null_not_symmetry: K2 changes the source coefficient (changes the answer).
- target_form: Q1, Q2 overlap with P2b's intervals and relative width in [0, 1e-6]; K1 interval containing exact 1/40.
- ceiling: Bowen-York chi_max ~ 0.928200 (numerical); j = 1000 gives chi ~ 0.92801.
- toy_model: J = 0 (u = 0) through the blind code.
- tolerances: from need.
- parameters: j = 300, 1000 (and a small j for K1)
- files: BAR-REP3.md, rep3/rep3.py
- command: memcap 8 python rep3/rep3.py
