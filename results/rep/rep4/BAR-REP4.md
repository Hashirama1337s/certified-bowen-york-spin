# BAR-REP4: blind replication of P2c at the highest spins

Written 2026-11-14 by the lead while P2c main runs, before any P2c result exists and before the blind seat starts. Fenced from p2c.py, BUILD-P2c.md, result_p2c.json*, controls_p2c.*, run_p2c.log, ceiling_p2c*, VERDICT-P2c.md, BAR-P2c.md, TARGETS-P2.md, seats/, P2C-SCOUT.md, this file, LOG.md, every VERDICT-*.md, result_p2b.json*, controls_p2b.*, and the gravity_engine package (the seat must not reuse the certified engine). The seat may read its own earlier blind code in rep3/ only if it was not written by the engine's author (rep3 was a blind seat: allowed).

instruments: rep4/rep4.py

## Object
Bowen-York spinning puncture (mu = 1, j = 4J): an independent rigorous enclosure of delta(j) = M_ADM - 1 at j = 3000 and j = 10000, by a method other than the engine's (exact-rational residual in s = a/(r+a), Chebyshev x Legendre, tensor Bernstein positivity).

## Lead's values
None exist yet at sealing time (P2c main running). Scoring compares to P2c's certified intervals at j = 3000 and j = 10000 once VERDICT-P2c is settled.

## Predictions
- R1: blind delta(10000) overlaps P2c's certified interval, relative width <= 1e-3.
- R2: blind delta(3000) overlaps P2c's certified interval, relative width <= 1e-3.

## Controls (any failure = VOID)
- K1: small-j delta/j^2 contains 1/40.
- K2: source 9/4 -> 9/4 (1 + 1e-3) moves blind delta(10000) off P2c's interval (checked at scoring).

## Verdict rule
REPLICATED if R1, R2 hold and K1, K2 pass. Disjoint = real discrepancy: stop and find the bug before any ceiling claim.

## Seal fields
- formula_ref: DESIGN.md:1-270
- calibrated_on: a669f700df
- depends_on: 3031217ac8, a669f700df
- null_not_symmetry: K2 changes the source coefficient (changes the answer).
- target_form: R1, R2 overlap with P2c's intervals and relative width in [0, 1e-3]; K1 interval containing exact 1/40.
- ceiling: Bowen-York chi_max ~ 0.928200 (numerical); j = 10000 gives chi ~ 0.9282.
- toy_model: J = 0 (u = 0) through the blind code.
- tolerances: from need; memory cap 8 GB; one process.
- parameters: j = 3000, 10000 (and a small j for K1)
- files: BAR-REP4.md, rep4/rep4.py
- command: memcap 8 python rep4/rep4.py
