# BAR-REP6: blind replication of P2e's top node j = 2^21 (the node the sixth-digit decision rests on)

Written 2026-11-14 by the lead after VERDICT-P2e / REFEREE-P2e and before the blind seat starts. Fenced from p2d.py, p2e.py, BUILD-P2d.md, BUILD-P2e.md, result_p2d.json*, result_p2e.json*, controls_p2d.*, controls_p2e.*, run_p2d.log, run_p2e.log, VERDICT-*.md, REFEREE-*.md, BAR-*.md, AMEND-*.md, TARGETS-*.md, DESIGN-P2d.md, LOG.md, ERRATA-PENDING.md, seats/, rep4/, rep5/, and gravity_engine/bowen_york_sqrt.py + its tests. The seat may read LEMMA-REP4.md (public mathematics: comparison lemma and tangent lemma TL, no target numbers) and gravity_engine/lichnerowicz_axisym.py.

instruments: rep6/rep6.py

## Object
Bowen-York spinning puncture at mu1 = 1, j = 2^21 (J = 2^19): (a) a rigorous enclosure of delta = M_ADM - 1; (b) a certified tangent super-solution psi_+ in the sense of LEMMA-REP4 Lemma TL (V1 from a certified sub-solution), giving the coefficient 2L_hi = 1 + 2 lim r psi_+; (c) the assembled tail bound chi_tail = (j/4) / (1 + delta_lo - max(2L_hi, 0))^2, valid for all j >= 2^21 by Lemma TL. Own code, own implementation; a formulation different from a sqrt radial map Y = 1/(1 + sqrt(r/a)) is preferred; if the seat must use a sqrt-type map for feasibility at this spin, it says so and its code is still written from the spec without reading P2e's module.

## Predictions
- R1: blind delta(2^21) overlaps P2e's certified delta(2^21) interval.
- R2: blind 2L_hi(2^21) <= 3 x P2e's 2L_hi(2^21) (REFEREE-P2e 3: the decision would flip only at ~3.1x).
- R3: blind chi_tail < 0.9282005 (the sixth-digit decision reproduced from blind certificates alone).

## Controls (any failure = VOID)
- K1: small-j delta/j^2 contains 1/40.
- KT: psi = 0 gives 2L = 1 and chi_tail = (j/4)/delta_lo^2 form (plain tail).
- KN: source 9/4 -> 9/4 (1 + 1e-3) moves blind delta(2^21) off P2e's interval.

## Verdict rule
REPLICATED if R1-R3 hold and K1, KT, KN pass. R1 disjoint = real discrepancy: stop, the P2e decision is suspended until resolved. If the seat cannot certify at j = 2^21 within limits, report INFEASIBLE (not a pass) with the largest j reached.

## Seal fields
- formula_ref: LEMMA-REP4.md:177-230
- calibrated_on: 560847f6a4, 9e0df55e7a
- depends_on: 560847f6a4, 9677804a51, 9e0df55e7a
- null_not_symmetry: KN changes the source coefficient; KT is an exact reduction.
- target_form: R1 overlap boolean; R2 ratio in [0, 3]; R3 chi_tail in [0, 0.9282005); K1 interval containing 1/40.
- ceiling: P2e certified window [0.9282004538, 0.9282004687]; numerical 0.928200.
- toy_model: J = 0 gives delta = 0 and 2L = 1.
- tolerances: from need; memory cap 8 GB; one process; chunks < 10 min.
- parameters: j = 2^21 (and a small j for K1)
- files: BAR-REP6.md, rep6/rep6.py
- command: memcap 8 python rep6/rep6.py
