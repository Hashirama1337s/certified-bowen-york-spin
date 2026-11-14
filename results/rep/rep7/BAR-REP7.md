# BAR-REP7: blind replication of B3T (spin of the extreme mu = 0 Bowen-York data)

Written 2026-11-14 by the lead after B3T (10a65e940a) settled and before the blind seat starts. Fenced from b3t.py, BUILD-B3T.md, result_b3t.json, controls_b3t.*, run_b3t.log, BAR-*.md, VERDICT-*.md, REFEREE-*.md, TARGETS-*.md, LOG.md, ERRATA-PENDING.md, RAID-WING1.md, result_p2*.json, scout/b3/, seats/, rep6/. The seat may read scout/B3-EXTREME.md sections on the equation and the comparison class ONLY if it does not read the probe values (Q4 table) - simpler: the brief gives the mathematics; the seat reads no BHSLICE scout files.

instruments: rep7/rep7.py

## Object
The extreme (mu = 0) Bowen-York data, J = 1: Delta u = -S (1 + u)^-7, S = 9 sin^2(theta)/(4 r^6), u -> 0 at infinity, u -> +infinity at 0 like r^-1/2. A rigorous two-sided enclosure of M_0 = 2 lim r u_0 and hence of chi(0) = 1/M_0^2, by the seat's own method; a formulation different from w = r^-1/2 G/(1 + sqrt r) with xi = tanh(sigma/8) is preferred (e.g. a different throat compactification, a sub-solution + super-solution pair instead of one super-solution with an integral lower bound, or a different positivity method).

## Predictions
- Q1: the blind chi(0) window overlaps B3T's certified window.
- Q2: blind window width <= 1e-8.

## Controls (any failure = VOID)
- RN: extreme Reissner-Nordstrom analogue (exponent 3, exact M = 1) gives a window containing 1.
- KN: source 9/4 -> 9/4 (1 + 1e-3) moves the blind window by more than its width.

## Verdict rule
REPLICATED if Q1, Q2 hold and RN, KN pass. Disjoint = stop and find the bug before any extreme-data claim.

## Seal fields
- formula_ref: LEMMA-REP4.md:177-230
- calibrated_on: 10a65e940a
- depends_on: 10a65e940a, 560847f6a4
- null_not_symmetry: KN changes the source; RN is exactly known.
- target_form: Q1 overlap boolean; Q2 width in [0, 1e-8].
- ceiling: P2e family window [0.9282004538 (down), 0.9282004687 (up)].
- toy_model: extreme RN (exact M = 1).
- tolerances: from need; 8 GB; one process; chunks < 10 min; about 60 min compute.
- parameters: mu = 0, J = 1
- files: BAR-REP7.md, rep7/rep7.py
- command: memcap 8 python rep7/rep7.py
