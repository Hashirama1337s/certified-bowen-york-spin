# BAR-P2d: a sharp certified spin ceiling for the Bowen-York spinning-puncture family (tangent lemma)

Written 2026-11-14 by the lead before p2d.py exists. Main quest, Level 2. Design: DESIGN-P2d.md (sealed f825f520) + AMEND-P2d-1 (ea950c756f) + AMEND-P2d-2. Gate: GATE-P2d.md (claim C OPEN; claim T OPEN as a method, neighbours DGC Lemma 3.1 and De Ceuster-Li 2609.10173). Referees: REFEREE-REP4 (M), REFEREE-TL (N). Engine: gravity_engine/bowen_york_sqrt.py module version 1 (ENGINE_VERSION 4, POCKETLAB 4864f55) + the tangent certificate added by the build seat. Lemma: LEMMA-REP4.md (Lemma C-REP4, C-REP4', TL; referee O must return CORRECT or CORRECT WITH FIXES before --main).

instruments: p2d.py

## Object
For the single, non-boosted, conformally flat, maximal Bowen-York spinning puncture (puncture BC, J = J_ADM), mu > 0: at nodes j_k = 2^k, k = 0..18, certify (a) the point enclosure delta(j_k) by the sqrt-map certifier and (b) a tangent super-solution psi_+ (TL-H) giving the certified slope bound 2L_hi; then by TL: cells [j_k, j_{k+1}] and the tail [j_18, oo) give chi bounds; [0, 1] by P1 (chi <= j/4). Ceiling c := max over these pieces, in exact arithmetic. Lower side X3 := max chi_lo over nodes.

## Predictions
- D1: every node's point and tangent certificates certify (19/19 each).
- D2: c <= 0.935 (tier A of DESIGN-P2d holds even if higher tiers degrade).
- D3: c - X3 <= 1e-4 (the certified window shrinks at least ~340x from P2c's 3.45e-2).
- D4: consistency with earlier certified results: (i) at every P2d node that is exactly a node of P1, P2a, P2b or P2c, the two certified delta intervals overlap; (ii) at every other P2d node j in [1, 11994.82], the new delta_hi >= sqrt(j/j_a) delta_lo(j_a), where j_a is the nearest earlier certified node below j (the monotone-cell lower bound); (iii) c <= P2c's c = 0.962234888831116 and X3 >= P2c's X3_max 0.928146570605.
- D5: K1 (delta/j^2 contains 1/40 at j = 2^-64) and K1T (2L slope consistent with 1 - (3/40) j^2 at small j: a CONSISTENCY control, AMEND-P2d-1 item 6) both hold.

## Controls (any failure = VOID; each its own command; GO via `lab.py go`)
- ENGINE: gravity_engine tests pass (including the bowen_york_sqrt reproduction of REP4 A at j = 3000, 10000).
- psi = 0 control: with psi_+ = 0 the instrument reproduces P2c-type plain cells (2L = 1) exactly (REFEREE-TL 1).
- N-TAN: a 1e-4 shift of the tangent fails the tangent-consistency test at node pairs where the expected gap (computed by the build seat and written into AMEND-P2d-3 before data) exceeds the test's resolution; pairs where it cannot bite are listed in advance.
- N1: wrong-ansatz point certificate at the top node -> FAIL; N4 near-threshold margin -> FAIL through the bound path.
- Structural asserts: Z ansatz structure, degree preconditions (N_Z = N, K >= 1), piece cover exact, tail and cell bracket 1 + delta_lo - tL > 0 in exact arithmetic (AMEND-P2d-1 items 2-4).
- CEIL: c recomputed independently (exact Fractions) from the stored node records; TAIL bracket checked independently.
- RESUME byte-identical; exit codes 0 / 3 / other; MEM peak per invocation; NaN/inf; `expected_false` lists by-design false flags ("<path> (<why>)").

## Verdict rule
P2d CERTIFIED if D1-D5 hold and every control behaves as stated. If c > 0.935 or c - X3 > 1e-4, report the certified c and window (they are results) and the binding piece; no retuning (AMEND-P2d-2 item 3). Any overlap failure in D4 = stop and investigate.

## Seal fields
- formula_ref: DESIGN-P2d.md:1-200
- calibrated_on: 3031217ac8, ad62cca292
- depends_on: 3031217ac8, ad62cca292, 1a9c1237d9
- null_not_symmetry: N-TAN shifts the tangent (changes the bound); N1/N4 perturb the sandwich; psi = 0 is an exact known reduction; CEIL/TAIL are independent recomputations.
- target_form: D1 booleans; D2 c in [0, 0.935]; D3 c - X3 in [0, 1e-4]; D4 overlap booleans and inequalities; D5 interval containment.
- ceiling: numerical estimate 0.928200 (Lovelace et al. 0805.4192 eq (5)); P2c certified window [0.928147, 0.962235]; Kerr bound 1.
- toy_model: J = 0 with nonzero W at CS = 0 must FAIL under a forced too-small margin; W = 0 REFUSED/TRIVIAL.
- tolerances: degrees from the build seat's schedule (fixed before main); REP4 fixed pieces by default (AMEND-P2d-2); one process; chunks <= 520 s.
- parameters: nodes j = 2^k, k = 0..18
- files: BAR-P2d.md, p2d.py
- command: memcap 8 python p2d.py --main result_p2d.json
