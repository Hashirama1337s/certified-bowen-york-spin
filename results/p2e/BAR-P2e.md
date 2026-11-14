# BAR-P2e: decide the sixth digit of the Bowen-York maximal spin

Written 2026-11-14 by the lead before p2e.py exists. Main quest, Level 2. Builds on P2d (HOLD 9677804a51; REFEREE-P2d CORRECT, wording fixed): certified sup chi in [0.9282003665, 0.9282009626]. REFEREE-P2d 4(d) (POST HOC, CHECK only): chi_lo increments k = 15..18 shrink by ~0.37 per doubling, suggesting a limit near 0.92820044; dropping eta_Z at node 2^18 alone gives a tail of 0.92820062, so more nodes and tighter top-node certificates are needed.

instruments: p2e.py

## Object
Same family, same certified route as P2d (sqrt-map point certificate + tangent certificate per node; tangent-lemma cells and tail; P1 on [0, 1]), engine gravity_engine/bowen_york_sqrt.py module v2 (no engine change unless a new module version is tested and committed before p2e.py is sealed). Nodes j = 2^k, k = 0..K_top with K_top in [20, 24] and per-band degrees chosen by the build seat from probes and FIXED in p2e.py before main.

## Predictions
- E1: every node's point and tangent certificates certify.
- E2: c - X3 < 1e-7.
- E3 (the decision): c < 0.9282005, i.e. the family's supremum rounds to 0.928200 at six decimals (Lovelace et al.'s 0.928200 is the correct rounding). Based on the POST HOC geometric extrapolation above (~0.92820044); stated before any P2e data.
- E4: consistency with P2d: c <= P2d's c and X3 >= P2d's X3; at every shared node the delta intervals overlap.

## Controls (any failure = VOID; each its own command; GO via `lab.py go`)
- ENGINE tests pass. psi = 0 reproduces the plain cell/tail exactly. N1/N4 and a tangent bump FAIL at the top node. N-TAN at the top pairs (expected gaps written in the build report; labelled as reproducibility, not blind - REFEREE-P2d 5).
- CEIL: independent exact-Fractions recomputation run by the instrument on the REAL result after main (not only a mock), plus TAIL.
- Structural asserts (Z ansatz, K >= 1, cover, bracket > 0); RESUME byte-identical; self-sha / GO row; exit codes 0 / 3 / other; MEM; NaN/inf; expected_false.

## Verdict rule
P2e CERTIFIED if E1, E2, E4 hold and every control behaves as stated; E3 is then scored as stated (HOLD if c < 0.9282005; FAILED if X3 >= 0.9282005 - which would also be a decision, the opposite rounding; UNDECIDED-FAILED if neither, report the certified window). No retuning after seeing values.

## Seal fields
- formula_ref: DESIGN-P2d.md:1-200
- calibrated_on: 9677804a51
- depends_on: 9677804a51, ad62cca292
- null_not_symmetry: N1/N4/tangent bump perturb certificates; psi = 0 exact reduction; CEIL/TAIL independent recomputations.
- target_form: E1 booleans; E2 c - X3 in [0, 1e-7); E3 c in [0, 0.9282005); E4 inequalities and overlap booleans.
- ceiling: P2d window [0.9282003665, 0.9282009626]; Lovelace et al. 0805.4192 eq (5) "eps_J <= 0.928200".
- toy_model: J = 0 refused/trivial; nonzero-W toy at CS = 0 fails under a forced margin.
- tolerances: per-band degrees fixed in p2e.py before main; REP4 fixed pieces unless AMEND-P2d-2 targeted refinement is used (recorded per node); one process; chunks <= 520 s.
- parameters: nodes 2^k, k = 0..K_top, K_top in [20, 24]
- files: BAR-P2e.md, p2e.py
- command: memcap 8 python p2e.py --main result_p2e.json
