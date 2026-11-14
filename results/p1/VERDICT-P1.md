# VERDICT-P1 (scored as sealed: claim 565736c4c4, BAR file e1c57f6005; instrument p1.py pinned + sealed 51253b7dc7; launch-check passed; main 1 invocation, exit 0, peak 0.44 GB)

PROCESS SLIP (disclosed): `lab controls-check 565736c4c4` returned NOT GO - p1.py has a checkpoint/resume path and its controls do not test it (doctrine #118, added the same morning after the K3 FAM1 checkpoint crashes). The lead's command chained the steps with ";" and launched anyway, and GO-P1 first claimed "controls-check GO" (corrected in GO-P1). The main ran in ONE invocation, so the untested resume path was never exercised; no reported number depends on it. Lesson: GO gates must be separate commands that stop on non-zero exit.

## Certified (Bowen-York spinning puncture, bare mass mu = 1, EVERY spin J in [0, 1/4] at once, 19 slabs in j = 4J)
For every J in [0, 1/4] an explicit sub/super-solution pair u- <= u <= u+ of the Hamiltonian constraint is certified on the whole compactified (s, cos theta) square, giving certified ADM masses M_ADM(J) = 1 + delta(J) and spins chi(J) = J / M_ADM^2:

| J | delta = M_ADM - 1 (certified ball) | rel. width | chi (certified ball) |
|---|---|---|---|
| 0.0025 | [2.49998740540817e-6 +/- 5.43e-21] | 1.19e-15 | [0.0024999875001098 +/- 4.74e-17] |
| 0.025 | [2.4987415806877e-4 +/- 5.04e-18] | 7.96e-15 | [0.024987510973319 +/- 2.22e-17] |
| 0.0625 | [1.5576057054402e-3 +/- 5.39e-17] | 4.84e-14 | [0.062305753244328 +/- 1.31e-16] |
| 0.125 | [6.172885617036e-3 +/- 8.62e-16] | 1.66e-13 | [0.12347095108294 +/- 2.01e-15] |
| 0.25 | [2.38361051533e-2 +/- 3.31e-14] | 5.96e-13 | [0.23849491691357 +/- 6.40e-15] |
Small-spin expansion, certified on the first slab: delta/j^2 -> [0.025000000000 +/- 2.13e-13] (contains the published 1/40) and (delta - j^2/40)/j^4 -> [-0.00125946969697 +/- 8.30e-15], containing -133/105600 (the j^4 coefficient derived by the design seat, found in no opened source; existence of the limit not proved, see BUILD-P1 weakness 5).

| | sealed bar | result | verdict |
|---|---|---|---|
| B1 | every slab closes | 19/19 "all slabs ok = True" | HOLD |
| B2 | certified delta overlaps [ref(1 - 1e-8), ref(1 + 1e-8)] at 5 J | all 5 overlap; agreement with the double references to ~1e-11 relative | HOLD |
| B3 | relative width <= 1e-6 | max 5.96e-13 (slab edges up to 6.9e-13) | HOLD |
| B4 | first-slab delta/j^2 contains 1/40 | [0.025000000000 +/- 2.13e-13] | HOLD |
Controls (controls_p1.json, all pass, 25 s): N1 (sub-solution slope +1e-3 rel) certified VIOLATION on slabs 0 and 18; N2 (u+ - 1e-6 j^2 s) certified VIOLATION; N3 stage-0 bracket (exact sympy identity, v1/s > 0 certified, sharp enclosures inside stage-0 brackets); known answers inside the enclosures; NaN/inf guards.

Verdict as sealed: P1 CERTIFIED (B1-B4 HOLD, N1-N3 as stated), with the GO-gate process slip disclosed above.
Not yet: nothing refereed (comparison lemma, puncture regularity, |P_l| <= rho^l bound, divided operator, box engine); box-split and resume paths never executed; no blind high-precision replication (the enclosures are tighter than any double reference, so a blind arb solve is the real check); M_ADM normalisation 1/(2 pi) inherited from DESIGN (literature misprint corrected there).
Main-quest status: Level 2 (real gravity) first certified result. Next: referee + blind replication, then larger spins toward the ~0.928 Bowen-York ceiling.
