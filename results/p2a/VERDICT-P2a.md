# VERDICT-P2a (scored as sealed: claim d62bfc0e0b, BAR file de758dbbf3, AMEND-P2a-1 9356a34707 (B4 = pointwise cell width, sealed before the main run); instrument p2a.py pinned + sealed 91b477f300; `lab.py go` GO (launch-check + controls-check, resume path tested); main 1 invocation, exit 0, peak 0.35 GB)

## Certified (Bowen-York spinning puncture, mu = 1, EVERY spin j = 4J in [1, 10.144])
123 point certificates (51/50 grid from j = 1 plus j = 2, 4, 6, 8, 10) and 122 monotone cells (delta(j)/sqrt(j) nondecreasing: corollary of Dain-Gabach Clement Lemma 3.1 + scaling; cell inequality Lean-checked, lean 036ea0f37b), no gaps. Exact-rational residual engine in P1's s-coordinate with Bernstein positivity.

| j | J | certified delta = M_ADM - 1 | rel. width | certified chi |
|---|---|---|---|---|
| 1 | 0.25 | [0.023836105153325914, 0.023836105153325914] (full: [0.02383610515332591435336432, 0.02383610515332591435337072]) | 2.68e-22 | 0.23849491691356 |
| 2 | 0.5 | [0.084685104039633339, 0.084685104039633340] | 5.36e-18 | 0.42497428536413 |
| 4 | 1 | [0.250423126410051656, 0.250423126410055940] | 1.71e-14 | [0.63956693845623, 0.63956693845624] |
| 6 | 1.5 | [0.423044107116482903, 0.423044107116483035] | 3.12e-16 | 0.74072079402632 |
| 8 | 2 | [0.586331612310199175, 0.586331612310201658] | 4.23e-15 | 0.79477104734360 |
| 10 | 2.5 | [0.738555349918561234, 0.738555349918583461] | 3.01e-14 | [0.82710941638178, 0.82710941638180] |
| 10.144 | 2.536 | [0.749129523043965742, 0.749129523043991110] | 3.39e-14 | [0.82892676839664, 0.82892676839666] |

| | sealed bar | result | verdict |
|---|---|---|---|
| B1 | every node certifies | 123/123 CERTIFIED (one Bernstein piece each) | HOLD |
| B2 | delta overlaps ref(1 +- 1e-9) at j = 2, 4, 6, 8, 10 (TARGETS T1) | all 5 inside their windows; agreement with the double references to ~1e-13 relative | HOLD |
| B3 | rel width <= 1e-10 at j = 2..10 | max 3.01e-14 (j = 10) | HOLD |
| B4 | every cell pointwise rel width <= 3e-2 (AMEND-P2a-1), no gaps | max 0.0282747 (first cell, set by the 51/50 grid ratio); no gaps = True | HOLD |
| B5 | j = 1 overlaps P1 and REP1 | inside P1's ball [0.0238361051533 +- 3.31e-14] and inside REP1's [0.0238361051533259143533505592, 0.0238361051533259143533844733] | HOLD |
Controls (controls_p2a.json, all_ok, resume_tested true): N1, N2 certified VIOLATION; N4 FAIL through the bound path at j = 2, 10; N-pos (corner and interior bump, with subdivision); N-cell; T6a exact v1 identity at j = 7/3 and 10 (319 terms); T0b REFUSED/TRIVIAL with S == 0, p(0) == 1/40; T0c; RESUME byte-identical (kill by own PID, budget interrupts, sha-change refusal); NaN/inf; --main refused without GO.

Verdict as sealed: P2a CERTIFIED. With P1 (J in [0, 1/4]) this gives chi enclosures for every j in [0, 10.144] (J up to 2.536); the largest certified chi is chi_lo(10.144) = 0.82892676839664, and the certified upper bound on chi over [1, 10.144] is 0.8339408723393343 (last cell [10, 10.144]). (Corrected after REFEREE-P2a fix 5b: an earlier line said "certified spins chi from 0 to 0.8289", which overclaimed - "every chi value in between is attained" needs continuity of j -> M_ADM, not cited.)
B4 disclosure (REFEREE-P2a fix 5a): under the whole-cell reading B4 would be FAILED (0.0385 > 0.03); the scored pointwise reading is the one TARGETS-P2 T3 (sealed b18e905ff1 before any code) already used ("relative width ~ e(j) x ln(ratio)", 0.0283 for cell [1, 1.02]), restored by AMEND-P2a-1.
REFEREE-P2a (code-level, CORRECT): R1-R5, Bernstein (2-D tensor), cells, controls verified; independent solver inside every certified interval at j = 2, 4, 6, 8, 10. For P2b, seal three more controls: exact Chebyshev reconstruction at random rational points, spherical-coordinate residual sign check with l > 2, Bernstein reconstruction at random points.
Not yet: nothing in p2a.py refereed (s-port of the exact residual incl. the tau*Q trick, Bernstein engine, pair checker, cell code); no blind replication of the j > 1 range; Bernstein subdivision never needed in real nodes (only in controls); results reproducible on this PC with single-thread linear algebra. Stage-0 lower bracket is weak at large j (T0c is a sanity check only). Main-quest next: P2b (j to 1000, chi ~ 0.928) and the engine extraction (gravity_engine).
