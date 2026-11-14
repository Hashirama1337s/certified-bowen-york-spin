# VERDICT-P2b (scored as sealed: claim eba765ebe0, BAR file d5e41fe748; instrument p2b.py pinned + sealed 571c21859b; engine gravity_engine/lichnerowicz_axisym.py (sha be006193...), first reuse; `lab.py go` GO (resume tested); main 3 invocations (runner convention: exit 0 + "time budget reached" = rerun; the lead's first loop stopped after chunk 1 and was restarted - no data affected), peak 2.36 GB in the last invocation; maximum 2.52 GB in invocation 2 - corrected per REFEREE-P2b 4)

## Certified (Bowen-York spinning puncture, mu = 1, EVERY spin j = 4J in [10, 1030.35])
100 point certificates (21/20 grid from 10 plus j = 30, 100, 300, 1000) and 99 monotone cells, no gaps; degrees (N, K) from (36, 10) to (112, 14).

| j | certified delta = M_ADM - 1 | rel. width |
|---|---|---|
| 10 | [0.738555349918368455, 0.738555349918776234] | 5.52e-13 |
| 30 | [1.879329556666385383, 1.879329557821432294] | 6.15e-10 |
| 100 | [4.202432956331064025, 4.202432958278665085] | 4.63e-10 |
| 300 | [7.993758667606936304, 7.993758727426259014] | 7.48e-9 |
| 1000 | [15.41318755062711644, 15.41318839322056450] | 5.47e-8 |
Physics outputs: X3 certified lower bound on the Bowen-York spin supremum: chi_lo(j = 1030.35) = 0.9280206464756176 (numerical ceiling 0.928200, Lovelace et al.); certified upper bound on chi over [10, 1030.35]: max cell chi_max = 0.9326738205130032 (cell [153.67, 161.36]; set by the 21/20 cell width, not by the point certificates).

| | sealed bar | result | verdict |
|---|---|---|---|
| B1 | every node certifies | 100/100 | HOLD |
| B2 | delta overlaps T5 windows ref(1 +- 1e-8) at j = 30, 100, 300, 1000 | all four T5 references (1.8793295572437, 4.2024329573050, 7.9937586975160, 15.4131879719) lie inside the certified intervals | HOLD |
| B3 | rel width <= 1e-5 every node | max 7.805e-8 (j = 1030.35) | HOLD |
| B4 | pointwise cell width <= 5e-2, no gaps | max 0.02404 (first cell [10, 10.5]); no gaps | HOLD |
| B5 | j = 10 overlaps P2a's delta(10) | contains P2a's interval [0.738555349918561234, 0.738555349918583461] | HOLD |
Controls (controls_p2b.json, all_ok, resume_tested): ENGINE tests 11 passed; N1 VIOLATION at j = 100, 1000; N4 FAIL through the bound path; REF Chebyshev reconstruction, spherical-coordinate residual signs with l >= 4 (j = 300, 1000), Bernstein reconstruction; T0 trivial; T0c; RESUME byte-identical; MEM 2.2 GB; NaN/inf; --main refused without GO.

Verdict as sealed: P2b CERTIFIED. With P1 + P2a + P2b: chi enclosures for every j in [0, 1030.35]; the largest certified lower bound is chi(1030.35) >= 0.92802064647; the certified chi(1030.35) lies in [0.92802064647, 0.92802078265] (CORRECTED per REFEREE-P2b 3). The cover "every spin in [10, 1030.35]" rests on Dain-Gabach Clement Lemma 3.1 + scaling (monotone cells). Not shown: that every chi value in [0, 0.92802] is attained (needs continuity of j -> delta, not cited).
Not yet: p2b.py runner and the engine's use beyond j = 10.144 not refereed; no blind replication in [10, 1000]; the spherical reference check's teeth are weak at large j (source errors below ~1e-7 do not flip its sign; BUILD-P2b); the engine fit bug for N > ~110 (P2C-SCOUT) did not stop any node here (CORRECTED per REFEREE-P2b 2: at N = 112 the arb refinement plateaus at ~1e-26 instead of 2^-120 and only the double residual is logged; harmless because the width is set by truncation in N - a float64-rounded fit gives the same width; damaged fits only widen or fail, never certify falsely) but should be fixed before P2c.

REFEREE-P2b (code-level, K): certificate path CORRECT at high degree (exact proof path; damaged-fit tests only widen; independent residual at j = 1000 agrees to 1e-68 S at 45 points; own solve delta(1000) = 15.4131879719239 inside). Before P2c: engine fit fix (midpoint/exact arb operators) and a distinct budget exit code.
