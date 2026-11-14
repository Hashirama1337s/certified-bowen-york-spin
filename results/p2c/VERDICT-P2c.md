# VERDICT-P2c (scored as sealed: claim 3031217ac8, BAR file e2da2198e4; instrument p2c.py pinned + sealed ca9e5cb837, sha 9d2452bd...; engine lichnerowicz_axisym v2 sha 3659ae21... unchanged; `lab.py go` GO; GO-P2c armed by the lead; main 3 invocations, exit codes 3, 3, 0; peak 5.01 GB (invocation 2); ceiling assembler run once on result_p2c.json sha 4b2dbf7b...)

## Certified: a uniform spin ceiling for the whole Bowen-York spinning-puncture family
For every single, non-boosted, conformally flat, maximal Bowen-York spinning puncture (puncture boundary conditions) with bare mass mu > 0 (mu = 0 excluded) and any angular momentum J = J_ADM (computer-assisted: exact rationals + outward arb; CORRECTED wording per REFEREE-P2c 6a, 6e):

  J / M_ADM^2 <= c = 0.96223488883111534543... (exact rational in ceiling_p2c.json; decimal rounded up)

The cover: [0, 1] by P1 (chi <= j/4 <= 1/4); [1, j1] by the monotone cells of P2a, P2b and P2c; [j1, oo) by the tail X1 at j1 = 11994.815017806193. Chain links P1 -> P2a -> P2b -> P2c all exact (chain_ok). It rests on existence and uniqueness of the solution for every mu > 0 (E2, cited) and J_ADM = J (cited) (CORRECTED per REFEREE-P2c 6b), the comparison lemma (refereed, P1), Dain-Gabach Clement Lemma 3.1 + scaling (delta/sqrt(j) nondecreasing; REFEREE-P2DESIGN; Lean-checked cell inequality 036ea0f37b) and the exact engine certificates. mu = 0 is not covered.

| piece | max chi bound (rounded up) |
|---|---|
| P1 range [0, 1] | 0.25 |
| P2a cells (worst [10, 10.144]) | 0.83394087233933423687 |
| P2b cells (worst [153.67, 161.36]) | 0.93267382051300310048 |
| P2c cells (worst [1030.35, 1287.93]) | 0.93989487162744096378 |
| X1 tail, j >= 11994.82 | **0.96223488883111534543 (binding)** |

New point certificates (14 nodes, (5/4) grid from P2b's last node, plus j = 3000, 10000; degrees N = 112-192, K = 12):

| j | certified delta = M_ADM - 1 (mu = 1) | rel. width |
|---|---|---|
| 3000 | [27.42520280038156340, 27.42726738308364238] | 7.53e-5 |
| 10000 | [50.89602943727101564, 50.90004503937838419] | 7.89e-5 |
| 11994.82 (top) | [55.82467866972199620, 55.85355114628200546] | 5.17e-4 |

(CORRECTED per REFEREE-P2c 4: upper ends were copied from the engine log line, which truncates (lichnerowicz_axisym.py:852 slices delta_hi_dec[:20]); now rounded up. The stored exact rationals were always right.)

Lower side: the largest certified lower bound on the family's spin is chi_lo(9595.85) = 0.92814657060579 (X3 at the top node: 0.92772309989694600). So the supremum of chi over the family (mu > 0) lies in [0.928146, 0.962235]; attainment is not shown; certified uncertainty c - X3(top) <= 0.03452 (c - X3_max <= 0.03409) (CORRECTED per REFEREE-P2c 4, 6c).

| | sealed bar | result | verdict |
|---|---|---|---|
| C1 | every new node certifies | 14/14 | HOLD |
| C2 | c < 1 | c = 0.962234888831115 | HOLD |
| C3 | c <= 0.98 | 0.962235 | HOLD |
| C4 | delta at j = 3000, 10000 overlaps T5 ref (1 +- 1e-6) | refs 27.42623509173 and 50.8980372 lie inside the certified intervals | HOLD |
| C5 | X3 and c - X3 reported | X3(top) 0.927723, X3_max 0.928147, c - X3 0.03451 | HOLD |

Controls (controls_p2c.json, all_ok, resume_tested; BUILD-P2c.md): ENGINE tests pass (28; the BAR says 17 - engine v3 added 11 K3 tests on 2026-11-14 without touching lichnerowicz_axisym; deviation in count only); CHAIN first-node string equal, modified P2b copy refused; N1 VIOLATION at j = 1e4; N4 FAIL through the bound path at the top degree; REF Chebyshev (S by a second route equal), spherical (6/6 signs, l >= 4 present), Bernstein 6/6; T0 trivial; TAIL X1 = mpmath-120; CEIL independent mpmath-120 recomputation, same c and binding piece (on a mock result before main - see Not yet); RESUME byte-identical, changed sha refused; MEM peak 4.58 GB; NaN/inf; --main refused without GO.

Verdict as sealed: P2c CERTIFIED.

Disclosures: (1) the build seat (not blind) read LOG's "c <= 0.98 (~0.975)" estimate while gathering context; the expected value it reported (0.962234888831...) came from its own probes and matches the run exactly. (2) Print slip in the run log: delta(10000)'s upper end printed ...418 here and ...419 in the build seat's probe; same exact rational expected - checked at referee step (display rounding, not a different certificate, to be confirmed).

REFEREE-P2c (L, hostile, strongest tier): ceiling argument CORRECT (tail and cell bounds re-derived); chain CORRECT (no gap in j, 237 records); independent exact-rational c CORRECT (identical rational, same binding piece, same per-piece maxima); disclosure (2) resolved: same certificate, display truncation (fixed above; P2a/P2b tables have the same pre-existing pattern, logged for the release errata); runner CORRECT WITH FIXES (self-sha check in --main, P2a->P2b first-node assert: next instrument); wording fixes applied. Not covered by L: R1-R5 logic (inherited from REFEREE-P2DESIGN/P2a/P2b); Lean 036ea0f37b not opened.

REP4 (blind, sealed ad62cca292 / f6c18cffaf) REPLICATED: independent rigorous delta(10000) = 50.898037238324703198... (rel width 2.4e-22) and delta(3000) = 27.42623509173264567580156... (4.9e-26), both inside P2c's intervals; K1 contains 1/40; K2 moves off; different radial map (sqrt), basis, comparison function and positivity method. Disclosure: the seat saw chi_max ~ 0.928 in its own earlier notes (implies ~3 digits). Its comparison lemma is its own derivation, not refereed.

Not yet: no claim that every chi value up to c is attained; the ceiling is set by the tail estimate (about 3.6 % above the numerical estimate 0.928200, Lovelace et al. 0805.4192 eq (5); CORRECTED per REFEREE-P2c 6d); tightening c needs a larger j1 or a sharper tail lemma.
