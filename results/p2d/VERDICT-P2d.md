# VERDICT-P2d (scored as sealed: claim 9677804a51, BAR file a703963e65; amendments AMEND-P2d-1 ea950c756f, AMEND-P2d-2, AMEND-P2d-3 sealed before data; instrument p2d.py pinned + sealed, sha f1fae733...; engine gravity_engine/bowen_york_sqrt.py module v2, sha 6f3a4099... (POCKETLAB 6a4a597); `lab.py go` GO; GO-P2d armed; main 2 invocations (exit 3, 0), peak 1.23 GB; result_p2d.json sha 5620dd7a...)

## Certified: the maximal spin of the Bowen-York family, first six digits 0.928200... (CORRECTED per REFEREE-P2d 4a: the correctly rounded sixth digit is not yet decided - the window straddles 0.9282005)
For every single, non-boosted, conformally flat, maximal Bowen-York spinning puncture (puncture boundary conditions, J = J_ADM) with bare mass mu > 0 (mu = 0 excluded), computer-assisted (exact rationals + outward arb):

  J / M_ADM^2 <= c = 0.92820096252204785788   (exact rational in result_p2d.json; decimal rounded up)

and the family attains spins at least X3 = 0.92820036654400608722 (certified chi_lo at j = 2^18; decimal rounded down). So

  sup over the family of J / M_ADM^2  lies in  [0.9282003665, 0.9282009626],   width c - X3 = 5.96e-7

(P2c's certified window was [0.928147, 0.962235], width 3.409e-2: about 57,000x narrower; CORRECTED per REFEREE-P2d 4b). Literature: Lovelace et al. 0805.4192 p.2 eq (5) print "εJ ≤ 0.928200" (εJ := J_ADM/E_ADM^2). Read literally, that printed bound is EXCEEDED by a certified family member (chi(2^18) >= 0.9282003665, by 3.7e-7); it is consistent only as a 6-decimal rounding of the limit, which requires sup < 0.9282005 - our window does not yet decide this (CORRECTED per REFEREE-P2d 4c; a correction finding, not a defect of theirs or ours).

Cover: [0, 1] by P1 (chi <= 1/4); [2^k, 2^(k+1)] for k = 0..17 by tangent-lemma cells from node 2^k; [2^18, oo) by the tangent-lemma tail from node 2^18 (binding). It rests on LEMMA-REP4 (Lemmas C-REP4, C-REP4', TL; REFEREE-LEMMA O CORRECT / CWF, fixes applied), REFEREE-TL (N), REFEREE-REP4 (M), E1-E3 cited (Brandt-Bruegmann gr-qc/9703066 p.2; Dain-Gabach Clement 0806.2180 pp.4-5), and the exact engine certificates. The tangent-lemma route does not use DGC Lemma 3.1 (REFEREE-TL 1). Supremum, not attainment: no claim that some member reaches c.

| piece | chi bound (rounded) |
|---|---|
| tail [262144, oo) (binding) | 0.928200962522 |
| cell [131072, 262144] | 0.928200428421 |
| cell [65536, 131072] | 0.928200349297 |
| P1 [0, 1] | 0.25 |

| | sealed bar | result | verdict |
|---|---|---|---|
| D1 | every node's point and tangent certificates certify | 19/19 points, 19/19 tangents | HOLD |
| D2 | c <= 0.935 | 0.928200962522 | HOLD |
| D3 | c - X3 <= 1e-4 | 5.96e-7 | HOLD |
| D4 | (i) shared nodes overlap P1/P2a/P2b/P2c; (ii) monotone-cell consistency; (iii) c <= P2c's c, X3 >= P2c's X3_max | (i) 5/5 overlap; (ii) 10/10; (iii) both hold exactly | HOLD |
| D5 | K1 contains 1/40; K1T consistent with 1 - (3/40) j^2 | both runs (192, 256 bits): 1/40 inside; K1T - 3/40 = -4.66e-17 | HOLD |

Controls (controls_p2d.json all_ok, resume_tested; BUILD-P2d.md): ENGINE 51 tests passed (incl. REP4 A reproduction, FIX F regression); psi = 0 reproduces P2c's X1 and plain cell exactly; N-TAN bites exactly at the pairs AMEND-P2d-3 predicted (14-15, 15-16, 16-17, 17-18) and nowhere it was predicted not to bite was it required; TC holds at all 18 pairs; N1, N4, N4-T, tangent bump all FAIL as they should; structural asserts 10/10 fire; RESUME byte-identical, changed-sha refused; CEIL independent Fractions recomputation identical (controls: mock grid; lead's own recomputation on the real result: identical rational, same binding piece, X3 identical - <SCRATCH> lead_ceil_p2d.py); TAIL identical; NaN/inf INVALID; --main refused without GO; self-sha check.

Verdict as sealed: P2d CERTIFIED.

Disclosures: the builder (not blind) predicted c = 0.92820096252205 and c - X3 = 5.96e-7 from its own rigorous probes at the scheduled settings before main; the run reproduced them (deterministic). The design seat's TARGETS-P2d (sealed) was not opened by the builder.

REFEREE-P2d (P, hostile): assembly CORRECT (cover, directions, safe sqrt(1/2) lower bound); tangent chain CORRECT (re-derived in (r, theta), independent finite-difference check 12/12 at node 10, E_- used throughout); own Fractions recomputation identical c, X3, binding piece; node 2^10 rerun from scratch: every exact field identical; wording fixes 4a-4c applied above. Disclosures (REFEREE-P2d 5): (1) AMEND-P2d-3 was sealed after controls had computed nodes 17-18 and after the builder's rigorous probes at k = 14..18, so the N-TAN "prediction" replays known deterministic numbers (reproducibility, not a blind prediction); (2) the degree schedule was fixed before main but chosen after rigorous probes of the same nodes (not blind; validity unaffected, BAR thresholds predate it); (3) the sealed CEIL control ran on a mock grid; the real-data CEIL checks (lead + referee P) are POST HOC.

REP5 (blind, s-coordinate, sealed 9e0df55e7a) REPLICATED: its own rigorous slope bounds at j = 64 and 4096 lie within 7e-10 and 3.5e-10 of P2d's 2L_hi (bar 1e-3), and its delta(64), delta(4096) lie inside P2d's intervals; an independent secant route (point certificates + convexity of M in mu at fixed J) also brackets the slope. Disclosures: the lead's BAR-REP5 S1 called 2L_hi a true-slope bound; it is TL's tangent coefficient (TL bounds M directly) - the seat caught the direction issue and built the correct one-sided bounds; the seat saw earlier P1-P2c summaries through a stray stdlib-shadowing inspect.py in the lead's <SCRATCH> (no slope values; renamed). REFEREE-P2d + REP5 complete the P2d chain.

Not yet: the eta_Z at the top node (5.0e-5) exceeds the slope itself (~3.4e-5) - it sets the window; a higher top degree would shrink it further (not done: degrees fixed before main). mu = 0 not covered.
