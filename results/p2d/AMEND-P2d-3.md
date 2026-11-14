# AMEND-P2d-3: expected N-TAN gaps at every node pair (build seat, 2026-11-14), written BEFORE any main run

Amends BAR-P2d control N-TAN and AMEND-P2d-1 item 7. NOT sealed by the build seat (the lead seals it). TARGETS-P2d.md was
not opened. Sources: the build seat's own probes only (<SCRATCH> p2dbuild/dprobe.py = double model, CHECK;
p2dbuild/tprobe.py = rigorous runs of the same engine code at the scheduled (N, K, prec) of p2d.py, k = 14..18).

## The test and its resolution
TC at pair (k, k+1), mu1 = 1 units at node j_k = 2^k, x = 1/sqrt(2):
    TC holds  <=>  P := 1 + delta_lo(k) - 2L_hi(k) + x 2L_hi(k)  <=  x (1 + delta_hi(k+1))   (exact; p2d.py tc_pair).
Gap g_k := x (1 + delta_hi(k+1) - 2L_hi(k)) - (1 + delta_lo(k) - 2L_hi(k)) (M units of node k).
N-TAN lowers 2L_hi(k) by s = 1e-4 (no clamping): (1 + delta_lo - t) rises by s and x (1 + delta_hi - t) by x s, so
the gap drops by exactly (1 - x) s = 2.9289e-5. N-TAN bites (TC becomes false) iff g_k < 2.9289e-5. The test's own
resolution is exact (rational arithmetic, no rounding), so only the gap size decides.
The certified gap is the model gap (second-order curvature of M(mu)) plus the certificate slack
(1 - x)(2L_hi - 2L) + (delta - delta_lo)(k) + x (delta_hi - delta)(k+1); the slack is negligible at k <= 13
(eta_Z <= ~6e-7, widths <= 1e-8) and is included exactly at k = 14..18 (rigorous probes).

## Expected gaps (g_k) and verdict per pair
| pair (k, k+1) | expected gap g_k | source | N-TAN |
|---|---|---|---|
| 0-1 | 1.64e-2 | rigorous mock run (24, 4) 1.6394e-2; model 1.6363e-2 | cannot bite |
| 1-2 | 3.18e-2 | model 3.1808e-2 | cannot bite |
| 2-3 | 3.31e-2 | model 3.3113e-2 | cannot bite |
| 3-4 | 2.22e-2 | model 2.2171e-2 | cannot bite |
| 4-5 | 1.25e-2 | model 1.2523e-2 | cannot bite |
| 5-6 | 6.78e-3 | model 6.7840e-3 | cannot bite |
| 6-7 | 3.66e-3 | model 3.6552e-3 | cannot bite |
| 7-8 | 1.97e-3 | model 1.9732e-3 | cannot bite |
| 8-9 | 1.07e-3 | model 1.0675e-3 | cannot bite |
| 9-10 | 5.78e-4 | model 5.7836e-4 | cannot bite |
| 10-11 | 3.14e-4 | model 3.1359e-4 | cannot bite |
| 11-12 | 1.70e-4 | model 1.7010e-4 | cannot bite |
| 12-13 | 9.23e-5 | model 9.2282e-5 | cannot bite |
| 13-14 | 5.0e-5 | model 5.0070e-5 (+ slack < 1e-6) | cannot bite (gap 1.7x the shift effect) |
| 14-15 | 2.754e-5 | rigorous 2.7542e-5 (model 2.7167e-5) | bites, margin 6 % (gap/effect 0.940) |
| 15-16 | 1.735e-5 | rigorous 1.7349e-5 (model 1.4741e-5) | bites (0.592) |
| 16-17 | 8.09e-6 | rigorous 8.0928e-6 (model 7.9985e-6) | bites (0.276) |
| 17-18 | 6.26e-6 | rigorous 6.2627e-6 (model 4.3400e-6) | bites (0.214) |

Model = double collocation of the point and tangent problems at (N, K) = (64, 10) for k <= 8 and (160, 16) for k >= 8
(coefficient tails <= 1e-14), gap computed from delta and 2L directly (CHECK, not a bound). Rigorous = engine
certificates (point + tangent) at exactly the degrees p2d.py's schedule uses for those nodes: k = 14, 15 at (110, 12, 192),
k = 16, 17, 18 at (140, 14, 224), REP4's 14 fixed pieces. These runs are deterministic (single-thread fit, exact
arithmetic), so the main run is expected to reproduce them.

## Pre-registered N-TAN statement
- Predicted to bite: pairs 14-15, 15-16, 16-17, 17-18 (p2d.py NTAN_EXPECTED_BITES; result key
  `ntan_bites_where_predicted`). Pair 14-15 has only a 6 % margin.
- Cannot bite (listed in advance, as BAR-P2d asks): pairs 0-1 through 13-14. There TC passes with 2L_hi lowered by 1e-4,
  because the curvature gap is larger than the shift effect.
- TC itself is expected to hold at all 18 pairs (gaps > 0).
- The controls already exercise N-TAN on the real top pair 17-18 (controls_p2d.json parts.node.NTAN_17_18).
