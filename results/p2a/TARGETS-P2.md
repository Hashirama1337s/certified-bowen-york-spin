# TARGETS-P2 (expected values for BHSLICE P2; double-precision mock, NOT proofs)

FENCE: blind seats must not read this file. All numbers here are CHECK level (double-precision collocation, this design
seat, 2026-11-14). Spread = the disagreement across resolutions, an indication only.
Normalisation: bare mass mu = 1, j = 4J, delta = M_ADM - 1, chi = J / M_ADM^2 = j / (4 (1 + delta)^2).
Sources: <LAB>/BHSLICE/p2design/probe_p2a.log (s-map, Ns/NL = 48/10, 64/12, 96/14), probe_throat_s.log (64/10,
96/12, 128/12), probe_throat_far160.log and probe_throat_far224.log (160/10, 224/10), cells.log, targets_calc.log,
exact_feas.log. The DESIGN.md (probe0) values agree where they overlap (j = 2, 4, 10, 30, 100, 300, 1000).

## T1. P2a sealed-probe targets (sharp points)

| j | J | delta = M_ADM - 1 | spread | chi | M_ADM |
|---|---|---|---|---|---|
| 2 | 0.5 | 8.46851040396335e-2 | 3e-15 | 0.4249742853641 | 1.0846851040396 |
| 4 | 1 | 0.250423126410053 | 1.3e-14 | 0.6395669384562 | 1.2504231264101 |
| 6 | 1.5 | 0.42304410711648 | 3.6e-14 | 0.7407207940263 | 1.4230441071165 |
| 8 | 2 | 0.58633161231020 | 4.4e-14 | 0.7947710473436 | 1.5863316123102 |
| 10 | 2.5 | 0.73855534991857 | 6.2e-14 | 0.8271094163818 | 1.7385553499186 |

Proposed B2 tolerance: overlap with [ref (1 - 1e-9), ref (1 + 1e-9)], about 1e4 times the double spread.
Engine feasibility (REP1's exact engine in y, non-target points, exact_feas.log): relative widths of delta at Ny = 56,
K = 10 were 2.3e-17 (j = 3), 3.6e-15 (j = 5), 2.8e-12 (j = 9), at 7 s per point. In s (lower degree needed) expect
equal or better at the same Ny. Proposed B3: relative width <= 1e-10 at all five points.
Non-target midpoints from that run (for a cross-check only): delta(3) = 0.1642320229426064, delta(5) =
0.3375078223848005, delta(9) = 0.6637837692593188.

## T2. Positivity margin (R3)
Not computed by this seat. Expected: u > 0 strictly inside and at the puncture (u(puncture) ~ delta, sec 8 oddity of
DESIGN-P2), u ~ delta / (2r) at infinity, so V = U/s is bounded away from 0. The build seat must report min V before the
seal.

## T3. Monotone-cell widths on [1, 10] (worst case over each grid; cells.log, targets_calc.log)

| grid | nodes | worst rel. width of delta (where) | rel. width near j = 10 | worst abs. width of chi |
|---|---|---|---|---|
| uniform h = 0.1 | 91 | 0.143 (cell [1, 1.1]) | ~0.015 | 0.0058 |
| uniform h = 0.05 | 181 | 0.071 (cell [1, 1.05]) | - | 0.0029 |
| geometric ratio 1.05 (+ targets) | 53 | 0.071 (cell [1, 1.05]) | 0.0047 | 0.0173 |
| geometric ratio 1.02 (+ targets) | 122 | 0.0283 (cell [1, 1.02]) | 0.0027 | 0.0069 |
| geometric ratio 1.01 (+ targets) | 237 | 0.0141 (cell [1, 1.01]) | 0.0020 | 0.0035 |

The relative width is ~ e(j) x ln(ratio), where e(j) = d ln(delta/sqrt j)/d ln j is 1.5 at small j and ~0.5 at j = 10.
Proposed P2a grid: geometric ratio 1.02 plus the targets (122 nodes). Proposed B4: every cell's rel. width <= 0.03.

## T4. Monotonicity sanity (double, 0.05 grid on [1, 10])
delta/sqrt(j) is nondecreasing on the grid: True. chi is increasing on the grid: True.

## T5. P2b and beyond (s-map; resolution spread in brackets)

| j | delta | chi | Chebyshev degree for tail < 1e-13 |
|---|---|---|---|
| 30 | 1.8793295572437 (3e-13) | 0.904645677694 | 37 |
| 100 | 4.2024329573050 (7e-13) | 0.923691663744 | 47 |
| 300 | 7.9937586975160 (2e-13) | 0.927211481646 | 64 |
| 1000 | 15.4131879719 (3e-12) | 0.92801313263 | 86 (97 at Ns = 96) |
| 3000 | 27.42623509173 (3e-12) | 0.9281594195 | 113 |
| 10000 | 50.8980372 (5e-9) | 0.9281926861 | 152 |

Deficit against Lovelace et al.'s 0.928200: 1.87e-4 (j = 1000), 4.06e-5 (3000), 7.31e-6 (10000). Local exponents 1.39
(1000-3000) and 1.42 (3000-10000), consistent with their fitted -1.4 (eq (50)), assuming their limit value.

## T6. Ceiling tail X1: chi(j) <= j1 / (4 delta(j1)^2) for all j >= j1 (with certified delta_- in place of delta)

| j1 | j1 / (4 delta^2) | meaning |
|---|---|---|
| 1000 | 1.0523 | useless (above 1) |
| 3000 | 0.9971 | just below the Kerr bound 1 |
| 10000 | 0.9650 | certified ceiling ~0.965 if a point certificate at j = 1e4 closes |

Cells on [1000, j1] stay far below these values (chi < 0.9282 there, T5).
