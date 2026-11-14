# REPORT-REP4 (blind replication seat, 2026-11-14) - saved by the lead from the seat's final message (the seat's own Write was refused by the harness); condensed, numbers verbatim

## Disclosures
- The brief gave K1's target 1/40, the mass normalisation (re-derived by the seat) and method hints; no target numbers for delta(3000), delta(10000).
- Opened: DOCTRINE.md; own earlier seat files rep3/REP3-NOTES.md, rep2/REP2-NOTES.md, rep3/rep3.py, head of rep3/rep3_result.json. No forbidden path opened, grepped or listed.
- Indirect view: REP3-NOTES quotes chi_max ~ 0.928, which via chi = j/(4(1+delta)^2) implies delta(3000) ~ 27.4, delta(10000) ~ 50.9 to ~3 digits. Not used in rep4.py (only hard-coded comparison: 1/40 for K1).

## Formulation (differs from the first seat in map, unknown, basis, comparison function, positivity)
- Radial map Y = 1/(1 + sqrt(r/a)), a = 1/2; D = s d/ds = -Y(1-Y) d/dY; (r^2 u_r)_r = (1/4)(D^2 u + 2Du).
- u = Y^2 W, W = P(Y) + Y(1-Y) sum_k q_k(Y) T_k(2v-1), v = cos^2 theta; Chebyshev in 2Y-1.
- Divided residual Gh[W] (E = 2a(1-Y)^2 Phi, CS = 18 a^3 j^2); Y^2 strictly superharmonic: Lin[1] = -3/2; sandwich u_+- = Y^2 (W_N +- eta), eta = (2/3) sup|Gh| (1 + 2^-20).
- Mass: delta = 2a W(0) (re-derived from the ADM surface integral); independent Gauss-law route as a cross-check.
- Residual bound: Q = E_N^7 Lin[W_N] + source, deg <= 8N+44 in Y, 8K in v (hand count); 14 Y-pieces; exact DCT on tensor Chebyshev-Gauss grids in arb; sup|Q| <= sum|c|; inf E_N by c00 - sum'|c|.
- Comparison lemma derived by the seat (Kato + weak maximum principle, H^1_loc cusp at the puncture), NOT refereed; existence/regularity of u assumed (standard, source not opened).

## Results (result_rep4.json; intersection of resolutions A = (140, 16, 224 bits) and B = (200, 24, 288 bits), A contains B)
| quantity | lower (down) | upper (up) | rel width |
|---|---|---|---|
| delta(3000) | 27.42623509173264567580156405725267 | 27.42623509173264567580156540633947 | 4.919e-26 |
| delta(10000) | 50.89803723832470319819716946884112 | 50.89803723832470319820952343563975 | 2.427e-22 |
| K2 delta(10000), source x(1+1e-3) | 50.91100673496594772349367188598665 | 50.91100673496594772350596742973728 | 2.415e-22 |
K1 (j = 2^-64): delta/j^2 in [0.02499999999999907290635182575185424, 0.02500000000000092704915746519705877] contains 1/40. Toy J = 0: delta = 0 exactly.

## Controls
Identity checks all True (DCT, angular identity, radial operator vs sympy, superharmonicity, source transform, Lin[1], arb vs exact, double vs arb); negative control (6W -> 5W) FAILS as it should; exact-rational sample at j = 1e4: Gh[W_+] max -2.39e-13 < 0, Gh[W_-] min +2.39e-13 > 0; nulls (1e-9 bump, K2 source vs honest source) break the signs; Gauss-law route agrees with the slope route to ~1e-28 (B, j = 1e4).

## Runs
Main: 4 resumable invocations (316.7, 310.8, 309.5, 87.4 s), peak 1.44 GB; JSON stage cache, no pickling.

## Weaknesses
Comparison lemma not refereed; degree bounds hand-counted (no extra-node check); rigour rests on python-flint 0.9.0 arb + DCT exactness; observed (unproved, unused) inversion symmetry Y -> 1-Y; superseded logs labelled; one relaunch from an unedited file stopped within seconds by own PIDs.
