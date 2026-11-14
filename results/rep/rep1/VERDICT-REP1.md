# VERDICT-REP1 (scored as sealed: claim e44dd02296, BAR file 06e406719a; instrument rep1/rep1.py pinned after delivery and re-run by the lead: rep1/lead_rerun.log reproduces the seat's numbers)

Fenced blind seat (fence 874a29db). Disclosed exposure: the Read tool returned all of DESIGN.md, so the seat saw DESIGN's double-precision reference table (delta at J = 0.0625, 0.25 to ~10 digits) and the small-spin series; none was an input to rep1.py. Lesson for future blind briefs: keep targets in a separate file, never in the spec file the seat must read.
Independent method: radial y = mu/(r + mu) (scale mu, not mu/2), Fourier cos(2k theta) angular basis, exact rational residual polynomial (fmpq_poly) with one rounding (sum of Chebyshev coefficients), sub/super pair w +- eta y.

| | sealed bar | blind result | P1 certified ball | verdict |
|---|---|---|---|---|
| Q1 | delta(0.25) overlaps P1's ball, rel width <= 1e-6 | [0.0238361051533259143533505592110, 0.0238361051533259143533844732458], rel width 1.42e-21 | [0.0238361051533 +/- 3.31e-14] | HOLD (inside P1's ball) |
| Q2 | delta(0.0625) overlaps P1's ball, rel width <= 1e-6 | [0.00155760570544021599044140515642, 0.00155760570544021599044140535156], rel width 1.25e-25 | [0.0015576057054402 +/- 5.39e-17] | HOLD (inside P1's ball) |
| K1 | delta/j^2 -> contains 1/40 | uniform on j in (0, 2^-16]: [0.02499999998201, 0.02500000001799] at two resolutions | - | PASS |
| K2 | source 9/4 (1 + 1e-3) moves delta(0.25) off P1's ball | shift [2.2758170122437972169e-5, 2.2758170122437972238e-5] >> P1 width 6.6e-14 | - | PASS |

Verdict as sealed: REPLICATED. The blind enclosures lie inside P1's certified balls and are 7 (J = 0.25) to 8 (J = 0.0625) orders tighter.
Certified now, by two independent codes: M_ADM(J = 0.25, mu = 1) in [1.02383610515332591435335055921, 1.02383610515332591435338447325] (21 significant digits); M_ADM(0.0625) to 26 digits.
Seat comments to fold into DESIGN-P2: the margin rule of DESIGN 2.3 is over-conservative (monotone source helps both sides); Lemma C really needs Phi+- >= kappa (1 + mu/2r) on both sides; exact rational residuals make the width truncation-limited (Chebyshev coefficients shrink ~3-5x per degree in y; nearest singularity y = 2).
