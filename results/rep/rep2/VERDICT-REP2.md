# VERDICT-REP2 (scored as sealed: claim 227cbcafc2, BAR file 0ab09c8b28; instrument rep2/rep2.py pinned after delivery and re-run by the lead: rep2/lead_rerun.log reproduces the seat's numbers incl. K2)

Fenced blind seat (fence d2ee9002). Disclosed exposure: DESIGN 1.3 series, DESIGN 2.2 stage-0 numbers at j = 0.5, 1, 4, REP1's j <= 1 enclosures, DESIGN-P2 hints (degrees, log-slope ~0.5 at j = 10); no value at j = 6 or 10 seen before its own run. Independent method: y = 1/(1 + r), angular v = cos^2 theta with T_k(2v - 1), exact rational residual, positivity of w_-/y by adaptive bisection with exact centred Taylor forms (not Bernstein), margin 64 S (1 + 2^-20).

| | sealed bar | blind result | P2a certified | verdict |
|---|---|---|---|---|
| Q1 | delta(10) overlaps P2a, rel width <= 1e-8 | [0.738555349918572347887568348426, 0.738555349918572347890115434849], 3.45e-21 | [0.738555349918561234, 0.738555349918583461] | HOLD (inside) |
| Q2 | delta(6) overlaps P2a, rel width <= 1e-8 | [0.423044107116482969357400750737, 0.423044107116482969357401280516], 1.25e-24 | [0.423044107116482903, 0.423044107116483035] | HOLD (inside) |
| K1 | small-j delta/j^2 contains 1/40 | uniform on (0, 2^-24]: [0.02499999999999967, 0.02500000000000033] at two resolutions | - | PASS |
| K2 | source 9/4 (1 + 1e-3) moves delta(10) off P2a | shift [3.672954695342064601705e-4, 3.672954695342064652807e-4] >> P2a width 2.2e-14 | - | PASS |
Verdict as sealed: REPLICATED. Also: j = 1 cross-check overlaps REP1; toy J = 0 exact (S = 0, refused as trivial); wrong-ansatz nulls inflate eta.
Certified by two independent codes: chi(10) in [0.827109416381791178798080260634, 0.827109416381791178800503789071]; M_ADM(10) = 1.738555349918572347888... (21 digits).
