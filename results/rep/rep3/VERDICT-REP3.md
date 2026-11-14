# VERDICT-REP3 (scored as sealed: claim a669f700df, BAR file faace3c1c4; instrument rep3/rep3.py pinned after delivery and fully re-run from empty stage cache by the lead: all certified values identical, only timing fields differ)

Fenced blind seat (fence 54136d00; gravity_engine fenced too). Disclosed exposure: DESIGN 1.3 series, the published ceiling ~0.928 (implies the rough size of delta at large j), DESIGN-P2 degree hints, REP2's j = 6, 10 values; no delta(300) or delta(1000) seen before its own run. Independent method: y = a/(r + a) with a = 9/20, Chebyshev in v = cos^2 theta, piecewise certificate on 21 y-pieces (discrete Chebyshev transform of Q per piece in arb, sup|Q| / inf D^7), Chebyshev-form positivity.

| | sealed bar | blind result | P2b certified | verdict |
|---|---|---|---|---|
| Q1 | delta(1000) overlaps P2b, rel width <= 1e-6 | [15.4131879719238393605151772769, 15.4131879719238420021732318176], 1.71e-16 | [15.41318755062711644, 15.41318839322056450] | HOLD (inside) |
| Q2 | delta(300) overlaps P2b, rel width <= 1e-6 | [7.99375869751659771889919847622, 7.99375869751659826345215360813], 6.81e-17 | [7.993758667606936304, 7.993758727426259014] | HOLD (inside) |
| K1 | small-j delta/j^2 contains 1/40 | [0.0249999999999999999999999813, 0.0250000000000000000000000187] on (0, 2^-40] | - | PASS |
| K2 | source 9/4 (1 + 1e-3) moves delta(1000) off P2b | shift [4.100615126796e-3, 4.100615126801e-3] >> P2b width 8.4e-7 | - | PASS |
Verdict as sealed: REPLICATED.
Certified by two independent codes: chi(1000) in [0.928013132630170991794, 0.928013132630171290517] and chi(300) in [0.927211481646314848568, 0.927211481646314960850].
