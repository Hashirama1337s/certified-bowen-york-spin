# BAR-B2V: certified-data tests of B2's theorems (monotone, convex, tangent-free tail) and Conjecture E (exponent nu = 1.382019136)

Written 2026-11-14 by the lead before b2_verdict.py is run on any real data (seat self-test only, on synthetic families). Source: scout/B2-PROOF.md section 6 (proof seat, unrefereed). Data: result_p2e.json (sha 2aa80d5f...), 22 nodes j = 2^k, k = 0..21.

instruments: scout/b2_verdict.py

## Predictions
- V1 (Theorem B, convexity): all 20 three-node chord inequalities chord_k > chord_{k+1} are decided TRUE; none decided reversed (undecided pairs reported, not failed).
- V2 (Theorem A, D1): all 21 chords lie in (0, 1) and chi strictly increases on all 21 consecutive pairs (none decided outside).
- V3 (D6, tangent-free tail): the chord tail T4 from nodes (2^20, 2^21) is < 0.9282005, and T4 - 0.92820046861667179726 lies in [0, 2e-8].
- V4 (Conjecture E): at the highest pair whose interval exponent has width < 2e-4: |e - 0.88201913609| < 2e-3; |e - 0.91421356237| > 0.02 (spherical rival); |e - 0.79642947556| > 0.05 (p = 6 null); and |e_top - 0.882019| < |e_first(k>=6) - 0.882019| (approach improves with j).

## Controls (seat self-test, b2_verdict_selftest.log sha fb328733...): positive m(x) = sqrt(x^2 + 1) passes T1-T3 and T5 -> 0.5; nonconvex null fails T2 at indices [0, 4] and T1_lt1 = False. Re-run on this machine before scoring.

## Verdict rule
Each of V1-V4 scored separately as HOLD / FAILED. V1-V3 failing would contradict a claimed theorem (stop and find out why); V4 failing kills Conjecture E's numeric form only.

## Seal fields
- formula_ref: scout/B2-PROOF.md:1-400
- calibrated_on: 560847f6a4
- depends_on: 560847f6a4
- null_not_symmetry: nonconvex synthetic family; exponent rivals 0.91421 and 0.79643
- target_form: V1/V2 counts; V3 T4 < 0.9282005 and gap in [0, 2e-8]; V4 inequalities as stated
- ceiling: P2e window; Lovelace et al. fitted exponent "-1.4"
- toy_model: sqrt(x^2 + 1) (extreme-RN-like) family
- tolerances: exact rationals + mpmath intervals at 80 digits
- parameters: P2e nodes k = 0..21
- files: BAR-B2V.md, scout/b2_verdict.py
- command: memcap 4 python scout/b2_verdict.py result_p2e.json
