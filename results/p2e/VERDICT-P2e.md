# VERDICT-P2e (scored as sealed: claim 560847f6a4, BAR file 535ea73d73; instrument p2e.py pinned + sealed, sha a25e6cc4...; engine gravity_engine/bowen_york_sqrt.py module v2, sha 6f3a4099... (unchanged since P2d); `lab.py go` GO; GO-P2e armed; main 4 invocations (exit 3, 3, 3, 0), peak 2.80 GB; result_p2e.json sha 2aa80d5f...)

## Certified: the rounded sixth digit of the supremum spin of the conformally flat Bowen-York puncture family (qualifiers below) is decided - it rounds to 0.928200
For every single, non-boosted, conformally flat, maximal Bowen-York spinning puncture (puncture BC, J = J_ADM) with bare mass mu > 0 (mu = 0 excluded), computer-assisted (exact rationals + outward arb):

  sup over the family of J / M_ADM^2  lies in  [0.92820045389621230212, 0.92820046861667179726]

(lower end rounded down, upper end rounded up; width c - X3 = 1.472e-8). Since c < 0.9282005, the supremum rounds to 0.928200 at six decimals: the numerical literature value 0.928200 (Lovelace et al. 0805.4192 eq (5); whether an earlier paper first printed it is to be confirmed by the release gate) agrees with the correctly rounded six-decimal value of the certified supremum. As noted in VERDICT-P2d, their printed "εJ ≤ 0.928200" read literally as a strict bound is exceeded by certified members (here chi(2^21) >= 0.9282004538, by 4.5e-7) (CORRECTED per REFEREE-P2e 5b). Window: P2c 3.409e-2 -> P2d 5.96e-7 -> P2e 1.47e-8 (about 2.3 million times narrower than P2c).

Cover: [0, 1] by P1 (chi <= 1/4); cells [2^k, 2^(k+1)], k = 0..20, by the tangent lemma from node 2^k; tail [2^21, oo) from node 2^21 (binding). Same lemma basis as P2d (LEMMA-REP4 with referee O fixes; REFEREE-TL; REFEREE-REP4; E1-E3 cited). Supremum, not attainment.

| piece | chi bound (rounded) |
|---|---|
| tail [2097152, oo) (binding) | 0.928200468617 |
| cell [1048576, 2097152] | 0.928200456416 (CORRECTED per REFEREE-P2e 5d: rounded up) |
| cell [524288, 1048576] | 0.928200454324 |

| | sealed bar | result | verdict |
|---|---|---|---|
| E1 | every node's point and tangent certificates certify | 22/22 + 22/22 | HOLD |
| E2 | c - X3 < 1e-7 | 1.472e-8 | HOLD |
| E3 | c < 0.9282005 (rounds to 0.928200) | c = 0.928200468617, 3.14e-8 below | HOLD |
| E4 | c <= P2d's c, X3 >= P2d's X3; shared nodes overlap | both inequalities hold; 19/19 shared nodes overlap (k = 0..17 identical to P2d, k = 18 nested) | HOLD |

Controls (controls_p2e.json all_ok, resume_tested; BUILD-P2e.md): ENGINE 51 tests; psi = 0 exact; structural asserts 10/10; N1, N4, N4-T, tangent bump FAIL at the top node as they should; TC holds at all 21 pairs; N-TAN bites at 14-15 ... 20-21 as tabulated in BUILD-P2e (reproducibility, not blind); RESUME byte-identical, changed sha refused; CEIL and TAIL independent recomputations run INSIDE the instrument on the real result: identical (REFEREE-P2d gap closed); the lead's own Fractions recomputation (<SCRATCH> lead_ceil_p2d.py) on result_p2e.json: identical c, X3, binding piece.

Verdict as sealed: P2e CERTIFIED; E3 HOLD - the rounded sixth digit is decided.

Disclosures: (1) the E3 prediction was stated before any P2e data but after P2d and REFEREE-P2d's post hoc extrapolation (~0.92820044); (2) the builder's degree schedule (K_TOP = 21; N raised to 180/220 for k >= 18) and N-TAN table were chosen after rigorous probes of the same nodes (not blind; validity unaffected; BAR predates them); the builder predicted c = 0.92820046861667179726 exactly and the deterministic run reproduced it; (3) the bar 0.9282005 lies 3.14e-8 above c and 4.61e-8 above X3; the decision does not need a tighter window; to flip it, delta_lo(2^21) would have to be ~260x its certified width too high or the slope ~3.1x its certified bound (REFEREE-P2e 3; CORRECTED wording per 5c); the PDE constants remain a common-mode risk shared with every earlier stage (guarded by K1 = 1/40 and agreement with the independent numerical literature).

REFEREE-P2e (Q, hostile): diff vs p2d.py CORRECT (no change to assembly, formulas, directions); node 2^19 rerun from scratch: 0 of 56 fields differ; independent Fractions c, X3, E3 identical, decision robust to every rounding direction; E4_mock override CORRECT (by design; add to expected_false next time); wording/display fixes applied above; LOG/ledger display slips -> ERRATA-PENDING.

Blind replication of node 2^21: REP6 (0405395ed8) FAILED as sealed (blind delta contained P2e's; blind tangent 2380x looser); REP6b (3e929bb723) REPLICATED the decision: blind chi_tail = 0.9282004720 < 0.9282005 from an independent formulation (quarter-power/log radial map, polyellipse-Cauchy positivity), 2L_hi 7.3e-6, delta nested inside P2e's.

Not yet: release gate to confirm who first printed 0.928200.
