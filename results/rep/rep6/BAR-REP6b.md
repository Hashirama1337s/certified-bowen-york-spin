# BAR-REP6b: tighter blind tangent at j = 2^21 (after REP6 FAILED on R2/R3)

Written 2026-11-14 by the lead after REP6 (0405395ed8) was settled FAILED: blind delta(2^21) contained P2e's interval (R1 HOLD) but the blind TL coefficient 2L_hi = 0.01407 was 2380x P2e's, so the blind tail 0.928235 did not reproduce the decision. Same blind seat, same fences (re-armed), same formulation (rep6), improved tangent certificate (the seat named the limits: odd residual from D(Rhat_N); V1 from u_-). The seat is told only that its 2L_hi was far too loose to decide the bar 0.9282005 - not P2e's values.

instruments: rep6/rep6b.py

## Predictions
- R2b: blind 2L_hi(2^21) <= 3 x P2e's 2L_hi(2^21).
- R3b: blind chi_tail(2^21) = (j/4)/(1 + delta_lo_blind - max(2L_hi_blind, 0))^2 < 0.9282005, with delta_lo_blind from REP6 (unchanged) or a re-certified one.

## Controls (any failure = VOID)
- KT (psi = 0 gives 2L = 1) through the machinery, not by hand.
- KZ: the improved certificate with eta_Z forced to (2/3) sup (1 - 2^-10) FAILS.
- The REP6 point certificate (delta) unchanged or re-certified and still overlapping P2e.

## Verdict rule
REPLICATED (decision) if R2b, R3b hold and controls pass. Otherwise FAILED again, reported as is; the P2e decision then stays single-engine on the tangent at the top node (disclosed in any release text).

## Seal fields
- formula_ref: LEMMA-REP4.md:177-230
- calibrated_on: 0405395ed8
- depends_on: 560847f6a4, 0405395ed8
- null_not_symmetry: KZ forces a too-small margin (must fail).
- target_form: R2b ratio in [0, 3]; R3b chi_tail in [0, 0.9282005).
- ceiling: P2e window [0.9282004538, 0.9282004687].
- toy_model: J = 0 gives 2L = 1.
- tolerances: from need; 8 GB; one process; chunks < 10 min; about 45 min compute.
- parameters: j = 2^21
- files: BAR-REP6b.md, rep6/rep6b.py
- command: memcap 8 python rep6/rep6b.py
