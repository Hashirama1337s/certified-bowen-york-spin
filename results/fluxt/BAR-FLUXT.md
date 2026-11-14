# BAR-FLUXT: the throat flux E(sigma) is constant on an independent 2-D solution

Written 2026-11-14 by the lead before fluxt.py exists (REFEREE-FLUX single next strike). Identity: scout/B5-DIGITS.md sec. 1 (REFEREE-FLUX: CORRECT / CORRECT WITH FIXES). The check uses the 2-D extreme-data solution of the B3T float solve (b3t.py, sealed; float stage only, re-run, not edited) - a solution computed with no knowledge of the flux identity.

instruments: fluxt.py

## Object
E(sigma) = Int_{S^2} [Omega_s^2/2 - |grad_S Omega|^2/2 - Omega^2/8 - (3/8) sin^2(theta) Omega^-6] dOmega, Omega = r^(1/2) Phi, evaluated at 10 slices sigma in [-6, 6] on the 2-D float solution (J = 1), by high-order quadrature in theta.

## Predictions
- F1: E(sigma) is constant across the 10 slices within 10x the solution's own residual-propagated tolerance (stated by the instrument before evaluation from the float residual sup).
- F2: the constant equals -pi M_0 with M_0 from the flux formula (2/3) Sum (2l+1) a_l^2 (scout/b5 probe value 1.03795639891365137036534431609), within the same tolerance.
- F3 (null): the same evaluation with the gradient coefficient 1/2 replaced by 1/4 DRIFTS across slices by more than 100x the tolerance.

## Controls (any failure = VOID)
- Extreme RN analogue (power 3): E constant = -pi q exactly (analytic Omega^2 = e^s + q).
- Quadrature convergence: doubling the theta nodes changes E by less than the tolerance.

## Verdict rule
FLUXT HOLD if F1, F2 hold and F3 drifts as stated. F1 failing = the identity or the solution is wrong: stop.

## Seal fields
- formula_ref: scout/B5-DIGITS.md:1-60
- calibrated_on: 10a65e940a
- depends_on: 10a65e940a, 133b99d241
- null_not_symmetry: F3 changes the flux density (must drift); RN exact.
- target_form: F1 spread in [0, 10 tol]; F2 |E/(-pi) - M_0| in [0, 10 tol]; F3 spread > 100 tol.
- ceiling: B3T window and REP7 window for chi(0).
- toy_model: extreme RN (exact).
- tolerances: from the float solution's residual; stated by the instrument before evaluation.
- parameters: 10 slices sigma in [-6, 6]
- files: BAR-FLUXT.md, fluxt.py
- command: memcap 8 python fluxt.py
