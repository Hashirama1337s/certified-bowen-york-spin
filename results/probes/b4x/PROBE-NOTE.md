# B4x PROBE-1: pre-stated, unsealed (listed so its failure is visible)

A double-precision sanity probe run during the no-horizon proof work on the extreme (mu = 0) data. It is not a
certificate and no claim of this repository rests on it.

- Status: **pre-stated, unsealed.** The bar was written in the probe script's header before the run; it was never
  hashed into the lab ledger, so it is not a sealed test and is not counted on the scorecard.
- Bar, as written in the header (verbatim):
  - P1: on the numerical extreme solution (J = 1, a = 8, N = 64, K = 8), Om0(sigma, c) = e^{sigma/2} + g G must be
    strictly increasing in sigma: every forward difference Om0(s + 0.05, c) - Om0(s, c) over s in [-25, 15],
    c in {0, 0.05, ..., 1} must be > +1e-12. FAIL if any difference < -1e-12.
  - P2: A0(sigma) = 2 pi int_{-1}^{1} Om0^4 dc must be strictly increasing on the same grid and must exceed
    A0_cyl = 2 pi int V^4 dc (V = G(-1, c)), with A0(sigma) -> A0_cyl as sigma -> -oo.
  - Null N1: Schwarzschild mu = 1, Om = e^{sigma/2} + (1/2) e^{-sigma/2}; the same check must report negative
    differences for sigma < log(1/2).
- Result (`b4x_mono_probe.log`): **P1 FAILED, P2 FAILED, N1 null OK.** Minimum forward difference -1.29e-10 at
  sigma = -17.35 (2318 entries below -1e-12); minimum dA -9.0e-9; min(A0 - A0_cyl) = -2.4e-7.
- Reading (by the proof seat; the hostile referee then relabelled the probe "pre-stated, unsealed" rather than
  "failed as sealed"): the 1e-12 bar sat far below the numerical solve's own residual (about 1.4e-7), so the probe
  could not test the monotonicity claim. POST HOC (`b4x_mono_posthoc.log`,
  labels only): at N = 104, K = 12 (residual 5.3e-12) the negative differences shrink to -7.6e-15 and move to
  sigma < -20.8; on sigma >= -10 every difference is positive (about 1.9e-8). A re-bar was proposed and dropped by the
  referee (it could not fail). The monotonicity used by the proof is proved, not probed.
- The probe scripts import the project's double-precision throat solver from the lab folder; they are not shipped
  (they certify nothing). The two logs above are the complete outputs.
