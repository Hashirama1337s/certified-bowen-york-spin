# FLUXT-REPORT: throat flux E(sigma) on the 2-D float solution (written by fluxt.py)

- Instrument: fluxt.py sha256 `957929cc1d4f3111f1f3d4ec33ffb0df7bed5e94f8465a9aca03871a10fa0c27`; b3t.py sha256 `0043585b666fd0f5b0afd2ab4f16dc98335600834943c6a8bd7fdee13aed679a` (imported, float stage).
- Verdict: **FLUXT HOLD**

## Tolerance (stated before any slice)
- sup|Q| (sampled) = 6.323830e-13; W = 4.0518922302 (doubled nodes 4.0518922302); tol = 2.562348e-12.
- Float solution: 2 g0 = 1.0379563989136513 (E(+oo) = -pi 2 g0 on this solution).

## Slices (nc = 64)

| sigma | E | E/(-pi) - M_0 | E (grad coef 1/4) |
|---|---|---|---|
| -6.0000 | -3.2608361975736457 | +6.661e-16 | -3.2537386589682384 |
| -4.6667 | -3.2608361975736453 | +4.441e-16 | -3.2537742348147964 |
| -3.3333 | -3.2608361975736457 | +6.661e-16 | -3.2539899626287805 |
| -2.0000 | -3.2608361975736444 | +2.220e-16 | -3.2551194634724565 |
| -0.6667 | -3.2608361975736431 | -2.220e-16 | -3.2584607481024435 |
| +0.6667 | -3.2608361975736457 | +6.661e-16 | -3.2606958354554720 |
| +2.0000 | -3.2608361975736466 | +8.882e-16 | -3.2608353277704221 |
| +3.3333 | -3.2608361975736408 | -8.882e-16 | -3.2608361957869580 |
| +4.6667 | -3.2608361975736417 | -6.661e-16 | -3.2608361975710500 |
| +6.0000 | -3.2608361975736435 | +0.000e+00 | -3.2608361975736408 |

## Scorecard

| test | value | bar | result |
|---|---|---|---|
| F1 spread | 5.773e-15 (0.002 tol) | <= 10 tol | HOLD |
| F2 max abs(E/(-pi) - M_0) | 8.882e-16 (0.000 tol) | <= 10 tol | HOLD |
| F3 null spread (1/4) | 7.098e-03 (2.770e+09 tol) | > 100 tol | DRIFTS |
| QUAD doubling | 7.994e-15 (0.003 tol) | < tol | ok |
| RN analytic | max dev 4.441e-16 | <= 1e-13 | ok |

## Diagnostic (not a bar control)
- RN through the ansatz path (b3t RN float solve N=64 K=0): max abs(E + pi) = 4.374e-13, its own tol 3.208e-11, within: True.

## Notes
- tol is a float tolerance from a sampled residual sup (not rigorous).
- On the float solution E(+oo) = -pi 2 g0 exactly, so F2 tests 2 g0 (and the slices) against the 1-D flux value M_0.
