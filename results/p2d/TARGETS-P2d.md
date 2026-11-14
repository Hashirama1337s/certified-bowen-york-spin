# TARGETS-P2d (design seat P2d, 2026-11-14) - SEAL THIS FILE; FENCE IT FROM EVERY BLIND SEAT

Reference values the P2d build must later hit. All are DOUBLE-PRECISION model values (CHECK, never proof) unless
marked "rigorous probe". Source logs: <LAB>/BHSLICE/seats/design-p2d/*.log (also fenced). The design file
DESIGN-P2d.md carries no digits beyond the public literature value 0.928200.

Normalisation: mu = 1, j = 4J, delta = M_ADM - 1, chi = j/(4(1+delta)^2), 2L = 1 + 2a Z(0) (a = 1/2), the tangent slope.
Identity check: 2L = 1 + delta - 2 j delta'(j) (finite differences agree to 6 digits up to j = 1e4; probe_a.log).

## T1. Node table, grid j = 2^k, k = 0..18 (eval_pow2.log; sqrt map, N = 160, K = 14, double)

| j | delta | chi | 2L | tangent bound on [j, 2j] (last row: tail [j, oo)) | plain P2c-style bound |
|---|---|---|---|---|---|
| 1 | 0.023836105153 | 0.2384949169 | 9.3280e-01 | 0.4437040654 | 0.467922 |
| 2 | 0.084685104040 | 0.4249742854 | 7.9316e-01 | 0.6881915284 | 0.797532 |
| 4 | 0.250423126410 | 0.6395669385 | 5.5252e-01 | 0.8438573585 | 1.090675 |
| 8 | 0.586331612310 | 0.7947710473 | 3.2493e-01 | 0.8994563772 | 1.195469 |
| 16 | 1.140176167096 | 0.8732949412 | 1.7834e-01 | 0.9175369019 | 1.172177 |
| 32 | 1.970505186728 | 0.9066284789 | 9.6372e-02 | 0.9241074121 | 1.115821 |
| 64 | 3.170604116827 | 0.9198606639 | 5.2041e-02 | 0.9266213732 | 1.064067 |
| 128 | 4.881738203658 | 0.9249933235 | 2.8150e-02 | 0.9275921024 | 1.024486 |
| 256 | 7.309164188346 | 0.9269691752 | 1.5249e-02 | 0.9279664957 | 0.995945 |
| 512 | 10.746126097422 | 0.9277279738 | 8.2672e-03 | 0.9281105831 | 0.975784 |
| 1024 | 15.608924396091 | 0.9280191731 | 4.4840e-03 | 0.9281659549 | 0.961636 |
| 2048 | 22.487152290092 | 0.9281309037 | 2.4326e-03 | 0.9281872166 | 0.951719 |
| 4096 | 32.215082251965 | 0.9281737721 | 1.3198e-03 | 0.9281953774 | 0.944762 |
| 8192 | 45.972803610439 | 0.9281902197 | 7.1612e-04 | 0.9281985090 | 0.939875 |
| 16384 | 65.429350109365 | 0.9281965302 | 3.8857e-04 | 0.9281997107 | 0.936436 |
| 32768 | 92.945165334584 | 0.9281989515 | 2.1084e-04 | 0.9282001717 | 0.934014 |
| 65536 | 131.858460450361 | 0.9281998804 | 1.1440e-04 | 0.9282003486 | 0.932306 |
| 131072 | 186.890200569630 | 0.9282002369 | 6.2075e-05 | 0.9282004165 | 0.931101 |
| 262144 | 264.716850307873 | 0.9282003736 | 3.3683e-05 | **tail 0.9282006089** | 0.935226 |

Model c for the 19-node plan = 0.9282006089 (binding: tail at 2^18). Model X3 (top-node chi) = 0.9282003736.
Tail values if the grid is cut earlier (T(j) = (j/4)/(1 + delta - 2L)^2, model): T(128) ~ 0.93386, T(512) ~ 0.92904,
T(1024) = 0.9285205 (rigorous probe below 0.928520461), T(4096) = 0.9282475, T(16384) ~ 0.928207.

## T2. Further model points (probe_a/b/c.log)

| j | delta | chi | 2L | T(j) |
|---|---|---|---|---|
| 10 | 0.738555349918572 | 0.8271094164 | 2.688804e-01 | 1.1574371 |
| 100 | 4.202432957304864 | 0.9236916637 | 3.502997e-02 | 0.9362576 |
| 300 | 7.993758697516597 | 0.9272114816 | 1.325483e-02 | 0.9299505 |
| 1000 | 15.413187971923843 | 0.9280131326 | 4.578832e-03 | 0.9285311 |
| 3000 | 27.426235091732647 | 0.9281594195 | 1.737061e-03 | 0.9282729 |
| 10000 | 50.898037238324690 | 0.9281926862 | 6.006117e-04 | 0.9282141704 |
| 1e5 | 163.115345137224978 | 0.9282001363 | 7.880707e-05 | 0.9282010277 |
| 1e6 | 517.978203197943685 | 0.9282004454 | 1.034059e-05 | 0.9282004824 |
| 1e7 | 1640.153166747590149 | 0.9282004582 | 1.356845e-06 | 0.9282004597 |

Model sup chi (j -> oo) ~ 0.92820046 (between chi(1e7) and T(1e7)); literature 0.928200 (Lovelace et al. 0805.4192 eq (5)).
REP4 rigorous delta(3000), delta(10000) agree with the T2 rows to all printed digits.

## T3. Small-j anchor (small_j.log)
(1 - 2L)/j^2 -> 3/40 (= 0.075): 0.07499999122 at j = 1e-3; delta/j^2 -> 1/40.

## T4. Rigorous feasibility probes (tl_cert_probe.py, imports rep4.py read-only; NOT instruments)

| j | N, K, prec | delta width (rel) | sup|H| bound | eta_Z | 2L_hi | chi_lo(j) | T(j) | time | peak |
|---|---|---|---|---|---|---|---|---|---|
| 1024 | 60, 10, 192 | 1.55e-7 | 3.07e-5 | 2.05e-5 | 4.504487e-3 | 0.928019038 | 0.928522887 | 22 s | 0.73 GB |
| 1024 | 90, 12, 192 | 1.02e-12 | 2.57e-9 | 1.71e-9 | 4.484005e-3 | 0.928019173 | 0.928520461 | 42 s | 0.92 GB |
| 4096 | 100, 12, 192 | 2.58e-11 | 2.05e-8 | 1.37e-8 | 1.319847e-3 | 0.928173772 | 0.928247541 | 48 s | 1.02 GB |
| 262144 | 128, 14, 192 | 5.04e-8 | 1.20e-3 | 7.99e-4 | 8.329e-4 | 0.928200327 | 0.928206240 | 82 s | 1.32 GB |

The 2^18 row is limited by the 14 fixed pieces (inf E^8 per piece; see cert_4096.log per-piece table), not by the
collocation residual (Z refinement reaches 2e-54 at the nodes, cert_2p18b.log). With adaptive pieces the build should
approach the T1 model row: c in [0.9282006, 0.9282063] for the 19-node plan, X3 >= 0.9282003.
Z_+ < 0 at all sample nodes (max about -0.996 at j = 1024, -0.999 at 2^18); G1 node minimum = 1.000000 (at Y = 0).

## T5. Bars the lead may seal (suggested)
- every node 2^k, k = 0..18: delta certificate AND tangent certificate AND G1 > 0 certified;
- c (exact rational) <= 0.9283 (tier C bar) ; stretch bar c <= 0.928201;
- X3 = max chi_lo >= 0.92820;
- delta at j = 1, 2, 4, 8 overlaps P2a's certified intervals; at j = 3000/10000 if added, overlaps P2c and REP4;
- tangent consistency (TC) holds at all 18 consecutive pairs; K1T: (1 - 2L_hi)/j^2 at j = 2^-64 within 1e-9 of 3/40.
