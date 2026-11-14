# CLAIMS: every claim with its status, and every number with its source

Status labels: **CERTIFIED** = computer-assisted enclosure (exact rationals + outward ball arithmetic) from a sealed run;
**PROVED** = written proof checked by at least one hostile referee seat; **CITED** = taken from an opened source (locator
in paper.md, References); **PROBE** = numerical value without an enclosure, or an enclosure from an unsealed run.
Seal ids are the lab ledger's claim hashes; referee seats are named by their letter (H ... V). Source files are the lab
file names; the repository copies them under `results/` and `code/`.

## Claims

| # | Claim (paper section) | Status | Seals | Referee seats, blind seats | Novelty gate |
|---|---|---|---|---|---|
| C1 | Sub/super-solution pair and δ, χ enclosures for every J in [0, 1/4] at μ = 1 (j in [0, 1]); "every J" at slope-route accuracy about 2e-11 relative (§3) | CERTIFIED | 565736c4c4 | referee H (CORRECT, fixes A, B, N4 folded into P2); blind REP1 e44dd02296 | OPEN (GATE-NOVELTY A1), wording "to our knowledge" |
| C2 | δ/j² → 1/40 on the first slab; (δ − j²/40)/j⁴ enclosure contains −133/105600 (§3) | CERTIFIED (enclosure); coefficient derived, re-derived by referee H and blind REP5 | 565736c4c4 | H; REP5 9e0df55e7a | J⁴ term OPEN, "not found in print" (GATE-NOVELTY B2); limit existence not proved |
| C3 | Point and cell certificates for every j in [1, 10.144] (§3) | CERTIFIED | d62bfc0e0b (+ AMEND-P2a-1 9356a34707) | referee J (CORRECT); blind REP2 227cbcafc2 | as C1 |
| C4 | Point and cell certificates for every j in [10, 1030.35] (§3) | CERTIFIED | eba765ebe0 | referee K (CORRECT); blind REP3 a669f700df | as C1 |
| C5 | Uniform ceiling J/M_ADM² ≤ 0.96223488883111534543 for every μ > 0; cells to j1 = 11994.815017806193245 (§3) | CERTIFIED | 3031217ac8 | referee L (CORRECT); blind REP4 ad62cca292 | OPEN (P2C-SCOUT, GATE-P2d claim C) |
| C6 | δ(j)/√j nondecreasing (monotone cells) | PROVED; CITED (corollary of DGC Lemma 3.1 + scaling, credit DGC) | - | referee I; Lean 036ea0f37b covers only the cell inequality | credit DGC |
| C7 | Comparison lemma C-REP4 (class K, Kato), C-REP4', Lemma TL (§2) | PROVED | - | referees M, N, O (CORRECT / CORRECT WITH FIXES) | TL: claim the use, not the inequality (GATE-P2d claim T) |
| C8 | sup over μ > 0 of J/M_ADM² in [0.9282004538, 0.9282004687]; rounds to 0.928200 (§4) | CERTIFIED | 9677804a51 (P2d), 560847f6a4 (P2e) | referees P, Q (CORRECT); blind REP5 9e0df55e7a (slopes), REP6 0405395ed8 FAILED, REP6b 3e929bb723 HOLD | certified ceiling OPEN (GATE-P2d) |
| C9 | Tangent-free chord tail ≤ 0.92820047282 < 0.9282005 on certified P2e data (§4) | CERTIFIED (given C12) | a580a8ad0f (B2V) | referee R (D6 direction) | method, not headline |
| C10 | Lovelace et al.'s "ε_J ≤ 0.928200" exceeded literally by certified members; correct as six-decimal rounding (§4, §10) | CERTIFIED (inequality) + CITED | 560847f6a4 | referees P, Q | earliest six-digit print among opened sources: LOPC 2008 |
| C11 | M(μ) nondecreasing, difference quotients in [0, 1], M(μ) − μ ≤ M_0 ≤ M(μ) (§5) | PROVED; credited to DGC Lemmas 3.1-3.2 ("we record") | - | referee R | PARTIAL, implicitly TAKEN (GATE-WING1 a, e) |
| C12 | Strict monotonicity, difference quotients in (0, 1), J ≠ 0 (§5) | PROVED | - | referee R | OPEN to floor (GATE-WING1 f); "to our knowledge" |
| C13 | Strict convexity of M(μ), J ≠ 0 (§5) | PROVED | - | referee R | OPEN (GATE-CONVEX F1); credit Laetsch 1975 (multiplier form, abstract only) |
| C14 | M ∈ C^1, M'(μ) = 1 + 2 lim r w in (0, 1) (exact tangent slope; not 2L) (§5) | PROVED | - | referee R (fix 4c applied) | not gated; supporting lemma |
| C15 | χ strictly increasing in j; supremum not attained; sup = χ(0) = J/M_0² of u_0 = lim u_μ, J_ADM = J at μ = 0 (§5) | PROVED | - | referees R, S | OPEN as theorem (GATE-WING1 f) |
| C16 | Conjecture E: M(μ) − M_0 ~ μ^{2ν}, ν = 1.38201913608948921114 (§5) | PROBE (heuristic); sealed numeric check V4 HOLD | a580a8ad0f | referee R ("exact" removed) | not gated |
| C17 | Comparison B3-5 at μ = 0 and two-sided M_0 bounds from one super-solution (§6) | PROVED | - | referee R (fixes applied), referee S | - |
| C18 | χ(0) in [0.92820045875619, 0.92820045878266]; M_0/√J in [1.0379563989038, 1.0379563989187] (§6) | CERTIFIED | 10a65e940a (B3T) | referee S (CORRECT, rerun byte-identical); blind REP7 133b99d241 | certified μ = 0 value OPEN; 3 digits at μ = 0 in print (HHO 2009) |
| C19 | E3_0: M_0 = (1/2π) ∫ SΦ_0^-7, u_0 = N[g_0], r²∂_r u_0 → −M_0/2 (Appendix A) | PROVED | - | referees R (item 1), U (item 3) | - |
| C20 | Throat flux identity and M_0 = (1/(6π)) ∫ (4|∇V|² + V²) = (2/3) Σ (2l+1) a_l² (§7) | PROVED | 0ce319a05a (FLUXT, float slice test) | referee U (CORRECT, two limit steps now written: Appendices A, B) | formula OPEN (GATE-FLUX F1); mechanism TAKEN (KMPS 1999; Marques; Corvino-Pollack); 1-D throat ODE TAKEN (HHO 2009) |
| C21 | Flux digits M_0 = 1.03795639891365137036534431609, χ(0) = 0.9282004587652052045098878 (§7) | PROBE | - | referee U (shooting reproduces 25 digits) | > 6 digits not in print |
| C22 | Lemma 8: for μ > 0 the throat is the only closed connected MOTS and minimal surface (immersed surfaces cover it); weakly outer-trapped surfaces lie in r ≤ μ/2; strict area minimum at r = μ/2 (§8.1) | PROVED | - | referees T (original route, CORRECT WITH FIXES), V (reflection route, CORRECT) | uniqueness OPEN (GATE-NOHORIZON N2); area minimum PARTIAL in print, now proved |
| C23 | Theorem 9: extreme data contain no closed minimal surface, MOTS or apparent horizon; every coordinate sphere strictly mean-convex (§8.2-8.4) | PROVED (two proofs; second blind) | - | referee V (CORRECT / CORRECT WITH FIXES); blind re-proof seat (REPROOF-NOHORIZON) | OPEN named problem (DGC Sec. 2); expectation in print (DLZ, HHO); sliding = Berestycki-Nirenberg; touching = Benko-Stavrov Allen, Reiris |
| C24 | Extreme area A_0(r) decreases strictly to A_0 = ∫V⁴; A_μ(μ/2) > A_0 for μ > 0 (§8.4) | PROVED | - | referees T (§4), V (§5) | limit value TAKEN (HHO eq. 26); monotone approach OPEN |
| C25 | c_BY = A_0/(8πJ) ≈ 1.1991; A_0/J = 30.13745891933985963534 (§8.4); sharpness (A_μ(μ/2) → A_0) not proved | PROBE | - | referee T (labelled probe); REFEREE-RELEASE | A ≥ 8π|J| TAKEN for μ > 0 throats (Dain-Reiris 2011; not the μ = 0 cylindrical-end data); c ≈ 1.199 TAKEN numerically (LOPC, HHO) |
| C26 | M_AH/a probe enclosures match all eight Cook-York values (§8.4) | PROBE (rigorous computation, unsealed) | - | referee T (j = 3 row reproduced) | numbers in print at 3-4 digits |
| C27 | Cook-York Table I Energy/a agrees with certified δ at j = 1, 10, 30, 300, 1000, 10000; one unit high at j = 100 (certified) and j = 3 (probe) (§3, §10) | CERTIFIED comparison + CITED | P2a, P2b, P2c seals | referee T (Table I digits read from the page image) | fairness note |
| C28 | Lousto et al. 2012 Table I BY row 0.9282 agrees with χ(10⁴) from REP4 (§3) | CERTIFIED comparison + CITED | ad62cca292 | - | spin of their row inferred |
| C29 | Repair of Gabach Clement (71)-(72) by a barrier (Appendix D) | PROVED | - | referee V (supplied the repair); blind re-proof seat flagged the related h-bound step | fairness note |
| C30 | Self-contained weak/Hopf/strong maximum principles for Δ − c, c ≥ 0 locally bounded (Appendix C) | PROVED (classical argument written out) | - | checked in REFEREE-RELEASE §4 | replaces the unopened textbook locator |
| C32 | For A = P = 0 the exact super-solution condition of the barrier √(1 + Q/r) − 1 is Q² ≥ 3|J|; GC's 9J in (28) is sufficient, not sharp (§10) | PROVED (derived; re-derived in REFEREE-RELEASE §1) | - | blind re-proof seat (check C5e); REFEREE-RELEASE | fairness note (GC eqs (28)-(31)) |
| C33 | De Ceuster-Li use "certified" for a 1e-10 floating-point residual gate, not a proof (§12) | CITED | - | - | GATE-NOVELTY A1; GATE.md (2609.10173 Sec. IV C) |
| C31 | C^1 throat convergence via W^{2,p} + Morrey (Appendix B) | PROVED modulo the standard interior L^p estimate (not opened) | - | referee U (item 4 supplied the route) | - |

## Errata applied in this text (from ERRATA-PENDING)

1. Item 1: all δ tables regenerated from the exact rationals with directed rounding (matches `results/DISPLAY-TABLES.md`).
2. Item 3: the text says Lean 036ea0f37b covers only the cell inequality.
3. Item 5: P2e window printed as [0.9282004538, 0.9282004687].
4. Item 6: "0.928200 first printed": earliest among opened sources is Lovelace et al. 2008 (DLT 2002 prints 0.928;
   no 0.928 was found in the opened Cook-York 1990 text).
5. Items 7-9: √r = ((1+ξ)/(1−ξ))²; the referee wording for the extreme spin; "width = integrated residual".
6. Item 2 (instrument self-sha and chain asserts) is for the next instrument version, not this text. Item 4 (re-run the
   novelty gate within a week of release) is open.

## Numbers and their sources

Direction: `down` = exact lower end rounded down; `up` = exact upper end rounded up; `ball` = arb ball as stored;
`exact` = exact rational or value quoted verbatim; `probe` = non-rigorous value quoted verbatim.

| Printed | Direction | Source file (field) |
|---|---|---|
| δ(J), χ(J) balls at J = 0.0025, 0.025, 0.0625, 0.125, 0.25 (Table 1) | ball | result_p1.json (points[].delta_str, chi_str) |
| [0.025000000000 ± 2.13e-13] | ball | result_p1.json (first_slab.delta_over_j2_at_0_str) |
| [−0.00125946969697 ± 8.30e-15]; −133/105600 | ball; exact | result_p1.json (first_slab.delta4_limit_str, delta4_contains_m133_105600) |
| 266/825 | exact | GATE-NOVELTY.md (conversion, checked: 4⁴ · 133/105600 = 266/825) |
| 2e-11 (slope-route accuracy) | exact | GATE-NOVELTY.md, scope note from REFEREE-P1 |
| Table 2, P2a rows (δ ends, χ ends, rel. widths) | down/up | result_p2a.json (nodes[].delta_lo, delta_hi, chi_lo, chi_hi, rel_width) |
| Table 2, P2b rows | down/up | result_p2b.json (nodes[]) |
| Table 3, P2c rows | down/up | result_p2c.json (nodes[]) |
| 0.9326738205130032; cell [153.67, 161.36] | up | result_p2b.json (max_cell_chi_max) |
| 0.96223488883111534543 | up | ceiling_p2c.json (c_dec_up) |
| 0.83394087233933423687, 0.93267382051300310048, 0.93989487162744096378 | up | VERDICT-P2c.md table (from ceiling_p2c.json per_piece_max) |
| 0.928146, 0.962235, 0.0340884 (P2c window, width rounded up) | down/up | ceiling_p2c.json (X3_max_over_nodes 0.9281465706057948783930773; c; width recomputed 0.034088318226) |
| Q² ≥ 3|J| | exact | scout/REPROOF-NOHORIZON.md (check C5e); REFEREE-RELEASE §1 |
| 1e-10 (De Ceuster-Li residual gate) | exact (print) | GATE-NOVELTY.md A1; src/2609.10173.txt Sec. IV C |
| 675/256 | exact | P2C-SCOUT.md (DGC barriers), via LOG.md 08:10 |
| 0.92820096252204785788, 0.92820036654400608722, 0.9282009626, 0.9282003665, 5.95979e-7 | up/down | result_p2d.json (c_dec_up, X3_dec_down, c_minus_X3_dec_up) |
| 0.92820046861667179726, 0.92820045389621230212, 0.9282004687, 0.9282004538, 1.47205e-8 | up/down | result_p2e.json (c, X3, c_minus_X3_dec_up) |
| 0.928200468617, 0.928200456416, 0.928200454324 | up | result_p2e.json (pieces[].value) |
| [750.5607143681, 750.5607144170] | down/up | result_p2e.json (nodes[21].point.delta_lo_exact, delta_hi_exact) |
| 5.911e-6 | up | result_p2e.json (nodes[21].tangent.twoL_hi_exact = 5.9107129690605e-06) |
| 3.1383e-8, 4.6103e-8 (margins to 0.9282005); 3.13e-8, 4.61e-8 | down | computed from result_p2e.json c, X3; VERDICT-P2e.md |
| 260 times, 3.1 times | exact | VERDICT-P2e.md disclosure (3), REFEREE-P2e §3 |
| about 2.3 million times | exact | REFEREE-P2e §5(a) (2,315,710) |
| 4.5e-7 (exceedance of 0.928200) | down | result_p2e.json X3 − 0.928200 = 4.53896e-7 |
| 0.0140665, 0.928235200531265, [750.5607142022, 750.5607145913] | up; as printed | rep6/RESULT.md (REP6 results (a)-(c)) |
| 2380 times | exact | VERDICT-P2e.md (REP6 paragraph) |
| [750.5607143858, 750.5607143991] | down/up | rep6/result_rep6b.json (delta_recertified.delta_lo_exact, delta_hi_exact) |
| 7.3143e-6 | up | rep6/result_rep6b.json (tangent.twoL_hi_exact) |
| 0.9282004721 | up | rep6/result_rep6b.json (chi_tail_exact; LOG.md prints the truncated 0.9282004720) |
| 1.8e-5 (margin of the cover) | exact | LOG.md 15:38 |
| s(64) = [0.0520407975655210747020518703815, 0.0520407975655210747020518836049] | down/up | rep5/RESULT.md §1 table (exact rationals in result_rep5.json) |
| 7e-10, 3.5e-10 | exact | VERDICT-P2d.md (REP5 paragraph) |
| [1.02383610515332591435335055921, 1.02383610515332591435338447325] | as printed | VERDICT-REP1.md |
| [0.827109416381791178798080260634, 0.827109416381791178800503789071] | as printed | VERDICT-REP2.md |
| [0.928013132630170991794, 0.928013132630171290517] | as printed | VERDICT-REP3.md |
| δ(3000), δ(10000) REP4 intervals | down/up | rep4/REPORT-REP4.md (results table; result_rep4.json) |
| [0.9281926861896467, 0.9281926861896468] (χ(10⁴)) | down/up | computed from rep4/REPORT-REP4.md δ(10000) ends, χ = j/(4(1+δ)²) |
| 0.92820047282 | up | scout/b2v_result.json (T4_tail_chi_up = 0.92820047281343580...) |
| 20/20, 21/21 | exact | scout/b2v_result.json (T2_decided, T3_decided) |
| 0.882059, width 5.1e-5 | exact (centre rounded) | scout/b2v_result.json (T5_exponents[12]) |
| 0.882019136, 0.914213562 | exact | BAR-B2V.md (V4 targets); scout/B2-PROOF.md table |
| λ_1 = 1.65997689251753810055, ν = 1.38201913608948921114 | probe | scout/B5-DIGITS.md §2b (K = 80); B2-PROOF.md table (1.659976892517514, 1.382019136089480) |
| [0.92820045875619, 0.92820045878266], width 2.647e-11 | down/up | result_b3t.json (window.chi_lo_dec, chi_hi_dec, width_dec) |
| [1.0379563989038, 1.0379563989187] | down/up | result_b3t.json (window.M_lo, M_up exact; rounded here) |
| 2.16576e-12, 2.51977e-12, 0.5183644946, 17 patches, N = 112, K = 12 | as stored | result_b3t.json (supQ_upper, eta_dec, inf_G_N_lower, patches, cfg) |
| 8.9e-16 (quadrature loss) | up | REFEREE-B3T §4; BUILD-B3T §4 |
| [1.037956398789003, 1.037956399038277], [0.928200458542311, 0.928200458988139], [0.928200458680027, 0.928200458850410] | down/up | rep7/result_rep7.json (M0_lo/hi, chi0_lo/hi, BY2 chi_lo/hi) |
| 5.773e-15, 2.562e-12, 8.882e-16, 7.098e-3, 2.77e9 | as computed (float) | result_fluxt.json (score.F1, tolerance.tol, F2, F3) |
| 1.03795639891365137036534431609 | probe | BAR-FLUXT.md / scout/B5-DIGITS.md §1b (flux_K30.log) |
| 0.9282004587652052045098878 | probe | REFEREE-FLUX §6 (chi0 from shooting); same leading digits in B5-DIGITS |
| 25 digits; about 180 digits | probe | REFEREE-FLUX §6 |
| A_0/J = 30.13745891933985963534 | probe | scout/B5-DIGITS.md §2b; LOG.md 16:59 |
| c_BY ≈ 1.1991 | probe | REFEREE-B4 §4 (1.19914); GATE-B4 normalisations (1.199133) |
| 2.033, 2.227, 3.034, 4.673, 8.071, 13.67, 24.69, 77.56 (M_AH/a); 2.048, 2.329, 3.477, 5.759, 10.41, 17.99, 32.83, 103.8 (Energy/a) | exact (print) | Cook-York 1990 Table I (lab copy cook1990.txt), read from the page image (scout/B4-HORIZON.md H1, REFEREE-B4 §3) |
| 10.404865-10.404866 (E/a at j = 100) | down/up | result_p2b.json (nodes j = 100), E/a = 2(1 + δ) |
| 2.32846 (E/a at j = 3) | probe | scout/B4-HORIZON.md §4 table |
| about 2% (Cook-York accuracy) | exact | LOG.md 16:09; scout/B4-HORIZON.md §4 |
| 0.9282, 0.9831, 0.8311 (Lousto et al. Table I) | exact (print) | scout/src/b4gate/1203.3223.txt Table I (GATE-B4) |
| 0.928 (DLT eq. 5; HHO); 0.928200, 0.9837 (LOPC eq. 5); 0.9837 (HHO) | exact (print) | src/gr-qc_0201062.txt l. 81; src/0805.4192.txt l. 122; scout/src/0908.1063.txt l. 690-694 |
| −1.4 (Lovelace exponent); 2ν = 2.764 | exact (print); probe | src/0805.4192.txt l. 673-691; GATE-FLUX D2 |
| 63J²/4 (C_V) | exact | REFEREE-WING1 §2 |
| j²/40, 1 − (3/40) j² | exact | REFEREE-WING1 §1, §4; REFEREE-TL §3 |
| 15 referee reports; 17 sealed tests (16 HOLD, 1 FAILED) | count | REFEREE-*.md files; LOG.md settle lines; SCORECARD.md |
| 123 nodes, 122 cells; 100 nodes, 99 cells; 14 nodes; 22 nodes | count | VERDICT-P2a, P2b, P2c, P2e |

Tokens allowed without a result-file source (definitions and exact constants): 9/4, 7, 18 (source coefficients);
1/2, 1/4, 3/4, 3/8, 1/8, 1/(6π), 2/3, 2l+1 (formula constants); 2^k, 2^21 = 2097152, 2^20 = 1048576, 2^19 = 524288
(node positions); 0.9282005 (sealed decision bar, BAR-P2e); 1e-3, 2e-3, 2e-4, 0.02 (sealed bar constants, BAR-B2V);
1.001 (null-control factor); 65 × 65, 256 bits (B3T configuration); 2026-11-14 (search floor of the gates).
