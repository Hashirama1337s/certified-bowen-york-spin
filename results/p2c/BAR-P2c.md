# BAR-P2c: a certified uniform spin ceiling for the whole Bowen-York spinning-puncture family

Written 2026-11-14 by the lead before p2c.py exists. Main quest, Level 2. Gate: P2C-SCOUT.md (OPEN to its floor: no published uniform bound sup chi <= c < 1 over the family; per-member chi < 1 is published; the ~0.928 ceiling is numerical). Engine: gravity_engine/lichnerowicz_axisym.py version 2 (commit 737cbbf).

instruments: p2c.py

## Object
For every mu > 0 and every J >= 0, chi = J / M_ADM^2 depends only on j = 4J/mu^2 (scale invariance; REFEREE-P2DESIGN item 4). Certify:
(i) point certificates (engine R1-R5) at every node of the grid j_{k+1} = (5/4) j_k starting at P2b's last node j_0 = 1030.3467644609451531 (exact rational string from result_p2b.json) up to the first node >= 10000, plus extra nodes j = 3000 and j = 10000;
(ii) K1 cells between neighbours (delta(j)/sqrt(j) nondecreasing; Lean-checked cell inequality lean 036ea0f37b; arb, outward) with per-cell chi_max;
(iii) the tail X1 at the top node j1: chi(j) <= j1 / (4 delta_-(j1)^2) for every j >= j1 (REFEREE-P2DESIGN item 4, CORRECT);
(iv) the ceiling c := max( max cell chi_max over P1 + P2a + P2b + P2c cells covering [0, j1], chi_max bound on P1's range [0, 1] (chi <= j/4 there), tail X1 at j1 ), computed in exact/arb arithmetic from the stored certified records (result_p1.json, result_p2a.json, result_p2b.json and the new nodes).

## Predictions
- C1: every new node certifies.
- C2: the certified ceiling satisfies c < 1 (the headline: every Bowen-York spinning puncture with mu > 0 has J / M_ADM^2 <= c < 1).
- C3: c <= 0.98 (lead estimate ~0.975: tail X1 at j1 = 1e4 ~0.97, new cells ~0.94, earlier cells <= 0.9327).
- C4: at j = 3000 and j = 10000 the certified delta overlaps TARGETS-P2 item T5's reference with window ref (1 +- 1e-6) (references: double precision; spread <= 5e-9).
- C5: the lower bound X3 = chi_lo at the top node is reported, and c - X3 (the certified uncertainty of the family's maximal spin) is reported.

## Controls (any failure = VOID; each its own command; GO via `lab.py go`)
- ENGINE: gravity_engine tests pass (17, incl. P2a and P2b reproduction).
- CHAIN: the first node string equals P2b's last node exactly; the ceiling script re-reads P1/P2a/P2b result files by sha256 and refuses on mismatch.
- N1 (u+ - 1e-6 j^2 s) at j = 10000 -> VIOLATION; N4 near-threshold (margin S/2 (1 - 2^-10)) at the top degree -> FAIL through the bound path (REFEREE-P2b item 5 suggestion).
- REF at j = 10000: exact Chebyshev reconstruction, spherical-coordinate residual signs (l >= 4), Bernstein reconstruction.
- TAIL: X1 recomputed independently in mpmath (120 digits) from the stored node record; CEIL: c recomputed independently from all stored cell records.
- RESUME and exit-code convention (0 finished, 3 rerun, other error); MEM peak per node reported; NaN/inf; `expected_false` lists any by-design false flag as "<path> (<why>)".

## Verdict rule
P2c CERTIFIED if C1-C5 hold and every control behaves as stated. If c >= 1, report the certified c (it is the result) and the binding piece; no retuning. Disjoint from a T5 window = investigate.

## Seal fields
- formula_ref: DESIGN-P2.md:68-88
- calibrated_on: eba765ebe0
- depends_on: 565736c4c4, d62bfc0e0b, eba765ebe0
- null_not_symmetry: N1, N4 perturb the sandwich; TAIL/CEIL are independent recomputations; CHAIN is an identity check.
- target_form: C1 boolean; C2 boolean (c < 1); C3 c in [0, 0.98]; C4 overlap with [ref(1-1e-6), ref(1+1e-6)]; C5 reported numbers.
- ceiling: Kerr chi < 1 (the bound being certified below it); numerical Bowen-York maximum 0.928200 (Lovelace et al. 0805.4192 eq (5)); certified lower bound so far 0.92802064647 (P2b).
- toy_model: J = 0 through the engine (refused/trivial).
- tolerances: degrees (N, K = 12) from need, PREC 160 with engine v2 midpoint fit operators; cache_ops False; one process; width target at j = 1e4 about 5e-3 (scout estimate), enough for C3 headroom.
- parameters: grid 5/4 from 1030.35 to >= 10000 plus j = 3000, 10000
- files: BAR-P2c.md, p2c.py
- command: memcap 8 python p2c.py --main result_p2c.json
