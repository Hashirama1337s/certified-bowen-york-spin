# BAR-P2b: certified Bowen-York spinning-puncture data for every spin j = 4J/mu^2 in [10, 1000]

Written 2026-11-14 by the lead before p2b.py exists. Main quest, Level 2. Engine reuse: gravity_engine/lichnerowicz_axisym.py (first reuse). Calibration probes (not scored; p2bprobe/probe_ladder.log, probe_cal.py output in LOG): j = 20, 40, 80, 150, 400, 1000 all certified; j = 1000 at (N, K) = (96, 30): 91 s, rel width 6.5e-6, peak 6.2 GB.

instruments: p2b.py

## Object
DESIGN-P2.md 4.3 with the engine module: point certificates (R1-R5) at every node of the grid j_{k+1} = (21/20) j_k from j_0 = 10 to the first node >= 1000, plus extra nodes j = 30, 100, 300, 1000; degrees (N, K) rising with j (chosen by the build seat from need, so that every node certifies with relative width <= 1e-6 where possible); K1 cells between neighbouring nodes (delta(j)/sqrt(j) nondecreasing; Lean-checked cell inequality lean 036ea0f37b) in arb with outward rounding, per-cell chi bounds (REFEREE-P2DESIGN fix 3a); K2 no-gap check. mu = 1, J = j/4. One process at a time (memory), resumable chunks <= 520 s.

## Predictions (expected values: TARGETS-P2.md item T5, sealed b18e905ff1)
- B1: every node certifies.
- B2: at j = 30, 100, 300, 1000 the certified delta overlaps [ref (1 - 1e-8), ref (1 + 1e-8)] with ref from TARGETS T5.
- B3: relative width of the certified delta <= 1e-5 at every node.
- B4: every K1 cell on the 21/20 grid covering [10, 1000] has POINTWISE relative delta width (as AMEND-P2a-1 defines it) <= 5e-2; no gaps.
- B5: the j = 10 node overlaps P2a's certified delta(10) = [0.738555349918561234, 0.738555349918583461] (VERDICT-P2a), so P1 + P2a + P2b chain over [0, 1000].

## Controls (any failure = VOID; each its own command; GO only via `lab.py go`)
- ENGINE: gravity_engine/tests pass on this machine before the run (reproduce P2a byte for byte).
- N1 (wrong ansatz u+ - 1e-6 j^2 s) at j = 100 and 1000 gives certified VIOLATION; N4 (eta' = 0.25 x needed margin) at j = 100 gives FAIL through the bound path.
- REF (REFEREE-P2a): exact Chebyshev reconstruction of Q at random rational points (j = 1000); spherical-coordinate residual sign check with l > 2 modes (mpmath) at j = 300 and 1000; Bernstein reconstruction at random rational points.
- T0: toy J = 0 refused/trivial (S == 0); T0c brackets and chi < 1 at every node.
- RESUME: interrupted + resumed run byte-identical (exact p/q strings; resume_tested true).
- MEM: peak memory per invocation under the 8 GB cap (memcap); report the peak.
- NaN/inf guards. Result JSON lists `expected_false` for any flag that must be false by design.

## Verdict rule
P2b CERTIFIED if B1-B5 hold and every control behaves as stated. Disjoint from a reference window = investigate, no retuning; a null that passes = VOID.

## Seal fields
- formula_ref: DESIGN-P2.md:155-180
- calibrated_on: d62bfc0e0b, 227cbcafc2
- depends_on: d62bfc0e0b, 227cbcafc2
- null_not_symmetry: N1, N4 perturb the sandwich; T0 exact known answer; REF checks are independent evaluations.
- target_form: B1 boolean; B2 overlap with [ref(1-1e-8), ref(1+1e-8)]; B3 [0, 1e-5]; B4 [0, 5e-2] and boolean; B5 boolean overlap.
- ceiling: Bowen-York chi_max ~ 0.928200 (numerical, Lovelace et al. 0805.4192 eq (5)); P2b reaches chi ~ 0.92801 at j = 1000 (TARGETS T5), deficit ~1.9e-4; Kerr chi < 1.
- toy_model: J = 0 (u = 0) through the engine.
- tolerances: degrees from need for B3 with headroom where feasible; memory cap 8 GB per process.
- parameters: grid 21/20 from j = 10 to >= 1000 (about 95 nodes) plus j = 30, 100, 300, 1000
- files: BAR-P2b.md, p2b.py
- command: memcap 8 python p2b.py --main result_p2b.json
