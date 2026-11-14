# AMEND-P2d-2: piece rule (build-seat finding, 2026-11-14), before any tangent code

Amends DESIGN-P2d (1a9c1237d9) after AMEND-P2d-1. The build seat found DESIGN-P2d's adaptive rule ((sup E / inf E)^8 <= 1.25 per piece, cap 128) infeasible for j >~ 1e3: the angular variation of E alone gives ~1.5 at Y = 1/2, needing 306 boxes at j = 1e3 and 590 at j = 1e4.

Decision (lead, technical default):
1. Default pieces = REP4's fixed 14 Y-pieces (the path backed by REP4's HOLD and REFEREE-REP4), for the point certificate and the tangent certificate alike. DESIGN-P2d's rigorous probe at 2^18 used these pieces.
2. Optional targeted refinement: after a fixed-piece pass, split (in Y or v, the larger pure-direction variation) only the boxes whose bound is within a factor 2 of the current maximum, repeating up to a cap of 128 boxes; every split cover passes check_cover (disjoint interiors, exact total area 1). Used only if the fixed pieces fail a node or miss a sealed width bar; every use is recorded per node.
3. No node is retried with different settings after seeing its value against TARGETS-P2d; degree and piece choices per node are fixed by the instrument's schedule before the main run.
