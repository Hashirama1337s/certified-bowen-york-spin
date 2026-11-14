# AMEND-P2a-1: which cell width B4 means (before the main run)

Written by the lead 2026-11-14, after the build seat's controls (which printed both widths at control nodes) and BEFORE any `--main` data. Parent claim d62bfc0e0b; bar file de758dbbf3; p2a.py unchanged.
B4's "relative delta width" of a K1 cell is the POINTWISE relative width of the cell enclosure at a given j in the cell, R(j_{k+1})/R(j_k) - 1 in REFEREE-P2DESIGN item 3 (constant across the cell; result key `rel_width`), as derived in DESIGN-P2 / REFEREE-P2DESIGN. It is not the spread of delta over the whole cell (key `rel_width_uniform`), which also includes delta's own growth across the cell; that number is reported but not scored.
Disclosure: the control run showed pointwise 2.83 % and whole-cell 3.85 % near j = 1. This amendment fixes the reading the design and referee already used; it is not chosen to pass.
Thresholds: no number changes (3e-2 stays). What counts as FAILED / VOID is unchanged.
