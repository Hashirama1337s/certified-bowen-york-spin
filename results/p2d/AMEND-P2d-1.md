# AMEND-P2d-1: fixes from REFEREE-TL (N) to DESIGN-P2d (sealed f825f520)

Written 2026-11-14 by the lead before any tangent-certificate code exists. Amends DESIGN-P2d; TARGETS-P2d unchanged.

1. Sub-solution side (REFEREE-TL 1): LEMMA-REP4.md states TL's sub side as "Delta u~ >= F(Phi') a.e. on {u~ > u}" (Kato on (u~ - u)^+, where Phi' > Phi_true > 0). TL-P (G1 >= 0, step 4.3.5) is dropped as a requirement; if computed it is a redundant check, not a gate.
2. Z ansatz (1): the module asserts that Z uses REP4's ansatz structure (angular part times Y(1-Y), Z(1, v) = P_Z(1)), N_Z = N and K >= 1, before the degree count (9N+50, 9K) is used.
3. Tail bracket (2): assert 1 + delta_lo - tL > 0 in exact arithmetic for the tail as well as for every cell.
4. Pieces (4i): assert adaptive pieces cover [0, 1] exactly with no gaps (rep4.py:463 pattern).
5. Rounding trap (4v): Z (and W) dyadic rounding uses REP4's coefficient-scaled grid (rep4.py:594-596), never an absolute 2^-(prec-40) grid; K1 and K1T run at 192/256 bits.
6. K1T wording (3): the small-spin tangent anchor 2L -> 1 - (3/40) j^2 is a CONSISTENCY control (tangent vs point problem; ratio 3 holds for any exponent), NOT an independent guard of the PDE constant 7. Status of TL8: derived by the first-order Gauss law on the tangent problem (REFEREE-TL 3), not the 2-line route.
7. N-TAN (4vi): its expected tangent-consistency gap at every node pair is computed and written into the BAR before sealing; pairs where it cannot bite are listed.
8. Wording (5): headline reports the certified window c - X3 from the run; never "within 1e-6 of the true supremum"; family qualifiers per REFEREE-P2c 6a; "a certified line below M on (0, mu1], close to the tangent".
9. Prior art: the GATE-P2d verdict (running) decides any "new lemma" wording.
