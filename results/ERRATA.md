# Errata (display corrections; the exact rationals were always right)

These items were found by the hostile referee rounds after the sealed runs. None changes a certified value or a
verdict: every certificate is stored as exact rationals in the result files, and those were always correct. The items
concern how some values were printed (decimal rounding direction), wording, and instrument hygiene for the next
version. The sealed files themselves (instruments, results, VERDICT files) are shipped unchanged, so the old prints
are still visible in them; the corrected displays are in `results/DISPLAY-TABLES.md` (regenerated from the exact
rationals with outward rounding by `tools/display_tables.py`).

1. **Upper ends printed rounded down (P2a, P2b, P2c logs).** The engine log line in
   `code/gravity_engine/lichnerowicz_axisym.py` (line 852) prints `delta_hi_dec[:20]`, which truncates upper ends
   (rounds them down in the last digit). VERDICT-P2c was fixed at the time; the tables in VERDICT-P2a and VERDICT-P2b
   copied the same truncated upper ends (example: P2a j = 1 upper printed 0.023836105153325914, below the exact
   0.0238361051533259143534...). Corrected display: `results/DISPLAY-TABLES.md`. The next engine version prints
   dec_ceil(delta_hi, k) / dec_floor(delta_lo, k) and also restores the arb working precision in a `finally` block.
2. **Chain checks to add (instrument hygiene).** The next instrument version adds a self-sha check in `--main` and
   asserts that each chain link's first node is >= the previous file's first node (referee note on P2c). The shipped
   instruments are the sealed versions and do not contain these checks.
3. **Scope of the Lean check (P2c).** The Lean file used for P2c covers only the cell inequality. DGC Lemma 3.1 and
   the scaling argument are cited, not Lean-checked. Release text must say so.
4. **Novelty gate to re-run.** The prior-art gate for the uniform ceiling (P2c) must be re-run (arXiv, INSPIRE full
   text, Zenodo) within a week of release.
5. **P2e lower end printed rounded the wrong way.** A lab log line printed P2e's lower end as 0.92820045390; rounded
   down it is 0.92820045389. The first ledger text gave the window as [0.9282004539, 0.9282004687] with both ends
   rounded the wrong way; the correct outward display is [0.9282004538, 0.9282004687] (ledger text corrected; see
   `results/ledger-bhslice.json`, seal 560847f6a4, "corrections"). Release text uses the exact-rational decimals.
6. **Attribution of the printed 0.928200.** Before release wording, confirm who first printed 0.928200 (Lovelace et
   al., arXiv:0805.4192, or earlier).
7. **A wrong power in a scout note (B3T).** A project note printed sqrt r = ((1+xi)/(1-xi))^4; the correct
   expression is ((1+xi)/(1-xi))^2 (the ^4 expression is r itself). `code/b3t.py` and the probes use the correct form.
8. **B3T wording for any release text.** Certified: sup over mu > 0 of J/M^2 equals chi(0) of the extreme (mu = 0)
   Bowen-York data (J_ADM = J; u_0 unique by Gabach Clement, arXiv:0911.0258, Lemma 4.3, with u_0 = lim u_mu in that
   class by DGC Theorem 2.1), and chi(0) lies in [0.92820045875619, 0.92820045878266]. Earlier prints: 0.928 at
   mu = 0 (Hannam et al. 2009, numerical); 0.928200 extrapolated (Lovelace et al. 2008, numerical). "Directly at
   mu = 0" describes the computation only. The field "referee R pending" inside `results/b3t/result_b3t.json` is
   stale: that referee reported "correct with fixes".
9. **B3T build-note wording.** A build note said "width ~ 2 eta"; correct is "width = integrated residual" (chi
   width 2.65e-11 against 2 eta = 5.0e-12).

Repository-level notes (added when this repository was assembled):

- `results/b2v/b2v_result.json` is the JSON object of the lab file without the shell residue that a redirect appended;
  see `results/b2v/NOTE.md`. Content unchanged and reproduced byte for byte.
- P2b ran on engine module version 1 (lichnerowicz_axisym.py sha256 be006193...). The shipped engine file is module
  version 2 (sha256 3659ae21...), the one P2c ran on; version 2 changed only the floating-point fit operators above
  N = 100 (proof path unchanged), and its own test suite reproduces the P2a and P2b nodes. A re-run of P2b from this
  repository therefore uses version 2 (`tools/reproduce.md`).
