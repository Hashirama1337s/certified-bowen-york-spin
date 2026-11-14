# Clean-room rebuild (2026-11-14)

Method: the whole repository was copied to a fresh folder outside the project; `python tools/check_pins.py` there
(26/26 pins and 117/117 manifest entries match), then the documented default run `python tools/reproduce.py --work <fresh empty folder>`
under an 8 GB memory cap, one process at a time. The engine was loaded from the work folder's copy of `code/gravity_engine`
(reported by `tools/run_sealed.py` on every engine-using run), not from the lab.

| Check | Run time | Result | sha256 of the output |
|---|---|---|---|
| pins | seconds | ALL MATCH | - |
| ceiling (P2c re-assembly) | < 1 s | IDENTICAL | a63509850d6b3ad2f26d0b66463f262446b8f48a67e143941d9b81e836ac9b8c |
| b2v (B2V on P2e nodes) | < 1 s | IDENTICAL | 7bbec769ed58b7482918247d04c5dc531abe7a779105f89179efd43ce644d6b8 |
| b3t (full certificate) | 0.2 min | IDENTICAL | 9a73e169c2e1c7811a6a7d9af084c0348a163d19be9b74eba750235e9ad28000 |
| fluxt (full run) | 0.3 min | IDENTICAL-EXCEPT-RUN-RECORD (`go`, `tolerance.solve_seconds`) | canonical b45731e1dc22171c2f64a941e58547f8938a9ceb908f6e1e6a2f53cca9dc88cf |
| p2e (full main run, 4 invocations) | 28.3 min | IDENTICAL | 2aa80d5fbda95e76a2e2f96c20244f95f647a8f4f674be7618825e0a1dbc8dde |

Verdict: CLEAN ROOM REPRODUCED (default checks). Not re-run: the optional checks (p2d, p2c, p2b, p2a, p1, rep4, rep7),
the control modes and the engine tests (see tools/reproduce.md). After the clean-room copy only text files changed
(LICENSE-CODE scope line, this report, MANIFEST.sha256); no code or result file changed.

## Re-check after the path patch (PATCHES.md)

After the lab paths in the shipped copies were replaced by reversible placeholders, the repository was copied to a
second fresh folder: `tools/check_pins.py` ALL MATCH (56 pins, 32 patched files restored to their original sha256,
124 manifest entries), then `python tools/reproduce.py --only ceiling b2v b3t fluxt` (the runner restores the sealed
bytes in its work folder first): ceiling IDENTICAL (a6350985...), b2v IDENTICAL (7bbec769...), b3t IDENTICAL
(9a73e169...), fluxt IDENTICAL-EXCEPT-RUN-RECORD (canonical b45731e1...); ALL CHECKS PASS. P2e (28 min) was not
repeated: none of its inputs, its instrument or the engine changed in their restored form, which check_pins verifies.
