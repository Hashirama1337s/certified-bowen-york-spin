# Reproduce: re-run the certificates

Requirements: Python 3.11 (tested with 3.11.9), python-flint 0.9.0, numpy 2.4.4, mpmath 1.3.0, scipy 1.17.1
(sympy 1.14.0 only for the engine tests). The lab ran every step under an 8 GB memory cap; the largest default step
peaked near 2.3 GB. Every step is a single process.

## One command

From the repository root:

```
python tools/check_pins.py        # every sealed file against its pinned sha256, every file against MANIFEST.sha256
python tools/reproduce.py         # the default checks below, about 30 minutes on one core
```

Expected last line: `[reproduce] ALL CHECKS PASS`. The runner works in `./work` (or `--work DIR`, which must be
empty); pass another interpreter with `--python "py -3.11"`. Subsets: `--only ceiling b2v b3t fluxt` (under a minute).
Optional checks are added by name (`python tools/reproduce.py p2d`) or all at once with `--all`.

## Default checks and expected outputs

| Check | What it re-runs | Command inside the work folder | Run time | Output | Expected sha256 of the output |
|---|---|---|---|---|---|
| pins | sha256 of all 56 pinned files (17 instruments, 2 engine modules, 7 result files, 23 BAR / AMEND / TARGETS files sealed in the ledger, 7 controls files read at GO; path-patched copies checked after restoring, PATCHES.md) | `tools/check_pins.py` | seconds | - | `tools/pins.json`, `tools/path_patch.json` |
| ceiling | P2c: the uniform ceiling c = 0.96223488883111534543 re-assembled from the exact rationals of P1, P2a, P2b, P2c (inputs refused unless their sha256 match the pins in `p2c.py`) | `p2c.py --ceiling ceiling_p2c.json --p2c-sha256 4b2dbf7b...` | seconds | `ceiling_p2c.json` | `a63509850d6b3ad2f26d0b66463f262446b8f48a67e143941d9b81e836ac9b8c` (byte-identical) |
| b2v | B2V: monotonicity, convexity, tangent-free tail and exponent tests on the certified P2e nodes | `scout/b2_verdict.py result_p2e.json` | seconds | stdout, stored as `b2v_result.json` (LF) | `7bbec769ed58b7482918247d04c5dc531abe7a779105f89179efd43ce644d6b8` (byte-identical) |
| b3t | B3T: the full certificate for the extreme data, chi(0) in [0.92820045875619, 0.92820045878266] | `b3t.py --main result_b3t.json` | ~20 s | `result_b3t.json` | `9a73e169c2e1c7811a6a7d9af084c0348a163d19be9b74eba750235e9ad28000` (byte-identical) |
| fluxt | FLUXT: the throat flux identity on an independent 2-D float solution | `fluxt.py` | ~20 s | `result_fluxt.json` | canonical `b45731e1dc22171c2f64a941e58547f8938a9ceb908f6e1e6a2f53cca9dc88cf` (run-record fields `go`, `tolerance.solve_seconds` removed; the raw file differs only there) |
| p2e | P2e: the full main run, 22 point + 22 tangent certificates, supremum in [0.9282004538, 0.9282004687] | `p2e.py --main result_p2e.json` (repeats while exit code 3) | ~27 min | `result_p2e.json` | `2aa80d5fbda95e76a2e2f96c20244f95f647a8f4f674be7618825e0a1dbc8dde` (byte-identical; this is also the pin in `BAR-B2V.md`) |

"Canonical" = sha256 of the JSON re-serialised with sorted keys and no spaces after removing the run-record fields
listed in `VOLATILE` in `tools/reproduce.py` (timings, the text of the launch row, recorded file locations).

## Optional checks (long; not part of the default run)

| Check | Command inside the work folder | Run time in the lab | Shipped output (sha256) |
|---|---|---|---|
| p2d | `p2d.py --main result_p2d.json` | ~15 min (2 invocations) | `results/p2d/result_p2d.json` (5620dd7a...) |
| p2c | `p2c.py --main result_p2c.json` | ~18 min (3 invocations) | `results/p2c/result_p2c.json` (4b2dbf7b...) |
| p2b | `p2b.py --main result_p2b.json` | ~25 min (3 invocations) | `results/p2b/result_p2b.json` (9f4fa988...) |
| p2a | `p2a.py --main result_p2a.json` | about an hour | `results/p2a/result_p2a.json` (6fcc3a94...) |
| p1 | `p1.py --main result_p1.json` | about an hour | `results/p1/result_p1.json` (927a72d7...) |
| rep4 | `rep/rep4/rep4.py` (stage cache; exit 3 = rerun) | ~30 min | `results/rep/rep4/result_rep4.json` |
| rep7 | `rep/rep7/rep7.py run BY`, `run BY2`, `run BYP`, `run RN`, `run KN`, then `assemble` | a few minutes per case | `results/rep/rep7/result_rep7.json` |

The optional checks were not re-run when this repository was assembled; their expected outputs are the shipped files.
Notes: P2b was certified on engine module version 1 (sha256 be006193...); the shipped engine is version 2, which
changed only the floating-point fit operators above N = 100, so a re-run of P2b uses version 2 and its recorded
engine hash will differ (`results/ERRATA.md`). P1 and P2a records contain per-slab / per-node timing fields, which
the comparison ignores. The REP1-REP3, REP5 and REP6/REP6b seats are shipped (`code/rep/`) with their results; their
re-runs need the seats' own intermediate state files (solutions and checkpoints), which are not shipped, so they are
not offered by the runner.

## What the runner does (and why each step is needed)

- **Flat work folder.** The instruments were written to sit next to their input and output files and to find each
  other by file name (`fluxt.py` imports `b3t.py`; `p2c.py`, `p2d.py`, `p2e.py` read earlier result files from
  their own folder). The runner copies `code/*.py`, `code/scout/`, `code/rep/` and the engine (`code/gravity_engine/`
  into `work/engine/`) and the input results from `results/` into one folder. The output being reproduced is never
  copied in.
- **Launch files.** Each `--main` refuses to start unless its launch file (`GO-P2e`, `GO-B3T`, `GO-FLUXT`, ...)
  exists. In the lab the lead wrote it after reading the controls; the runner writes it.
- **Self-sha gate.** `b3t.py`, `fluxt.py`, `p2d.py` and `p2e.py` also refuse unless a ledger file contains a `go` row
  naming the instrument with exactly its own sha256. The runner writes a one-file local ledger
  (`work/labdata/ledger.jsonl`) with one such row per gated instrument, using the pinned sha256, and points
  `LAB_DATA` at it. A modified instrument therefore still refuses to run.
- **Engine location.** Five instruments hold one string literal naming the lab's engine folder (`ENGINE_DIR = '...'`
  in p2b, p2c, p2d, p2e; `sys.path.insert(0, '...')` in rep5). Editing that line would break the pinned sha256, so
  the files are not edited: `tools/run_sealed.py` reads the source, points only that literal at `work/engine`,
  compiles it under the instrument's own file name and runs it; everything that hashes the file still sees the sealed
  bytes. It prints where `gravity_engine` was actually loaded from. `p2e.py` additionally checks the engine file
  against its pinned sha256 (6f3a4099...) before it starts.
- **Resume.** The P-series mains stop after a time budget (about 9 minutes) with exit code 3 and keep their state in
  `<output>.ckpt`; running the same command again continues. The runner repeats until exit code 0.
- **Compare.** Byte-for-byte first; if the bytes differ, field by field after removing the run-record fields.

## Running one instrument by hand

```
mkdir run && cd run
cp ../code/b3t.py . && python ../tools/pathpatch.py restore . && echo go > GO-B3T   # restore the sealed bytes first (PATCHES.md)
mkdir labdata && python -c "import json,hashlib;print(json.dumps(dict(kind='go',seal_id='10a65e940a',instruments=[dict(path='run/b3t.py',sha256=hashlib.sha256(open('b3t.py','rb').read()).hexdigest())])))" > labdata/ledger.jsonl
LAB_DATA=labdata python b3t.py --main result_b3t.json
sha256sum result_b3t.json     # 9a73e169...
```

For an instrument with an engine line, prefix the command with `python ../tools/run_sealed.py --engine ENGINE_DIR`
where `ENGINE_DIR` is a folder containing `gravity_engine/` (for example a copy of `code/`).

## Controls and engine tests

- The control modes (`--controls OUT.json`) start child processes of themselves and write into the lab's per-seat
  scratch folders; they are not offered by the runner. Their recorded outputs are shipped (`results/*/controls_*.json`,
  each with its all-pass flag), and the lab's GO step read them before every main run.
- `code/gravity_engine/tests/` are the engine's own tests as they stand in the lab. They read the lab's result files
  by absolute lab paths, so they do not run from this repository as shipped; they are included as the record of what
  the engine was tested against (P2a and P2b nodes byte for byte, REP4 resolution, the tangent certificate).
