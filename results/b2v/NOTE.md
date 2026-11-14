# Note on b2v_result.json

`b2v_result.json` is the output of `code/scout/b2_verdict.py` (pinned sha256 2be6e953...) run on
`results/p2e/result_p2e.json` (sha256 2aa80d5f...), as named in `BAR-B2V.md`. The instrument prints one JSON object.

The lab's copy of this file (sha256 3a6a3c75528917879021f20d8a44488c0c2147f1502af8262dbaaa6eefd9ef3c) was written
through a shell redirect, which added a byte-order mark, CRLF line ends and, after the JSON object, the shell's
record of the memory-cap wrapper's exit line. The copy here keeps the JSON object only, with LF line ends
(sha256 7bbec769ed58b7482918247d04c5dc531abe7a779105f89179efd43ce644d6b8). Its content is unchanged: a fresh run of
the instrument on the shipped P2e result reproduces it byte for byte (`python tools/reproduce.py --only b2v`).
