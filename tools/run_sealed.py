#!/usr/bin/env python
"""run_sealed.py -- run a sealed instrument byte-for-byte unchanged, with its engine import pointed at this repository.

Some instruments (code/p2b.py, p2c.py, p2d.py, p2e.py, rep/rep5/rep5.py) were written to import the shared engine
from the lab's own folder: they hold one string literal naming that folder (the line `ENGINE_DIR = '...'`, or for rep5
the line `sys.path.insert(0, '...')`). Editing that line would change the file's sha256, which the lab ledger pinned
before the certified run and which the instruments check themselves (self-sha gate). So the file on disk is never
edited. Instead this launcher reads the source, replaces only that one string value with the folder given by
--engine (it must contain gravity_engine/), compiles the result under the instrument's own file name and runs it as
__main__ (the value it replaces may be the lab folder or, in a shipped path-patched copy, its placeholder; see
PATCHES.md). Everything that reads the file itself (sha256 of __file__, the GO gate, recorded instrument hashes) still
sees the untouched sealed bytes. The substitution is printed on stderr. Instruments without such a line run as is.

Usage:  python tools/run_sealed.py --engine DIR INSTRUMENT.py [instrument arguments ...]
"""
import os
import re
import sys

PATTERNS = (re.compile(r"^(ENGINE_DIR\s*=\s*)'([^']*)'", re.M),
            re.compile(r"^(sys\.path\.insert\(0,\s*)'([^']*)'(\))", re.M))


def main(argv):
    if len(argv) < 3 or argv[0] != '--engine':
        sys.exit(__doc__)
    engine = os.path.abspath(argv[1]).replace('\\', '/')
    if not os.path.isfile(os.path.join(engine, 'gravity_engine', '__init__.py')):
        sys.exit('run_sealed: %s has no gravity_engine/__init__.py' % engine)
    path = os.path.abspath(argv[2])
    with open(path, encoding='utf-8') as f:
        src = f.read()
    n = 0
    for pat in PATTERNS:
        def sub(m):
            return m.group(1) + repr(engine) + (m.group(3) if m.lastindex >= 3 else '')
        src, k = pat.subn(sub, src)
        n += k
    if n > 1:
        sys.exit('run_sealed: more than one engine-path line in %s; refusing' % path)
    print('[run_sealed] %s: %s' % (os.path.basename(path),
                                    'engine path line pointed at %s (file on disk unchanged)' % engine if n
                                    else 'no engine path line; run unchanged'), file=sys.stderr, flush=True)
    sys.argv = [path] + argv[3:]
    sys.path.insert(0, os.path.dirname(path))
    code = compile(src, path, 'exec')
    g = {'__name__': '__main__', '__file__': path, '__builtins__': __builtins__}
    try:
        exec(code, g)
    finally:
        ge = sys.modules.get('gravity_engine')
        if ge is not None:
            print('[run_sealed] gravity_engine loaded from %s' % os.path.dirname(os.path.abspath(ge.__file__)),
                  file=sys.stderr, flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
