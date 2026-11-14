#!/usr/bin/env python
"""reproduce.py -- re-run the certificates of this repository from a clean work folder and compare with results/.

Default (about 30 minutes on one core): pins, P2c ceiling re-assembly, B2V, B3T, FLUXT, P2e.
Optional (pass their names, or --all): p2d, p2c, p2b, p2a, p1, rep4, rep7 (see tools/reproduce.md for run times).

Usage (from the repository root):
    python tools/reproduce.py                       # default checks
    python tools/reproduce.py --only ceiling b3t    # a subset
    python tools/reproduce.py p2d                   # default checks plus an optional one
    python tools/reproduce.py --all                 # everything listed above (several hours)
Options: --work DIR (default ./work; must not exist or be empty), --python "CMD" (default: this interpreter).

What it does: copies the sealed instruments and the engine from code/ and the input results from results/ into a flat
work folder (the layout the instruments were written for), restores the sealed bytes of the path-patched copies
(PATCHES.md, tools/pathpatch.py), writes the launch files the instruments demand (GO-*) and
a local ledger file with one 'go' row per gated instrument naming its pinned sha256 (LAB_DATA points at it), runs each
check through tools/run_sealed.py (the files are never edited; see that file), repeats a run while it exits with 3
(resume after its time budget), and compares every output with the shipped file: first byte for byte (sha256), then,
if the bytes differ, field by field after removing the run-record fields named in VOLATILE below. A check passes when
the sha256 matches, or when the only differences are in those fields (reported as IDENTICAL-EXCEPT-RUN-RECORD).
Last line: [reproduce] ALL CHECKS PASS (or the list of failures). Exit code 0 on success.
"""
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOOLS = os.path.join(ROOT, 'tools')
sys.path.insert(0, TOOLS)
import pathpatch  # noqa: E402

# fields that record the run rather than the mathematics (timings, launch-row text, file locations)
VOLATILE = {'go', 'started', 'finished', 'time', 'timing', 'secs', 'seconds', 'total_seconds_last_invocation',
            'solve_seconds', 'file', 'at'}

# gated instruments: (file, seal id of the sealed claim it was launched under)
GATED = {'b3t.py': '10a65e940a', 'fluxt.py': '0ce319a05a', 'p2d.py': '9677804a51', 'p2e.py': '560847f6a4'}
GO_FILES = {'p1.py': 'GO-P1', 'p2a.py': 'GO-P2a', 'p2b.py': 'GO-P2b', 'p2c.py': 'GO-P2c', 'p2d.py': 'GO-P2d',
            'p2e.py': 'GO-P2e', 'b3t.py': 'GO-B3T', 'fluxt.py': 'GO-FLUXT'}

INPUTS = {  # shipped result files copied into the work folder (inputs of later stages; pinned by sha256)
    'result_p1.json': 'results/p1/result_p1.json', 'result_p2a.json': 'results/p2a/result_p2a.json',
    'result_p2b.json': 'results/p2b/result_p2b.json', 'result_p2c.json': 'results/p2c/result_p2c.json',
    'ceiling_p2c.json': 'results/p2c/ceiling_p2c.json', 'result_p2d.json': 'results/p2d/result_p2d.json',
    'result_p2e.json': 'results/p2e/result_p2e.json'}

# name: (instrument relative to the work folder, arguments, output file in the work folder, shipped file,
#        output files that must be absent before the run, expected minutes, default?)
CHECKS = [
    ('ceiling', 'p2c.py', ['--ceiling', 'ceiling_p2c.json', '--p2c-sha256',
                           '4b2dbf7b95742ce309a361fa693074779fb4f895f1f8afaecff5d874e60c33e7'],'ceiling_p2c.json', 'results/p2c/ceiling_p2c.json',
     ['ceiling_p2c.json'], 0.2, True),
    ('b2v', 'scout/b2_verdict.py', ['../result_p2e.json'], 'b2v_result.json', 'results/b2v/b2v_result.json', [], 0.1, True),
    ('b3t', 'b3t.py', ['--main', 'result_b3t.json'], 'result_b3t.json', 'results/b3t/result_b3t.json',
     ['result_b3t.json'], 0.5, True),
    ('fluxt', 'fluxt.py', [], 'result_fluxt.json', 'results/fluxt/result_fluxt.json', ['result_fluxt.json'], 0.5, True),
    ('p2e', 'p2e.py', ['--main', 'result_p2e.json'], 'result_p2e.json', 'results/p2e/result_p2e.json',
     ['result_p2e.json'], 27, True),
    ('p2d', 'p2d.py', ['--main', 'result_p2d.json'], 'result_p2d.json', 'results/p2d/result_p2d.json',
     ['result_p2d.json'], 15, False),
    ('p2c', 'p2c.py', ['--main', 'result_p2c.json'], 'result_p2c.json', 'results/p2c/result_p2c.json',
     ['result_p2c.json'], 18, False),
    ('p2b', 'p2b.py', ['--main', 'result_p2b.json'], 'result_p2b.json', 'results/p2b/result_p2b.json',
     ['result_p2b.json'], 25, False),
    ('p2a', 'p2a.py', ['--main', 'result_p2a.json'], 'result_p2a.json', 'results/p2a/result_p2a.json',
     ['result_p2a.json'], 60, False),
    ('p1', 'p1.py', ['--main', 'result_p1.json'], 'result_p1.json', 'results/p1/result_p1.json',
     ['result_p1.json'], 60, False),
    ('rep4', 'rep/rep4/rep4.py', [], 'rep/rep4/result_rep4.json', 'results/rep/rep4/result_rep4.json', [], 30, False),
    ('rep7', 'rep/rep7/rep7.py', ['run BY', 'run BY2', 'run BYP', 'run RN', 'run KN', 'assemble'], 'rep/rep7/result_rep7.json',
     'results/rep/rep7/result_rep7.json', [], 120, False),
]


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def strip(x, dropped, path=''):
    if isinstance(x, dict):
        out = {}
        for k, v in x.items():
            if k in VOLATILE or k.endswith('_seconds') or k.endswith('time_s') or k.startswith('secs'):
                dropped.add(path + k)
                continue
            out[k] = strip(v, dropped, path + k + '.')
        return out
    if isinstance(x, list):
        return [strip(v, dropped, path + '[].') for v in x]
    return x


def canon_sha(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def diff(a, b, path='', out=None, limit=8):
    out = [] if out is None else out
    if len(out) >= limit:
        return out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(path + k + ' (present in one file only)')
            else:
                diff(a[k], b[k], path + k + '.', out, limit)
    elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        for i, (u, v) in enumerate(zip(a, b)):
            diff(u, v, path + '[%d].' % i, out, limit)
    elif a != b:
        out.append(path.rstrip('.'))
    return out


def compare(produced, shipped):
    """shipped files are compared in their restored (original) form; see PATCHES.md"""
    pb = open(produced, 'rb').read()
    sb = pathpatch.restore_bytes(open(shipped, 'rb').read())
    hs, hp = hashlib.sha256(sb).hexdigest(), hashlib.sha256(pb).hexdigest()
    if hs == hp:
        return 'IDENTICAL', hp, hs, []
    try:
        A, B = json.loads(pb.decode('utf-8-sig')), json.loads(sb.decode('utf-8-sig'))
    except ValueError as e:
        return 'DIFFERENT (not JSON: %s)' % e, hp, hs, []
    dropped = set()
    a, b = strip(A, dropped), strip(B, dropped)
    if a == b:
        return 'IDENTICAL-EXCEPT-RUN-RECORD', canon_sha(a), canon_sha(b), sorted(dropped)
    return 'DIFFERENT', canon_sha(a), canon_sha(b), diff(a, b)


def setup(work):
    os.makedirs(work, exist_ok=True)
    if os.listdir(work):
        sys.exit('[reproduce] work folder %s is not empty' % work)
    code = os.path.join(ROOT, 'code')
    for name in os.listdir(code):
        if name.endswith('.py'):
            shutil.copy2(os.path.join(code, name), work)
    shutil.copytree(os.path.join(code, 'scout'), os.path.join(work, 'scout'))
    shutil.copytree(os.path.join(code, 'rep'), os.path.join(work, 'rep'))
    shutil.copytree(os.path.join(code, 'gravity_engine'), os.path.join(work, 'engine', 'gravity_engine'),
                    ignore=shutil.ignore_patterns('__pycache__'))
    for dst, src in INPUTS.items():
        shutil.copy2(os.path.join(ROOT, src), os.path.join(work, dst))
    # restore the sealed bytes of every path-patched copy (PATCHES.md), so pins and self-sha gates see the originals
    for d, _, fs in os.walk(work):
        for f in fs:
            pathpatch.restore_file(os.path.join(d, f))
    for f, g in GO_FILES.items():
        with open(os.path.join(work, g), 'w', encoding='utf-8') as fh:
            fh.write('written by tools/reproduce.py (the lab wrote this file after reading the controls)\n')
    lab = os.path.join(work, 'labdata')
    os.makedirs(lab)
    with open(os.path.join(lab, 'ledger.jsonl'), 'w', encoding='utf-8', newline='\n') as fh:
        for f, seal in GATED.items():
            p = os.path.join(work, f)
            fh.write(json.dumps(dict(kind='go', seal_id=seal, at='reproduce.py',
                                     instruments=[dict(path=p.replace('\\', '/'), sha256=sha(p))])) + '\n')


def run_check(work, py, chk, log):
    name, inst, args, out, shipped, absent, minutes, _ = chk
    for a in absent:
        p = os.path.join(work, a)
        for q in (p, p + '.ckpt', p + '.runinfo.json'):
            if os.path.exists(q):
                os.remove(q)
    env = dict(os.environ, LAB_DATA=os.path.join(work, 'labdata'), PYTHONDONTWRITEBYTECODE='1')
    ipath = os.path.join(work, inst)
    base = py + [os.path.join(TOOLS, 'run_sealed.py'), '--engine', os.path.join(work, 'engine'), ipath]
    cwd = os.path.dirname(ipath)
    t0 = time.time()
    runs = [a.split() for a in args] if name == 'rep7' else [args]
    for argv in runs:
        for attempt in range(200):
            if name == 'b2v':     # the instrument prints its JSON; store it with LF line ends (as shipped)
                r = subprocess.run(base + argv, cwd=cwd, env=env, stdout=subprocess.PIPE)
                rc = r.returncode
                with open(os.path.join(work, out), 'wb') as fo:
                    fo.write(r.stdout.replace(b'\r\n', b'\n'))
            else:
                rc = subprocess.run(base + argv, cwd=cwd, env=env).returncode
            if rc != 3:
                break
            log('  %s: exit 3 (time budget used), resuming' % name)
        if rc != 0:
            return dict(name=name, status='RUN FAILED (exit %d)' % rc, minutes=(time.time() - t0) / 60)
    status, hp, hs, info = compare(os.path.join(work, out), os.path.join(ROOT, shipped))
    return dict(name=name, status=status, produced=hp, expected=hs, info=info, minutes=(time.time() - t0) / 60,
                output=out, shipped=shipped)


def main(argv):
    work, py, only, extra, every = os.path.join(os.getcwd(), 'work'), [sys.executable], None, [], False
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == '--work':
            work = os.path.abspath(argv[i + 1]); i += 2; continue
        if a == '--python':
            py = shlex.split(argv[i + 1]); i += 2; continue
        if a == '--only':
            only = []
            i += 1
            while i < len(argv) and not argv[i].startswith('--'):
                only.append(argv[i]); i += 1
            continue
        if a == '--all':
            every = True; i += 1; continue
        extra.append(a); i += 1
    names = [c[0] for c in CHECKS]
    for n in (only or []) + extra:
        if n not in names:
            sys.exit('[reproduce] unknown check %s (known: %s)' % (n, ', '.join(names)))
    todo = [c for c in CHECKS if (c[0] in only if only is not None else (every or c[7] or c[0] in extra))]

    def log(s):
        print('[reproduce] ' + s, flush=True)

    rc = subprocess.run(py + [os.path.join(TOOLS, 'check_pins.py')]).returncode
    if rc != 0:
        log('pins FAILED: a shipped file differs from its pinned sha256; stopping')
        return 1
    log('pins: every sealed file matches its pinned sha256')
    setup(work)
    log('work folder: %s' % work)
    results = []
    for chk in todo:
        log('running %s (about %s min)' % (chk[0], chk[6]))
        r = run_check(work, py, chk, log)
        results.append(r)
        log('%s: %s in %.1f min' % (r['name'], r['status'], r['minutes']))
        if 'produced' in r:
            log('  sha256 produced %s' % r['produced'])
            log('  sha256 expected %s%s' % (r['expected'], '' if r['status'] == 'IDENTICAL'
                                             else ' (canonical JSON without run-record fields)'))
            if r['info']:
                log('  %s: %s' % ('run-record fields ignored' if r['status'].startswith('IDENTICAL') else 'differs at',
                                  ', '.join(r['info'][:12])))
    with open(os.path.join(work, 'reproduce_report.json'), 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=1)
    bad = [r['name'] for r in results if not r['status'].startswith('IDENTICAL')]
    log('ALL CHECKS PASS' if not bad else 'FAILED: ' + ', '.join(bad))
    return 0 if not bad else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
