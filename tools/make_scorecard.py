#!/usr/bin/env python
"""make_scorecard.py -- build scorecard.svg and SCORECARD.md from the verdict files in results/.

Every sealed test of this project is listed. The test-level outcome (HOLD / FAILED) is read from the lab ledger export
results/ledger-bhslice.json. The per-prediction verdicts are read from the scoring table of each VERDICT file (rows whose
first cell is a prediction id such as B1, D3, E4 and whose last cell starts with HOLD / FAILED / PASS); for tests whose
record has no VERDICT file, from the scored text the lab ledger recorded for the test, from the bar fields of the result
file (B3T), or from the scorecard table the instrument wrote (FLUXT). FAILED rows are never dropped or moved. One
pre-stated but unsealed probe (B4x PROBE-1) is listed in its own section and is not counted among the sealed tests.
Run from anywhere:  python tools/make_scorecard.py
"""
import json
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, 'results')

# (section, [(test id, seal, source kind, source file under results/, one-line description)])
TESTS = [
    ('CERTIFICATE RUNS (INTERVAL ARITHMETIC; EVERY SPIN IN THE RANGE)', [
        ('P1', '565736c4c4', 'verdict', 'p1/VERDICT-P1.md', 'Small spins J in [0, 1/4] (j in [0, 1])'),
        ('P2a', 'd62bfc0e0b', 'verdict', 'p2a/VERDICT-P2a.md', 'Every j in [1, 10]'),
        ('P2b', 'eba765ebe0', 'verdict', 'p2b/VERDICT-P2b.md', 'Every j in [10, 1000] (shared engine)'),
        ('P2c', '3031217ac8', 'verdict', 'p2c/VERDICT-P2c.md', 'Uniform ceiling J/M^2 <= c < 1, c <= 0.98'),
        ('P2d', '9677804a51', 'verdict', 'p2d/VERDICT-P2d.md', 'Sharp ceiling via the tangent lemma'),
        ('P2e', '560847f6a4', 'verdict', 'p2e/VERDICT-P2e.md', 'Decide the sixth digit (rounds to 0.928200)'),
        ('B3T', '10a65e940a', 'b3t', 'b3t/result_b3t.json', 'Spin of the extreme (mu = 0) data, on the throat'),
    ]),
    ('BLIND REPLICATIONS (INDEPENDENT SEATS, OWN FORMULATIONS)', [
        ('REP1', 'e44dd02296', 'verdict', 'rep/rep1/VERDICT-REP1.md', 'Blind delta(0.0625), delta(0.25) vs P1'),
        ('REP2', '227cbcafc2', 'verdict', 'rep/rep2/VERDICT-REP2.md', 'Blind delta(6), delta(10) vs P2a'),
        ('REP3', 'a669f700df', 'verdict', 'rep/rep3/VERDICT-REP3.md', 'Blind delta(300), delta(1000) vs P2b'),
        ('REP4', 'ad62cca292', 'ledger', 'rep/rep4/REPORT-REP4.md', 'Blind delta(3000), delta(10000) vs P2c'),
        ('REP5', '9e0df55e7a', 'ledger', 'rep/rep5/RESULT.md', 'Blind tangent (slope) bounds vs P2d'),
        ('REP6', '0405395ed8', 'ledger', 'rep/rep6/RESULT.md', 'Blind top node 2^21 and tail vs P2e'),
        ('REP6b', '3e929bb723', 'ledger', 'rep/rep6/RESULT.md', 'Tighter blind tangent at 2^21 (re-test)'),
        ('REP7', '133b99d241', 'ledger', 'rep/rep7/RESULT.md', 'Blind extreme-data spin vs B3T'),
    ]),
    ('TESTS ON THE CERTIFIED DATA', [
        ('B2V', 'a580a8ad0f', 'ledger', 'b2v/b2v_result.json', 'Monotone, convex, tangent-free tail, exponent'),
        ('FLUXT', '0ce319a05a', 'fluxt', 'fluxt/FLUXT-REPORT.md', 'Throat flux identity, slice by slice'),
    ]),
]
UNSEALED = ('PRE-STATED BUT UNSEALED PROBE (NOT COUNTED ABOVE)', [
    ('B4x P-1', 'unsealed', 'probes/b4x/b4x_mono_probe.log',
     'Monotone throat profile (double precision)',
     [('FAILED', 'P1 FAILED'), ('FAILED', 'P2 FAILED'), ('UNSEALED', 'pre-stated, unsealed')],
     'P1 FAILED; P2 FAILED as written. Bar stated in the script header before the run, never hashed in the ledger. '
     'The 1e-12 bar sat below the solve\'s own residual (~1e-7), so the probe did not test the claim; post hoc the '
     'negative differences shrink with resolution (results/probes/b4x/b4x_mono_posthoc.log). No re-bar was run.'),
])

ID_RE = re.compile(r'^\s*([A-Z]{1,2}\d+)\b')
LEDGER_ITEM = re.compile(r'\b([A-Z]\d+) (HOLD|FAILED)\b')


def ledger():
    d = json.load(open(os.path.join(R, 'ledger-bhslice.json'), encoding='utf-8'))
    return {t['id']: t for t in d['sealed_claims']}


def predictions(path):
    rows = []
    for line in open(os.path.join(R, path), encoding='utf-8'):
        if not line.startswith('|'):
            continue
        cells = [c.strip() for c in line.strip().strip('|').split('|')]
        m = ID_RE.match(cells[0])
        if not m:
            continue
        last = [c for c in cells if c][-1]
        v = re.match(r'(HOLD|FAILED|PASS|VOID)', last)
        if v:
            rows.append((m.group(1), v.group(1)))
    return rows


def from_ledger(t):
    txt = t.get('measured_corrected') or t.get('measured') or ''
    seen, rows = set(), []
    for k, v in LEDGER_ITEM.findall(txt):
        if k not in seen:
            seen.add(k)
            rows.append((k, v))
    return rows


def from_b3t(path):
    bar = json.load(open(os.path.join(R, path), encoding='utf-8'))['bar']
    keys = [('T1', 'T1_certificate'), ('T2', 'T2_width_le_1e9'), ('T3', 'T3_inside_P2e'), ('T4', 'T4_within_1e9_of_probe')]
    return [(k, 'HOLD' if bar[f] is True else 'FAILED') for k, f in keys]


def from_fluxt(path):
    rows = []
    for line in open(os.path.join(R, path), encoding='utf-8'):
        m = re.match(r'\|\s*(F\d)\b.*\|\s*(HOLD|FAILED|DRIFTS|FAILED \(no drift\))\s*\|\s*$', line)
        if m:
            rows.append((m.group(1), {'DRIFTS': 'PASS'}.get(m.group(2), m.group(2).split()[0])))
    return rows


def compress(ids):
    out, i = [], 0
    while i < len(ids):
        j = i
        while (j + 1 < len(ids) and ids[j + 1][0] == ids[i][0] and ids[j + 1][1:].isdigit() and ids[j][1:].isdigit()
               and int(ids[j + 1][1:]) == int(ids[j][1:]) + 1):
            j += 1
        out.append(ids[i] if j - i < 2 else ids[i] + '-' + ids[j])
        if j - i == 1:
            out.append(ids[j])
        i = j + 1
    return ' '.join(out)


FONT = 'Segoe UI, Helvetica, Arial, sans-serif'
INK, MUTED, ACCENT, BAND, LINE, BORDER = '#18192b', '#555b74', '#3a3fa8', '#eef0fa', '#eceef5', '#d9dbe8'
CHIP = {'HOLD': ('#dff3e8', '#1f7a4d'), 'PASS': ('#dff3e8', '#1f7a4d'), 'FAILED': ('#fbe3e1', '#b3322a'),
        'VOID': ('#fbf0d4', '#8a6100'), 'UNSEALED': ('#fbf0d4', '#8a6100'), 'REPLICATED': ('#dff3e8', '#1f7a4d')}
NEON = {'#ffffff': '#1c2026', '#d9dbe8': '#2f353d', '#18192b': '#e8edf2', '#555b74': '#9aa4ae', '#eef0fa': '#262b31',
        '#3a3fa8': '#39ff14', '#eceef5': '#22272e', '#dff3e8': '#12301c', '#1f7a4d': '#39ff14', '#fbe3e1': '#3a1a1a',
        '#b3322a': '#ff6b6b', '#fbf0d4': '#33290f', '#8a6100': '#ffc94d'}


def tw(s, per=7.4):
    return round(len(s) * per)


def esc(s):
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def build():
    led = ledger()
    sections, md_rows, counts = [], [], {'HOLD': 0, 'FAILED': 0}
    for head, tests in TESTS:
        rows = []
        for tid, seal, kind, src, desc in tests:
            t = led.get(seal)
            if not t or t.get('outcome') not in counts:
                sys.exit('ledger outcome missing for %s %s' % (tid, seal))
            outcome = t['outcome']
            counts[outcome] += 1
            preds = {'verdict': lambda: predictions(src), 'ledger': lambda: from_ledger(t),
                     'b3t': lambda: from_b3t(src), 'fluxt': lambda: from_fluxt(src)}[kind]()
            chips = []
            for key in ('HOLD', 'PASS', 'FAILED'):
                ids = [p for p, v in preds if v == key]
                if ids:
                    chips.append((key, compress(ids) + ' ' + key))
            if preds:
                derived = 'FAILED' if any(v == 'FAILED' for _, v in preds) else 'HOLD'
                if derived != outcome:
                    sys.exit('verdict source and ledger disagree for %s %s' % (tid, seal))
                detail = '; '.join('%s %s' % pv for pv in preds)
            else:
                chips.append(('REPLICATED', 'REPLICATED'))
                detail = 'not itemised; recorded score: ' + (t.get('measured_corrected') or t['measured']).split(';')[0][:150]
            if t.get('corrections'):
                detail += ' (recorded text later corrected: %s)' % t['corrections'][-1]['reason'][:120]
            rows.append((tid + ' ' + seal, desc, chips, outcome))
            md_rows.append((head, tid, seal, desc, outcome, detail, src))
        sections.append((head, rows))
    uhead, utests = UNSEALED
    urows = [(tid + ' ' + seal, desc, chips, None) for tid, seal, src, desc, chips, _ in utests]
    sections.append((uhead, urows))
    n = counts['HOLD'] + counts['FAILED']
    title = 'Scorecard: %d sealed tests (%d HOLD, %d FAILED) + 1 unsealed probe' % (n, counts['HOLD'], counts['FAILED'])
    sub = 'Sealed = prediction and pass bar hashed before the run. Failures are shown as failures. Details: SCORECARD.md.'
    W, y, body = 900, 66, []
    for head, rows in sections:
        body.append(f'<rect x="10" y="{y}" width="{W - 20}" height="22" fill="{BAND}"/>')
        body.append(f'<text x="20" y="{y + 16}" font-size="11.5" font-weight="600" fill="{ACCENT}" letter-spacing="0.5">{esc(head)}</text>')
        y += 22
        for tid, what, chips, outcome in rows:
            body.append(f'<text x="20" y="{y + 22}" font-size="12.5" font-family="Consolas, Menlo, monospace" fill="{INK}">{esc(tid)}</text>')
            body.append(f'<text x="170" y="{y + 22}" font-size="12.5" fill="{INK}">{esc(what)}</text>')
            x = 530.0
            for key, label in chips:
                fill, fg = CHIP[key]
                t = tw(label)
                body.append(f'<rect x="{x}" y="{y + 8}" width="{t + 20}" height="20" rx="10" fill="{fill}"/>')
                body.append(f'<text x="{x + 10}" y="{y + 22}" font-size="11.5" font-weight="600" fill="{fg}" '
                            f'textLength="{t}" lengthAdjust="spacingAndGlyphs">{esc(label)}</text>')
                x += t + 26
            y += 30
            body.append(f'<line x1="10" y1="{y + 1}" x2="{W - 10}" y2="{y + 1}" stroke="{LINE}"/>')
        y += 8
    H = y + 6
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" font-family="{FONT}" '
           f'role="img" aria-label="{esc(title)}">\n<title>{esc(title)}</title>\n'
           f'<rect width="{W}" height="{H}" rx="10" fill="#ffffff" stroke="{BORDER}"/>\n'
           f'<text x="20" y="34" font-size="18" font-weight="600" fill="{INK}">{esc(title)}</text>\n'
           f'<text x="20" y="54" font-size="12" fill="{MUTED}">{esc(sub)}</text>\n' + '\n'.join(body) + '\n</svg>\n')
    pat = re.compile('|'.join(re.escape(k) for k in NEON), re.I)
    svg = pat.sub(lambda m: NEON[m.group(0).lower()], svg)
    with open(os.path.join(ROOT, 'scorecard.svg'), 'w', encoding='utf-8', newline='\n') as f:
        f.write(svg)
    chip = {'HOLD': '![HOLD](brand/chips/hold.svg)', 'FAILED': '![FAILED](brand/chips/failed.svg)'}
    md = ['# Scorecard', '', title + '.', '',
          'Generated by `tools/make_scorecard.py` from the verdict files in `results/` and the ledger export',
          '`results/ledger-bhslice.json`. A sealed test is one whose predictions and pass bars were hashed before its run.', '',
          '| Section | Test | Seal | What was predicted | Verdict | Per prediction (as scored) | Source |',
          '|---|---|---|---|---|---|---|']
    for head, tid, seal, desc, outcome, detail, src in md_rows:
        sec = head.split(' (')[0].capitalize()
        md.append(f'| {sec} | {tid} | `{seal}` | {desc} | {chip[outcome]} | {detail} | [results/{src}](results/{src}) |')
    for tid, seal, src, desc, chips, detail in utests:
        md.append(f'| Unsealed probe | {tid} | (none) | {desc} | pre-stated, unsealed (not counted) | {detail} | '
                  f'[results/{src}](results/{src}) |')
    md += ['', 'Notes:',
           '- REP6 is FAILED as sealed: its blind point enclosure at j = 2^21 contained P2e\'s (R1 HOLD), but its blind',
           '  tangent certificate was 2380x looser than P2e\'s (R2 FAILED) and its blind tail bound 0.928235 missed the',
           '  0.9282005 threshold (R3 FAILED). No contradiction (all one-sided bounds are consistent). The same seat then ran',
           '  REP6b under a new seal with a higher-resolution sub-solution: tail 0.9282004721 < 0.9282005 (HOLD).',
           '- P2d and P2e: the recorded result text was corrected after referee rounds (display directions, wording); the',
           '  outcomes did not change and the exact rationals in the result files were always the certified values',
           '  (results/ERRATA.md).',
           '- REP4, REP5, REP6b, REP7 were scored by the lead as REPLICATED without itemised prediction ids; the recorded',
           '  score text is in results/ledger-bhslice.json and the seat reports in results/rep/.',
           '- FLUXT F3 is a null that must drift (> 100 tol); PASS means it drifted as required.',
           '- B4x PROBE-1 was a double-precision probe whose bar was written in the script header but never sealed; it is',
           '  listed so that its failure is visible. It certifies nothing and no claim rests on it.']
    with open(os.path.join(ROOT, 'SCORECARD.md'), 'w', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(md) + '\n')
    print(title)


if __name__ == '__main__':
    build()
