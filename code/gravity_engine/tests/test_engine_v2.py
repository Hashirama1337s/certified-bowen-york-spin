"""Tests for engine version 2 of gravity_engine.lichnerowicz_axisym (fit-only changes; the proof path is unchanged).

(c)+(b) memory and high degree, in a fresh child process (its peak is the node's alone): j = 8000 at (N, K) = (152, 12),
        PREC 160, fit_ops 'mid', cache_ops False -> CERTIFIED, no NaN in the arb refinement, refinement reaches its
        2^-(PREC-40) stop, caches empty afterwards, peak private commit < 4 GB. Width/time/peak are printed (-s).
(a)     j = 1000 at (112, 14) with the default fit_ops ('auto' -> 'mid' at N = 112): the arb refinement reaches its
        2^-120 stop, and the certified interval lies INSIDE the sealed P2b interval
        [15.41318755062711644..., 15.41318839322056450...] with width <= the sealed width. (Two certified enclosures
        of the same number; the v2 one is nested in the sealed one. "Contains the sealed interval" AND "width <= sealed"
        together would force equality, so the test checks nesting + width, plus a non-empty intersection.)
legacy  fit_ops='legacy' reproduces the sealed P2b j = 1000 node record byte for byte (and still plateaus, as logged).
units   ENGINE_VERSION, 'auto' resolution, zero radii of the 'mid' operator, agreement of 'mid' and 'legacy' entries,
        cache_ops / clear_caches behaviour.
The v1 tests (test_lichnerowicz_axisym.py, P2a byte-for-byte at N = 40) run unchanged with the default 'auto'.
Run under memcap:  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 -B -m pytest <LAB>/POCKETLAB/gravity_engine/tests -q -p no:cacheprovider
"""
import sys
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # before any numpy import (deterministic fit)
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
sys.path.insert(0, '<LAB>/POCKETLAB')

import json
import subprocess
from fractions import Fraction

import pytest

from gravity_engine import lichnerowicz_axisym as la

RESULT_P2B = '<LAB>/BHSLICE/result_p2b.json'
CHILD = os.path.join(os.path.dirname(os.path.abspath(__file__)), '_v2_node_child.py')


def _frac(pq):
    a, b = pq.split('/')
    return Fraction(int(a), int(b))


@pytest.fixture(scope='module')
def p2b_1000():
    with open(RESULT_P2B) as f:
        d = json.load(f)
    return next(n for n in d['nodes'] if n['j'] == '1000/1')


# ------------------------------------------------------------------ (c) + (b): j = 8000 at (152, 12), fresh process
def test_j8000_N152_prec160_memory():
    r = subprocess.run([sys.executable, '-B', CHILD, '8000', '152', '12', 'mid', '0'],
                       capture_output=True, text=True, timeout=560)
    line = [x for x in r.stdout.splitlines() if x.startswith('RESULT ')]
    assert r.returncode == 0 and line, r.stderr[-2000:]
    out = json.loads(line[-1][7:])
    print('\n[v2] j=8000 (152,12):', json.dumps(out))
    assert out['prec'] == 160 and out['fit_ops'] == 'mid'
    assert out['status'] == 'CERTIFIED'
    assert not out['any_nan']
    assert out['reached_stop']
    assert out['caches_empty']
    assert out['peak_commit_gb'] is not None and out['peak_commit_gb'] < 4.0


# ------------------------------------------------------------------ (a): j = 1000 at (112, 14), default fit_ops
def test_j1000_N112_default_converges_and_nests(p2b_1000):
    P = la.BowenYorkSpinningPuncture(cache_ops=False)
    assert P.fit_ops_mode(112) == 'mid'
    rec, _ = la.certify_node(P, Fraction(1000), N=112, K=14, log=lambda m: None)
    fit = rec['_fit']
    print('\n[v2] j=1000 (112,14) arb_hist:', fit['arb_hist'], 'interval', rec['delta_lo_dec'], rec['delta_hi_dec'],
          'rel', rec['rel_width'])
    assert rec['status'] == 'CERTIFIED'
    assert fit['fit_ops'] == 'mid' and fit['reached_stop'] and fit['arb_hist'][-1] < 2.0 ** -120
    lo, hi = _frac(rec['delta_lo']), _frac(rec['delta_hi'])
    slo, shi = _frac(p2b_1000['delta_lo']), _frac(p2b_1000['delta_hi'])
    assert max(lo, slo) <= min(hi, shi)            # the two certified enclosures intersect
    assert slo <= lo and hi <= shi                  # v2 interval nested in the sealed P2b interval
    assert hi - lo <= shi - slo
    assert not P._fops and not P._aops


# ------------------------------------------------------------------ legacy: sealed P2b j = 1000 byte for byte
def test_legacy_reproduces_p2b_j1000(p2b_1000):
    P = la.BowenYorkSpinningPuncture(fit_ops='legacy', cache_ops=False)
    rec, _ = la.certify_node(P, Fraction(1000), N=112, K=14, log=lambda m: None)
    mine = {k: v for k, v in rec.items() if not k.startswith('_')}
    assert json.dumps(mine, sort_keys=True) == json.dumps(p2b_1000, sort_keys=True)
    assert rec['_fit']['fit_ops'] == 'legacy' and not rec['_fit']['reached_stop']   # the logged v1 plateau


# ------------------------------------------------------------------ units
def test_version_and_auto_resolution():
    assert la.ENGINE_VERSION == 2 and la.FIT_OPS == 'auto' and la.LEGACY_MAX_N == 100
    P = la.BowenYorkSpinningPuncture()
    assert P.cache_ops is True and P.fit_ops is None
    assert P.fit_ops_mode(40) == 'legacy' and P.fit_ops_mode(100) == 'legacy' and P.fit_ops_mode(101) == 'mid'
    assert la.BowenYorkSpinningPuncture(fit_ops='mid').fit_ops_mode(40) == 'mid'
    assert la.BowenYorkSpinningPuncture(fit_ops='legacy').fit_ops_mode(160) == 'legacy'
    with pytest.raises(ValueError):
        la.BowenYorkSpinningPuncture(fit_ops='bogus').fit_ops_mode(40)


def test_mid_operator_exact_and_close_to_legacy():
    N, K = 30, 4
    Pm = la.BowenYorkSpinningPuncture(fit_ops='mid')
    Pl = la.BowenYorkSpinningPuncture(fit_ops='legacy')
    Lm, Um, omm, srcm = Pm.arb_ops(N, K)
    Ll, Ul, oml, srcl = Pl.arb_ops(N, K)
    n = (N + 1) * (K + 1)
    assert Lm.nrows() == Lm.ncols() == n
    for M in (Lm, Um):
        assert all(M[i, k].rad() == 0 for i in range(n) for k in range(n))
    assert all(v.rad() == 0 for v in list(omm) + list(srcm))
    worst = 0.0
    for A, B in ((Lm, Ll), (Um, Ul)):
        for i in range(n):
            for k in range(n):
                d = abs(float((A[i, k] - B[i, k]).mid()))
                worst = max(worst, d / max(1.0, abs(float(B[i, k].mid()))))
    assert worst < 1e-30
    # cache key stays (N, K, PREC) (v1 subclasses check it); switching the construction rebuilds the entry
    assert (N, K, la.PREC) in Pm._aops and Pm._aops[(N, K, la.PREC)][4] == 'mid'
    Pm.fit_ops = 'legacy'
    Pm.arb_ops(N, K)
    assert Pm._aops[(N, K, la.PREC)][4] == 'legacy'


def test_cache_ops_and_clear_caches():
    P = la.BowenYorkSpinningPuncture()                   # v1 behaviour: caches kept per (N, K)
    P.float_ops(20, 3)
    P.arb_ops(20, 3)
    P.float_ops(22, 3)
    P.arb_ops(22, 3)
    assert len(P._fops) == 2 and len(P._aops) == 2
    P.clear_caches()
    assert not P._fops and not P._aops
    Q = la.BowenYorkSpinningPuncture(cache_ops=False)    # only the current degree is held
    Q.float_ops(20, 3)
    Q.arb_ops(20, 3)
    Q.float_ops(22, 3)
    Q.arb_ops(22, 3)
    assert list(Q._fops) == [(22, 3)] and list(Q._aops) == [(22, 3, la.PREC)]
    rec, _ = la.certify_node(Q, Fraction(2), log=lambda m: None)   # fit_u frees the operators before the proof
    assert rec['status'] == 'CERTIFIED' and not Q._fops and not Q._aops
