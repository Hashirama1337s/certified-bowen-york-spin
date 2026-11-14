"""Tests for gravity_engine.bowen_york_sqrt (rebuilt from BHSLICE rep4/rep4.py, REP4 seal ad62cca292 HOLD (replication);
REFEREE-REP4 fixes; comparison lemma <LAB>/BHSLICE/LEMMA-REP4.md).

(a) REP4 resolution A (N, K, prec) = (140, 16, 224), REP4's 14 fixed pieces (adaptive=False, the unchanged arithmetic
    path): delta(3000) and delta(10000) give REP4's exact rationals (eta, W0, both endpoints) IDENTICALLY, and the
    interval contains REP4's final (A cap B) certified interval.
(b) K1: j = 2^-64, delta/j^2 contains 1/40; at (40, 4, 192) with 4 start pieces the fixed-piece run equals REP4's K1_A
    exactly, and the adaptive run (DESIGN-P2d rule) is certified and nested inside it.
(c) toy J = 0 (CS = 0) with a NONZERO dyadic W (known residual Gh = Lin[W], computed exactly): the forced margin
    eta = (2/3) * (exact sampled max |Lin W|) * (1 - 2^-10) must FAIL; the honest certificate's interval must contain the
    true delta = 0; the full-degree residual (8N+44, 8K) has coefficients beyond REP4's deleted CS = 0 branch degree
    (N+2, K) (the bug the toy now catches), and none beyond (8N+44, 8K) on a larger grid (DEG spot check). W = 0 gives
    REFUSED/TRIVIAL.
(d) null: the j = 10000 A ansatz bumped by 1e-9 (relative, P_0 coefficient) breaks the exact-rational sample signs with
    the honest eta, while the unbumped ansatz keeps them.
(e) overlap with the s-engine's sealed P2a intervals (<LAB>/BHSLICE/result_p2a.json) at j = 2 and j = 8
    ((64, 10, 192), REP4 fixed pieces; DESIGN-P2d band j <= 64).
Also: identity controls, outward display rounding, ctx.prec restored (doctrine #115), adaptive cap -> FAIL at j = 1e4.
Total about 2.5 min, peak under 1 GB. Run under memcap:
  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 -B -m pytest <LAB>/POCKETLAB/gravity_engine/tests/test_bowen_york_sqrt.py -q -p no:cacheprovider
"""
import sys
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # before any numpy import (deterministic fit)
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
sys.path.insert(0, '<LAB>/POCKETLAB')

import json
import random
from fractions import Fraction

import pytest
from flint import fmpq, arb, ctx

from gravity_engine import bowen_york_sqrt as by

REP4_RESULT = '<LAB>/BHSLICE/rep4/result_rep4.json'
P2A_RESULT = '<LAB>/BHSLICE/result_p2a.json'
Q = by.parse_q
QUIET = None


@pytest.fixture(scope='module')
def rep4():
    with open(REP4_RESULT) as f:
        return json.load(f)


@pytest.fixture(scope='module')
def nodes_A():
    """REP4 resolution A, fixed pieces; keep the exact ansatz (about 85 s for both)."""
    out = {}
    for j in (3000, 10000):
        out[j] = by.certify_point(j, *by.RES_A, n_pieces=8, adaptive=False, keep=True, log=QUIET)
    return out


# ------------------------------------------------------------------ (a) REP4 resolution A
@pytest.mark.parametrize('j', [3000, 10000])
def test_a_reproduce_rep4_A_identically(nodes_A, rep4, j):
    rec, (lo, hi), c, wq, CS = nodes_A[j]
    ref = rep4['runs_j%d' % j][0]
    assert ref['label'] == 'j%d_A' % j and (ref['N'], ref['K'], ref['prec_bits']) == by.RES_A
    assert rec['status'] == 'CERTIFIED' and rec['n_pieces_used'] == 14
    # unchanged arithmetic path => identical exact rationals
    assert rec['eta_exact'] == ref['eta_exact']
    assert rec['W0_exact'] == ref['W0_exact']
    assert rec['delta_lo_exact'] == ref['delta_lo_exact']
    assert rec['delta_hi_exact'] == ref['delta_hi_exact']
    # and it contains REP4's final certified interval (A cap B)
    fin = rep4['delta_j%d' % j]
    assert lo <= Q(fin['delta_lo_exact']) and Q(fin['delta_hi_exact']) <= hi
    # displayed ends are outward: floor <= lo, ceil >= hi
    assert Fraction(rec['delta_lo']) <= Fraction(int(lo.p), int(lo.q))
    assert Fraction(rec['delta_hi']) >= Fraction(int(hi.p), int(hi.q))


# ------------------------------------------------------------------ (b) K1
def test_b_K1_contains_one_fortieth(rep4):
    j = fmpq(1, 2 ** 64)
    ref = rep4['K1']['runs'][0]
    assert (ref['N'], ref['K'], ref['prec_bits'], ref['n_pieces']) == (40, 4, 192, 4)
    rec, (lo, hi), c, wq, CS = by.certify_point(j, 40, 4, 192, n_pieces=4, adaptive=False, keep=True, log=QUIET)
    assert rec['status'] == 'CERTIFIED'
    assert rec['delta_lo_exact'] == ref['delta_lo_exact'] and rec['delta_hi_exact'] == ref['delta_hi_exact']
    t = j * j
    assert lo / t <= fmpq(1, 40) <= hi / t
    recA, (loA, hiA), *_ = by.certify_point(j, 40, 4, 192, n_pieces=4, adaptive=True, keep=True, log=QUIET)
    assert recA['status'] == 'CERTIFIED' and recA['adaptive'] and recA['n_pieces_used'] > 10
    assert loA / t <= fmpq(1, 40) <= hiA / t
    assert lo <= loA and hiA <= hi                     # the adaptive enclosure is nested (narrower) here


# ------------------------------------------------------------------ (c) toy J = 0 with nonzero W
def _toy_W(N, K, seed=11):
    rng = random.Random(seed)
    return [fmpq(rng.randint(-2 ** 12, 2 ** 12), 2 ** 24) for _ in range((N + 1) * (K + 1))]


def test_c_toy_nonzero_W_at_CS0_must_fail():
    N, K, prec = 6, 2, 192
    wq = _toy_W(N, K)
    # known residual: at CS = 0, Gh[W] = Lin[W] exactly; its max over a rational grid is a lower bound of sup|Gh|
    Ex = by.ExactEval(wq, N, K)
    G_lo = max(abs(Ex.at(fmpq(i, 40), fmpq(l, 8))[1]) for i in range(41) for l in range(9))
    assert G_lo > 0
    honest = by.certify(wq, N, K, fmpq(0), n_pieces=2, prec=prec, adaptive=False)
    assert honest['sup_Gh_le'] >= G_lo                 # the rigorous bound is an upper bound of the known residual
    # forced margin below the known residual: the certificate must FAIL (strict check)
    forced = (2 * G_lo / 3) * (1 - fmpq(1, 2 ** 10))
    bad = by.certify(wq, N, K, fmpq(0), n_pieces=2, prec=prec, adaptive=False, eta=forced)
    assert bad['status'] == 'FAIL' and not bad['strict']
    # the honest certificate holds and must enclose the true delta = 0 (u = 0 solves J = 0)
    assert honest['status'] == 'CERTIFIED'
    a = by.A_MAP
    assert 2 * a * (honest['W0'] - honest['eta']) <= 0 <= 2 * a * (honest['W0'] + honest['eta'])
    assert honest['W0'] != 0                           # a genuinely nonzero ansatz at the read-off point
    # degree: with nonzero W at CS = 0 the residual Q = E^7 Lin[W] has coefficients beyond REP4's deleted branch
    # degree (N+2, K), and none beyond (8N+44, 8K) on a larger grid
    with by._prec(prec):
        ans = by.Ansatz(wq, N, K)
        dY, dV = by.residual_degrees(N, K)
        NQ, MQ = dY + 1 + 9, dV + 1 + 3
        ys = by.nodes_on(fmpq(0), fmpq(1), NQ)
        vs = by.nodes_on(fmpq(0), fmpq(1), MQ)
        W, LW = ans.F.grid(ys, vs)
        a2 = 2 * arb(by.A_MAP)
        vals = []
        for i in range(NQ):
            for l in range(MQ):
                E = a2 * (1 - ys[i]) ** 2 * (1 + ys[i] ** 2 * W[i, l]) + ys[i] ** 2
                vals.append(E ** 7 * LW[i, l])
        from flint import arb_mat
        C = by.cheb_coeffs(arb_mat(NQ, MQ, vals), NQ, MQ)
        beyond_true = [C[i, l] for i in range(NQ) for l in range(MQ) if i > dY or l > dV]
        assert all(x.contains(0) for x in beyond_true)
        beyond_legacy = max(abs(float(C[i, l].mid())) for i in range(NQ) for l in range(MQ)
                            if (i > N + 2 or l > K) and i <= dY and l <= dV)
        assert beyond_legacy > 1e-6


def test_c_toy_zero_W_refused():
    N, K = 6, 2
    wq = [fmpq(0)] * ((N + 1) * (K + 1))
    c = by.certify(wq, N, K, fmpq(0), n_pieces=2, prec=128, adaptive=False)
    assert c['status'] == 'REFUSED/TRIVIAL' and c['sup_Gh_le'] == 0 and not c['certified']
    rec = by.certify_point(0, 24, 3, 128, n_pieces=2, log=QUIET)
    assert rec['status'] == 'REFUSED/TRIVIAL' and rec['delta_lo_exact'] is None


# ------------------------------------------------------------------ (d) null: 1e-9 bump breaks the signs
def test_d_null_bump_breaks_signs(nodes_A):
    rec, (lo, hi), c, wq, CS = nodes_A[10000]
    N, K, _ = by.RES_A
    honest = by.sample_signs(wq, N, K, CS, c['eta'])
    assert honest['signs_ok']
    wb = list(wq)
    wb[0] += wq[0] * fmpq(1, 10 ** 9)
    bumped = by.sample_signs(wb, N, K, CS, c['eta'])
    assert not bumped['signs_ok']


# ------------------------------------------------------------------ (e) overlap with the s-engine (P2a)
@pytest.mark.parametrize('j', [2, 8])
def test_e_overlap_with_P2a(j):
    with open(P2A_RESULT) as f:
        ref = next(n for n in json.load(f)['nodes'] if n['j'] == '%d/1' % j)
    assert ref['status'] == 'CERTIFIED'
    rec = by.certify_point(j, 64, 10, 192, n_pieces=8, adaptive=False, log=QUIET)
    assert rec['status'] == 'CERTIFIED'
    lo, hi = Q(rec['delta_lo_exact']), Q(rec['delta_hi_exact'])
    plo, phi = Q(ref['delta_lo']), Q(ref['delta_hi'])
    assert max(lo, plo) <= min(hi, phi)
    print('j=%d: sqrt-map [%s, %s] vs P2a [%s, %s]; nested in P2a: %s' % (
        j, rec['delta_lo'], rec['delta_hi'], by.dec_floor(plo, 34), by.dec_ceil(phi, 34), plo <= lo and hi <= phi))


# ------------------------------------------------------------------ extras
def test_identity_controls():
    assert by.identity_controls()['all_ok']


def test_display_rounding_outward():
    rng = random.Random(7)
    for _ in range(300):
        q = fmpq(rng.randint(-10 ** 30, 10 ** 30), rng.randint(1, 10 ** rng.randint(1, 60)))
        if q == 0:
            continue
        F = Fraction(int(q.p), int(q.q))
        nd = rng.randint(1, 30)
        assert Fraction(by.dec_floor(q, nd)) <= F <= Fraction(by.dec_ceil(q, nd))
        assert Fraction(by.sci_down(q, 4)) <= F <= Fraction(by.sci_up(q, 4))
    assert by.sci_up(fmpq(99999, 10 ** 4), 4) == '1.000e+01'


def test_precision_restored():
    ctx.prec = 77
    try:
        by.certify(_toy_W(4, 1), 4, 1, by.CS_of(fmpq(1)), n_pieces=2, prec=160, adaptive=False)
        by.identity_controls()
        assert ctx.prec == 77
    finally:
        ctx.prec = 53


def test_adaptive_cap_fails_at_large_j():
    """the DESIGN-P2d rule needs ~590 boxes at j = 1e4 (N = 60): with the cap of 128 the certificate must FAIL, not
    silently fall back."""
    wq, rf, hist, reached, CS = by.solve_point(fmpq(10000), 60, 8, 192)
    c = by.certify(wq, 60, 8, CS, prec=192, adaptive=True)
    assert c['status'] == 'FAIL' and 'more than 128 boxes' in c['fail_reason']


def test_ansatz_shape_assert():
    with pytest.raises(AssertionError):
        by.Ansatz([fmpq(0)] * 5, 2, 1)
    with pytest.raises(AssertionError):
        by.check_cover([(fmpq(0), fmpq(1, 2), fmpq(0), fmpq(1))])
