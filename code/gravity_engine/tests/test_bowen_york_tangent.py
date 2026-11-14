"""Tests for the tangent certificate of gravity_engine.bowen_york_sqrt (module version 2; DESIGN-P2d sec. 2 and 4.3,
AMEND-P2d-1, AMEND-P2d-2; Lemma TL in <LAB>/BHSLICE/LEMMA-REP4.md sec. 6) and for REFEREE-LEMMA FIX F.

(1) K1T: j = 2^-64 at (N, K) = (64, 8), 8 start pieces, 192 and 256 bits: the point certificate contains 1/40 and the
    tangent certificate gives (1 - 2L_hi)/j^2 within 1e-9 of 3/40 (a CONSISTENCY control of tangent vs point problem,
    AMEND-P2d-1 item 6; DESIGN-P2d 4.4 tolerance), and <= 3/40 + 1e-9.
(2) N-7: the same with the tangent problem built with 6 in place of 7 (fit and bound) must miss 3/40 (K1T fails).
(3) psi = 0: 2L_hi = 1 exactly, eta_Z = 0, CERTIFIED; the TL tail with 2L = 1 equals P2c's X1 = j/(4 delta_lo^2) and the
    TL cell with 2L = 1 equals P2c's plain cell j_b/(4(1 + sqrt(j_b/j_a) delta_lo)^2) exactly (ratio-4 cell, s = 1/2).
(4) tangent nulls that must FAIL: forced margin eta_Z = (2/3) sup (1 - 2^-10) (near threshold), and Z_N bumped by
    1e-6 max|Z| with the honest eta_Z.
(5) structural asserts: K = 0, wrong Z shape, cell/tail bracket <= 0, s_lo not a lower bound of sqrt(j_a/j_b).
(6) NaN / inf / huge tangent fit -> INVALID.
(7) degree spot check: Q_L of a real fit on a larger grid has no coefficient beyond (9N+50, 9K) (balls contain 0).
(8) FIX F regression: N = 2, K = 0, W = 0, j = 1: sup_Gh_le must be >= 9/32 (exact value of |Gh| at Y = 1/2, v = 0);
    before the patch the module gave 0.2552 (v-degree 8K = 0 under-sampled (1-v)).
(9) ctx.prec restored.
About 1 min, peak under 0.5 GB:
  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 -B -m pytest <LAB>/POCKETLAB/gravity_engine/tests/test_bowen_york_tangent.py -q -p no:cacheprovider
"""
import sys
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
sys.path.insert(0, '<LAB>/POCKETLAB')

import pytest
from flint import fmpq, arb, arb_mat, ctx

from gravity_engine import bowen_york_sqrt as by

QUIET = None
J_K1 = fmpq(1, 2 ** 64)


def _node(j, N, K, prec, n_pieces=8):
    rec, (lo, hi), c, wq, CS = by.certify_point(j, N, K, prec, n_pieces=n_pieces, keep=True, log=QUIET)
    assert rec['status'] == 'CERTIFIED'
    return rec, lo, hi, c, wq, CS


@pytest.fixture(scope='module')
def node8():
    return _node(fmpq(8), 40, 6, 192)


# ------------------------------------------------------------------ (1) K1T, (2) N-7
@pytest.mark.parametrize('prec', [192, 256])
def test_K1T_consistent_with_3_over_40(prec):
    rec, lo, hi, c, wq, CS = _node(J_K1, 64, 8, prec)
    t2 = J_K1 * J_K1
    assert lo / t2 <= fmpq(1, 40) <= hi / t2                      # K1
    T = by.certify_tangent(wq, 64, 8, CS, c['eta'], prec=prec)
    assert T['status'] == 'CERTIFIED' and T['eta_Z'] > 0
    x = (1 - T['twoL_hi']) / t2
    assert abs(x - fmpq(3, 40)) <= fmpq(1, 10 ** 9)               # K1T
    assert x <= fmpq(3, 40) + fmpq(1, 10 ** 9)
    # the scaled dyadic grid keeps Z (~ (3/40) 2^-128) alive (REFEREE-TL sec. 4 (v))
    assert T['Z0'] != 0 and T['refine_target_reached']


def test_N7_wrong_tangent_misses_3_over_40():
    rec, lo, hi, c, wq, CS = _node(J_K1, 64, 8, 192)
    T = by.certify_tangent(wq, 64, 8, CS, c['eta'], prec=192, factor=6)
    x = (1 - T['twoL_hi']) / (J_K1 * J_K1)
    assert T['status'] == 'CERTIFIED'                             # certified for the WRONG problem ...
    assert abs(x - fmpq(3, 40)) > fmpq(1, 1000)                   # ... and K1T catches it (6/7 * 3/40 = 9/140)


# ------------------------------------------------------------------ (3) psi = 0
def test_psi_zero_gives_2L_one_and_P2c_cells(node8):
    rec, lo, hi, c, wq, CS = node8
    T = by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=192, psi_zero=True)
    assert T['status'] == 'CERTIFIED' and T['eta_Z'] == 0 and T['twoL_hi'] == 1 and T['Z0'] == 0
    j = fmpq(8)
    assert by.tl_tail_bound(j, lo, T['twoL_hi']) == j / (4 * lo ** 2)                         # P2c X1
    jb = 4 * j
    assert by.tl_cell_bound(j, jb, fmpq(1, 2), lo, T['twoL_hi']) == jb / (4 * (1 + 2 * lo) ** 2)   # P2c plain cell
    # an honest tangent beats psi = 0 (2L < 1 at j = 8) and its tail is smaller
    Th = by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=192)
    assert Th['status'] == 'CERTIFIED' and Th['twoL_hi'] < 1
    assert by.tl_tail_bound(j, lo, Th['twoL_hi']) < by.tl_tail_bound(j, lo, T['twoL_hi'])


# ------------------------------------------------------------------ (4) nulls
def test_tangent_nulls_fail(node8):
    rec, lo, hi, c, wq, CS = node8
    T = by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=192)
    assert T['status'] == 'CERTIFIED'
    forced = fmpq(2, 3) * T['sup_H_le'] * (1 - fmpq(1, 2 ** 10))
    F = by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=192, eta_Z=forced)
    assert F['status'] == 'FAIL' and not F['strict'] and F['sup_H_le'] == T['sup_H_le']
    zq = list(T['zq'])
    mx = max(abs(x) for x in zq)
    zq[0] = zq[0] + mx / 10 ** 6
    B = by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=192, zq_override=zq, eta_Z=T['eta_Z'])
    assert B['status'] == 'FAIL' and B['sup_H_le'] > T['sup_H_le']


# ------------------------------------------------------------------ (5) structural asserts
def test_structural_asserts(node8):
    rec, lo, hi, c, wq, CS = node8
    with pytest.raises(AssertionError):
        by.tangent_degrees(10, 0)
    w0 = [fmpq(0)] * 11
    with pytest.raises(AssertionError):
        by.certify_tangent(w0, 10, 0, CS, fmpq(0), prec=128)                 # K = 0 refused
    with pytest.raises(AssertionError):
        by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=128, zq_override=[fmpq(0)] * 5)   # wrong Z shape
    j = fmpq(8)
    with pytest.raises(AssertionError):
        by.tl_tail_bound(j, lo, 1 + lo)                                       # bracket 0
    with pytest.raises(AssertionError):
        by.tl_cell_bound(j, 2 * j, fmpq(7, 10), lo, 10 * (1 + lo))            # bracket < 0
    with pytest.raises(AssertionError):
        by.tl_cell_bound(j, 2 * j, fmpq(71, 100), lo, fmpq(1, 2))             # 0.71 > sqrt(1/2)
    assert by.tl_cell_bound(j, 2 * j, fmpq(70, 100), lo, fmpq(1, 2)) > 0
    assert by.tangent_degrees(10, 1) == (140, 9)


# ------------------------------------------------------------------ (6) NaN / inf
def test_nan_inf_invalid(node8):
    rec, lo, hi, c, wq, CS = node8
    n = len(wq)
    for bad in (float('nan'), float('inf'), -float('inf'), 1e300):
        z = [0.0] * n
        z[3] = bad
        T = by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=128, z_override=z)
        assert T['status'] == 'INVALID' and not T['certified'] and T['twoL_hi'] is None


# ------------------------------------------------------------------ (7) degree spot check
def test_tangent_degree_spot_check():
    N, K, prec = 6, 2, 192
    rec, lo, hi, c, wq, CS = _node(fmpq(8), N, K, prec, n_pieces=2)
    T = by.certify_tangent(wq, N, K, CS, c['eta'], prec=prec, n_pieces=2)
    dY, dV = by.tangent_degrees(N, K)
    with by._prec(prec):
        wm = list(wq)
        wm[0] -= c['eta']
        aW, aZ = by.Ansatz(wm, N, K), by.Ansatz(T['zq'], N, K)
        for (y0, y1) in ((fmpq(0), fmpq(1)), (fmpq(1, 8), fmpq(1, 4))):
            NQ, MQ = dY + 1 + 25, dV + 1 + 7
            ys, vs = by.nodes_on(y0, y1, NQ), by.nodes_on(fmpq(0), fmpq(1), MQ)
            Wg, _ = aW.F.grid(ys, vs)
            Zg, LZ = aZ.F.grid(ys, vs)
            a2, cs = 2 * arb(by.A_MAP), arb(CS)
            vals = []
            for i in range(NQ):
                y = ys[i]
                c2 = a2 * (1 - y) ** 2
                for l in range(MQ):
                    E = c2 * (1 + y * y * Wg[i, l]) + y * y
                    vals.append(E ** 8 * LZ[i, l] - 7 * cs * (1 - vs[l]) * y ** 7 * (1 - y) ** 5 * (1 + c2 * Zg[i, l]))
            C = by.cheb_coeffs(arb_mat(NQ, MQ, vals), NQ, MQ)
            assert all(C[i, l].contains(0) for i in range(NQ) for l in range(MQ) if i > dY or l > dV)
            assert max(abs(float(C[i, l].mid())) for i in range(dY - 3, dY + 1) for l in range(MQ)) > 0


# ------------------------------------------------------------------ (8) FIX F regression
def test_fixF_K0_source_v_degree():
    N, K = 2, 0
    CS = by.CS_of(fmpq(1))
    assert by.residual_degrees(N, K) == (60, 1)
    c = by.certify([fmpq(0)] * (N + 1), N, K, CS, n_pieces=2, prec=128, adaptive=False)
    y = fmpq(1, 2)
    exact = CS * y ** 5 * (1 - y) ** 5 / ((1 - y) ** 2 + y * y) ** 7         # |Gh[0]| at Y = 1/2, v = 0
    assert exact == fmpq(9, 32)
    assert c['sup_Gh_le'] >= exact                                             # was 0.2552 before FIX F
    s = by.sample_signs([fmpq(0)] * (N + 1), N, K, CS, c['eta'], n_y=40, n_v=4)
    assert s['signs_ok'] if c['status'] == 'CERTIFIED' else True


# ------------------------------------------------------------------ (9) precision
def test_precision_restored(node8):
    rec, lo, hi, c, wq, CS = node8
    ctx.prec = 77
    try:
        by.certify_tangent(wq, 40, 6, CS, c['eta'], prec=160, psi_zero=True)
        assert ctx.prec == 77
    finally:
        ctx.prec = 53
