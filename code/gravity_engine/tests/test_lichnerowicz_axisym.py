"""Tests for gravity_engine.lichnerowicz_axisym (seal d62bfc0e0b; source instrument <LAB>/BHSLICE/p2a.py).

(a) reproduction: j = 2 and j = 10 certified with the module give the exact rational delta endpoints of
    <LAB>/BHSLICE/result_p2a.json byte for byte (same degrees); also the whole node records, and the K1/K2 cells
    rebuilt from the stored nodes equal the stored result file.
(b) toy J = 0: REFUSED/TRIVIAL with S == 0 and p(0) == 1/40.
(c) the three extra controls of REFEREE-P2a: exact Chebyshev reconstruction of Q at random rational points; residual
    signs from the spherical-coordinate Laplacian with l > 2 modes (mpmath, no use of the engine's operator);
    Bernstein reconstruction at random rational points.
(d) null: u+ - 1e-6 j^2 s (U-scale: U_N + (eta - 1e-6) s) gives a certified VIOLATION.
Run under memcap:  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 -B -m pytest <LAB>/POCKETLAB/gravity_engine/tests -q -p no:cacheprovider
"""
import sys
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # before any numpy import (deterministic fit)
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
sys.path.insert(0, '<LAB>/POCKETLAB')

import json
import math
import random
from fractions import Fraction

import pytest
import mpmath as mp
from flint import fmpq, fmpq_poly

from gravity_engine import lichnerowicz_axisym as la

RESULT = '<LAB>/BHSLICE/result_p2a.json'
BY = la.BOWEN_YORK


@pytest.fixture(scope='module')
def stored():
    with open(RESULT) as f:
        return json.load(f)


@pytest.fixture(scope='module')
def nodes():
    """certify j = 2 and j = 10 once with the module; keep the exact ansatz objects."""
    out = {}
    for j in (2, 10):
        rec, fn = la.certify_bowen_york(Fraction(j), keep=True, log=lambda m: None)
        out[j] = (rec, fn)
    return out


# ------------------------------------------------------------------ (a) reproduction
@pytest.mark.parametrize('j', [2, 10])
def test_reproduce_delta_endpoints(nodes, stored, j):
    rec, _ = nodes[j]
    ref = next(n for n in stored['nodes'] if n['j'] == '%d/1' % j)
    assert rec['status'] == 'CERTIFIED'
    assert rec['degrees'] == ref['degrees']
    assert rec['delta_lo'] == ref['delta_lo']
    assert rec['delta_hi'] == ref['delta_hi']
    # whole record (S, eta, Bernstein result, point checks, chi, T0c), timing excluded
    mine = {k: v for k, v in rec.items() if not k.startswith('_')}
    assert json.dumps(mine, sort_keys=True) == json.dumps(ref, sort_keys=True)


def test_reproduce_cells_and_cover(stored):
    """K1 cells and K2 cover rebuilt by the module from the stored node records equal the stored result exactly."""
    res = BY.build_result(stored['nodes'], stored['instrument_sha256'])
    assert res['all_nodes_certified'] and res['no_gaps'] and len(res['cells']) == 122
    assert json.dumps(res, sort_keys=True) == json.dumps(stored, sort_keys=True)


# ------------------------------------------------------------------ (b) toy J = 0
def test_toy_zero_refused():
    rec, fn0 = la.certify_bowen_york(Fraction(0), N=40, K=8, keep=True)
    assert rec['status'] == 'REFUSED/TRIVIAL'
    assert la.qparse(rec['S']) == 0
    assert fn0.p0 == fmpq(1, 40)


# ------------------------------------------------------------------ (c1) Chebyshev reconstruction of Q
def _cheb_int(n, x):
    """T_k(x) * q^k as python ints for x = p/q, k = 0..n (independent recurrence)."""
    p, q = int(x.p), int(x.q)
    a = [1, p]                     # a_k = T_k(x) q^k
    for k in range(1, n):
        a.append(2 * p * a[k] - q * q * a[k - 1])
    return a[:n + 1], q


@pytest.mark.parametrize('j,npts', [(2, 2), (10, 3)])
def test_chebyshev_reconstruction(nodes, j, npts):
    rec, fn = nodes[j]
    S, Q = fn.residual_S(want_Q=True)
    assert la.qstr(S) == rec['S']
    W = fn.W
    C, sc = la.chebmat(Q, W)
    rows, cols = C.nrows(), C.ncols()
    Ci = [[int(C[i, k]) for k in range(cols)] for i in range(rows)]
    # S == sum |cheb(tau Q)| / tau, recomputed from the Chebyshev matrix of tau Q
    Ct, sct = la.chebmat(Q * fn.tau, W)
    tot = sum(abs(int(Ct[i, k])) for i in range(Ct.nrows()) for k in range(Ct.ncols()))
    assert fmpq(tot) * sct / fn.tau == S
    mono = la.unpack(Q, W)
    rng = random.Random(1000 + j)
    for _ in range(npts):
        s0 = fmpq(rng.randint(0, 997), 997)
        c0 = fmpq(rng.randint(-991, 991), 991)
        # Chebyshev series, exact integer arithmetic
        As, qs = _cheb_int(rows - 1, 2 * s0 - 1)
        Ac, qc = _cheb_int(cols - 1, c0)
        tot_num = 0
        for k in range(cols):
            inner = sum(Ci[i][k] * As[i] * qs ** (rows - 1 - i) for i in range(rows) if Ci[i][k])
            tot_num += inner * Ac[k] * qc ** (cols - 1 - k)
        series = fmpq(tot_num) / (fmpq(qs) ** (rows - 1) * fmpq(qc) ** (cols - 1)) * sc
        # (i) the packed monomial form
        v_mono = sum((v * s0 ** i * c0 ** k for (i, k), v in mono.items()), fmpq(0))
        # (ii) straight from the ansatz terms, no packing: Q = D^7 Lh U + src, u = tau U
        u = sum((f(s0) * g(c0) for f, g in fn.Ut), fmpq(0))
        Lu = sum((f(s0) * g(c0) for f, g in fn.Lt), fmpq(0))
        D = 1 + (1 - s0) * u
        src = fmpq(9, 4) * (1 - c0 ** 2) * s0 ** 2 * (1 - s0) ** 2
        v_direct = D ** 7 * Lu / fn.tau + src
        assert series == v_mono == v_direct
        assert abs(series) <= S


# ------------------------------------------------------------------ (c2) residual signs from the spherical Laplacian
def _mpq(q):
    return mp.mpf(int(q.p)) / int(q.q)


@pytest.mark.parametrize('j', [2, 10])
def test_spherical_residual_signs(nodes, j):
    """E[u] = u_rr + 2u_r/r + (u_thth + cot u_th)/r^2 + (9/4) J^2 sin^2 r^-6 Phi^-7, Phi = 1 + 1/(2r) + u, built from the
    dyadic coefficients with T_n = cos(n acos X) and mpmath's Legendre; must be < 0 for u+ and > 0 for u-."""
    rec, fn = nodes[j]
    K = len(fn.qcs)
    assert K >= 2 and any(any(c != 0 for c in fn.qcs[k]) for k in range(1, K))   # modes l >= 4 present
    mp.mp.dps = 60
    pcm = [_mpq(c) for c in fn.pc]
    qcm = [[_mpq(c) for c in q] for q in fn.qcs]
    taum = mp.mpf(j) ** 2
    etam = _mpq(la.qparse(rec['eta']))
    Sm = _mpq(la.qparse(rec['S']))
    J = mp.mpf(j) / 4

    def cheb(cs, X):
        t = mp.acos(X)
        return mp.fsum(c * mp.cos(n * t) for n, c in enumerate(cs))

    def uN(r, th):                 # u-scale fitted polynomial in (r, theta)
        s = 1 / (2 * r + 1)
        X = 2 * s - 1
        c = mp.cos(th)
        val = s * cheb(pcm, X)
        for k, q in enumerate(qcm, 1):
            val += s ** 2 * (1 - s) * cheb(q, X) * mp.legendre(2 * k, c)
        return val

    rng = random.Random(77 + j)
    pts = []
    for _ in range(6):
        r = mp.mpf(rng.choice([1, 2, 3, 5, 7])) / rng.choice([1, 2, 4, 8, 16]) * rng.choice([1, 1, 3, 10])
        th = mp.pi * rng.randint(1, 15) / 32
        pts.append((r, th))
    teeth = 0
    for r, th in pts:
        ur = mp.diff(uN, (r, th), (1, 0))
        urr = mp.diff(uN, (r, th), (2, 0))
        ut = mp.diff(uN, (r, th), (0, 1))
        utt = mp.diff(uN, (r, th), (0, 2))
        lapN = urr + 2 * ur / r + (utt + mp.cot(th) * ut) / r ** 2
        u0 = uN(r, th)
        s = 1 / (2 * r + 1)
        lap_s = 8 / (2 * r + 1) ** 3 + 2 * (-2 / (2 * r + 1) ** 2) / r      # Laplacian of the radial s(r)
        src = mp.mpf(9) / 4 * J ** 2 * mp.sin(th) ** 2 * r ** -6

        def E(sgn, srcfac=1):
            u = u0 + sgn * taum * etam * s
            return lapN + sgn * taum * etam * lap_s + srcfac * src * (1 + 1 / (2 * r) + u) ** -7

        Ep, Em, E0 = E(+1), E(-1), E(0)
        assert Ep < 0 < Em, (float(r), float(th), float(Ep), float(Em))
        fac = 4 * s ** 4 / (1 - s) * taum                    # E = fac * Gh (REFEREE-P2a item 1)
        Q0 = E0 / fac * (1 + (1 - s) * u0) ** 7
        assert abs(Q0) <= Sm
        # teeth: a source wrong by 1 part in 1e12 must break the super-solution sign (at least at one point)
        teeth += bool(E(+1, 1 + mp.mpf(10) ** -12) >= 0)
    assert teeth >= 1


# ------------------------------------------------------------------ (c3) Bernstein reconstruction
def _bern_check(terms, a, b, rng, npts=5):
    B, n, m = la.bernstein_matrix(terms, a, b)
    for _ in range(npts):
        t = fmpq(rng.randint(0, 97), 97)
        w = fmpq(rng.randint(0, 89), 89)
        bt = [math.comb(n, i) * t ** i * (1 - t) ** (n - i) for i in range(n + 1)]
        bw = [math.comb(m, k) * w ** k * (1 - w) ** (m - k) for k in range(m + 1)]
        lhs = sum((B[i, k] * bt[i] * bw[k] for i in range(n + 1) for k in range(m + 1)), fmpq(0))
        s0, c0 = a + (b - a) * t, 2 * w - 1
        rhs = sum((f(s0) * g(c0) for f, g in terms), fmpq(0))
        assert lhs == rhs
    return B, n, m


def test_bernstein_reconstruction(nodes):
    rec, fn = nodes[10]
    eta = la.qparse(rec['eta'])
    T = fn.V_terms(eta)
    rng = random.Random(4242)
    B, n, m = _bern_check(T, fmpq(0), fmpq(1), rng)
    assert m == 2 * len(fn.qcs)
    bmin = min(B[i, k] for i in range(n + 1) for k in range(m + 1))
    assert la.qstr(bmin) == rec['positivity']['bern_min']     # the one-piece minimum the certificate used
    _bern_check(T, fmpq(1, 4), fmpq(3, 8), rng)
    R = [(fmpq_poly([fmpq(rng.randint(-9, 9), 7) for _ in range(6)]),
          fmpq_poly([fmpq(rng.randint(-9, 9), 5) for _ in range(5)]))]
    _bern_check(R, fmpq(1, 3), fmpq(1, 2), rng)


# ------------------------------------------------------------------ (d) null
@pytest.mark.parametrize('j', [2, 10])
def test_null_shifted_super_violation(nodes, j):
    rec, fn = nodes[j]
    eta = la.qparse(rec['eta'])
    r = la.certify_pair(fn, eta, e_sup=eta - fmpq(1, 10 ** 6), log=lambda m: None, tag='null j=%d' % j)
    assert r['status'] == 'VIOLATION'
