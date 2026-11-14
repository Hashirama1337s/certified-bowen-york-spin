"""BHSLICE P2a instrument: certified Bowen-York spinning puncture, mu = 1, J = j/4, j in [1, 10].

Spec: DESIGN-P2.md sections 1-4.1, 5 (nodes R1-R5, K1, K2) with every fix of REFEREE-P2DESIGN.md
(3a per-cell chi, 3b K1 in arb with outward rounding, 3c grid j_{k+1} = (51/50) j_k, 6a nonzero-j v1 identity,
6b T0b expected REFUSED/TRIVIAL, S > 0 assertion) and REFEREE-P1.md (N4 near-threshold null).

Problem (P1's s-coordinate, BUILD-P1 sec 1): s = a/(r+a), a = 1/2, c = cos(theta), tau = j^2, u = tau U,
   Gh[U] = Lh U + src * D^-7,   src = (9/4)(1-c^2) s^2 (1-s)^2,   D = 1 + (1-s) tau U,
   Lh U  = (1-s) U_0'' - 2 U_0' + sum_l [ (1-s) U_l'' - 2 U_l' - l(l+1) q_l ] P_l(c),
   U = s p(s) + sum_{l=2,4,..,2K} s^2 (1-s) q_l(s) P_l(c),  p, q_l Chebyshev polynomials in X = 2s - 1.
Certificate per node (all exact rationals, python-flint fmpq_poly / fmpz_mat / fmpq_mat):
 R1  Q := D_N^7 Lh U_N + src exactly (computed as tau Q = D_N^7 Lh u_N + tau src, so the big tau never enters D^7).
 R2  sup |Q| <= S := sum |Chebyshev coefficients of Q| on [0,1]_s x [-1,1]_c.
 R4  S > 0 asserted (else REFUSED/TRIVIAL); eta := S (1 + 2^-20)/2 rounded UP to a dyadic; U+- = U_N +- eta s.
 R3  V_- := U_-/s > 0 on the closed square by exact Bernstein coefficients on s-pieces (halve failing piece, depth 8).
     Then D_- >= 1, D_N >= 1, |Gh[U_N]| <= |Q| <= S, Gh[U+] <= S - 2 eta < 0 < -S + 2 eta <= Gh[U-] (Lh(s) = -2).
 R5  slope route: delta = 2 lim r u = tau p_true(0) in tau [p_N(0) - eta, p_N(0) + eta] (exact rationals).
     chi = j / (4 (1 + delta)^2) (decreasing in delta, exact rationals).
 Extra (belt and braces, not needed by the proof): exact sign check of Gh[U+-] at 27 rational points (a certified
 wrong sign there is a VIOLATION), stage-0 / RN / chi < 1 brackets (T0c).
K1 cells (DGC Lemma 3.1 + scaling: delta(j)/sqrt(j) nondecreasing), in arb with outward rounding:
   delta(j) in [sqrt(j/j_a) delta_-(j_a), sqrt(j/j_b) delta_+(j_b)] for j in [j_a, j_b];
   chi_max = j_b / (4 (1 + sqrt(j_b/j_a) delta_-(j_a))^2),  chi_min = j_a / (4 (1 + sqrt(j_a/j_b) delta_+(j_b))^2).

Usage (every run under memcap):
  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 p2a.py --probe [--j 10 --N 48 --K 10]
  ... p2a.py --controls controls_p2a.json      (controls; writes controls_p2a.log)
  ... p2a.py --main result_p2a.json            (refuses without GO-P2a; resumable; <= 520 s per invocation)
No pickle anywhere: checkpoints store exact rationals as "p/q" strings.
"""
import sys
sys.dont_write_bytecode = True
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # single-thread LAPACK => deterministic fits
    os.environ[_v] = '1'
import json, time, math, hashlib, argparse, subprocess, shutil
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat, fmpq, fmpz, fmpq_poly, fmpz_mat, fmpq_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-P2a')
TIME_BUDGET = 520.0
PREC = 160                # arb refinement precision (bits)
RBITS = 128               # fitted coefficients rounded to multiples of 2^-RBITS (exact dyadics)
COEF_CAP = 1.0e3          # guard: any |fitted coefficient| above this (or non-finite) => INVALID
BERN_DEPTH = 8            # max halvings of an s-piece in the Bernstein positivity check
EXTRA_NODES = (2, 4, 6, 8, 10)


def degrees(j):
    """(N, K): Chebyshev degree of p, q_l in X = 2s - 1, and Legendre modes l = 2..2K. Justified in BUILD-P2a.md."""
    jf = float(j)
    if jf <= 4.0:
        return 40, 8
    return 40, 10


# ================================================================== grid (exact rationals)
def grid_nodes():
    nodes = []
    j = Fraction(1)
    while True:
        nodes.append(j)
        if j >= 10:
            break
        j = j * Fraction(51, 50)
    for e in EXTRA_NODES:
        if Fraction(e) not in nodes:
            nodes.append(Fraction(e))
    return sorted(nodes)


def fq(x):
    """Fraction/int -> fmpq"""
    x = Fraction(x)
    return fmpq(int(x.numerator), int(x.denominator))


def qstr(q):
    return "%d/%d" % (int(q.p), int(q.q))


def qparse(s):
    a, b = s.split('/')
    return fmpq(int(a), int(b))


def qfloat(q):
    return float(Fraction(int(q.p), int(q.q)))


def qfloat_up(q):
    f = qfloat(q)
    return float(np.nextafter(np.nextafter(f, np.inf), np.inf))


# ================================================================== double / arb collocation (fit only; never a proof)
def _cheb_all(x, N, one):
    T = [one, x]
    d1 = [0 * one, one]
    d2 = [0 * one, 0 * one]
    for n in range(1, N):
        T.append(2 * x * T[n] - T[n - 1])
        d1.append(2 * T[n] + 2 * x * d1[n] - d1[n - 1])
        d2.append(4 * d1[n] + 2 * x * d2[n] - d2[n - 1])
    return T[:N + 1], d1[:N + 1], d2[:N + 1]


def _leg_all(c, L, one):
    P = [one, c]
    for l in range(1, L):
        P.append(((2 * l + 1) * c * P[l] - l * P[l - 1]) / (l + 1))
    return P[:L + 1]


def build_ops(N, K, mode):
    """Square collocation: s_i = (1 + cos(pi i/N))/2 (i = 0..N), c_m = cos(pi m/(2K)) (m = 0..K).
    Unknown vector: block 0 = p_n (n = 0..N), block k = q_{2k,n}. Returns rows of Lh and of U, and (1-s), src."""
    if mode == 'f':
        one, pi, cosf = 1.0, math.pi, math.cos
    else:
        one, pi, cosf = arb(1), arb.pi(), (lambda t: t.cos())
    ss = [(one + cosf(pi * i / N)) / 2 for i in range(N + 1)]
    ss[0] = one * 1
    ss[N] = 0 * one
    cs = [cosf(pi * m / (2 * K)) for m in range(K + 1)]
    cs[0] = one * 1
    cs[K] = 0 * one
    legs = [_leg_all(c, 2 * K, one) for c in cs]
    rows_L, rows_U, om, src = [], [], [], []
    for s in ss:
        X = 2 * s - one
        T, d1, d2 = _cheb_all(X, N, one)
        T1 = [2 * v for v in d1]
        T2 = [4 * v for v in d2]
        o = one - s
        f1 = [T[n] + s * T1[n] for n in range(N + 1)]
        f2 = [2 * T1[n] + s * T2[n] for n in range(N + 1)]
        Rf = [o * f2[n] - 2 * f1[n] for n in range(N + 1)]
        fv = [s * T[n] for n in range(N + 1)]
        h = s * s - s * s * s
        h1 = 2 * s - 3 * s * s
        h2 = 2 * one - 6 * s
        g1 = [h1 * T[n] + h * T1[n] for n in range(N + 1)]
        g2 = [h2 * T[n] + 2 * h1 * T1[n] + h * T2[n] for n in range(N + 1)]
        Rg = [o * g2[n] - 2 * g1[n] for n in range(N + 1)]
        gv = [h * T[n] for n in range(N + 1)]
        for mi, c in enumerate(cs):
            P = legs[mi]
            rl = list(Rf)
            ru = list(fv)
            for k in range(1, K + 1):
                l = 2 * k
                rl += [(Rg[n] - l * (l + 1) * T[n]) * P[l] for n in range(N + 1)]
                ru += [gv[n] * P[l] for n in range(N + 1)]
            rows_L.append(rl)
            rows_U.append(ru)
            om.append(o)
            src.append(9 * (one - c * c) * s * s * o * o / 4)
    return rows_L, rows_U, om, src


_FOPS, _AOPS = {}, {}


def float_ops(N, K):
    if (N, K) not in _FOPS:
        L, U, om, src = build_ops(N, K, 'f')
        _FOPS[(N, K)] = (np.array(L), np.array(U), np.array(om), np.array(src))
    return _FOPS[(N, K)]


def arb_ops(N, K):
    if (N, K, PREC) not in _AOPS:
        ctx.prec = PREC
        L, U, om, src = build_ops(N, K, 'a')
        _AOPS[(N, K, PREC)] = (arb_mat(L), arb_mat(U), om, src)
    return _AOPS[(N, K, PREC)]


def newton_float(N, K, tau, w0, itmax=40):
    L, U, om, src = float_ops(N, K)
    w = w0.copy()
    for it in range(itmax):
        D = 1 + om * (U @ w)
        F = L @ w + tau * src * D ** -7
        J = L - (7 * tau * src * om * D ** -8)[:, None] * U
        dw = np.linalg.solve(J, F)
        w = w - dw
        if np.max(np.abs(dw)) < 1e-15 * max(1.0, np.max(np.abs(w))):
            break
    D = 1 + om * (U @ w)
    return w, float(np.max(np.abs(L @ w + tau * src * D ** -7)))


def fit_u(j, N, K, log=None):
    """u = tau U (u-scale) coefficients: double Newton with a fixed continuation in j (1, 2, 4, ..., j),
    then iterative refinement with residuals in arb at PREC bits. Deterministic given (j, N, K)."""
    jf = float(j)
    n = (N + 1) * (K + 1)
    w = np.zeros(n)
    sched = []
    x = min(1.0, jf)
    while x < jf:
        sched.append(x)
        x *= 2
    sched.append(jf)
    for jj in sched:
        w, resf = newton_float(N, K, jj * jj, w)
    if not np.all(np.isfinite(w)):
        return None, dict(float_residual=float('nan'), arb_hist=[])
    ctx.prec = PREC
    La, Ua, om, src = arb_ops(N, K)
    Lf, Uf, omf, srcf = float_ops(N, K)
    tau_q = fq(j) ** 2
    tau_a = arb(tau_q)
    tf = jf * jf
    wa = arb_mat([[arb(float(v))] for v in w])
    hist = []
    for it in range(14):
        Lw = La * wa
        Uw = Ua * wa
        F = [Lw[i, 0] + tau_a * src[i] / (1 + om[i] * Uw[i, 0]) ** 7 for i in range(n)]
        Ff = np.array([float(f.mid()) for f in F])
        res = float(np.max(np.abs(Ff)))
        hist.append(res)
        if res < 2.0 ** (-(PREC - 40)):
            break
        wd = np.array([float(wa[i, 0].mid()) for i in range(n)])
        D = 1 + omf * (Uf @ wd)
        J = Lf - (7 * tf * srcf * omf * D ** -8)[:, None] * Uf
        dw = np.linalg.solve(J, Ff)
        wa = arb_mat([[(wa[i, 0] - arb(float(dw[i]))).mid()] for i in range(n)])
    return [wa[i, 0].mid() for i in range(n)], dict(float_residual=resf, arb_hist=hist)


def to_dyadic(a, bits=RBITS):
    """exact arb midpoint -> nearest multiple of 2^-bits (exact fmpq)."""
    m, e = a.man_exp()
    m, e = int(m), int(e)
    if e >= -bits:
        return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))
    sh = -bits - e
    mm = (m + (1 << (sh - 1))) >> sh
    return fmpq(mm, 2 ** bits)


def arb_to_fmpq(a):
    m, e = a.man_exp()
    m, e = int(m), int(e)
    return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))


def split_coefs(vec, N, K):
    pc = list(vec[:N + 1])
    qcs = [list(vec[(N + 1) * k:(N + 1) * (k + 1)]) for k in range(1, K + 1)]
    return pc, qcs


# ================================================================== exact polynomial engine (R1, R2)
SP = fmpq_poly([0, 1])
OMS = fmpq_poly([1, -1])
ONE = fmpq_poly([1])
HP = fmpq_poly([0, 0, 1, -1])            # s^2 (1 - s)
SRC_S = fmpq(9, 4) * SP ** 2 * OMS ** 2  # (9/4) s^2 (1-s)^2
SRC_C = fmpq_poly([1, 0, -1])           # 1 - c^2
_TY, _LEG = [], []


def cheb_s(N):
    """T_n(2s - 1) as exact fmpq_poly in s, n = 0..N."""
    if not _TY:
        _TY.extend([ONE, fmpq_poly([-1, 2])])
    while len(_TY) <= N:
        n = len(_TY) - 1
        _TY.append(2 * fmpq_poly([-1, 2]) * _TY[n] - _TY[n - 1])
    return _TY[:N + 1]


def leg_c(L):
    if not _LEG:
        _LEG.extend([ONE, fmpq_poly([0, 1])])
    while len(_LEG) <= L:
        l = len(_LEG) - 1
        _LEG.append((fmpq(2 * l + 1) * fmpq_poly([0, 1]) * _LEG[l] - fmpq(l) * _LEG[l - 1]) / fmpq(l + 1))
    return _LEG[:L + 1]


def Rop(f):
    return OMS * f.derivative().derivative() - 2 * f.derivative()


def build_terms(pc, qcs):
    """U (or u) and Lh U as lists of (poly in s, poly in c); p(s) and q_l(s); p(0)."""
    N = len(pc) - 1
    Ty = cheb_s(N)
    Pl = leg_c(2 * len(qcs))
    pp = fmpq_poly([0])
    for n in range(N + 1):
        if pc[n] != 0:
            pp += pc[n] * Ty[n]
    f0 = SP * pp
    Ut = [(f0, ONE)]
    Lt = [(Rop(f0), ONE)]
    qs = []
    for k, qc in enumerate(qcs, 1):
        l = 2 * k
        q = fmpq_poly([0])
        for n in range(len(qc)):
            if qc[n] != 0:
                q += qc[n] * Ty[n]
        qs.append(q)
        if q == 0:
            continue
        g = HP * q
        Ut.append((g, Pl[l]))
        Lt.append((Rop(g) - fmpq(l * (l + 1)) * q, Pl[l]))
    p0 = fmpq(0)
    for n in range(N + 1):
        p0 += pc[n] if n % 2 == 0 else -pc[n]
    return Ut, Lt, pp, qs, p0


def pack(terms, W):
    """Kronecker packing of sum f(s) g(c) into one fmpq_poly (index = i W + k for s^i c^k)."""
    terms = [(f, g) for f, g in terms if f != 0 and g != 0]
    if not terms:
        return fmpq_poly([0])
    ds = max(f.degree() for f, g in terms)
    co = [fmpq(0)] * ((ds + 1) * W)
    for f, g in terms:
        gc = g.coeffs()
        assert len(gc) <= W
        for i, fi in enumerate(f.coeffs()):
            if fi == 0:
                continue
            base = i * W
            for k, gk in enumerate(gc):
                if gk != 0:
                    co[base + k] += fi * gk
    return fmpq_poly(co)


def unpack(P, W):
    out = {}
    for idx, v in enumerate(P.coeffs()):
        if v != 0:
            out[(idx // W, idx % W)] = v
    return out


_AX = {}


def Ax_int(d):
    """integer A with x^m = 2^-d sum_k A[k,m] T_k(x), k, m = 0..d."""
    if d not in _AX:
        A = [[0] * (d + 1) for _ in range(d + 1)]
        for m in range(d + 1):
            sc = 1 << (d - m)
            for l in range(m + 1):
                A[abs(m - 2 * l)][m] += sc * math.comb(m, l)
        _AX[d] = A
    return _AX[d]


_AY = {}


def Ay_mat(dy):
    """integer matrix: s^i = 2^-(2 dy) sum_k Ay[k,i] T_k(2s-1)."""
    if dy not in _AY:
        Ax = Ax_int(dy)
        B = fmpz_mat(dy + 1, dy + 1, [((1 << (dy - i)) * math.comb(i, m) if m <= i else 0)
                                       for m in range(dy + 1) for i in range(dy + 1)])
        _AY[dy] = fmpz_mat(dy + 1, dy + 1, [Ax[k][m] for k in range(dy + 1) for m in range(dy + 1)]) * B
    return _AY[dy]


def chebmat(P, W):
    """Exact Chebyshev coefficients of the packed P(s, c) on [0,1] x [-1,1]:
    coefficient of T_k(2s-1) T_l(c) = C[k,l] * scale. Returns (C fmpz_mat, scale fmpq) or (None, 0) for P = 0."""
    if P == 0:
        return None, fmpq(0)
    num = P.numer().coeffs()
    den = int(P.denom())
    Ln = len(num)
    dy = (Ln - 1) // W
    rows = dy + 1
    flat = [0] * (rows * W)
    for idx, v in enumerate(num):
        flat[idx] = int(v)
    dc = 0
    for i in range(rows):
        for jj in range(W - 1, -1, -1):
            if flat[i * W + jj] != 0:
                dc = max(dc, jj)
                break
    M = fmpz_mat(rows, dc + 1, [flat[i * W + jj] for i in range(rows) for jj in range(dc + 1)])
    Ac = Ax_int(dc)
    AcT = fmpz_mat(dc + 1, dc + 1, [Ac[l][jj] for jj in range(dc + 1) for l in range(dc + 1)])
    C = Ay_mat(dy) * M * AcT
    return C, fmpq(1, den * (1 << (2 * dy + dc)))


def abs_sum_combined(CA, scA, CB, scB, lam):
    """exact sum over (k,l) of |CA[k,l] scA + lam CB[k,l] scB| (CB small)."""
    tot = 0
    if CA is not None:
        for i in range(CA.nrows()):
            for jj in range(CA.ncols()):
                v = int(CA[i, jj])
                tot += -v if v < 0 else v
    S = fmpq(tot) * scA if CA is not None else fmpq(0)
    if CB is not None and lam != 0:
        for i in range(CB.nrows()):
            for jj in range(CB.ncols()):
                b = int(CB[i, jj])
                if b == 0:
                    continue
                a = int(CA[i, jj]) if (CA is not None and i < CA.nrows() and jj < CA.ncols()) else 0
                val = fmpq(a) * scA + lam * fmpq(b) * scB
                S += (val if val >= 0 else -val) - (fmpq(abs(a)) * scA if CA is not None else fmpq(0))
    return S


class ExactFn:
    """U_N given by Chebyshev/Legendre coefficients. For tau > 0 the coefficients are in u-scale (u = tau U, exact
    dyadics from the fit); for tau == 0 they are the U-scale coefficients themselves (stage-0 toy model)."""

    def __init__(self, pc, qcs, tau):
        self.pc, self.qcs, self.tau = pc, qcs, tau
        self.K = len(qcs)
        self.W = 16 * self.K + 3
        self.Ut, self.Lt, self.pp, self.qs, self.p0 = build_terms(pc, qcs)

    def p0_U(self):
        return self.p0 / self.tau if self.tau != 0 else self.p0

    def residual_S(self, want_Q=False):
        """R1 + R2: exact S = sum |Chebyshev coefficients of Q|, Q = D^7 Lh U + src (U-scale)."""
        W = self.W
        srcP = pack([(SRC_S, SRC_C)], W)
        CB, scB = chebmat(srcP, W)
        if self.tau == 0:
            Q = pack(self.Lt, W) + srcP
            C, sc = chebmat(Q, W)
            S = abs_sum_combined(C, sc, None, 0, 0)
            return (S, Q) if want_Q else (S, None)
        DP = pack([(ONE, ONE)] + [(OMS * f, g) for f, g in self.Ut], W)
        A = DP ** 7 * pack(self.Lt, W)                # = tau (Q - src)
        CA, scA = chebmat(A, W)
        S_tau = abs_sum_combined(CA, scA, CB, scB, self.tau)   # = sum |cheb(tau Q)|
        S = S_tau / self.tau
        if want_Q:
            return S, A / self.tau + srcP
        return S, None

    # exact point evaluation of tau * Gh[U_N + sgn e s] (same sign as Gh)
    def tauGh_at(self, s0, c0, e_signed):
        u = fmpq(0)
        Lu = fmpq(0)
        for f, g in self.Ut:
            u += f(s0) * g(c0)
        for f, g in self.Lt:
            Lu += f(s0) * g(c0)
        if self.tau == 0:
            return Lu - 2 * e_signed + SRC_S(s0) * SRC_C(c0)
        u += self.tau * e_signed * s0
        Lu += -2 * self.tau * e_signed
        D = 1 + (1 - s0) * u
        if D <= 0:
            return None
        return Lu + self.tau * SRC_S(s0) * SRC_C(c0) / D ** 7

    def V_terms(self, shift_U):
        """terms of tau * (U/s - shift_U) for U = U_N (u-scale polys): p(s) - tau shift + s(1-s) sum q_l P_l."""
        t = self.tau if self.tau != 0 else fmpq(1)
        Pl = leg_c(2 * self.K)
        terms = [(self.pp - t * shift_U, ONE)]
        for k, q in enumerate(self.qs, 1):
            if q != 0:
                terms.append((SP * OMS * q, Pl[2 * k]))
        return terms


PT_S = [fmpq(i, 8) for i in range(9)]
PT_C = [fmpq(0), fmpq(1, 2), fmpq(1)]


def point_checks(fn, e_signed, want_sign):
    """exact sign of Gh[U_N + e_signed s] at 27 rational points; want_sign -1 (super) or +1 (sub)."""
    for s0 in PT_S:
        for c0 in PT_C:
            v = fn.tauGh_at(s0, c0, e_signed)
            if v is None:
                return dict(status='VIOLATION', at=[qstr(s0), qstr(c0)], why='D <= 0 at a point')
            if (want_sign < 0 and v >= 0) or (want_sign > 0 and v <= 0):
                return dict(status='VIOLATION', at=[qstr(s0), qstr(c0)], value_float=qfloat(v))
    return dict(status='OK', points=len(PT_S) * len(PT_C))


# ================================================================== R3: exact Bernstein positivity
_BT = {}


def bern_T(n):
    """fmpq_mat T with T[j,i] = C(j,i)/C(n,i) (monomial -> Bernstein on [0,1])."""
    if n not in _BT:
        _BT[n] = fmpq_mat(n + 1, n + 1, [fmpq(math.comb(j, i), math.comb(n, i)) if i <= j else fmpq(0)
                                         for j in range(n + 1) for i in range(n + 1)])
    return _BT[n]


_GW = {}


def c_to_w(m):
    """G[k, mm]: c^k = sum_mm G[k,mm] w^mm with c = 2w - 1."""
    if m not in _GW:
        _GW[m] = fmpq_mat(m + 1, m + 1, [fmpq(math.comb(k, mm) * 2 ** mm * (-1) ** (k - mm)) if mm <= k else fmpq(0)
                                         for k in range(m + 1) for mm in range(m + 1)])
    return _GW[m]


def bernstein_positive(terms, max_depth=BERN_DEPTH):
    """certify V(s,c) = sum f(s) g(c) > 0 on [0,1] x [-1,1]. Corner Bernstein coefficients are exact point values:
    a corner <= 0 is a certified VIOLATION. Otherwise halve the s-piece until all coefficients > 0, or FAIL at depth."""
    m = max(g.degree() for f, g in terms)
    m = max(m, 0)
    Vc = [fmpq_poly([0]) for _ in range(m + 1)]
    for f, g in terms:
        for k, gk in enumerate(g.coeffs()):
            if gk != 0:
                Vc[k] += gk * f
    n = max(max(p.degree() for p in Vc), 0)
    G = c_to_w(m)
    Tc = bern_T(m).transpose()
    Ts = bern_T(n)
    queue = [(fmpq(0), fmpq(1), 0)]
    pieces, gmin = 0, None
    while queue:
        a, b, d = queue.pop(0)
        sub = fmpq_poly([a, b - a])
        rows = []
        cols = [Vc[k](sub).coeffs() if Vc[k] != 0 else [] for k in range(m + 1)]
        ent = []
        for i in range(n + 1):
            for k in range(m + 1):
                ent.append(cols[k][i] if i < len(cols[k]) else fmpq(0))
        Mt = fmpq_mat(n + 1, m + 1, ent)
        B = Ts * (Mt * G) * Tc
        corners = [B[0, 0], B[n, 0], B[0, m], B[n, m]]
        bmin = min(B[i, k] for i in range(n + 1) for k in range(m + 1))
        if min(corners) <= 0:
            return dict(status='VIOLATION', piece=[qstr(a), qstr(b)], corner_min=qfloat(min(corners)), pieces=pieces)
        if bmin > 0:
            pieces += 1
            gmin = bmin if gmin is None or bmin < gmin else gmin
            continue
        if d >= max_depth:
            return dict(status='FAIL', piece=[qstr(a), qstr(b)], bern_min=qfloat(bmin), pieces=pieces)
        mid = (a + b) / 2
        queue += [(a, mid, d + 1), (mid, b, d + 1)]
    return dict(status='PASS', pieces=pieces, bern_min=qstr(gmin), bern_min_float=qfloat(gmin))


def V_point_min(fn, shift_U, ns=64, nc=8):
    """exact minimum of tau V over a rational grid (used only to build the N-pos control)."""
    terms = fn.V_terms(shift_U)
    best = None
    for i in range(ns + 1):
        s0 = fmpq(i, ns)
        for k in range(nc + 1):
            c0 = fmpq(k, nc)
            v = fmpq(0)
            for f, g in terms:
                v += f(s0) * g(c0)
            if best is None or v < best:
                best = v
    return best


# ================================================================== decimal output (outward)
def dec_floor(q, nd=25):
    if q == 0:
        return "0"
    e = math.floor(math.log10(abs(qfloat(q)))) if abs(qfloat(q)) > 0 else 0
    scale = nd - 1 - e
    v = q * fmpq(10) ** scale if scale >= 0 else q / fmpq(10) ** (-scale)
    iv = int(v.p) // int(v.q)
    sgn = "-" if iv < 0 else ""
    s = str(abs(iv))
    if scale <= 0:
        return sgn + s + "0" * (-scale)
    if len(s) <= scale:
        s = "0" * (scale - len(s) + 1) + s
    return sgn + s[:-scale] + "." + s[-scale:]


def dec_ceil(q, nd=25):
    s = dec_floor(-q, nd)
    return s[1:] if s.startswith("-") else ("-" + s if s != "0" else s)


def round_up_dyadic(q, bits=64):
    if q <= 0:
        return q
    e = int(q.p).bit_length() - int(q.q).bit_length()
    sh = bits - e
    v = q * fmpq(2) ** sh if sh >= 0 else q / fmpq(2) ** (-sh)
    c = -((-int(v.p)) // int(v.q))
    r = fmpq(c) / fmpq(2) ** sh if sh >= 0 else fmpq(c) * fmpq(2) ** (-sh)
    assert r >= q
    return r


# ================================================================== stage-0 toy model (T0b, T0c, T6a)
def v1_cheb(N, K):
    """v1 = s p + s^2 (1-s) q2 P2 (GNPP A1, B1 in P1's s): exact Chebyshev coefficients in X = 2s - 1."""
    s = SP
    om = OMS
    p = (om ** 4 + 5 * om ** 3 * s + 10 * om ** 2 * s ** 2 + 5 * om * s ** 3 + s ** 4) / fmpq(40)
    q2 = -(om * s) / fmpq(20)

    def to_cheb(poly):
        # substitute s = (X+1)/2, then X^m = 2^-m sum_l C(m,l) T_|m-2l|
        PX = poly(fmpq_poly([fmpq(1, 2), fmpq(1, 2)]))
        co = PX.coeffs()
        out = [fmpq(0)] * (len(co) + 1)
        for m_, a in enumerate(co):
            for l in range(m_ + 1):
                out[abs(m_ - 2 * l)] += a * fmpq(math.comb(m_, l), 2 ** m_)
        return out

    pc = to_cheb(p) + [fmpq(0)] * (N + 1)
    q2c = to_cheb(q2) + [fmpq(0)] * (N + 1)
    pc = pc[:N + 1]
    qcs = [q2c[:N + 1]] + [[fmpq(0)] * (N + 1) for _ in range(K - 1)]
    return pc, qcs


def v1_box(sb, cb, one):
    """arb enclosure of v1 on a box (s ball, c ball)."""
    s, c = sb, cb
    om = one - s
    p = (om ** 4 + 5 * om ** 3 * s + 10 * om ** 2 * s ** 2 + 5 * om * s ** 3 + s ** 4) / 40
    q2 = -(om * s) / 20
    return s * p + s * s * om * q2 * (3 * c * c - 1) / 2


def stage0_lower(j, nt=160, nc=16):
    """certified lower bound of the stage-0 bracket I(j^2 v1) <= delta (lower Riemann sum in arb; v1 = super-solution,
    integrand decreasing in U).  delta >= (9 tau/4) (1/2) int_0^1 int_-1^1 t^2 (1-t)^3 (1-c^2) (1+(1-t) tau v1)^-7."""
    ctx.prec = 128
    tau = arb(fq(j)) ** 2
    tot = arb(0)
    for i in range(nt):
        ta, tb = arb(i) / nt, arb(i + 1) / nt
        tball = (ta + tb) / 2 + arb(0, 1) * (tb - ta) / 2
        for k in range(nc):
            ca, cb = arb(k) / nc, arb(k + 1) / nc
            cball = (ca + cb) / 2 + arb(0, 1) * (cb - ca) / 2
            Dmax = (1 + (1 - tball) * tau * v1_box(tball, cball, arb(1))).upper()
            low = ta * ta * (1 - tb) ** 3 * (1 - cb * cb) / Dmax ** 7
            tot += low.lower() * (tb - ta) * (cb - ca)
    # c in [-1, 1] = 2 x [0, 1]; prefactor (9 tau / 4) (1/2)
    val = (arb(9) / 4 * tau * tot).lower()
    return arb_to_fmpq(val)


# ================================================================== one node certificate
def certify_node(j, N=None, K=None, log=print, keep=False, coef_override=None):
    """Point certificate at the exact rational j (R1-R5 + extras). Returns a JSON-able record (deterministic)."""
    j = Fraction(j)
    if N is None:
        N, K = degrees(j)
    tau = fq(j) ** 2
    t0 = time.time()
    if j == 0:
        # T0a/T0b: u = 0 is the exact solution. The j -> 0 limit problem in U-scale is linear with exact solution v1;
        # it is fed to the exact engine, which must give S == 0 and therefore REFUSE (never a PASS of the sign checks).
        pc0, qcs0 = v1_cheb(N, K)
        fn0 = ExactFn(pc0, qcs0, fmpq(0))
        S0, _ = fn0.residual_S()
        rec = dict(j='0/1', j_dec='0', degrees=dict(N=N, K=K), S=qstr(S0), S_float=qfloat(S0), p0_U=qstr(fn0.p0),
                   status='REFUSED/TRIVIAL' if not S0 > 0 else 'ERROR: S > 0 for the exact stage-0 solution',
                   why='j = 0: u = 0 exact (delta = 0); U = v1 gives S == 0, a non-strict certificate')
        return rec, (fn0 if keep else None)
    if coef_override is None:
        wa, diag = fit_u(j, N, K)
        if wa is None:
            return dict(j=qstr(fq(j)), status='INVALID', why='non-finite double fit'), None
        wf = [float(a) for a in wa]
    else:
        wa, diag = None, dict(float_residual=None, arb_hist=[])
        wf = list(coef_override)
    # guard (F7): non-finite or absurd coefficients => INVALID, never PASS
    if not all(math.isfinite(v) for v in wf) or max(abs(v) for v in wf) > COEF_CAP:
        return dict(j=qstr(fq(j)), status='INVALID', why='non-finite or absurd coefficient (> %g)' % COEF_CAP), None
    if coef_override is None:
        uq = [to_dyadic(a) for a in wa]
    else:
        uq = [fmpq(Fraction(v).numerator, Fraction(v).denominator) for v in wf]
    t1 = time.time()
    pc, qcs = split_coefs(uq, N, K)
    fn = ExactFn(pc, qcs, tau)
    S, _ = fn.residual_S()
    t2 = time.time()
    rec = dict(j=qstr(fq(j)), j_dec=dec_floor(fq(j), 20), degrees=dict(N=N, K=K, prec_bits=PREC, round_bits=RBITS),
               S=qstr(S), S_float=qfloat_up(S))
    if not S > 0:
        rec.update(status='REFUSED/TRIVIAL', why='S <= 0 (non-strict certificate)')
        return rec, (fn if keep else None)
    eta = round_up_dyadic(S * (1 + fmpq(1, 2 ** 20)) / 2)
    assert eta >= S * (1 + fmpq(1, 2 ** 20)) / 2
    rec['eta'] = qstr(eta)
    rec['eta_float'] = qfloat_up(eta)
    # R3 positivity of V_- = U_-/s (exact Bernstein)
    pos = bernstein_positive(fn.V_terms(eta))
    t3 = time.time()
    rec['positivity'] = pos
    # extras: exact sign checks at 27 rational points
    pc_sup = point_checks(fn, eta, -1)
    pc_sub = point_checks(fn, -eta, +1)
    rec['point_checks'] = dict(sup=pc_sup, sub=pc_sub)
    margin_ok = (2 * eta > S) and S > 0
    ok = pos['status'] == 'PASS' and margin_ok and pc_sup['status'] == 'OK' and pc_sub['status'] == 'OK'
    if not ok:
        rec['status'] = ('VIOLATION' if 'VIOLATION' in (pos['status'], pc_sup['status'], pc_sub['status']) else 'FAIL')
        log(f"  node j={float(j):.6f}: {rec['status']} pos={pos['status']} sup={pc_sup['status']} sub={pc_sub['status']}")
        return rec, (fn if keep else None)
    # R5 slope route (exact rationals): delta = tau p_U(0) +- tau eta = p_u(0) +- tau eta
    p0u = fn.p0 if tau != 0 else fn.p0 * tau
    dlo = p0u - tau * eta
    dhi = p0u + tau * eta
    jq = fq(j)
    chi_lo = jq / (4 * (1 + dhi) ** 2)
    chi_hi = jq / (4 * (1 + dlo) ** 2)
    rel = (dhi - dlo) / dlo
    # T0c brackets: stage-0 (I(tau v1) <= delta <= tau/40), RN ((1+delta)^2 <= 1 + 3J), chi < 1
    s0lo = stage0_lower(j)
    t0c = dict(stage0_lower=dec_floor(s0lo, 20), stage0_upper=dec_ceil(tau / 40, 20),
               inside_stage0=bool(s0lo <= dlo and dhi <= tau / 40),
               rn_ok=bool((1 + dhi) ** 2 <= 1 + 3 * jq / 4), chi_lt_1=bool(chi_hi < 1))
    t0c['ok'] = bool(t0c['inside_stage0'] and t0c['rn_ok'] and t0c['chi_lt_1'])
    rec.update(status='CERTIFIED', delta_lo=qstr(dlo), delta_hi=qstr(dhi), delta_lo_dec=dec_floor(dlo), delta_hi_dec=dec_ceil(dhi),
               chi_lo=qstr(chi_lo), chi_hi=qstr(chi_hi), chi_lo_dec=dec_floor(chi_lo), chi_hi_dec=dec_ceil(chi_hi),
               rel_width=qfloat_up(rel), chi_rel_width=qfloat_up((chi_hi - chi_lo) / chi_lo), t0c=t0c)
    t4 = time.time()
    log(f"  node j={float(j):.6f} N={N} K={K}: CERTIFIED delta in [{rec['delta_lo_dec'][:20]}, {rec['delta_hi_dec'][:20]}] "
        f"rel {rec['rel_width']:.2e} eta {rec['eta_float']:.2e} chi in [{rec['chi_lo_dec'][:16]}, {rec['chi_hi_dec'][:16]}] "
        f"Bernstein pieces {pos['pieces']} min {pos['bern_min_float']:.3e}; T0c {t0c['ok']} "
        f"(fit {t1-t0:.1f}s, S {t2-t1:.1f}s, Bern {t3-t2:.1f}s, rest {t4-t3:.1f}s; float res {diag['float_residual']})")
    rec['_timing'] = dict(fit=t1 - t0, S=t2 - t1, bern=t3 - t2, rest=t4 - t3)
    return rec, (fn if keep else None)


def certify_pair(fn, eta, base_sub=None, e_sup=None, e_sub=None, log=print, tag=''):
    """Generic pair check used by the nulls: super = base_sup + e_sup s, sub = base_sub - e_sub s (U-scale e's).
    Order: exact point signs (VIOLATION) -> R3 Bernstein on V_- -> bound path (2 e > S(base), S > 0)."""
    e_sup = eta if e_sup is None else e_sup
    e_sub = eta if e_sub is None else e_sub
    fsub = fn if base_sub is None else base_sub
    S_sup, _ = fn.residual_S()
    S_sub = S_sup if fsub is fn else fsub.residual_S()[0]
    path = []
    ps = point_checks(fn, e_sup, -1)
    pb = point_checks(fsub, -e_sub, +1)
    path.append(f"point checks: sup {ps['status']}, sub {pb['status']}")
    if ps['status'] != 'OK' or pb['status'] != 'OK':
        log(f"  {tag}: VIOLATION at a rational point ({'; '.join(path)}) sup={ps} sub={pb}")
        return dict(status='VIOLATION', path=path, sup=ps, sub=pb)
    pos = bernstein_positive(fsub.V_terms(e_sub))
    path.append(f"Bernstein V_-: {pos['status']}")
    if pos['status'] != 'PASS':
        log(f"  {tag}: {pos['status']} in positivity ({'; '.join(path)})")
        return dict(status=pos['status'], path=path, pos=pos)
    if fsub is not fn:   # super base must be >= 0 too (D >= 1)
        pos2 = bernstein_positive(fn.V_terms(fmpq(0)))
        path.append(f"Bernstein base_sup/s: {pos2['status']}")
        if pos2['status'] != 'PASS':
            return dict(status=pos2['status'], path=path, pos=pos2)
    bsup = bool(S_sup > 0 and 2 * e_sup > S_sup)
    bsub = bool(S_sub > 0 and 2 * e_sub > S_sub)
    path.append(f"bound path: sup 2e = {qfloat(2*e_sup):.3e} vs S = {qfloat(S_sup):.3e} -> {'ok' if bsup else 'FAILS'}; "
                f"sub 2e = {qfloat(2*e_sub):.3e} vs S = {qfloat(S_sub):.3e} -> {'ok' if bsub else 'FAILS'}")
    st = 'PASS' if (bsup and bsub) else 'FAIL'
    log(f"  {tag}: {st} ({'; '.join(path)})")
    return dict(status=st, path=path)


# ================================================================== K1 cells and K2 coverage (arb, outward)
def cell_record(na, nb):
    ctx.prec = 256
    ja, jb = qparse(na['j']), qparse(nb['j'])
    dla, dhb = qparse(na['delta_lo']), qparse(nb['delta_hi'])
    A_ja, A_jb, A_dla, A_dhb = arb(ja), arb(jb), arb(dla), arb(dhb)
    R_lo = (A_dla / A_ja.sqrt()).lower()
    R_hi = (A_dhb / A_jb.sqrt()).upper()
    relw = ((A_ja / A_jb).sqrt() * A_dhb / A_dla - 1).upper()
    chi_max = (A_jb / (4 * (1 + (A_jb / A_ja).sqrt() * A_dla) ** 2)).upper()
    chi_min = (A_ja / (4 * (1 + (A_ja / A_jb).sqrt() * A_dhb) ** 2)).lower()
    nonempty = bool(R_hi >= R_lo)
    relu = (dhb - dla) / dla
    return dict(j_lo=na['j'], j_hi=nb['j'], j_lo_dec=na['j_dec'], j_hi_dec=nb['j_dec'],
                delta_lo=qstr(dla), delta_hi=qstr(dhb), delta_lo_dec=dec_floor(dla), delta_hi_dec=dec_ceil(dhb),
                R_lo=qstr(arb_to_fmpq(R_lo)), R_hi=qstr(arb_to_fmpq(R_hi)),
                chi_max=qstr(arb_to_fmpq(chi_max)), chi_min=qstr(arb_to_fmpq(chi_min)),
                chi_max_dec=dec_ceil(arb_to_fmpq(chi_max), 16), chi_min_dec=dec_floor(arb_to_fmpq(chi_min), 16),
                rel_width=qfloat_up(arb_to_fmpq(relw)), rel_width_uniform=qfloat_up(relu), nonempty=nonempty)


def build_result(node_recs, sha):
    nodes = sorted(node_recs, key=lambda r: Fraction(*map(int, r['j'].split('/'))))
    clean = [{k: v for k, v in r.items() if not k.startswith('_')} for r in nodes]
    all_cert = all(r['status'] == 'CERTIFIED' for r in clean)
    cells = []
    if all_cert:
        for a, b in zip(clean[:-1], clean[1:]):
            cells.append(cell_record(a, b))
    js = [Fraction(*map(int, r['j'].split('/'))) for r in clean]
    no_gaps = bool(all_cert and js and js[0] == 1 and js[-1] >= 10 and all(x < y for x, y in zip(js[:-1], js[1:]))
                   and all(cells[i]['j_hi'] == cells[i + 1]['j_lo'] for i in range(len(cells) - 1))
                   and len(cells) == len(clean) - 1 and all(c['nonempty'] for c in cells))
    return dict(instrument_sha256=sha, object='Bowen-York spinning puncture, mu = 1, J = j/4, delta = M_ADM - 1, chi = J/M_ADM^2',
                grid='j_{k+1} = (51/50) j_k from 1 to the first node >= 10, plus j = 2, 4, 6, 8, 10',
                nodes=clean, cells=cells, all_nodes_certified=all_cert, no_gaps=no_gaps,
                all_t0c_ok=all(r.get('t0c', {}).get('ok', False) for r in clean),
                max_cell_rel_width=max((c['rel_width'] for c in cells), default=None),
                max_cell_rel_width_uniform=max((c['rel_width_uniform'] for c in cells), default=None))


# ================================================================== resumable grid runner (main and RESUME control)
def sha_self():
    return hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()


def atomic_json(path, obj, sort=True):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=sort)
    os.replace(tmp, path)


def run_grid(out, nodes, budget, log):
    """Exit codes: 0 = finished or budget reached (rerun to resume), 3 = checkpoint from another instrument sha."""
    T0 = time.time()
    sha = sha_self()
    ck = out + '.ckpt'
    spec = [qstr(fq(j)) for j in nodes]
    if os.path.exists(ck):
        st = json.load(open(ck))
        if st.get('instrument_sha256') != sha:
            log(f"REFUSED: checkpoint {ck} was written by instrument sha {st.get('instrument_sha256')}, this is {sha}.")
            return 3
        if st.get('nodes_spec') != spec:
            log("REFUSED: checkpoint node list differs from this run's node list.")
            return 3
    else:
        st = dict(instrument_sha256=sha, nodes_spec=spec, done={}, timing={}, invocations=0)
    st['invocations'] += 1
    atomic_json(ck, st)
    worst = 10.0
    did = 0
    for j in nodes:
        key = qstr(fq(j))
        if key in st['done']:
            continue
        if did >= 1 and time.time() - T0 + 1.5 * worst > budget - 20:
            log(f"time budget reached: {len(st['done'])}/{len(nodes)} nodes done; rerun the same command to resume")
            return 0
        t = time.time()
        rec, _ = certify_node(j, log=log)
        dt = time.time() - t
        worst = max(worst, dt)
        st['timing'][key] = rec.pop('_timing', None)
        st['timing'][key + '_total'] = dt
        st['done'][key] = rec
        atomic_json(ck, st)
        did += 1
    res = build_result(list(st['done'].values()), sha)
    atomic_json(out, res)
    log(f"grid finished: {len(nodes)} nodes, all certified = {res['all_nodes_certified']}, no gaps = {res['no_gaps']}, "
        f"max cell rel width = {res['max_cell_rel_width']}, invocations = {st['invocations']}")
    return 0


# ================================================================== controls
def run_controls(out, log):
    T0 = time.time()
    sha = sha_self()
    R = dict(instrument_sha256=sha, started=time.strftime('%Y-%m-%d %H:%M:%S'))
    oks = {}

    # --- T0b: j = 0, U_N = v1 -> REFUSED/TRIVIAL with S == 0, p(0) == 1/40 exactly
    pc, qcs = v1_cheb(40, 8)
    r0, fn0 = certify_node(Fraction(0), 40, 8, log=log, keep=True)
    S0 = qparse(r0['S'])
    st0 = r0['status']
    R['T0b'] = dict(status=st0, S=qstr(S0), S_is_zero=bool(S0 == 0), p0=qstr(fn0.p0), p0_is_1_40=bool(fn0.p0 == fmpq(1, 40)))
    oks['T0b'] = bool(st0 == 'REFUSED/TRIVIAL' and S0 == 0 and fn0.p0 == fmpq(1, 40))
    log(f"T0b j=0, U=v1: status {st0}, S == 0: {S0 == 0}, p(0) == 1/40: {fn0.p0 == fmpq(1, 40)}")

    # --- T6a: nonzero rational j, U_N = v1 -> Q == src (1 - (1 + (1-s) j^2 v1)^7) term for term (sympy, independent)
    import sympy as sp
    t6 = {}
    for jt in (Fraction(7, 3), Fraction(10)):
        tau = fq(jt) ** 2
        fnv = ExactFn([c * tau for c in pc], [[c * tau for c in q] for q in qcs], tau)   # u-scale = tau * v1
        Sv, Qe = fnv.residual_S(want_Q=True)
        s_, c_ = sp.symbols('s c')
        pS = ((1 - s_) ** 4 + 5 * (1 - s_) ** 3 * s_ + 10 * (1 - s_) ** 2 * s_ ** 2 + 5 * (1 - s_) * s_ ** 3 + s_ ** 4) / 40
        q2S = -(1 - s_) * s_ / 20
        v1S = s_ * pS + s_ ** 2 * (1 - s_) * q2S * (3 * c_ ** 2 - 1) / 2
        jr = sp.Rational(jt.numerator, jt.denominator)
        base = sp.Poly(sp.expand(1 + (1 - s_) * jr ** 2 * v1S), s_, c_)
        rhs = sp.Poly(sp.Rational(9, 4) * (1 - c_ ** 2) * s_ ** 2 * (1 - s_) ** 2, s_, c_) * (sp.Poly(1, s_, c_) - base ** 7)
        rd = {(int(m[0]), int(m[1])): sp.Rational(v) for m, v in rhs.terms()}
        ed = unpack(Qe, fnv.W)
        same = set(rd.keys()) == set(ed.keys()) and all(fmpq(int(rd[k].p), int(rd[k].q)) == ed[k] for k in rd)
        t6[str(jt)] = dict(identity_exact=bool(same), terms=len(rd), S_positive=bool(Sv > 0))
        log(f"T6a j={jt}: engine Q == (9/4)(1-c^2)s^2(1-s)^2(1-(1+(1-s)j^2 v1)^7) term for term: {same} ({len(rd)} terms); S > 0: {Sv > 0}")
    R['T6a'] = t6
    oks['T6a'] = all(v['identity_exact'] for v in t6.values())

    # --- NaN / inf / 1e300 guards
    g = {}
    N, K = degrees(2)
    nco = (N + 1) * (K + 1)
    wa, _ = fit_u(Fraction(2), N, K)
    base_f = [float(a) for a in wa]
    for bad in (float('nan'), float('inf'), 1e300):
        cf = list(base_f)
        cf[N + 1 + 3] = bad
        r, _ = certify_node(Fraction(2), N, K, log=log, coef_override=cf)
        g[str(bad)] = r['status']
    R['nan_inf_guard'] = g
    oks['nan_inf_guard'] = all(v == 'INVALID' for v in g.values())
    log(f"guard: {g}")

    # --- point certificates at control nodes
    g36, g37, g116 = Fraction(51, 50) ** 36, Fraction(51, 50) ** 37, Fraction(51, 50) ** 116
    cnodes = [Fraction(1), Fraction(2), g36, Fraction(41, 20), g37, g116, Fraction(199, 20), Fraction(10)]
    recs, fns = {}, {}
    for jc in cnodes:
        t = time.time()
        rec, fnc = certify_node(jc, log=log, keep=True)
        rec['_secs'] = time.time() - t
        recs[jc], fns[jc] = rec, fnc
    R['control_nodes'] = [{k: v for k, v in recs[jc].items() if k not in ('_timing',)} for jc in cnodes]
    oks['control_nodes_certified'] = all(recs[jc]['status'] == 'CERTIFIED' for jc in cnodes)
    R['node_seconds'] = {str(float(jc)): recs[jc]['_secs'] for jc in cnodes}

    # --- T0c: inside stage-0 bracket, RN, chi < 1 (all control nodes)
    oks['T0c'] = all(recs[jc].get('t0c', {}).get('ok', False) for jc in cnodes)
    R['T0c'] = {str(float(jc)): recs[jc].get('t0c') for jc in cnodes}
    log(f"T0c (stage-0 bracket, RN, chi < 1) at all control nodes: {oks['T0c']}")

    # --- N1, N2, N4 at j = 2 and j = 10
    n1, n2, n4 = {}, {}, {}
    for jc in (Fraction(2), Fraction(10)):
        fn = fns[jc]
        eta = qparse(recs[jc]['eta'])
        tau = fn.tau
        # N1: sub-solution slope raised by 1e-3 relative: base_sub = U_N + 1e-3 p_N(0) s  (T_0 coefficient of p)
        pc1 = list(fn.pc)
        pc1[0] = pc1[0] + fn.p0 / 1000
        fsub = ExactFn(pc1, fn.qcs, tau)
        n1[str(jc)] = certify_pair(fn, eta, base_sub=fsub, log=log, tag=f"N1 j={jc}")
        # N2: U_+ - 1e-6 s
        n2[str(jc)] = certify_pair(fn, eta, e_sup=eta - fmpq(1, 10 ** 6), log=log, tag=f"N2 j={jc}")
        # N4: U_+ = U_N + eta' s, eta' = 0.25 x (S/2): a quarter of the needed margin
        S = qparse(recs[jc]['S'])
        n4[str(jc)] = certify_pair(fn, eta, e_sup=S / 8, log=log, tag=f"N4 j={jc}")
    R['N1'], R['N2'], R['N4'] = n1, n2, n4
    oks['N1'] = all(v['status'] == 'VIOLATION' for v in n1.values())
    oks['N2'] = all(v['status'] == 'VIOLATION' for v in n2.values())
    oks['N4'] = all(v['status'] == 'FAIL' and any('bound path' in p and 'FAILS' in p for p in v['path'])
                    and v['path'][0].startswith('point checks: sup OK') for v in n4.values())
    # positive version of the pair checker (same code path, true margins) must PASS
    R['pair_positive'] = {str(jc): certify_pair(fns[jc], qparse(recs[jc]['eta']), log=log, tag=f"pair+ j={jc}")['status']
                          for jc in (Fraction(2), Fraction(10))}
    oks['pair_positive'] = all(v == 'PASS' for v in R['pair_positive'].values())

    # --- N-pos at j = 10: shift V_- below its certified Bernstein minimum -> not PASS; shift up -> PASS
    fn = fns[Fraction(10)]
    eta = qparse(recs[Fraction(10)]['eta'])
    tau = fn.tau
    bmin_u = qparse(recs[Fraction(10)]['positivity']['bern_min'])      # tau-scaled
    mpt_u = V_point_min(fn, eta)                                         # exact min of tau V_- on a 65 x 9 grid
    down = eta + mpt_u / tau + fmpq(1, 10 ** 6)                          # U-scale shift: V_- - (m_pt + 1e-6)
    rdn = bernstein_positive(fn.V_terms(down))
    rup = bernstein_positive(fn.V_terms(eta - fmpq(1, 10 ** 6)))
    R['N_pos'] = dict(bern_min_U=qfloat(bmin_u / tau), point_min_U=qfloat(mpt_u / tau),
                      shift_down_is_below_bern_min=bool(mpt_u / tau + fmpq(1, 10 ** 6) > bmin_u / tau),
                      down=rdn, up=rup)
    log(f"N-pos j=10: Bernstein min {qfloat(bmin_u / tau):.6e}, grid point min {qfloat(mpt_u / tau):.6e}; "
        f"shifted down -> {rdn['status']}; shifted up -> {rup['status']}")
    # interior variant (corners stay positive, so the subdivision path is exercised):
    # V'' = V_- - beta 4 s (1-s). Down: beta = min_c V_-(1/2, c) + 1e-6 (negative at s = 1/2) -> never PASS.
    # Up: beta = 0.95 x grid-min of V_-/(4 s (1-s)) over interior grid points -> PASS, expected after subdivision.
    base_terms = fn.V_terms(eta)
    bump = fmpq_poly([0, 4, -4])

    def tV(s0, c0):
        return sum((f(s0) * g(c0) for f, g in base_terms), fmpq(0))
    vmid = min(tV(fmpq(1, 2), fmpq(k, 8)) for k in range(9))
    ratio = min(tV(fmpq(i, 64), fmpq(k, 8)) / (4 * fmpq(i, 64) * (1 - fmpq(i, 64))) for i in range(1, 64) for k in range(9))
    rdn2 = bernstein_positive(base_terms + [(-(vmid + tau * fmpq(1, 10 ** 6)) * bump, ONE)])
    for fac in (fmpq(95, 100), fmpq(99, 100), fmpq(999, 1000), fmpq(9999, 10000)):
        rup2 = bernstein_positive(base_terms + [(-(fac * ratio) * bump, ONE)])
        rup2['factor'] = qstr(fac)
        if rup2['status'] != 'PASS' or rup2['pieces'] > 1:
            break
    R['N_pos']['interior_down'] = rdn2
    R['N_pos']['interior_up'] = rup2
    oks['N_pos'] = bool(rdn['status'] != 'PASS' and rup['status'] == 'PASS' and R['N_pos']['shift_down_is_below_bern_min']
                        and rdn2['status'] != 'PASS' and rup2['status'] == 'PASS' and rup2['pieces'] > 1)
    log(f"N-pos interior j=10: bump down -> {rdn2['status']} ({rdn2.get('pieces')} pieces passed before stop); "
        f"bump up (factor {rup2['factor']}) -> {rup2['status']} ({rup2.get('pieces')} pieces)")

    # --- N-cell: correct cells contain the mid-cell certificates; swapped cells are empty or exclude them
    ncell = {}
    ctx.prec = 256
    for (ja, jm, jb) in ((g36, Fraction(41, 20), g37), (g116, Fraction(199, 20), Fraction(10))):
        ra, rm, rb = recs[ja], recs[jm], recs[jb]
        cr = cell_record(ra, rb)
        Ajm = arb(fq(jm))
        lo_ok = arb(qparse(cr['R_lo'])) * Ajm.sqrt()
        hi_ok = arb(qparse(cr['R_hi'])) * Ajm.sqrt()
        mlo, mhi = arb(qparse(rm['delta_lo'])), arb(qparse(rm['delta_hi']))
        contains = bool(lo_ok.lower() <= mlo and mhi <= hi_ok.upper())
        # swapped orientation: lower from j_b's sub-solution, upper from j_a's super-solution
        sw_lo = (arb(qparse(rb['delta_lo'])) * (Ajm / arb(fq(jb))).sqrt()).lower()
        sw_hi = (arb(qparse(ra['delta_hi'])) * (Ajm / arb(fq(ja))).sqrt()).upper()
        sw_empty = bool(sw_lo > sw_hi)
        sw_excl = bool(mhi < sw_lo or mlo > sw_hi)
        ncell[str(jm)] = dict(correct_contains_mid=contains, swapped_empty=sw_empty, swapped_excludes_mid=sw_excl,
                              cell_rel_width=cr['rel_width'], chi_min=cr['chi_min_dec'], chi_max=cr['chi_max_dec'])
        log(f"N-cell mid j={jm}: correct cell contains mid certificate {contains}; swapped empty {sw_empty}, excludes {sw_excl}; "
            f"cell rel width {cr['rel_width']:.4e}")
    R['N_cell'] = ncell
    oks['N_cell'] = all(v['correct_contains_mid'] and (v['swapped_empty'] or v['swapped_excludes_mid']) for v in ncell.values())

    # --- RESUME: uninterrupted vs killed+resumed vs budget-interrupted, byte-identical; sha change refused
    scratch = os.path.join(HERE, 'seats', 'build-p2a', 'resume')
    if os.path.isdir(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)
    me = [sys.executable, os.path.abspath(__file__)]
    oa, ob, oc, od = (os.path.join(scratch, x) for x in ('A.json', 'B.json', 'C.json', 'D.json'))
    t = time.time()
    rA = subprocess.run(me + ['--_grid', oa, '--_n', '6'], capture_output=True, text=True)
    tA = time.time() - t
    pB = subprocess.Popen(me + ['--_grid', ob, '--_n', '6'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    killed_at = None
    while pB.poll() is None:
        try:
            nd = len(json.load(open(ob + '.ckpt'))['done'])
        except Exception:
            nd = 0
        if nd >= 3:
            pB.kill()          # our own child, by PID
            pB.wait()
            killed_at = nd
            break
        time.sleep(0.1)
    rB2 = subprocess.run(me + ['--_grid', ob, '--_n', '6'], capture_output=True, text=True)
    invC = 0
    while not os.path.exists(oc) and invC < 10:
        subprocess.run(me + ['--_grid', oc, '--_n', '6', '--_budget', '1'], capture_output=True, text=True)
        invC += 1
    bA = open(oa, 'rb').read() if os.path.exists(oa) else b'A'
    bB = open(ob, 'rb').read() if os.path.exists(ob) else b'B'
    bC = open(oc, 'rb').read() if os.path.exists(oc) else b'C'
    stD = json.load(open(oa + '.ckpt'))
    stD['instrument_sha256'] = '0' * 64
    stD['done'].pop(list(stD['done'].keys())[-1])
    atomic_json(od + '.ckpt', stD)
    rD = subprocess.run(me + ['--_grid', od, '--_n', '6'], capture_output=True, text=True)
    invB = json.load(open(ob + '.ckpt'))['invocations']
    R['RESUME'] = dict(uninterrupted_rc=rA.returncode, uninterrupted_seconds=tA, killed_after_nodes=killed_at,
                       resumed_rc=rB2.returncode, B_invocations=invB, C_invocations=json.load(open(oc + '.ckpt'))['invocations'],
                       A_eq_B=bool(bA == bB), A_eq_C=bool(bA == bC), sha_changed_rc=rD.returncode,
                       sha_changed_no_output=not os.path.exists(od), A_sha256=hashlib.sha256(bA).hexdigest())
    oks['RESUME'] = bool(rA.returncode == 0 and killed_at is not None and killed_at < 6 and invB >= 2 and bA == bB and bA == bC
                         and R['RESUME']['C_invocations'] >= 2 and rD.returncode == 3 and not os.path.exists(od))
    log(f"RESUME: A==B (killed by PID after {killed_at} nodes, resumed) {bA == bB}; A==C (budget-interrupted, "
        f"{R['RESUME']['C_invocations']} invocations) {bA == bC}; changed-sha checkpoint refused rc={rD.returncode}")

    # --- --main refuses without GO-P2a
    if os.path.exists(GO_FILE):
        R['main_refusal'] = 'not exercised: GO-P2a exists (never start --main from the controls)'
        oks['main_refusal'] = True
    else:
        rM = subprocess.run(me + ['--main', os.path.join(scratch, 'M.json')], capture_output=True, text=True)
        R['main_refusal'] = dict(rc=rM.returncode, out=rM.stdout.strip()[-200:])
        oks['main_refusal'] = bool(rM.returncode == 2 and not os.path.exists(os.path.join(scratch, 'M.json')))
    log(f"--main without GO-P2a: {R['main_refusal']}")

    # --- j = 1 node for the lead's comparison with P1 and REP1
    r1 = recs[Fraction(1)]
    R['j1_delta'] = dict(lo=r1.get('delta_lo_dec'), hi=r1.get('delta_hi_dec'), rel_width=r1.get('rel_width'))
    R['oks'] = oks
    R['resume_tested'] = bool(oks['RESUME'])
    R['all_ok'] = all(oks.values())
    R['seconds'] = time.time() - T0
    for rr in R['control_nodes']:
        rr.pop('_secs', None)
    atomic_json(out, R, sort=False)
    log(f"controls done in {R['seconds']:.1f}s; all_ok = {R['all_ok']}; {oks}; instrument sha256 = {sha}")
    return 0 if R['all_ok'] else 1


def run_probe(args, log):
    js = [Fraction(args.j)] if args.j else [Fraction(1), Fraction(2), Fraction(10)]
    for j in js:
        N, K = (args.N, args.K) if args.N else degrees(j)
        t = time.time()
        rec, _ = certify_node(j, N, K, log=log)
        log(f"probe j={j} N={N} K={K}: {rec['status']} rel {rec.get('rel_width')} S {rec.get('S_float')} "
            f"total {time.time()-t:.1f}s timing {rec.get('_timing')}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--probe', action='store_true')
    ap.add_argument('--j')
    ap.add_argument('--N', type=int)
    ap.add_argument('--K', type=int)
    ap.add_argument('--controls')
    ap.add_argument('--main')
    ap.add_argument('--_grid')
    ap.add_argument('--_n', type=int)
    ap.add_argument('--_budget', type=float, default=TIME_BUDGET)
    args = ap.parse_args()
    logf = None

    def log(msg):
        print(msg, flush=True)
        if logf:
            logf.write(msg + '\n')
            logf.flush()

    if args.controls:
        logf = open(os.path.splitext(args.controls)[0] + '.log', 'w')
        sys.exit(run_controls(args.controls, log))
    elif args.main:
        if not os.path.exists(GO_FILE):
            print("REFUSED: <LAB>/BHSLICE/GO-P2a does not exist (the lead arms it with lab.py go after sealing).")
            sys.exit(2)
        logf = open(os.path.splitext(args.main)[0] + '.log', 'a')
        log(f"--- main invocation {time.strftime('%Y-%m-%d %H:%M:%S')} sha256 {sha_self()}")
        sys.exit(run_grid(args.main, grid_nodes(), TIME_BUDGET, log))
    elif args._grid:
        nodes = grid_nodes()[:args._n] if args._n else grid_nodes()
        sys.exit(run_grid(args._grid, nodes, args._budget, log))
    elif args.probe:
        run_probe(args, log)
