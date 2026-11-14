"""bowen_york_sqrt: certified sqrt-map sandwich for the Bowen-York spinning puncture (module version 2, package v4).

In plain words: for one spin value it proves that the true solution of the constraint lies between two explicit
functions, and so gives an exact rational interval for delta = M_ADM - 1 (bare mass mu = 1). It is a second,
independent certifier of the same PDE as lichnerowicz_axisym.py, with a different radial map, unknown, basis, margin
function and residual bound. At j = 1e4 it is many orders of magnitude sharper at lower cost (rep4/probe_A10000.log).

Object:  Delta u = -(9/4) J^2 sin^2(theta) / (r^6 Phi^7),  Phi = 1 + 1/(2r) + u  (mu = 1),  u -> 0 at infinity,
         j := 4J,  delta := M_ADM - 1 = 2 lim r u.
Method (rep4.py, re-derived by REFEREE-REP4 sec. 2):
  * Y = 1/(1 + sqrt(r/a)), a = 1/2 (Y = 0 is infinity, Y = 1 the puncture); v = cos^2(theta), tau_k(v) = T_k(2v - 1);
  * u = Y^2 W,  W = P(Y) + Y(1-Y) sum_{k=1..K} q_k(Y) tau_k(v),  P, q_k Chebyshev in 2Y - 1 of degree N;
  * r^2 (Delta u + S Phi^-7) = (1-Y) Y^3 Gh[W],
        Gh[W] = -(1/4)[6W + (6Y-3) W_Y - Y(1-Y) W_YY] + sum_k q_k A_k(v) + CS (1-v) Y^5 (1-Y)^5 / E^7,
        E = 2a(1-Y)^2 (1 + Y^2 W) + Y^2 = 2a(1-Y)^2 Phi,  CS = 18 a^3 j^2,  A_k = 4v(1-v) tau_k'' + (2-6v) tau_k';
  * Lin[1] = -3/2 (Y^2 strictly superharmonic), so W_N +- eta is a super/sub pair when (3/2) eta > sup|Gh[W_N]| and
    E_- = E_N - 2a(1-Y)^2 Y^2 eta > 0;
  * sup|Gh[W_N]| <= max over Y-pieces of sup|Q| / inf E_N^7, Q = E_N^7 Lin[W_N] + CS(1-v)Y^5(1-Y)^5, a polynomial of
    degree (8N+44, 8K), sampled in arb on a tensor Chebyshev-Gauss grid with one more node than the degree per piece, so
    the discrete Chebyshev transform gives its coefficients exactly (balls); sup|Q| <= sum|c|; inf E_N by the exact-degree
    form c00 - sum'|c| with bisection;
  * comparison: LEMMA-REP4.md (H^1_loc / Kato class; the sandwich has a r^-1/2 gradient cusp at the puncture, no point
    mass, Laplacian in L^1_loc), giving delta in [2a(W_N(0) - eta), 2a(W_N(0) + eta)].

Provenance and entry rule (gravity_engine/__init__.py, doctrine #121):
  * Rebuilt (re-typed, not copied blind) from <LAB>/BHSLICE/rep4/rep4.py, sha256 SOURCE_INSTRUMENT_SHA256 (sealed;
    read-only). REP4 is sealed claim ad62cca292 and its verdict is HOLD (REPLICATED; BHSLICE LOG.md 10:40, VERDICT-P2c
    REP4 paragraph). That HOLD is a blind-REPLICATION verdict: REP4 reproduced P2c's delta(3000), delta(10000) inside
    P2c's intervals. It is the HOLD verdict backing this code. Referee: REFEREE-REP4 (referee M; method sound, five
    fixes). Comparison lemma: <LAB>/BHSLICE/LEMMA-REP4.md (written for this module, NOT yet refereed).
  * Design: <LAB>/BHSLICE/DESIGN-P2d.md sec. 3 and 4.3 (module plan and fix list).
  * Backs: the REP4 certified nodes delta(3000), delta(10000) (resolution A reproduced exactly, tests (a)), K1 (j = 2^-64,
    delta/j^2 contains 1/40). DESIGN-P2d nodes C-REP, DEG, TOY, LOG, WORD are the fixes folded in here.
Fixes relative to rep4.py (REFEREE-REP4):
  1. residual degree is ALWAYS (8N+44, 8K); the CS = 0 branch (rep4.py:410-411, degree (N+2, K)) is deleted (FIX 4);
  2. displayed numbers: upper ends rounded UP (dec_ceil / sci_up), lower ends DOWN (dec_floor / sci_down); nothing is
     truncated (REFEREE-REP4 sec. 4; REFEREE-P2c item 4);
  3. python-flint ctx.prec is set only inside `_prec(...)` blocks that restore the caller's precision in `finally`
     (doctrine #115); importing or calling this module leaves the global precision unchanged;
  4. adaptive pieces (DESIGN-P2d sec. 4.3 item 4): a box is halved while (sup E_N / inf E_N)^8 > split_ratio
     (default 5/4), with sup/inf from the exact-degree Chebyshev form of E_N (c00 +- sum'|c|); cap max_pieces (128)
     boxes, then the certificate FAILS. DEVIATION from DESIGN-P2d's Y-pieces: boxes are split in Y OR in v (whichever
     pure-direction Chebyshev variation is larger), because at Y = 1/2 the v-variation of E alone gives
     (max_v E / min_v E)^8 ~ 1.5 > 1.25 for j >~ 1e3 (double CHECK), so Y-only splitting never meets the rule.
     adaptive=False gives REP4's 14 fixed Y-pieces x [0, 1] and its arithmetic path exactly;
  5. arb iterative refinement of the double fit (rep4's refine: residual in arb at the collocation nodes, correction with
     the double LU of the Newton Jacobian, midpoints kept), now with its history and stop flag in the record;
  6. a status word: 'CERTIFIED' | 'FAIL' | 'REFUSED/TRIVIAL' (sup bound exactly 0: nothing to certify, e.g. J = 0 with
     W = 0). The REP4 toy (W = 0, CS = 0) could not fail; the replacement toy (nonzero W at CS = 0) is in the tests;
  7. asserts (REFEREE-TL sec. 4): adaptive pieces cover [0, 1] with no gap or overlap; every coefficient vector enters
     through `Ansatz`, which asserts the (N+1)(K+1) shape of the class of LEMMA-REP4 sec. 4 (structural, not certified);
     the fit is rounded to a dyadic grid scaled to the coefficient size (rep4.py:594-596), never an absolute 2^-(prec-40).
Wording (REFEREE-REP4 FIX 5): an "independent second certifier of the same PDE", not an independent replication of the
physics. Common-mode risk = the PDE and its constants; guarded only by K1 = 1/40 and the literature chi_max ~ 0.928.
  8. (v2, REFEREE-LEMMA FIX F) residual v-degree max(8K, 1): at K = 0 the source's (1-v) has v-degree 1 and the old
     count 8K = 0 under-sampled Q (a false bound: N = 2, K = 0, W = 0, j = 1 gave 0.2552 < 9/32). K >= 1 values unchanged,
     so every stored certificate (REP4 A, K1) is reproduced identically (tests).
Version 2 (2026-11-14) adds section 7, the TANGENT certificate of DESIGN-P2d sec. 2 / 4.3 as amended by AMEND-P2d-1 and
AMEND-P2d-2 (Lemma TL: LEMMA-REP4.md sec. 6): certify_tangent (TL-H for psi_+ = Y^2 (Z_N + eta_Z), bound with E_-, arb
refinement of Z, scaled dyadic rounding, structural asserts K >= 1, N_Z = N, Z ansatz, exact piece cover; 2L_hi exact),
tl_cell_bound / tl_tail_bound (exact rationals, brackets asserted > 0), psi = 0 reduction (2L = 1). Tests:
tests/test_bowen_york_tangent.py. The tangent certificate backs no sealed node yet (BHSLICE P2d, claim 9677804a51, is the
first run that uses it).
NOT in this module: mu = 0; adaptive pieces for the tangent certificate (AMEND-P2d-2: REP4's fixed pieces).
Determinism: the fit is double precision LAPACK; byte-identical certificates need the same thread count as the run being
reproduced (the module does not touch thread variables). Certificates are valid either way. No pickle anywhere: exact
rationals are stored as "p/q" strings.
"""
import math
import time
from contextlib import contextmanager
from fractions import Fraction

import numpy as np
from scipy.linalg import lu_factor, lu_solve
from flint import arb, arb_mat, fmpq, fmpq_poly, ctx

MODULE_VERSION = 2
SEAL_ID = 'ad62cca292'                        # REP4 sealed claim (HOLD, replication)
SOURCE_INSTRUMENT_SHA256 = 'eb4f8cc460b417ef0ec4789002b38ae385dbcfaf54c5f7f5e0f526bc54b1d10a'   # rep4/rep4.py
DESIGN_FILE = '<LAB>/BHSLICE/DESIGN-P2d.md'
LEMMA_FILE = '<LAB>/BHSLICE/LEMMA-REP4.md'
REFEREE_FILE = '<LAB>/BHSLICE/REFEREE-REP4.md'
CERTIFIED_NODES = ('REP4 delta(3000)', 'REP4 delta(10000)', 'REP4 K1')

A_MAP = fmpq(1, 2)
RES_A = (140, 16, 224)                         # REP4 resolution A (N, K, prec bits)
RES_B = (200, 24, 288)                         # REP4 resolution B
SPLIT_RATIO = fmpq(5, 4)                       # adaptive pieces: split while (sup E / inf E)^8 > this
MAX_PIECES = 128
# default False (REP4's fixed pieces, the HOLD-backed arithmetic path). The DESIGN-P2d rule (adaptive=True) needs more
# than MAX_PIECES boxes for j >~ 1e3 (measured at (N, K) = (60, 8): 306 boxes at j = 1e3, 590 at j = 1e4), so with the
# design's cap it FAILS there; see README CHANGELOG v4.
ADAPTIVE_DEFAULT = False


@contextmanager
def _prec(bits):
    """run a block at `bits` of python-flint precision; the caller's precision is restored in `finally` (#115)."""
    old = ctx.prec
    ctx.prec = bits
    try:
        yield
    finally:
        ctx.prec = old


def CS_of(j_q, eps=fmpq(0)):
    """CS = 18 a^3 j^2 (1 + eps); eps != 0 only for the K2-type perturbation control."""
    return 18 * A_MAP ** 3 * j_q * j_q * (1 + eps)


def qf(q):
    return float(Fraction(int(q.p), int(q.q)))


def qstr(q):
    return "%s/%s" % (q.p, q.q)


def parse_q(s):
    a, b = s.split('/')
    return fmpq(int(a), int(b))


def _q_of(m, e):
    return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))


def arb_upper_q(x):
    """exact rational upper bound of an arb ball (midpoint + radius, both exact dyadics)."""
    m, e = x.mid().man_exp()
    r, er = x.rad().mid().man_exp()
    return _q_of(int(m), int(e)) + _q_of(int(r), int(er))


def arb_lower_q(x):
    return -arb_upper_q(-x)


# ===================================================================================== display (rounded outward)
def _ilog10(q):
    """floor(log10(q)) for an exact rational q > 0."""
    p, d = int(q.p), int(q.q)
    e = math.floor(math.log10(p) - math.log10(d))
    while fmpq(10) ** (e + 1) <= q if e + 1 >= 0 else q * fmpq(10) ** (-(e + 1)) >= 1:
        e += 1
    while (q < fmpq(10) ** e) if e >= 0 else (q * fmpq(10) ** (-e) < 1):
        e -= 1
    return e


def _scale10(q, k):
    return q * fmpq(10) ** k if k >= 0 else q / fmpq(10) ** (-k)


def _floor_q(v):
    return int(v.p) // int(v.q)


def dec_floor(q, nd):
    """decimal string <= q with nd significant digits (a lower end)."""
    if q == 0:
        return "0"
    e = _ilog10(abs(q))
    sc = nd - 1 - e
    iv = _floor_q(_scale10(q, sc))
    sg = "-" if iv < 0 else ""
    s = str(abs(iv))
    if sc <= 0:
        return sg + s + "0" * (-sc)
    if len(s) <= sc:
        s = "0" * (sc - len(s) + 1) + s
    return sg + s[:-sc] + "." + s[-sc:]


def dec_ceil(q, nd):
    """decimal string >= q with nd significant digits (an upper end)."""
    s = dec_floor(-q, nd)
    return s[1:] if s.startswith("-") else ("-" + s if s != "0" else s)


def sci_up(q, nd=4):
    """scientific string >= q (upper bounds in logs and piece records)."""
    return _sci(q, nd, True)


def sci_down(q, nd=4):
    """scientific string <= q (lower bounds)."""
    return _sci(q, nd, False)


def _sci(q, nd, up):
    if q == 0:
        return "0"
    neg = q < 0
    x = -q if neg else q
    e = _ilog10(x)
    v = _scale10(x, nd - 1 - e)
    ceil_mag = (up and not neg) or (not up and neg)
    iv = -_floor_q(-v) if ceil_mag else _floor_q(v)
    if iv >= 10 ** nd:                        # 9.999.. rounded up to 10.00..
        iv //= 10
        e += 1
    s = str(iv)
    return ("-" if neg else "") + s[0] + ("." + s[1:] if len(s) > 1 else "") + "e%+03d" % e


# ===================================================================================== 1. double collocation (fit only)
def cheb_cf(x, N):
    ph = np.arccos(x)
    n = np.arange(N + 1)[None, :]
    P = ph[:, None]
    T = np.cos(n * P)
    s = np.sin(P)
    T1 = n * np.sin(n * P) / s
    T2 = (x[:, None] * T1 - n ** 2 * T) / (1 - x[:, None] ** 2)
    return T, T1, T2


def gauss01(N):
    i = np.arange(N)
    return (1 + np.cos(np.pi * (i + 0.5) / N)) / 2


class Disc:
    """Gauss collocation of Gh on (N+1) x (K+1) nodes; unknowns = Chebyshev coefficients of P, q_1..q_K."""

    def __init__(self, N, K):
        self.N, self.K = N, K
        a = qf(A_MAP)
        Y = gauss01(N + 1)
        v = gauss01(K + 1)
        T, T1, T2 = cheb_cf(2 * Y - 1, N)
        TY, TYY = 2 * T1, 4 * T2
        y = Y[:, None]

        def LR(f, f1, f2):
            return -0.25 * (6 * f + (6 * y - 3) * f1 - y * (1 - y) * f2)
        RF = LR(T, TY, TYY)
        h = y - y * y
        h1 = 1 - 2 * y
        h2 = -2.0
        g = h * T
        g1 = h1 * T + h * TY
        g2 = h2 * T + 2 * h1 * TY + h * TYY
        RG = LR(g, g1, g2)
        Tv, Tv1, Tv2 = cheb_cf(2 * v - 1, K)
        Tv = Tv[:, 1:]
        Tv1 = 2 * Tv1[:, 1:]
        Tv2 = 4 * Tv2[:, 1:]
        vv = v[:, None]
        AK = 4 * vv * (1 - vv) * Tv2 + (2 - 6 * vv) * Tv1
        n1 = N + 1
        L = np.zeros((n1 * (K + 1), n1 * (K + 1)))
        U = np.zeros_like(L)
        for jv in range(K + 1):
            rows = np.arange(n1) * (K + 1) + jv
            L[rows, :n1] = RF
            U[rows, :n1] = T
            for k in range(1, K + 1):
                L[rows, k * n1:(k + 1) * n1] = RG * Tv[jv, k - 1] + T * AK[jv, k - 1]
                U[rows, k * n1:(k + 1) * n1] = g * Tv[jv, k - 1]
        self.L, self.U = L, U
        self.Y = np.repeat(Y, K + 1)
        self.V = np.tile(v, n1)
        self.c2 = 2 * a * (1 - self.Y) ** 2

    def F(self, w, CS):
        W = self.U @ w
        E = self.c2 * (1 + self.Y ** 2 * W) + self.Y ** 2
        src = CS * (1 - self.V) * self.Y ** 5 * (1 - self.Y) ** 5 * E ** -7
        return self.L @ w + src, E, src

    def jac(self, w, CS):
        R, E, src = self.F(w, CS)
        return self.L - (7 * src / E * self.c2 * self.Y ** 2)[:, None] * self.U

    def newton(self, CS, w, iters=60):
        for it in range(iters):
            R, E, s = self.F(w, CS)
            dw = np.linalg.solve(self.jac(w, CS), R)
            w = w - dw
            if np.max(np.abs(dw)) <= 1e-15 * max(1e-300, np.max(np.abs(w))):
                break
        R, E, s = self.F(w, CS)
        return w, float(np.max(np.abs(R))), it + 1


def prolong(w, N, K, N2, K2):
    n1, m1 = N + 1, N2 + 1
    out = np.zeros(m1 * (K2 + 1))
    for k in range(min(K, K2) + 1):
        nn = min(n1, m1)
        out[k * m1:k * m1 + nn] = w[k * n1:k * n1 + nn]
    return out


def continuation(j_target, CS_scale, N=60, K=8, ratio=1.15, j_start=0.5):
    """coarse continuation in j (CS = CS_scale j^2) at (N, K); returns the coarse solution at j_target."""
    disc = Disc(N, K)
    w = np.zeros(disc.L.shape[0])
    x = min(j_start, j_target)
    while x < j_target:
        w, res, it = disc.newton(CS_scale * x * x, w)
        x *= ratio
    w, res, it = disc.newton(CS_scale * j_target ** 2, w)
    return w, (N, K), res


# ===================================================================================== 2. arb fields and refinement
def cheb_deriv(c):
    N = len(c) - 1
    if N == 0:
        return [c[0] * 0]
    d = [c[0] * 0] * (N + 2)
    for k in range(N - 1, -1, -1):
        d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1]
    d[0] = d[0] / 2
    return d[:N]


def ydiff(c):
    """Chebyshev coefficients of d/dY (variable 2Y - 1)."""
    return [2 * x for x in cheb_deriv(c)]


def pad(c, n):
    return list(c) + [c[0] * 0] * (n - len(c))


def tau_polys(K):
    """tau_k(v) = T_k(2v-1) and A_k = 4v(1-v) tau_k'' + (2-6v) tau_k' as exact fmpq_polys, k = 0..K."""
    x = fmpq_poly([-1, 2])
    V = fmpq_poly([0, 1])
    T = [fmpq_poly([1]), x]
    for n in range(1, K):
        T.append(2 * x * T[n] - T[n - 1])
    T = T[:K + 1]
    A = [4 * V * (1 - V) * t.derivative().derivative() + (2 - 6 * V) * t.derivative() for t in T]
    return T, A


def horner(poly, x):
    cs = [arb(c) for c in poly.coeffs()]
    if not cs:
        return arb(0)
    r = cs[-1]
    for c in reversed(cs[:-1]):
        r = r * x + c
    return r


def LR_arb(y, f, f1, f2):
    return -(6 * f + (6 * y - 3) * f1 - y * (1 - y) * f2) / 4


class Fields:
    """W = P + Y(1-Y) sum q_k tau_k from Chebyshev coefficient lists P, q[1..K]; arb grids of W and Lin[W]."""

    def __init__(self, p, q):
        self.p = list(p)
        self.q = [list(x) for x in q]
        self.N = len(p) - 1
        self.K = len(q)
        n1 = self.N + 1
        cols = []
        for c in [self.p] + self.q:
            d1 = ydiff(c)
            d2 = ydiff(d1)
            cols += [pad(c, n1), pad(d1, n1), pad(d2, n1)]
        self.cols = cols

    def radial(self, ys):
        n1 = self.N + 1
        rows = []
        for y in ys:
            ph = (2 * y - 1).acos()
            rows.append([(ph * n).cos() for n in range(n1)])
        TY = arb_mat(rows)
        C = arb_mat([[arb(col[n]) for col in self.cols] for n in range(n1)])
        P = TY * C
        Wr = []
        Lr = []
        K = self.K
        for i, y in enumerate(ys):
            p, p1, p2 = P[i, 0], P[i, 1], P[i, 2]
            wr = [p]
            lr = [LR_arb(y, p, p1, p2)]
            h = y - y * y
            h1 = 1 - 2 * y
            qs = []
            for k in range(1, K + 1):
                q, q1, q2 = P[i, 3 * k], P[i, 3 * k + 1], P[i, 3 * k + 2]
                g = h * q
                g1 = h1 * q + h * q1
                g2 = -2 * q + 2 * h1 * q1 + h * q2
                wr.append(g)
                lr.append(LR_arb(y, g, g1, g2))
                qs.append(q)
            Wr.append(wr)
            Lr.append(lr + qs)
        return arb_mat(Wr), arb_mat(Lr)

    def angular(self, vs):
        T, A = tau_polys(self.K)
        Aw = [[arb(1)] * len(vs)] + [[horner(T[k], v) for v in vs] for k in range(1, self.K + 1)]
        Al = [[arb(1)] * len(vs)] + [[horner(T[k], v) for v in vs] for k in range(1, self.K + 1)] + \
             [[horner(A[k], v) for v in vs] for k in range(1, self.K + 1)]
        return arb_mat(Aw), arb_mat(Al)

    def grid(self, ys, vs):
        Wr, Lr = self.radial(ys)
        Aw, Al = self.angular(vs)
        return Wr * Aw, Lr * Al


def split(w, N, K):
    n1 = N + 1
    return w[:n1], [w[k * n1:(k + 1) * n1] for k in range(1, K + 1)]


def arb_nodes01(N):
    pi = arb.pi()
    return [(1 + (pi * fmpq(2 * i + 1, 2 * N)).cos()) / 2 for i in range(N)]


def refine(disc, wf, CS_q, prec, maxit=25):
    """arb iterative refinement of the double Newton fit: residual Gh in arb at the collocation nodes, correction from
    the double LU of the Jacobian, keep midpoints. Stops at 2^-(prec-40) relative residual or on a plateau.
    Returns (coefficients as exact arb midpoints, history, target_reached). Fit only: the proof does not trust it."""
    N, K = disc.N, disc.K
    with _prec(prec):
        ys = arb_nodes01(N + 1)
        vs = arb_nodes01(K + 1)
        cs = arb(CS_q)
        a2 = 2 * arb(A_MAP)
        lu = lu_factor(disc.jac(wf, qf(CS_q)))
        w = [arb(float(x)) for x in wf]
        hist = []
        target = 2.0 ** (-(prec - 40))
        for it in range(maxit):
            p, q = split(w, N, K)
            W, LW = Fields(p, q).grid(ys, vs)
            R = []
            for i in range(N + 1):
                y = ys[i]
                om = 1 - y
                c2 = a2 * om * om
                y2 = y * y
                s0 = cs * y2 * y2 * y * om ** 5
                for jv in range(K + 1):
                    E = c2 * (1 + y2 * W[i, jv]) + y2
                    R.append(LW[i, jv] + s0 * (1 - vs[jv]) / E ** 7)
            Rf = np.array([float(r.mid()) for r in R])
            scale = max(1e-300, max(abs(float(x.mid())) for x in w))
            res = float(np.max(np.abs(Rf))) / scale
            hist.append(res)
            if res < target or (len(hist) > 3 and res > 0.5 * hist[-2]):
                break
            dw = lu_solve(lu, Rf)
            w = [(w[i] - arb(float(dw[i]))).mid() for i in range(len(w))]
        return w, hist, bool(hist[-1] < target)


def to_dyadic(x, bits):
    """arb midpoint -> nearest multiple of 2^-bits (exact fmpq)."""
    m, e = x.mid().man_exp()
    m = int(m)
    e = int(e)
    if e >= -bits:
        return _q_of(m, e)
    sh = -bits - e
    return fmpq((m + (1 << (sh - 1))) >> sh, 2 ** bits)


# ===================================================================================== 3. certificate
_COS = {}


def dct_mat(N):
    """N x N discrete Chebyshev transform for N Chebyshev-Gauss nodes (exact for degree <= N-1), cached per prec."""
    key = (N, ctx.prec)
    if key in _COS:
        return _COS[key]
    pi = arb.pi()
    tab = [(pi * fmpq(m, 2 * N)).cos() for m in range(4 * N)]
    rows = []
    for k in range(N):
        sc = (1 if k == 0 else 2)
        rows.append([tab[(k * (2 * i + 1)) % (4 * N)] * sc / N for i in range(N)])
    M = arb_mat(rows)
    _COS[key] = M
    return M


def clear_caches():
    _COS.clear()


def nodes_on(lo, hi, N):
    pi = arb.pi()
    mid = arb((lo + hi) / 2)
    half = arb((hi - lo) / 2)
    return [mid + half * (pi * fmpq(2 * i + 1, 2 * N)).cos() for i in range(N)]


def cheb_coeffs(vals, Ny_, Mv_):
    return dct_mat(Ny_) * vals * dct_mat(Mv_).transpose()


def abs_sum(C, skip00=False):
    s = arb(0)
    for i in range(C.nrows()):
        for j in range(C.ncols()):
            if skip00 and i == 0 and j == 0:
                continue
            s += abs(C[i, j])
    return s


def chebform_lower(vals, Ny_, Mv_):
    C = cheb_coeffs(vals, Ny_, Mv_)
    return C[0, 0] - abs_sum(C, skip00=True)


class Ansatz:
    """the exact dyadic ansatz W_N = P(Y) + Y(1-Y) sum_k q_k(Y) tau_k(v); W0 = W(Y = 0) = P(0) exactly.
    Structural class assert (LEMMA-REP4 sec. 4 and sec. 6, REFEREE-TL sec. 4 (iii)): every coefficient vector the proof
    sees enters here and is read in exactly this form (angular part times Y(1-Y), so W(1, v) = P(1)); a future tangent
    Z_+ must come through this class with the same (N, K)."""

    def __init__(self, wq, N, K):
        assert N >= 0 and K >= 0 and len(wq) == (N + 1) * (K + 1), 'ansatz shape must be (N+1)(K+1) coefficients'
        assert all(isinstance(x, type(fmpq(0))) for x in wq), 'ansatz coefficients must be exact fmpq'
        self.wq = wq
        self.N = N
        self.K = K
        p, q = split(wq, N, K)
        self.F = Fields(p, q)
        self.W0 = sum((wq[n] * (1 if n % 2 == 0 else -1) for n in range(N + 1)), fmpq(0))


def E_grid(ans, ylo, yhi, vlo, vhi):
    N = ans.N + 7                          # E has Y-degree N+6, v-degree K
    M = ans.K + 1
    ys = nodes_on(ylo, yhi, N)
    vs = nodes_on(vlo, vhi, M)
    W, _ = ans.F.grid(ys, vs)
    a2 = 2 * arb(A_MAP)
    vals = arb_mat(N, M, [a2 * (1 - ys[i]) ** 2 * (1 + ys[i] ** 2 * W[i, j]) + ys[i] ** 2
                          for i in range(N) for j in range(M)])
    return vals, N, M


def E_lower(ans, ylo, yhi, max_boxes=400, rel=fmpq(9, 10), vlo=fmpq(0), vhi=fmpq(1)):
    """rigorous lower bound of E_N on [ylo, yhi] x [vlo, vhi]; bisect until each box's form bound is > 0 and
    >= rel * (its smallest node value). Returns (lower bound or None, boxes)."""
    stack = [(ylo, yhi, vlo, vhi, 0)]
    nbox = 0
    lbmin = None
    while stack:
        a0, a1, b0, b1, dep = stack.pop()
        nbox += 1
        vals, N, M = E_grid(ans, a0, a1, b0, b1)
        lb = chebform_lower(vals, N, M)
        mn = min(vals[i, j].mid() for i in range(N) for j in range(M))
        if lb > 0 and lb >= arb(rel) * mn:
            lq = arb_lower_q(lb)
            lbmin = lq if lbmin is None or lq < lbmin else lbmin
            continue
        if nbox > max_boxes or dep > 16:
            if lb > 0:                     # still a rigorous bound, only looser
                lq = arb_lower_q(lb)
                lbmin = lq if lbmin is None or lq < lbmin else lbmin
                continue
            return None, nbox
        if (a1 - a0) >= (b1 - b0) / 4:
            m = (a0 + a1) / 2
            stack += [(a0, m, b0, b1, dep + 1), (m, a1, b0, b1, dep + 1)]
        else:
            m = (b0 + b1) / 2
            stack += [(a0, a1, b0, m, dep + 1), (a0, a1, m, b1, dep + 1)]
    return lbmin, nbox


def residual_degrees(N, K):
    """Y- and v-degree of Q = E^7 Lin[W] + CS (1-v) Y^5 (1-Y)^5: W (N+2, K), Lin[W] (N+2, K), E (N+6, K),
    E^7 Lin (8N+44, 8K); the source term has (10, 1). The same for every CS, including CS = 0 (REFEREE-REP4 FIX 4).
    v-degree max(8K, 1) (REFEREE-LEMMA FIX F, module v2): at K = 0 the source's (1-v) still has v-degree 1; with 8K = 0
    the single v-node under-sampled Q and certified a false bound (N = 2, K = 0, W = 0, j = 1: 0.2552 < 9/32)."""
    return 8 * N + 44, max(8 * K, 1)


def piece_bound(ans, ylo, yhi, CS_q, vlo=fmpq(0), vhi=fmpq(1)):
    """rigorous sup|Gh[W_N]| <= sup|Q| / inf E_N^7 on the box [ylo, yhi] x [vlo, vhi]. Returns ((bound, Elow) or None, info)."""
    N, K = ans.N, ans.K
    dQy, dQv = residual_degrees(N, K)
    NQ = dQy + 1
    MQ = dQv + 1
    ys = nodes_on(ylo, yhi, NQ)
    vs = nodes_on(vlo, vhi, MQ)
    W, LW = ans.F.grid(ys, vs)
    cs = arb(CS_q)
    a2 = 2 * arb(A_MAP)
    flat = []
    for i in range(NQ):
        y = ys[i]
        om = 1 - y
        y2 = y * y
        c2 = a2 * om * om
        s0 = cs * y2 * y2 * y * om ** 5
        for jv in range(MQ):
            E = c2 * (1 + y2 * W[i, jv]) + y2
            E2 = E * E
            E7 = E2 * E2 * E2 * E
            flat.append(E7 * LW[i, jv] + s0 * (1 - vs[jv]))
    Qv = arb_mat(NQ, MQ, flat)
    SQ = arb_upper_q(abs_sum(cheb_coeffs(Qv, NQ, MQ)))
    Elow, nb = E_lower(ans, ylo, yhi, vlo=vlo, vhi=vhi)
    if Elow is None or Elow <= 0:
        return None, dict(ylo=str(ylo), yhi=str(yhi), vlo=str(vlo), vhi=str(vhi), sup_Q_le=sci_up(SQ), Elow="FAIL")
    bound = SQ / Elow ** 7
    return (bound, Elow), dict(ylo=str(ylo), yhi=str(yhi), vlo=str(vlo), vhi=str(vhi), sup_Q_le=sci_up(SQ), Elow_ge=sci_down(Elow, 6),
                               E_boxes=nb, bound_le=sci_up(bound))


def round_up_dyadic(q, bits=40):
    """smallest multiple of 2^(e) >= q with e = floor(log2 q) - bits (exact)."""
    if q <= 0:
        return fmpq(0)
    e = math.floor(math.log2(qf(q))) - bits
    sc = fmpq(2) ** (-e)
    v = q * sc
    r = fmpq(-((-v.p) // v.q)) / sc
    assert r >= q
    return r


def make_pieces(n):
    """n uniform Y-pieces, the end pieces split geometrically 3 more times toward Y = 0 and Y = 1 (REP4 layout)."""
    h = fmpq(1, n)
    left = [h / 2 ** m for m in range(3, 0, -1)]
    right = [1 - h / 2 ** m for m in range(1, 4)]
    mid = [fmpq(i, n) for i in range(1, n)]
    return [fmpq(0)] + left + mid + right + [fmpq(1)]


def E_ratio_up(ans, ylo, yhi, vlo=fmpq(0), vhi=fmpq(1)):
    """from one exact-degree Chebyshev form of E_N on the box [ylo, yhi] x [vlo, vhi]: returns (ratio, sY, sV) with
    ratio = (c00 + s)/(c00 - s) >= sup E / inf E (None if c00 - s <= 0), s = sum'|c|, sY = sum_{i>=1} |c_i0| (pure Y
    variation), sV = sum_{k>=1} |c_0k| (pure v variation). Used only to choose boxes (never as a proof input)."""
    vals, N, M = E_grid(ans, ylo, yhi, vlo, vhi)
    C = cheb_coeffs(vals, N, M)
    s = abs_sum(C, skip00=True)
    sY = sum((abs(C[i, 0]) for i in range(1, C.nrows())), arb(0))
    sV = sum((abs(C[0, k]) for k in range(1, C.ncols())), arb(0))
    lo = C[0, 0] - s
    rt = arb_upper_q((C[0, 0] + s) / lo) if lo > 0 else None
    return rt, float(sY.mid()), float(sV.mid())


def adaptive_boxes(ans, n_pieces=8, split_ratio=SPLIT_RATIO, max_pieces=MAX_PIECES):
    """start from the Y-pieces make_pieces(n_pieces) x [0, 1]; halve a box while (sup E / inf E)^8 > split_ratio
    (DESIGN-P2d sec. 4.3 item 4), along Y when its pure-Y Chebyshev variation sY >= the pure-v variation sV, else along
    v. (Y-only splitting cannot meet the rule for j >~ 1e3: at Y = 1/2 the v-variation alone gives (max/min)^8 ~ 1.5.)
    Returns (boxes, None) or (None, reason) when more than max_pieces boxes would be needed."""
    work = make_pieces(n_pieces)
    done = []
    stack = [(work[i], work[i + 1], fmpq(0), fmpq(1)) for i in range(len(work) - 2, -1, -1)]
    while stack:
        y0, y1, v0, v1 = stack.pop()
        rt, sY, sV = E_ratio_up(ans, y0, y1, v0, v1)
        if rt is not None and rt ** 8 <= split_ratio:
            done.append((y0, y1, v0, v1))
            continue
        if len(done) + len(stack) + 2 > max_pieces:
            return None, "more than %d boxes needed" % max_pieces
        if sY >= sV:
            m = (y0 + y1) / 2
            stack += [(m, y1, v0, v1), (y0, m, v0, v1)]
        else:
            m = (v0 + v1) / 2
            stack += [(y0, y1, m, v1), (y0, y1, v0, m)]
    check_cover(done)
    done.sort(key=lambda b: (qf(b[0]), qf(b[2])))
    return done, None


def check_cover(boxes):
    """assert that closed boxes inside [0,1]^2 with pairwise disjoint interiors have total area exactly 1, i.e. they
    cover the square with no gap (rep4.py:463 pattern; REFEREE-TL sec. 4 (i)). Exact rational arithmetic."""
    assert boxes, 'no boxes'
    for (y0, y1, v0, v1) in boxes:
        assert 0 <= y0 < y1 <= 1 and 0 <= v0 < v1 <= 1, 'box outside the unit square or empty'
    for i in range(len(boxes)):
        a = boxes[i]
        for b in boxes[i + 1:]:
            assert not (a[0] < b[1] and b[0] < a[1] and a[2] < b[3] and b[2] < a[3]), 'overlapping boxes'
    area = sum(((y1 - y0) * (v1 - v0) for (y0, y1, v0, v1) in boxes), fmpq(0))
    assert area == 1, 'boxes do not cover [0,1]^2 (area %s)' % area
    return True


def certify(wq, N, K, CS_q, n_pieces=8, prec=256, adaptive=ADAPTIVE_DEFAULT, split_ratio=SPLIT_RATIO, max_pieces=MAX_PIECES,
            eta=None):
    """certificate for the exact dyadic ansatz wq (Chebyshev coefficients of P, q_1..q_K, degree N).
    adaptive=False: REP4's 14 Y-pieces x [0, 1] (arithmetic path of rep4.py); True: adaptive_boxes.
    eta: None = honest margin up((2/3) sup (1 + 2^-20)); an exact rational here is for CONTROLS only (forced margin).
    Returns a dict with exact fmpq values and status 'CERTIFIED' | 'FAIL' | 'REFUSED/TRIVIAL'."""
    tm = time.time()
    with _prec(prec):
        ans = Ansatz(wq, N, K)
        pieces = []
        sup = fmpq(0)
        Emin = None
        ok = True
        why = None
        if adaptive:
            boxes, why = adaptive_boxes(ans, n_pieces, split_ratio, max_pieces)
            if boxes is None:
                ok = False
                boxes = []
        else:
            Ed = make_pieces(n_pieces)
            assert Ed[0] == 0 and Ed[-1] == 1 and all(Ed[i] < Ed[i + 1] for i in range(len(Ed) - 1))
            boxes = [(Ed[i], Ed[i + 1], fmpq(0), fmpq(1)) for i in range(len(Ed) - 1)]
        if boxes:
            check_cover(boxes)
        for (y0, y1, v0, v1) in boxes:
            bd, info = piece_bound(ans, y0, y1, CS_q, v0, v1)
            pieces.append(info)
            if bd is None:
                ok = False
                why = why or "E_N lower bound failed on a box"
                continue
            sup = max(sup, bd[0])
            Emin = bd[1] if Emin is None or bd[1] < Emin else Emin
        if eta is None:
            eta = round_up_dyadic(sup * (1 + fmpq(1, 2 ** 20)) * fmpq(2, 3))
        strict = bool(ok and eta > 0 and fmpq(3, 2) * eta > sup)
        # E_- = E_N - 2a (1-Y)^2 Y^2 eta >= Emin - (2a/16) eta   ((1-Y) Y <= 1/4)
        Em_low = (Emin - 2 * A_MAP * eta / 16) if Emin is not None else None
        Em_ok = bool(Em_low is not None and Em_low > 0)
        if ok and sup == 0:
            status = 'REFUSED/TRIVIAL'
        elif strict and Em_ok:
            status = 'CERTIFIED'
        else:
            status = 'FAIL'
    return dict(sup_Gh_le=sup, eta=eta, W0=ans.W0, strict=strict, E_N_min_lower=Emin, E_minus_lower=Em_low,
                E_minus_positive=Em_ok, certified=(status == 'CERTIFIED'), status=status, fail_reason=why,
                pieces=pieces, n_pieces_used=len(pieces), all_pieces_ok=ok, adaptive=bool(adaptive),
                residual_degrees=residual_degrees(N, K), secs=round(time.time() - tm, 1))


# ===================================================================================== 4. exact-rational evaluation
def cheb_eval_q(c, X):
    b1 = fmpq(0)
    b2 = fmpq(0)
    for ck in reversed(c[1:]):
        b1, b2 = 2 * X * b1 - b2 + ck, b1
    return X * b1 - b2 + c[0]


class ExactEval:
    """W and Lin[W] (no source) at exact rational (Y, v); independent of the arb engine."""

    def __init__(self, wq, N, K):
        p, q = split(wq, N, K)
        self.P = [p, ydiff(p), ydiff(ydiff(p))]
        self.Q = [[c, ydiff(c), ydiff(ydiff(c))] for c in q]
        self.T, self.A = tau_polys(K)
        self.K = K

    def at(self, y, v):
        X = 2 * y - 1
        W, W1, W2 = [cheb_eval_q(c, X) for c in self.P]
        h = y - y * y
        h1 = 1 - 2 * y
        angA = fmpq(0)
        for k in range(1, self.K + 1):
            q, q1, q2 = [cheb_eval_q(c, X) for c in self.Q[k - 1]]
            tk = self.T[k](v)
            W += h * q * tk
            W1 += (h1 * q + h * q1) * tk
            W2 += (-2 * q + 2 * h1 * q1 + h * q2) * tk
            angA += q * self.A[k](v)
        L = -(6 * W + (6 * y - 3) * W1 - y * (1 - y) * W2) / 4 + angA
        return W, L


def sample_signs(wq, N, K, CS_q, eta, n_y=24, n_v=6):
    """exact rational Gh[W_N +- eta] on a (n_y+1) x (n_v+1) grid (a sample, not a proof): super < 0 < sub expected."""
    Ex = ExactEval(wq, N, K)
    a2 = 2 * A_MAP
    mx_plus = None
    mn_minus = None
    mx_abs = fmpq(0)
    for i in range(n_y + 1):
        y = fmpq(i, n_y)
        for l in range(n_v + 1):
            v = fmpq(l, n_v)
            W, L = Ex.at(y, v)
            src = CS_q * (1 - v) * y ** 5 * (1 - y) ** 5
            c2 = a2 * (1 - y) ** 2
            EN = c2 * (1 + y * y * W) + y * y
            Gn = L + src / EN ** 7
            Gp = L - fmpq(3, 2) * eta + src / (EN + c2 * y * y * eta) ** 7
            Gm = L + fmpq(3, 2) * eta + src / (EN - c2 * y * y * eta) ** 7
            mx_plus = Gp if mx_plus is None or Gp > mx_plus else mx_plus
            mn_minus = Gm if mn_minus is None or Gm < mn_minus else mn_minus
            mx_abs = max(mx_abs, abs(Gn))
    return dict(max_Gh_plus=mx_plus, min_Gh_minus=mn_minus, max_abs_Gh_N=mx_abs, signs_ok=bool(mx_plus < 0 < mn_minus))


def gauss_law_delta_double(wq, N, K, CS):
    """CHECK only (double, not rigorous): delta = 2 * 2a CS int_0^1 Y^5 (1-Y)^7 <(1-v) E^-7>_c dY (Gauss's law)."""
    a = qf(A_MAP)
    w = np.array([qf(x) for x in wq])
    n1 = N + 1
    xg, wg = np.polynomial.legendre.leggauss(600)
    yq = (xg + 1) / 2
    wy = wg / 2
    wc = wg / 2
    vq = xg ** 2
    T, _, _ = cheb_cf(2 * yq - 1, N)
    Tv, _, _ = cheb_cf(2 * vq - 1, K)
    p = T @ w[:n1]
    qv = np.stack([T @ w[k * n1:(k + 1) * n1] for k in range(1, K + 1)]) if K > 0 else None
    tot = 0.0
    for iy in range(len(yq)):
        y = yq[iy]
        W = p[iy] + (y - y * y) * (qv[:, iy] @ Tv[:, 1:].T if K > 0 else 0.0)
        E = 2 * a * (1 - y) ** 2 * (1 + y * y * W) + y * y
        tot += wy[iy] * y ** 5 * (1 - y) ** 7 * np.sum(wc * (1 - vq) * E ** -7)
    return 2 * 2 * a * CS * tot


# ===================================================================================== 5. drivers
def solve_point(j_q, N, K, prec, eps=fmpq(0)):
    """fit: continuation at (60, 8), Newton at (N, K), arb refinement, dyadic rounding.
    Returns (wq exact dyadics, fine Newton residual, refinement history, target reached, CS_q)."""
    CS_q = CS_of(j_q, eps)
    scale = qf(CS_q) / qf(j_q) ** 2 if j_q != 0 else qf(CS_of(fmpq(1), eps))
    disc = Disc(N, K)
    if j_q == 0:
        wf, rf, it = disc.newton(0.0, np.zeros(disc.L.shape[0]))
    else:
        wc, res, rc = continuation(qf(j_q), scale)
        w0 = prolong(wc, res[0], res[1], N, K)
        wf, rf, it = disc.newton(qf(CS_q), w0)
    if not np.all(np.isfinite(wf)):           # v2: NaN/inf guard before any arb work (fit only; values unchanged)
        raise ValueError('non-finite double fit at j = %s' % qstr(j_q))
    wa, hist, reached = refine(disc, wf, CS_q, prec)
    mx = max(abs(float(x.mid())) for x in wa)
    if mx == 0:
        return [fmpq(0)] * len(wa), rf, hist, reached, CS_q
    e_sc = math.ceil(math.log2(mx))       # dyadic grid relative to the coefficient scale (K1: |W| ~ j^2)
    with _prec(prec):
        wq = [to_dyadic(x, prec - 32 - e_sc) for x in wa]
    return wq, rf, hist, reached, CS_q


def _log_default(msg):
    print(msg, flush=True)


def certify_point(j, N=RES_A[0], K=RES_A[1], prec=RES_A[2], n_pieces=8, adaptive=ADAPTIVE_DEFAULT, eps=fmpq(0),
                  split_ratio=SPLIT_RATIO, max_pieces=MAX_PIECES, keep=False, log=_log_default):
    """certify delta(j) (mu = 1, j = 4J; j an int, Fraction or fmpq, exact). Returns a record (strings and plain
    values only); with keep=True returns (record, (lo, hi) fmpq, certify-dict, wq, CS_q)."""
    j_q = j if isinstance(j, fmpq) else fmpq(int(Fraction(j).numerator), int(Fraction(j).denominator))
    t0 = time.time()
    wq, rf, hist, reached, CS_q = solve_point(j_q, N, K, prec, eps)
    c = certify(wq, N, K, CS_q, n_pieces=n_pieces, prec=prec, adaptive=adaptive, split_ratio=split_ratio,
                max_pieces=max_pieces)
    a = A_MAP
    lo = 2 * a * (c['W0'] - c['eta'])
    hi = 2 * a * (c['W0'] + c['eta'])
    n1 = N + 1
    sc = max([abs(qf(x)) for x in wq[:n1]] + [1e-300])
    tail_p = max(abs(qf(x)) for x in wq[n1 - 4:n1]) / sc
    tail_q = max(abs(qf(wq[k * n1 + n])) for k in range(1, K + 1) for n in range(n1 - 4, n1)) / sc if K else 0.0
    rec = dict(j=qstr(j_q), source_eps=qstr(eps), N=N, K=K, prec_bits=prec, n_pieces_start=n_pieces,
               adaptive=bool(adaptive), split_ratio=qstr(split_ratio), n_pieces_used=c['n_pieces_used'],
               residual_degrees=list(c['residual_degrees']),
               rel_tail_P_last4_approx="%.2e" % tail_p, rel_tail_q_last4_approx="%.2e" % tail_q,   # diagnostics, not bounds
               fine_newton_residual="%.2e" % rf, refine_history=["%.1e" % h for h in hist],
               refine_target_reached=reached,
               sup_Gh_N_le=sci_up(c['sup_Gh_le']), sup_Gh_N_le_exact=qstr(c['sup_Gh_le']),
               eta_exact=qstr(c['eta']), W0_exact=qstr(c['W0']),
               E_N_min_ge=sci_down(c['E_N_min_lower'], 6) if c['E_N_min_lower'] is not None else None,
               E_minus_min_ge=sci_down(c['E_minus_lower'], 6) if c['E_minus_lower'] is not None else None,
               E_minus_positive=c['E_minus_positive'], strict=c['strict'], status=c['status'],
               certified=c['certified'], fail_reason=c['fail_reason'], pieces=c['pieces'],
               delta_lo=None, delta_hi=None, delta_lo_exact=None, delta_hi_exact=None, rel_width_le=None,
               gauss_law_delta_double_CHECK=repr(float(gauss_law_delta_double(wq, N, K, qf(CS_q)))),
               secs_cert=c['secs'], secs_total=round(time.time() - t0, 1),
               engine=dict(module='gravity_engine/bowen_york_sqrt.py', module_version=MODULE_VERSION, seal=SEAL_ID))
    if c['certified']:                    # an interval is reported only when the certificate holds
        rec.update(delta_lo=dec_floor(lo, 34), delta_hi=dec_ceil(hi, 34), delta_lo_exact=qstr(lo),
                   delta_hi_exact=qstr(hi), rel_width_le=sci_up((hi - lo) / lo, 3) if lo > 0 else None)
        shown = "delta in [%s, %s] rel.w<=%s" % (dec_floor(lo, 24), dec_ceil(hi, 24), rec['rel_width_le'])
    else:
        shown = "no interval (%s)" % (c['fail_reason'] or 'not strict / E_- not positive / trivial')
    if log:
        log("bowen_york_sqrt j=%s (N,K,prec)=(%d,%d,%d) boxes=%d%s: sup|Gh|<=%s eta=%s E_N>=%s | %s %s (%.1fs)"
            % (qstr(j_q), N, K, prec, c['n_pieces_used'], ' adaptive' if adaptive else '', rec['sup_Gh_N_le'],
               sci_up(c['eta']), rec['E_N_min_ge'], shown, c['status'], rec['secs_total']))
    if keep:
        return rec, (lo, hi), c, wq, CS_q
    return rec


# ===================================================================================== 6. identity controls
def identity_controls():
    """exact/symbolic identities behind the method (rep4 controls (a)-(g)); every flag must be True."""
    import random
    import sympy as sp
    out = {}
    rnd = random.Random(4)
    with _prec(192):
        # (a) DCT recovers known local Chebyshev coefficients (degree < node count)
        Nn, Mn = 12, 8
        c = [[fmpq(rnd.randint(-50, 50), rnd.randint(1, 9)) for l in range(6)] for k in range(10)]
        pi = arb.pi()
        vals = []
        for i in range(Nn):
            th = pi * fmpq(2 * i + 1, 2 * Nn)
            for m in range(Mn):
                ph = pi * fmpq(2 * m + 1, 2 * Mn)
                vals.append(sum((arb(c[k][l]) * (th * k).cos() * (ph * l).cos() for k in range(10) for l in range(6)),
                                arb(0)))
        C = cheb_coeffs(arb_mat(Nn, Mn, vals), Nn, Mn)
        out['dct_recovers_known_coefficients'] = bool(all(C[k, l].contains(c[k][l] if (k < 10 and l < 6) else 0)
                                                          for k in range(Nn) for l in range(Mn)))
        # (b) angular operator: A_k(c^2) = [(1-c^2) g']' with g = tau_k(c^2) = T_2k(c)
        ok = True
        Cc = fmpq_poly([0, 1])
        C2 = Cc * Cc
        T, A = tau_polys(8)
        Tc = [fmpq_poly([1]), Cc]
        for n in range(1, 16):
            Tc.append(2 * Cc * Tc[n] - Tc[n - 1])
        for k in range(1, 9):
            g = T[k](C2)
            ok = ok and (((1 - C2) * g.derivative()).derivative() == A[k](C2)) and (g == Tc[2 * k])
        out['angular_operator_identity'] = bool(ok)
        # (c) radial operator (sympy): (r^2 u_r)_r = (1-Y) Y^3 LR[W] for u = Y^2 W(Y); and Lin of Y^2 is -3/2
        r, Ys = sp.symbols('r Y', positive=True)
        a = sp.Rational(1, 2)
        Wp = sp.Rational(3, 7) - 2 * Ys + sp.Rational(5, 3) * Ys ** 3 + Ys ** 5 / 11
        Yr = 1 / (1 + sp.sqrt(r / a))
        u = (Ys ** 2 * Wp).subs(Ys, Yr)
        lhs = sp.diff(r ** 2 * sp.diff(u, r), r)
        LRW = -(6 * Wp + (6 * Ys - 3) * sp.diff(Wp, Ys) - Ys * (1 - Ys) * sp.diff(Wp, Ys, 2)) / 4
        rhs = ((1 - Ys) * Ys ** 3 * LRW).subs(Ys, Yr)
        good = True
        for sv in (sp.Rational(1, 3), sp.Rational(5, 2), sp.Rational(7, 9), sp.Integer(4)):
            good = good and sp.simplify((lhs - rhs).subs(r, a * sv ** 2)) == 0
        g2 = Yr ** 2
        lg = sp.diff(r ** 2 * sp.diff(g2, r), r)
        for sv in (sp.Rational(1, 3), sp.Rational(5, 2)):
            good = good and sp.simplify((lg + sp.Rational(3, 2) * (1 - Yr) * Yr ** 3).subs(r, a * sv ** 2)) == 0
        out['radial_operator_symbolic'] = bool(good)
        # (d) source transformation and r Y^2 = a (1-Y)^2, exact rationals
        F = Fraction
        aF = F(1, 2)
        good = True
        for (Y, v, W, J) in ((F(1, 3), F(1, 5), F(2, 7), F(5, 2)), (F(4, 5), F(0), F(1, 9), F(2500)),
                             (F(1, 50), F(1, 2), F(3), F(750)), (F(9, 10), F(3, 4), F(-1, 5), F(1, 3))):
            rr = aF * ((1 - Y) / Y) ** 2
            Phi = 1 + 1 / (2 * rr) + Y * Y * W
            E = 2 * aF * (1 - Y) ** 2 * (1 + Y * Y * W) + Y * Y
            jj = 4 * J
            CS = 18 * aF ** 3 * jj * jj
            good = good and (E == 2 * aF * (1 - Y) ** 2 * Phi)
            good = good and (rr * rr * F(9, 4) * J * J * (1 - v) / (rr ** 6 * Phi ** 7)
                             == (1 - Y) * Y ** 3 * CS * (1 - v) * Y ** 5 * (1 - Y) ** 5 / E ** 7)
            good = good and (rr * Y * Y == aF * (1 - Y) ** 2)
        out['source_transformation_and_rY2'] = bool(good)
        # (e) Lin[1] = -3/2 through the exact evaluator
        N, K = 6, 2
        wq = [fmpq(0)] * ((N + 1) * (K + 1))
        wq[0] = fmpq(1)
        Ex = ExactEval(wq, N, K)
        out['Lin_of_1_is_minus_3_2'] = bool(all(Ex.at(fmpq(i, 7), fmpq(l, 3))[1] == fmpq(-3, 2)
                                                for i in range(8) for l in range(4)))
        # (f) arb field evaluator contains the exact rational evaluation
        N, K = 24, 3
        wq = [fmpq(rnd.randint(-2 ** 20, 2 ** 20), 2 ** 24) for _ in range((N + 1) * (K + 1))]
        Ex = ExactEval(wq, N, K)
        An = Ansatz(wq, N, K)
        good = True
        for (y, v) in ((fmpq(1, 3), fmpq(2, 7)), (fmpq(19, 20), fmpq(1, 2)), (fmpq(1, 100), fmpq(1))):
            we, Le = Ex.at(y, v)
            W, LW = An.F.grid([arb(y)], [arb(v)])
            good = good and W[0, 0].contains(we) and LW[0, 0].contains(Le)
        out['arb_evaluator_contains_exact'] = bool(good)
        # (g) double operator agrees with the arb evaluator
        d = Disc(N, K)
        wf = np.array([qf(x) for x in wq])
        LWd = d.L @ wf
        mx = 0.0
        for idx in range(0, len(LWd), 7):
            W, LW = An.F.grid([arb(d.Y[idx])], [arb(d.V[idx])])
            mx = max(mx, abs(float(LW[0, 0].mid()) - LWd[idx]))
        out['double_operator_matches_arb'] = bool(mx < 1e-9)
    out['all_ok'] = all(out.values())
    return out


# ===================================================================================== 7. tangent certificate (v2)
# DESIGN-P2d sec. 2 (Lemma TL, LEMMA-REP4.md sec. 6), sec. 4.3 items 3-4, as amended by AMEND-P2d-1 (REFEREE-TL fixes) and
# AMEND-P2d-2 (REP4's fixed pieces by default). psi_+ = Y^2 Z_+, Z_+ = Z_N + eta_Z, must satisfy (TL-H) <=> H[Z_+] <= 0 on
# the square, with (REP4 variables, mu1 = 1, W_- = W_N - eta the certified sub-solution's ansatz, E_- = E(W_-))
#     H[Z] = Lin[Z] - 7 Y^2 src_- (1 + 2a(1-Y)^2 Z) / E_-,   src_- = CS (1-v) Y^5 (1-Y)^5 / E_-^7,
#     E_-^8 H[Z] = Q_L := E_-^8 Lin[Z] - 7 CS (1-v) Y^7 (1-Y)^5 (1 + 2a(1-Y)^2 Z),   degree (9N+50, 9K) for K >= 1.
# TL6: Lin[1] = -3/2 and the potential term only drops for a constant shift eta_Z >= 0, so
#     H[Z_N + eta_Z] <= sup|H[Z_N]| - (3/2) eta_Z < 0   when (3/2) eta_Z > sup|H[Z_N]| >= sup|Q_L| / inf E_-^8 per piece.
# The bound uses inf E_- (from the W_- ansatz), never E_N (REFEREE-TL sec. 4 (b)). Then (TL4)
#     2L_hi = 1 + 2a (Z_N(0) + eta_Z)   exactly, and for mu' in (0, 1]:  M(mu') >= 1 + delta_lo - (1 - mu') 2L_hi.
TANGENT_FACTOR = 7          # the 7 of d/dPhi Phi^-7; any other value is for the N-7 control only (wrong tangent problem)


def _need(cond, msg):
    """structural assert that survives python -O (AMEND-P2d-1 items 2-4)."""
    if not cond:
        raise AssertionError(msg)


def tangent_degrees(N, K):
    """Y- and v-degree of Q_L: E_-^8 (8N+48, 8K) times Lin[Z] (N+2, K) gives (9N+50, 9K); the source term
    Y^7 (1-Y)^5 (1-v)(1 + 2a(1-Y)^2 Z) has (N+16, K+1). K >= 1 is asserted (AMEND-P2d-1 item 2, REFEREE-TL sec. 4 (iv));
    the v-degree is written max(9K, K+1) so that the count stays right even if the assert were ever relaxed (FIX F)."""
    _need(N >= 0 and K >= 1, 'tangent degree count (9N+50, 9K) needs N >= 0 and K >= 1 (at K = 0 the v-degree is 1)')
    return max(9 * N + 50, N + 16), max(9 * K, K + 1)


def tangent_operator(disc, w, CS, factor=TANGENT_FACTOR):
    """double collocation of the tangent problem at the ansatz w (W_-): A z = b with A = L - (f src/E 2a(1-Y)^2 Y^2) U,
    b = f src Y^2 / E. For f = 7, A is REP4's Newton Jacobian disc.jac(w, CS) (DESIGN-P2d TL7). Fit only."""
    R, E, src = disc.F(w, CS)
    A = disc.L - (factor * src / E * disc.c2 * disc.Y ** 2)[:, None] * disc.U
    return A, factor * src * disc.Y ** 2 / E


def refine_tangent(disc, wm, zf, CS_q, prec, factor=TANGENT_FACTOR, maxit=25):
    """arb iterative refinement of the double tangent fit (REQUIRED, DESIGN-P2d sec. 4.3 item 3): residual H[Z] in arb at
    the collocation nodes with the exact W_- (wm, fmpq), correction from the double LU of the tangent operator, midpoints
    kept. Stops at 2^-(prec-40) relative to max|z| or on a plateau. Fit only: the proof does not trust it."""
    N, K = disc.N, disc.K
    with _prec(prec):
        ys = arb_nodes01(N + 1)
        vs = arb_nodes01(K + 1)
        cs = arb(CS_q)
        a2 = 2 * arb(A_MAP)
        Wn, _ = Ansatz(wm, N, K).F.grid(ys, vs)
        A, b = tangent_operator(disc, np.array([qf(x) for x in wm]), qf(CS_q), factor)
        lu = lu_factor(A)
        z = [arb(float(x)) for x in zf]
        hist = []
        target = 2.0 ** (-(prec - 40))
        for it in range(maxit):
            p, q = split(z, N, K)
            Zn, LZn = Fields(p, q).grid(ys, vs)
            R = []
            for i in range(N + 1):
                y = ys[i]
                om = 1 - y
                c2 = a2 * om * om
                y2 = y * y
                s0 = factor * cs * y2 * y2 * y2 * y * om ** 5
                for jv in range(K + 1):
                    E = c2 * (1 + y2 * Wn[i, jv]) + y2
                    E2 = E * E
                    E4 = E2 * E2
                    R.append(LZn[i, jv] - s0 * (1 - vs[jv]) * (1 + c2 * Zn[i, jv]) / (E4 * E4))
            Rf = np.array([float(r.mid()) for r in R])
            scale = max(1e-300, max(abs(float(x.mid())) for x in z))
            res = float(np.max(np.abs(Rf))) / scale
            hist.append(res)
            if res < target or (len(hist) > 3 and res > 0.5 * hist[-2]):
                break
            dz = lu_solve(lu, Rf)
            z = [(z[i] - arb(float(dz[i]))).mid() for i in range(len(z))]
        return z, hist, bool(hist[-1] < target)


def _puncture_value_independent_of_v(zq, N, K):
    """exact check of the class structure (LEMMA-REP4 sec. 6 note): Z(1, v) = P_Z(1) for v = 0, 1/2, 1 (the angular part
    carries Y(1-Y)). Structural (true for every Ansatz); kept as an executable statement of AMEND-P2d-1 item 2."""
    Ex = ExactEval(zq, N, K)
    vals = [Ex.at(fmpq(1), v)[0] for v in (fmpq(0), fmpq(1, 2), fmpq(1))]
    return bool(vals[0] == vals[1] == vals[2] == cheb_eval_q(split(zq, N, K)[0], fmpq(1)))


def tangent_piece_bound(ansW, ansZ, ylo, yhi, CS_q, factor=TANGENT_FACTOR, vlo=fmpq(0), vhi=fmpq(1)):
    """rigorous sup|H[Z_N]| <= sup|Q_L| / inf E_-^8 on [ylo, yhi] x [vlo, vhi]; ansW is the W_- ansatz, ansZ the Z_N
    ansatz. Q_L is sampled on (dY+1) x (dV+1) Chebyshev-Gauss nodes (exact DCT for its degree), sup|Q_L| <= sum|c|,
    inf E_- by E_lower on the W_- ansatz. Returns ((bound, E_-low) or None, info)."""
    N, K = ansZ.N, ansZ.K
    dY, dV = tangent_degrees(N, K)
    NQ = dY + 1
    MQ = dV + 1
    ys = nodes_on(ylo, yhi, NQ)
    vs = nodes_on(vlo, vhi, MQ)
    Wg, _ = ansW.F.grid(ys, vs)
    Zg, LZ = ansZ.F.grid(ys, vs)
    cs = arb(CS_q)
    a2 = 2 * arb(A_MAP)
    flat = []
    for i in range(NQ):
        y = ys[i]
        om = 1 - y
        y2 = y * y
        c2 = a2 * om * om
        s7 = factor * cs * y2 * y2 * y2 * y * om ** 5
        for jv in range(MQ):
            E = c2 * (1 + y2 * Wg[i, jv]) + y2
            E2 = E * E
            E4 = E2 * E2
            flat.append(E4 * E4 * LZ[i, jv] - s7 * (1 - vs[jv]) * (1 + c2 * Zg[i, jv]))
    SQ = arb_upper_q(abs_sum(cheb_coeffs(arb_mat(NQ, MQ, flat), NQ, MQ)))
    Elow, nb = E_lower(ansW, ylo, yhi, vlo=vlo, vhi=vhi)
    base = dict(ylo=str(ylo), yhi=str(yhi), vlo=str(vlo), vhi=str(vhi), sup_QL_le=sci_up(SQ))
    if Elow is None or Elow <= 0:
        return None, dict(base, E_minus_low="FAIL")
    bound = SQ / Elow ** 8
    return (bound, Elow), dict(base, E_minus_low_ge=sci_down(Elow, 6), E_boxes=nb, bound_le=sci_up(bound))


def _tangent_boxes(n_pieces):
    Ed = make_pieces(n_pieces)
    _need(Ed[0] == 0 and Ed[-1] == 1 and all(Ed[i] < Ed[i + 1] for i in range(len(Ed) - 1)), 'Y-pieces not ordered')
    boxes = [(Ed[i], Ed[i + 1], fmpq(0), fmpq(1)) for i in range(len(Ed) - 1)]
    check_cover(boxes)                     # AMEND-P2d-1 item 4: exact cover of [0,1]^2, no gap, no overlap
    return boxes


def certify_tangent(wq, N, K, CS_q, eta, prec=256, n_pieces=8, eta_Z=None, factor=TANGENT_FACTOR, psi_zero=False,
                    z_override=None, zq_override=None):
    """tangent certificate (TL-H) at one node, for the point certificate's exact ansatz wq (W_N) and margin eta.
    W_- = W_N - eta (constant Chebyshev coefficient; W = P + ..., T_0 = 1, no c0/2 convention).
    Fit: double tangent solve at W_-, arb refinement, dyadic rounding on a grid scaled to max|z| (AMEND-P2d-1 item 5).
    Certificate: per piece (REP4's 14 fixed Y-pieces x [0,1], AMEND-P2d-2) sup|H[Z_N]| <= sup|Q_L| / inf E_-^8;
    eta_Z = up((2/3) sup (1 + 2^-20)); strict iff (3/2) eta_Z > sup; 2L_hi = 1 + 2a (Z_N(0) + eta_Z) exactly.
    psi_zero=True: psi_+ = 0, for which TL-H holds identically (H[0] = -f Y^2 src_-/E_- <= 0 once E_- > 0, CS >= 0):
    only E_- > 0 per piece is certified and 2L_hi = 1 exactly (REFEREE-TL sec. 1; P2c's plain cell and X1).
    Controls only: eta_Z (forced margin), factor != 7 (N-7: wrong tangent problem, fit and bound alike), z_override
    (double fit replaced; non-finite or |z| > 1e100 gives status INVALID), zq_override (exact Z_N replaced).
    The caller must pass the eta of a CERTIFIED point certificate of wq (TL needs u_- = Y^2 W_- a certified sub-solution).
    Returns a dict of exact fmpq values; status 'CERTIFIED' | 'FAIL' | 'REFUSED/TRIVIAL' | 'INVALID'."""
    tm = time.time()
    dY, dV = tangent_degrees(N, K)                       # asserts K >= 1
    _need(len(wq) == (N + 1) * (K + 1), 'W ansatz shape')
    _need(eta >= 0 and CS_q >= 0, 'eta >= 0 and CS >= 0 required')
    base = dict(N=N, K=K, prec_bits=prec, factor=factor, tangent_degrees=(dY, dV), psi_zero=bool(psi_zero))
    with _prec(prec):
        wm = list(wq)
        wm[0] = wm[0] - eta
        ansW = Ansatz(wm, N, K)
        hist, reached = [], None
        if psi_zero:
            zq = [fmpq(0)] * len(wq)
        elif zq_override is not None:
            zq = list(zq_override)
        else:
            disc = Disc(N, K)
            if z_override is None:
                A, b = tangent_operator(disc, np.array([qf(x) for x in wm]), qf(CS_q), factor)
                zf = np.linalg.solve(A, b)
            else:
                zf = np.asarray(z_override, dtype=float)
            if zf.shape != (len(wq),) or not np.all(np.isfinite(zf)) or float(np.max(np.abs(zf))) > 1e100:
                return dict(base, status='INVALID', certified=False, fail_reason='non-finite or huge tangent fit',
                            twoL_hi=None, eta_Z=None, Z0=None, sup_H_le=None, E_minus_min_lower=None, strict=False,
                            pieces=[], n_pieces_used=0, refine_history=[], refine_target_reached=None, zq=None,
                            secs=round(time.time() - tm, 1))
            za, hist, reached = refine_tangent(disc, wm, zf, CS_q, prec, factor)
            mx = max(abs(float(x.mid())) for x in za)
            if mx == 0:
                zq = [fmpq(0)] * len(za)
            else:
                e_sc = math.ceil(math.log2(mx))  # scaled grid (rep4.py:594-596): never an absolute 2^-(prec-40)
                zq = [to_dyadic(x, prec - 32 - e_sc) for x in za]
        ansZ = Ansatz(zq, N, K)                          # the class of LEMMA-REP4 sec. 4 / sec. 6 (shape asserted)
        _need(ansZ.N == ansW.N == N and ansZ.K == ansW.K == K and len(zq) == len(wm), 'N_Z = N and K_Z = K required')
        _need(_puncture_value_independent_of_v(zq, N, K), 'Z(1, v) must equal P_Z(1) (ansatz structure)')
        boxes = _tangent_boxes(n_pieces)
        pieces = []
        ok = True
        why = None
        sup = fmpq(0)
        Emin = None
        for (y0, y1, v0, v1) in boxes:
            if psi_zero:
                El, nb = E_lower(ansW, y0, y1, vlo=v0, vhi=v1)
                info = dict(ylo=str(y0), yhi=str(y1), E_minus_low_ge=(sci_down(El, 6) if El is not None else 'FAIL'))
                bd = None if (El is None or El <= 0) else (fmpq(0), El)
            else:
                bd, info = tangent_piece_bound(ansW, ansZ, y0, y1, CS_q, factor, v0, v1)
            pieces.append(info)
            if bd is None:
                ok = False
                why = why or 'E_- lower bound failed on a box'
                continue
            sup = max(sup, bd[0])
            Emin = bd[1] if Emin is None or bd[1] < Emin else Emin
        if psi_zero:
            eta_Z = fmpq(0)
            strict = ok                                  # H[0] <= 0 exactly once E_- > 0 on every piece
        else:
            if eta_Z is None:
                eta_Z = round_up_dyadic(sup * (1 + fmpq(1, 2 ** 20)) * fmpq(2, 3))
            strict = bool(ok and eta_Z > 0 and fmpq(3, 2) * eta_Z > sup)
        if not psi_zero and ok and sup == 0:
            status = 'REFUSED/TRIVIAL'
        elif strict:
            status = 'CERTIFIED'
        else:
            status = 'FAIL'
            why = why or 'not strict: (3/2) eta_Z <= sup|H|'
        twoL_hi = 1 + 2 * A_MAP * (ansZ.W0 + eta_Z)
    return dict(base, status=status, certified=(status == 'CERTIFIED'), fail_reason=(None if status == 'CERTIFIED' else why),
                sup_H_le=sup, eta_Z=eta_Z, Z0=ansZ.W0, twoL_hi=twoL_hi, E_minus_min_lower=Emin, strict=strict,
                pieces=pieces, n_pieces_used=len(pieces), refine_history=hist, refine_target_reached=reached,
                zq=zq, secs=round(time.time() - tm, 1))


def tangent_record(t):
    """JSON-safe summary of a certify_tangent dict (exact values as "p/q" strings; display rounded outward)."""
    def q(x):
        return None if x is None else qstr(x)
    return dict(status=t['status'], certified=t['certified'], fail_reason=t.get('fail_reason'), N=t['N'], K=t['K'],
                prec_bits=t['prec_bits'], factor=t['factor'], psi_zero=t['psi_zero'],
                tangent_degrees=list(t['tangent_degrees']),
                twoL_hi_exact=q(t['twoL_hi']), eta_Z_exact=q(t['eta_Z']), Z0_exact=q(t['Z0']),
                sup_H_le_exact=q(t['sup_H_le']),
                sup_H_le=(sci_up(t['sup_H_le']) if t['sup_H_le'] is not None else None),
                twoL_hi_up=(sci_up(t['twoL_hi'], 16) if t['twoL_hi'] is not None else None),
                eta_Z_up=(sci_up(t['eta_Z']) if t['eta_Z'] is not None else None),
                E_minus_min_ge=(sci_down(t['E_minus_min_lower'], 6) if t.get('E_minus_min_lower') is not None else None),
                strict=t.get('strict'), n_pieces_used=t.get('n_pieces_used'), pieces=t.get('pieces'),
                refine_history=["%.1e" % h for h in (t.get('refine_history') or [])],
                refine_target_reached=t.get('refine_target_reached'))


# ---- cell and tail bounds from TL (DESIGN-P2d sec. 2, TL4 + TL5; exact rationals; brackets asserted > 0)
def tl_cell_bound(j_a, j_b, s_lo, delta_lo, twoL_hi):
    """chi(j) <= (j_a/4) / (1 + delta_lo - tL (1 - s_lo))^2 for every j in [j_a, j_b], tL = max(2L_hi, 0), where s_lo is
    an exact rational LOWER bound of sqrt(j_a/j_b) (s_lo^2 j_b <= j_a is checked exactly): a smaller s_lo only enlarges the
    bound. Raises AssertionError unless the bracket is > 0 (AMEND-P2d-1 item 3)."""
    _need(0 < j_a < j_b, 'cell needs 0 < j_a < j_b')
    _need(0 < s_lo <= 1 and s_lo * s_lo * j_b <= j_a, 's_lo must be a rational lower bound of sqrt(j_a/j_b)')
    tL = twoL_hi if twoL_hi > 0 else 0 * twoL_hi
    br = 1 + delta_lo - tL * (1 - s_lo)
    _need(br > 0, 'cell bracket 1 + delta_lo - tL (1 - s_lo) must be > 0')
    return (j_a / 4) / (br * br)


def tl_tail_bound(j_a, delta_lo, twoL_hi):
    """chi(j) <= (j_a/4) / (1 + delta_lo - tL)^2 for every j >= j_a (TL4 with mu' -> 0, TL5). Raises AssertionError
    unless the bracket 1 + delta_lo - tL is > 0 (AMEND-P2d-1 item 3, REFEREE-TL sec. 2)."""
    _need(j_a > 0, 'tail needs j_a > 0')
    tL = twoL_hi if twoL_hi > 0 else 0 * twoL_hi
    br = 1 + delta_lo - tL
    _need(br > 0, 'tail bracket 1 + delta_lo - tL must be > 0')
    return (j_a / 4) / (br * br)
