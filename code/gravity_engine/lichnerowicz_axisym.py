"""Exact-rational certificates for axisymmetric Lichnerowicz (Hamiltonian-constraint) equations on the compactified square.

In plain words: give it a puncture-type constraint equation  Delta u + (source) Phi^-7 = 0  in one of the lab's
compactified coordinates, and it returns a PROOF (exact rational arithmetic, no floating point in any inequality) that the
true solution lies between two explicit polynomials, together with an exact rational interval for the mass term.

Provenance (doctrine #121; gravity_engine/__init__.py):
- DESIGN file: <LAB>/BHSLICE/DESIGN-P2.md (node table sec 1; `produces: ... exact-rational s-engine with Bernstein
  positivity`).
- Extracted from the sealed instrument <LAB>/BHSLICE/p2a.py (sha256 feee7d6f273fbfe97581dca0d224c0e387b6cdf31f9f10ae6f
  949543de2dbf04, read-only, never edited). The mathematics and the order of every exact operation are identical to p2a.py;
  only the interfaces are generalised (problem object: exact source numerator, collocation source, coordinate/basis, the
  exponent of D, the scale tau(j)).
- Certificate seal: claim d62bfc0e0b (BAR de758dbbf3, AMEND-P2a-1 9356a34707, instrument seal 91b477f300).
  VERDICT-P2a.md: HOLD (123/123 nodes CERTIFIED on j in [1, 10.144]). REFEREE-P2a.md (code-level referee J): CORRECT
  (R1-R5, Bernstein, cells, controls). VERDICT-REP2.md (blind seat, claim 227cbcafc2): REPLICATED at j = 6 and j = 10.
- CERTIFIED nodes of DESIGN-P2 this module backs:
  R1 exact residual Q (tau*Q trick), R2 abs-Chebyshev bound S, R3 exact Bernstein positivity on the 2-D square,
  R4 margin eta and strict sub/super signs, R5 slope route for the mass term, P1 point certificates,
  K1 monotone cell enclosures (arb, outward), K2 gap-free cover.
- Status: CERTIFIED for the Bowen-York spinning puncture with mu = 1, j in [1, 10.144] (the sealed P2a run); the
  reproduction tests in tests/test_lichnerowicz_axisym.py re-derive j = 2 and j = 10 byte for byte and the three extra
  controls REFEREE-P2a asked for. Any NEW problem object is a new claim: its source numerator, the sign hypotheses
  listed below and its mass normalisation must be checked (and refereed) before a certificate built with it is trusted.

The certificate (DESIGN-P2 nodes; p2a.py docstring):
  coordinate  s = a/(r + a), a = 1/2, c = cos(theta); tau = tau(j) > 0 the source scale; u = tau U.
  equation    Gh[U] = Lh U + src(s, c) D^-P,  D = 1 + (1 - s) tau U  (P = 7 for Lichnerowicz),
              Lh U = (1-s) U_0'' - 2 U_0' + sum_l [(1-s) U_l'' - 2 U_l' - l(l+1) q_l] P_l(c)
              (the divided flat Laplacian in s; referee J re-derived it from spherical coordinates:
               Delta u + tau src-term = (4 s^4/(1-s)) tau Gh[U]).
  ansatz      U = s p(s) + sum_{l=2,4,..,2K} s^2 (1-s) q_l(s) P_l(c), p and q_l Chebyshev in X = 2s - 1.
  R1  Q := D_N^P Lh U_N + src, exact (computed as tau Q = D_N^P Lh u_N + tau src: tau never enters D^P).
  R2  sup |Q| <= S := sum |Chebyshev coefficients of Q| on [0,1]_s x [-1,1]_c (exact rational).
  R4  S > 0 asserted, else REFUSED/TRIVIAL; eta := S (1 + 2^-20)/2 rounded UP to a 64-bit dyadic; U+- = U_N +- eta s.
  R3  V_- := U_-/s > 0 on the closed square by exact tensor Bernstein coefficients (halve the failing s-piece, depth 8;
      a corner coefficient <= 0 is an exact point value and therefore a certified VIOLATION).
      Then D_- >= 1, D_N >= 1, |Gh[U_N]| <= |Q| <= S, and since Lh(s) = -2 and the source term decreases in D:
      Gh[U+] <= S - 2 eta < 0 < -S + 2 eta <= Gh[U-]; the comparison lemma (C1, P1 refereed) gives u_- <= u <= u_+.
  R5  slope route: delta = 2 lim r u = p_u(0) +- tau eta (exact rationals; a = 1/2 makes 2 r s -> 1).
Hypotheses a new problem object must satisfy for R4 to hold as written: src >= 0 on the square, the nonlinearity is
src * D^-P with P > 0 (decreasing in D), and the operator/coordinate are the ones above (Lh(s) = -2 constant). The
module asserts only what it can check exactly (src terms polynomial, packing width); the rest is the caller's claim.

What this module does NOT do: the extreme limit mu = 0 (needs a weighted comparison lemma, DESIGN-P2 Y2), other
spacetimes or other coordinates (a new coordinate needs its own divided operator, Lh(s) and slope normalisation), the
integral route R6 and spin slabs S1 (never built), and it proves nothing about j outside the nodes it is run on except
through the K1 cells, which rest on the Bowen-York-specific monotonicity delta(j)/sqrt(j) nondecreasing (DGC Lemma 3.1
+ scaling).

Determinism: the fit (not part of the proof) is bit-reproducible only with single-thread LAPACK. This module sets the
thread variables before importing numpy; if numpy was imported earlier in the process, `NUMPY_PREIMPORTED` is True
and byte-identical reproduction of stored certificates is not guaranteed (the certificates themselves stay valid).
No pickle anywhere: exact rationals are stored as "p/q" strings.

Engine version 2 (ENGINE_VERSION; CHANGELOG in gravity_engine/README.md): fit-only options. `fit_ops` picks the arb
fit-operator construction: 'legacy' = v1 balls, 'mid' = exact midpoint entries whose radii stay at working precision,
'auto' = legacy for N <= LEGACY_MAX_N, else mid. `cache_ops=False` and `clear_caches()` free the fit operators between
degrees. The proof path is unchanged.
"""
import sys
import os

NUMPY_PREIMPORTED = 'numpy' in sys.modules
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # single-thread LAPACK => deterministic fits
    os.environ[_v] = '1'
import math
import time
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat, fmpq, fmpz_mat, fmpq_poly, fmpq_mat, ctx

SEAL_ID = 'd62bfc0e0b'
SOURCE_INSTRUMENT_SHA256 = 'feee7d6f273fbfe97581dca0d224c0e387b6cdf31f9f10ae6f949543de2dbf04'
DESIGN_FILE = '<LAB>/BHSLICE/DESIGN-P2.md'
CERTIFIED_NODES = ('R1', 'R2', 'R3', 'R4', 'R5', 'P1', 'K1', 'K2')

PREC = 160                # arb refinement precision (bits)
RBITS = 128               # fitted coefficients rounded to multiples of 2^-RBITS (exact dyadics)
COEF_CAP = 1.0e3          # guard: any |fitted coefficient| above this (or non-finite) => INVALID
BERN_DEPTH = 8            # max halvings of an s-piece in the Bernstein positivity check

# engine version 2 (2026-11-14): fit-only changes, see gravity_engine/README.md CHANGELOG. The proof path (exact Q, S,
# eta, Bernstein, point checks, slope route, cells) is byte-identical to version 1.
ENGINE_VERSION = 2
FIT_OPS = 'auto'          # default arb fit-operator construction: 'auto' = 'legacy' for N <= LEGACY_MAX_N, else 'mid'
LEGACY_MAX_N = 100        # 'auto' keeps the v1 ball construction up to here (P2a, N = 40, reproduces byte for byte)


# ================================================================== small exact helpers
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


# ================================================================== exact polynomial building blocks
SP = fmpq_poly([0, 1])                   # s
OMS = fmpq_poly([1, -1])                 # 1 - s
ONE = fmpq_poly([1])
HP = fmpq_poly([0, 0, 1, -1])            # s^2 (1 - s)
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
    """Legendre P_l(c) as exact fmpq_poly, l = 0..L."""
    if not _LEG:
        _LEG.extend([ONE, fmpq_poly([0, 1])])
    while len(_LEG) <= L:
        l = len(_LEG) - 1
        _LEG.append((fmpq(2 * l + 1) * fmpq_poly([0, 1]) * _LEG[l] - fmpq(l) * _LEG[l - 1]) / fmpq(l + 1))
    return _LEG[:L + 1]


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
    """packed poly -> {(i, k): coefficient of s^i c^k}."""
    out = {}
    for idx, v in enumerate(P.coeffs()):
        if v != 0:
            out[(idx // W, idx % W)] = v
    return out


# ================================================================== R2: exact Chebyshev coefficients and abs sums
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


# ================================================================== R3: exact Bernstein positivity on the 2-D square
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


def _bern_setup(terms):
    m = max(g.degree() for f, g in terms)
    m = max(m, 0)
    Vc = [fmpq_poly([0]) for _ in range(m + 1)]
    for f, g in terms:
        for k, gk in enumerate(g.coeffs()):
            if gk != 0:
                Vc[k] += gk * f
    n = max(max(p.degree() for p in Vc), 0)
    return Vc, n, m, c_to_w(m), bern_T(m).transpose(), bern_T(n)


def _bern_piece(setup, a, b):
    Vc, n, m, G, Tc, Ts = setup
    sub = fmpq_poly([a, b - a])
    cols = [Vc[k](sub).coeffs() if Vc[k] != 0 else [] for k in range(m + 1)]
    ent = []
    for i in range(n + 1):
        for k in range(m + 1):
            ent.append(cols[k][i] if i < len(cols[k]) else fmpq(0))
    Mt = fmpq_mat(n + 1, m + 1, ent)
    return Ts * (Mt * G) * Tc


def bernstein_matrix(terms, a=fmpq(0), b=fmpq(1)):
    """Exact tensor Bernstein coefficients B[i,k] of V(s,c) = sum f(s) g(c) on [a,b]_s x [-1,1]_c:
    V(a + (b-a) t, 2w - 1) = sum_{i,k} B[i,k] C(n,i) t^i (1-t)^(n-i) C(m,k) w^k (1-w)^(m-k). Returns (B, n, m)."""
    setup = _bern_setup(terms)
    return _bern_piece(setup, a, b), setup[1], setup[2]


def bernstein_positive(terms, max_depth=BERN_DEPTH):
    """certify V(s,c) = sum f(s) g(c) > 0 on [0,1] x [-1,1]. Corner Bernstein coefficients are exact point values:
    a corner <= 0 is a certified VIOLATION. Otherwise halve the s-piece until all coefficients > 0, or FAIL at depth."""
    setup = _bern_setup(terms)
    n, m = setup[1], setup[2]
    queue = [(fmpq(0), fmpq(1), 0)]
    pieces, gmin = 0, None
    while queue:
        a, b, d = queue.pop(0)
        B = _bern_piece(setup, a, b)
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


# ================================================================== coordinate + basis
def _cheb_all(x, N, one, mid=False):
    """T_n, T_n', T_n'' at x (n = 0..N) by the three-term recurrences. mid=True (arb only): keep only the midpoint after
    every step, so interval radii never grow (the ball recurrence grows radii roughly like 3^N; engine version 2)."""
    T = [one, x]
    d1 = [0 * one, one]
    d2 = [0 * one, 0 * one]
    for n in range(1, N):
        t = 2 * x * T[n] - T[n - 1]
        e1 = 2 * T[n] + 2 * x * d1[n] - d1[n - 1]
        e2 = 4 * d1[n] + 2 * x * d2[n] - d2[n - 1]
        if mid:
            t, e1, e2 = t.mid(), e1.mid(), e2.mid()
        T.append(t)
        d1.append(e1)
        d2.append(e2)
    return T[:N + 1], d1[:N + 1], d2[:N + 1]


def _leg_all(c, L, one, mid=False):
    P = [one, c]
    for l in range(1, L):
        v = ((2 * l + 1) * c * P[l] - l * P[l - 1]) / (l + 1)
        P.append(v.mid() if mid else v)
    return P[:L + 1]


class SCompactBasis:
    """P1's compactified coordinate s = a/(r + a), a = 1/2, with the even-Legendre axisymmetric basis
    U = s p(s) + sum_{l=2..2K} s^2 (1-s) q_l(s) P_l(c) (p, q_l Chebyshev in X = 2s - 1).
    Divided operator Lh (BUILD-P1 sec 1; re-derived by REFEREE-P2a item 1): Rop(f) = (1-s) f'' - 2 f', angular -l(l+1) q_l.
    D = 1 + (1-s) u (so Phi = D/(1-s) and Phi >= 1 + 1/(2r) when u >= 0); Lh(s) = -2; delta = 2 lim r u = p_u(0)."""
    name = 's = a/(r+a), a = 1/2; even Legendre in c = cos(theta)'
    a = Fraction(1, 2)
    lh_of_s = -2                     # Rop(s): used in the strict-sign step R4 and in exact point checks
    D_factor = OMS                   # D = 1 + D_factor * u
    radial0 = SP                     # U_0 = s p(s)
    radial_l = HP                    # U_l = s^2 (1-s) q_l(s)
    radial_l_over_s = SP * OMS       # U_l / s

    @staticmethod
    def Rop(f):
        return OMS * f.derivative().derivative() - 2 * f.derivative()

    @staticmethod
    def angular_eig(l):
        return l * (l + 1)

    def build_ops(self, N, K, mode, src_eval):
        """Square collocation: s_i = (1 + cos(pi i/N))/2 (i = 0..N), c_m = cos(pi m/(2K)) (m = 0..K).
        Unknown vector: block 0 = p_n (n = 0..N), block k = q_{2k,n}. Returns rows of Lh and of U, (1-s), src.
        mode 'f': float64; 'a': arb balls (legacy, engine v1: radii grow ~3^N, harmful for N >~ 110 at 160 bits);
        'm': arb midpoints (engine v2): collocation points are rounded to exact dyadics (their midpoints), every
        recurrence step keeps only its midpoint and every returned entry is an exact arb (radius 0). The 'm' operator
        is the collocation operator at those dyadic points up to rounding at working precision; it only steers the
        fit, never the proof."""
        mid = (mode == 'm')
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
        if mid:
            ss = [v.mid() for v in ss]
            cs = [v.mid() for v in cs]
        legs = [_leg_all(c, 2 * K, one, mid) for c in cs]
        rows_L, rows_U, om, src = [], [], [], []
        for s in ss:
            X = 2 * s - one
            T, d1, d2 = _cheb_all(X, N, one, mid)
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
                sv = src_eval(s, c, one, o)
                if mid:
                    rl = [v.mid() for v in rl]
                    ru = [v.mid() for v in ru]
                    sv = sv.mid()
                rows_L.append(rl)
                rows_U.append(ru)
                om.append(o.mid() if mid else o)
                src.append(sv)
        return rows_L, rows_U, om, src

    def build_terms(self, pc, qcs):
        """U (or u) and Lh U as lists of (poly in s, poly in c); p(s) and q_l(s); p(0)."""
        N = len(pc) - 1
        Ty = cheb_s(N)
        Pl = leg_c(2 * len(qcs))
        pp = fmpq_poly([0])
        for n in range(N + 1):
            if pc[n] != 0:
                pp += pc[n] * Ty[n]
        f0 = self.radial0 * pp
        Ut = [(f0, ONE)]
        Lt = [(self.Rop(f0), ONE)]
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
            g = self.radial_l * q
            Ut.append((g, Pl[l]))
            Lt.append((self.Rop(g) - fmpq(self.angular_eig(l)) * q, Pl[l]))
        p0 = fmpq(0)
        for n in range(N + 1):
            p0 += pc[n] if n % 2 == 0 else -pc[n]
        return Ut, Lt, pp, qs, p0


# ================================================================== problem objects
class AxisymLichnerowiczProblem:
    """Base problem: Gh[U] = Lh U + src D^-power on the compactified square, u = tau(j) U.
    Subclasses give: src_terms (exact source numerator as [(f(s), g(c))], U-scale), src_eval (the same source for the
    collocation fit, float or arb; it only steers the fit, the proof uses src_terms), tau(j), degrees(j), and optional
    hooks observables(...) / brackets(...) / certify_zero(...). See the module docstring for the sign hypotheses."""
    name = 'abstract'
    power = 7
    basis = SCompactBasis()
    point_s = [fmpq(i, 8) for i in range(9)]
    point_c = [fmpq(0), fmpq(1, 2), fmpq(1)]
    src_terms = []

    # engine v2 options (class defaults; set per problem object):
    #   fit_ops   None -> module FIT_OPS ('auto'); 'legacy' (v1 arb balls), 'mid' (v2 arb midpoints), 'auto'
    #             ('legacy' for N <= LEGACY_MAX_N, else 'mid'). Fit only: the proof path never reads these operators.
    #   cache_ops True (v1 behaviour): fit operators stay cached per (N, K) on the problem object.
    #             False: only the current (N, K) is held during a fit, and fit_u frees it before the exact proof part.
    fit_ops = None
    cache_ops = True

    def __init__(self, fit_ops=None, cache_ops=True):
        self._fops, self._aops = {}, {}
        if fit_ops is not None:
            self.fit_ops = fit_ops
        if not cache_ops:
            self.cache_ops = False

    def clear_caches(self):
        """free all cached fit operators (float and arb) of this problem object."""
        self._fops.clear()
        self._aops.clear()

    def fit_ops_mode(self, N):
        """resolved arb fit-operator construction for radial degree N: 'legacy' or 'mid'."""
        m = self.fit_ops if self.fit_ops is not None else FIT_OPS
        if m == 'auto':
            m = 'legacy' if N <= LEGACY_MAX_N else 'mid'
        if m not in ('legacy', 'mid'):
            raise ValueError("fit_ops must be 'legacy', 'mid' or 'auto', got %r" % (m,))
        return m

    def src_eval(self, s, c, one, o):
        raise NotImplementedError

    def tau(self, j):
        raise NotImplementedError

    def degrees(self, j):
        raise NotImplementedError

    def width(self, K):
        """packing width: c-degree of D^P Lh u is 2K (P + 1); +3 keeps every product carry-free (16K + 3 at P = 7)."""
        return 2 * K * (self.power + 1) + 3

    def src_at(self, s0, c0):
        return sum((f(s0) * g(c0) for f, g in self.src_terms), fmpq(0))

    def observables(self, jq, tau, dlo, dhi):
        return {}

    def brackets(self, j, jq, tau, dlo, dhi, chi):
        return None

    def certify_zero(self, N, K, keep):
        return dict(j='0/1', status='REFUSED/TRIVIAL', why='j = 0 has no source'), None

    # ---- fit operators (cached per problem object; cache_ops=False keeps only the current degree)
    def float_ops(self, N, K):
        if (N, K) not in self._fops:
            if not self.cache_ops:
                self._fops.clear()
            L, U, om, src = self.basis.build_ops(N, K, 'f', self.src_eval)
            L, U = np.array(L), np.array(U)
            self._fops[(N, K)] = (L, U, np.array(om), np.array(src))
        return self._fops[(N, K)]

    def arb_ops(self, N, K):
        """(Lh, U) as arb_mat and (1-s), src lists at the collocation points, at PREC bits. The construction is
        fit_ops_mode(N); the cache key stays (N, K, PREC) as in engine v1 and the entry remembers its construction."""
        mode = self.fit_ops_mode(N)
        key = (N, K, PREC)
        if key in self._aops and self._aops[key][4] != mode:
            del self._aops[key]
        if key not in self._aops:
            if not self.cache_ops:
                self._aops.clear()
            ctx.prec = PREC
            L, U, om, src = self.basis.build_ops(N, K, 'a' if mode == 'legacy' else 'm', self.src_eval)
            La = arb_mat(L)
            del L
            Ua = arb_mat(U)
            del U
            self._aops[key] = (La, Ua, om, src, mode)
        return self._aops[key][:4]


# ================================================================== fit (never a proof)
def newton_float(problem, N, K, tau, w0, itmax=40):
    L, U, om, src = problem.float_ops(N, K)
    P = problem.power
    w = w0.copy()
    for it in range(itmax):
        D = 1 + om * (U @ w)
        F = L @ w + tau * src * D ** -P
        J = L - (P * tau * src * om * D ** -(P + 1))[:, None] * U
        dw = np.linalg.solve(J, F)
        w = w - dw
        if np.max(np.abs(dw)) < 1e-15 * max(1.0, np.max(np.abs(w))):
            break
    D = 1 + om * (U @ w)
    return w, float(np.max(np.abs(L @ w + tau * src * D ** -P)))


def fit_u(problem, j, N, K):
    """u = tau U (u-scale) coefficients: double Newton with a fixed continuation in j (1, 2, 4, ..., j),
    then iterative refinement with residuals in arb at PREC bits. Deterministic given (j, N, K).
    (The continuation feeds float(j)^2 to the double solver as in p2a.py; tau(j) = j^2 is assumed there.)"""
    jf = float(j)
    P = problem.power
    n = (N + 1) * (K + 1)
    w = np.zeros(n)
    sched = []
    x = min(1.0, jf)
    while x < jf:
        sched.append(x)
        x *= 2
    sched.append(jf)
    for jj in sched:
        w, resf = newton_float(problem, N, K, jj * jj, w)
    if not np.all(np.isfinite(w)):
        if not problem.cache_ops:
            problem.clear_caches()
        return None, dict(float_residual=float('nan'), arb_hist=[])
    ctx.prec = PREC
    La, Ua, om, src = problem.arb_ops(N, K)
    Lf, Uf, omf, srcf = problem.float_ops(N, K)
    tau_q = problem.tau(j)
    tau_a = arb(tau_q)
    tf = jf * jf
    wa = arb_mat([[arb(float(v))] for v in w])
    hist = []
    for it in range(14):
        Lw = La * wa
        Uw = Ua * wa
        F = [Lw[i, 0] + tau_a * src[i] / (1 + om[i] * Uw[i, 0]) ** P for i in range(n)]
        Ff = np.array([float(f.mid()) for f in F])
        res = float(np.max(np.abs(Ff)))
        hist.append(res)
        if res < 2.0 ** (-(PREC - 40)):
            break
        wd = np.array([float(wa[i, 0].mid()) for i in range(n)])
        D = 1 + omf * (Uf @ wd)
        J = Lf - (P * tf * srcf * omf * D ** -(P + 1))[:, None] * Uf
        dw = np.linalg.solve(J, Ff)
        wa = arb_mat([[(wa[i, 0] - arb(float(dw[i]))).mid()] for i in range(n)])
    out = [wa[i, 0].mid() for i in range(n)]
    diag = dict(float_residual=resf, arb_hist=hist, fit_ops=problem.fit_ops_mode(N), prec_bits=PREC,
                stop=2.0 ** (-(PREC - 40)), reached_stop=bool(hist and hist[-1] < 2.0 ** (-(PREC - 40))))
    if not problem.cache_ops:
        del La, Ua, om, src, Lf, Uf, omf, srcf
        problem.clear_caches()
    return out, diag


def split_coefs(vec, N, K):
    pc = list(vec[:N + 1])
    qcs = [list(vec[(N + 1) * k:(N + 1) * (k + 1)]) for k in range(1, K + 1)]
    return pc, qcs


# ================================================================== R1 + R2 on one exact ansatz
class ExactFn:
    """U_N given by Chebyshev/Legendre coefficients. For tau > 0 the coefficients are in u-scale (u = tau U, exact
    dyadics from the fit); for tau == 0 they are the U-scale coefficients themselves (linear limit problem)."""

    def __init__(self, problem, pc, qcs, tau):
        self.problem = problem
        self.basis = problem.basis
        self.P = problem.power
        self.pc, self.qcs, self.tau = pc, qcs, tau
        self.K = len(qcs)
        self.W = problem.width(self.K)
        self.Ut, self.Lt, self.pp, self.qs, self.p0 = self.basis.build_terms(pc, qcs)

    def p0_U(self):
        return self.p0 / self.tau if self.tau != 0 else self.p0

    def residual_S(self, want_Q=False):
        """R1 + R2: exact S = sum |Chebyshev coefficients of Q|, Q = D^P Lh U + src (U-scale).
        With want_Q also returns the packed exact Q (width self.W)."""
        W = self.W
        srcP = pack(self.problem.src_terms, W)
        CB, scB = chebmat(srcP, W)
        if self.tau == 0:
            Q = pack(self.Lt, W) + srcP
            C, sc = chebmat(Q, W)
            S = abs_sum_combined(C, sc, None, 0, 0)
            return (S, Q) if want_Q else (S, None)
        DP = pack([(ONE, ONE)] + [(self.basis.D_factor * f, g) for f, g in self.Ut], W)
        A = DP ** self.P * pack(self.Lt, W)            # = tau (Q - src)
        CA, scA = chebmat(A, W)
        S_tau = abs_sum_combined(CA, scA, CB, scB, self.tau)   # = sum |cheb(tau Q)|
        S = S_tau / self.tau
        if want_Q:
            return S, A / self.tau + srcP
        return S, None

    def tauGh_at(self, s0, c0, e_signed):
        """exact point value of tau * Gh[U_N + e_signed s] (same sign as Gh); None if D <= 0 there."""
        lhs = self.basis.lh_of_s
        u = fmpq(0)
        Lu = fmpq(0)
        for f, g in self.Ut:
            u += f(s0) * g(c0)
        for f, g in self.Lt:
            Lu += f(s0) * g(c0)
        if self.tau == 0:
            return Lu + lhs * e_signed + self.problem.src_at(s0, c0)
        u += self.tau * e_signed * s0
        Lu += lhs * self.tau * e_signed
        D = 1 + self.basis.D_factor(s0) * u
        if D <= 0:
            return None
        return Lu + self.tau * self.problem.src_at(s0, c0) / D ** self.P

    def V_terms(self, shift_U):
        """terms of tau * (U/s - shift_U) for U = U_N (u-scale polys): p(s) - tau shift + s(1-s) sum q_l P_l."""
        t = self.tau if self.tau != 0 else fmpq(1)
        Pl = leg_c(2 * self.K)
        terms = [(self.pp - t * shift_U, ONE)]
        for k, q in enumerate(self.qs, 1):
            if q != 0:
                terms.append((self.basis.radial_l_over_s * q, Pl[2 * k]))
        return terms


def point_checks(fn, e_signed, want_sign):
    """exact sign of Gh[U_N + e_signed s] at the problem's rational points; want_sign -1 (super) or +1 (sub)."""
    ps, pcs = fn.problem.point_s, fn.problem.point_c
    for s0 in ps:
        for c0 in pcs:
            v = fn.tauGh_at(s0, c0, e_signed)
            if v is None:
                return dict(status='VIOLATION', at=[qstr(s0), qstr(c0)], why='D <= 0 at a point')
            if (want_sign < 0 and v >= 0) or (want_sign > 0 and v <= 0):
                return dict(status='VIOLATION', at=[qstr(s0), qstr(c0)], value_float=qfloat(v))
    return dict(status='OK', points=len(ps) * len(pcs))


def V_point_min(fn, shift_U, ns=64, nc=8):
    """exact minimum of tau V over a rational grid (used only to build positivity controls)."""
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


# ================================================================== one node certificate (P1)
def certify_node(problem, j, N=None, K=None, log=print, keep=False, coef_override=None):
    """Point certificate at the exact rational j (R1-R5 + extras). Returns (JSON-able record, ExactFn or None).
    Record keys for the Bowen-York problem are exactly those of p2a.py's certify_node."""
    j = Fraction(j)
    if N is None:
        N, K = problem.degrees(j)
    tau = problem.tau(j)
    t0 = time.time()
    if j == 0:
        return problem.certify_zero(N, K, keep)
    if coef_override is None:
        wa, diag = fit_u(problem, j, N, K)
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
    fn = ExactFn(problem, pc, qcs, tau)
    S, _ = fn.residual_S()
    t2 = time.time()
    rec = dict(j=qstr(fq(j)), j_dec=dec_floor(fq(j), 20), degrees=dict(N=N, K=K, prec_bits=PREC, round_bits=RBITS),
               S=qstr(S), S_float=qfloat_up(S))
    rec['_fit'] = diag            # fit diagnostics (not part of the proof; '_' keys are dropped from stored results)
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
    # extras: exact sign checks at rational points
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
    obs = problem.observables(jq, tau, dlo, dhi)
    rel = (dhi - dlo) / dlo
    br = problem.brackets(j, jq, tau, dlo, dhi, obs)
    rec.update(status='CERTIFIED', delta_lo=qstr(dlo), delta_hi=qstr(dhi), delta_lo_dec=dec_floor(dlo),
               delta_hi_dec=dec_ceil(dhi), rel_width=qfloat_up(rel))
    for k, v in obs.items():
        rec[k] = v
    if br is not None:
        rec['t0c'] = br
    t4 = time.time()
    log(f"  node j={float(j):.6f} N={N} K={K}: CERTIFIED delta in [{rec['delta_lo_dec'][:20]}, {rec['delta_hi_dec'][:20]}] "
        f"rel {rec['rel_width']:.2e} eta {rec['eta_float']:.2e} Bernstein pieces {pos['pieces']} "
        f"min {pos['bern_min_float']:.3e}; brackets {None if br is None else br.get('ok')} "
        f"(fit {t1-t0:.1f}s, S {t2-t1:.1f}s, Bern {t3-t2:.1f}s, rest {t4-t3:.1f}s; float res {diag['float_residual']})")
    rec['_timing'] = dict(fit=t1 - t0, S=t2 - t1, bern=t3 - t2, rest=t4 - t3)
    return rec, (fn if keep else None)


def certify_pair(fn, eta, base_sub=None, e_sup=None, e_sub=None, log=print, tag=''):
    """Generic pair check used by nulls: super = base_sup + e_sup s, sub = base_sub - e_sub s (U-scale e's).
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


# ================================================================== Bowen-York spinning puncture (the certified problem)
class BowenYorkSpinningPuncture(AxisymLichnerowiczProblem):
    """Bowen-York spinning puncture, bare mass mu = 1, angular momentum J = j/4 (DESIGN-P2 E1-E4):
    Delta u = -(9/4) J^2 sin^2(theta) r^-6 Phi^-7, Phi = 1 + 1/(2r) + u.  In s: tau = j^2,
    src = (9/4)(1-c^2) s^2 (1-s)^2, delta = M_ADM - 1 = 2 lim r u, chi = J/M_ADM^2 = j/(4 (1 + delta)^2).
    Degrees (N, K) = (40, 8) for j <= 4, (40, 10) above (BUILD-P2a sec 3)."""
    name = 'Bowen-York spinning puncture, mu = 1, J = j/4'
    object_line = 'Bowen-York spinning puncture, mu = 1, J = j/4, delta = M_ADM - 1, chi = J/M_ADM^2'
    SRC_S = fmpq(9, 4) * SP ** 2 * OMS ** 2  # (9/4) s^2 (1-s)^2
    SRC_C = fmpq_poly([1, 0, -1])           # 1 - c^2
    src_terms = [(SRC_S, SRC_C)]
    EXTRA_NODES = (2, 4, 6, 8, 10)

    def src_eval(self, s, c, one, o):
        return 9 * (one - c * c) * s * s * o * o / 4

    def tau(self, j):
        return fq(j) ** 2

    def degrees(self, j):
        jf = float(j)
        if jf <= 4.0:
            return 40, 8
        return 40, 10

    def observables(self, jq, tau, dlo, dhi):
        chi_lo = jq / (4 * (1 + dhi) ** 2)
        chi_hi = jq / (4 * (1 + dlo) ** 2)
        return dict(chi_lo=qstr(chi_lo), chi_hi=qstr(chi_hi), chi_lo_dec=dec_floor(chi_lo), chi_hi_dec=dec_ceil(chi_hi),
                    chi_rel_width=qfloat_up((chi_hi - chi_lo) / chi_lo))

    def brackets(self, j, jq, tau, dlo, dhi, obs):
        """T0c: stage-0 (I(tau v1) <= delta <= tau/40), RN ((1+delta)^2 <= 1 + 3J), chi < 1. Sanity only, not proof."""
        chi_hi = qparse(obs['chi_hi'])
        s0lo = self.stage0_lower(j)
        t0c = dict(stage0_lower=dec_floor(s0lo, 20), stage0_upper=dec_ceil(tau / 40, 20),
                   inside_stage0=bool(s0lo <= dlo and dhi <= tau / 40),
                   rn_ok=bool((1 + dhi) ** 2 <= 1 + 3 * jq / 4), chi_lt_1=bool(chi_hi < 1))
        t0c['ok'] = bool(t0c['inside_stage0'] and t0c['rn_ok'] and t0c['chi_lt_1'])
        return t0c

    # ---- toy model: the order-j^2 closed form v1 (GNPP A1, B1 in P1's s)
    @staticmethod
    def v1_cheb(N, K):
        """v1 = s p + s^2 (1-s) q2 P2: exact Chebyshev coefficients in X = 2s - 1 (U-scale)."""
        s = SP
        om = OMS
        p = (om ** 4 + 5 * om ** 3 * s + 10 * om ** 2 * s ** 2 + 5 * om * s ** 3 + s ** 4) / fmpq(40)
        q2 = -(om * s) / fmpq(20)

        def to_cheb(poly):
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

    @staticmethod
    def v1_box(sb, cb, one):
        """arb enclosure of v1 on a box (s ball, c ball)."""
        s, c = sb, cb
        om = one - s
        p = (om ** 4 + 5 * om ** 3 * s + 10 * om ** 2 * s ** 2 + 5 * om * s ** 3 + s ** 4) / 40
        q2 = -(om * s) / 20
        return s * p + s * s * om * q2 * (3 * c * c - 1) / 2

    def stage0_lower(self, j, nt=160, nc=16):
        """certified lower bound of the stage-0 bracket I(j^2 v1) <= delta (lower Riemann sum in arb)."""
        ctx.prec = 128
        tau = arb(fq(j)) ** 2
        tot = arb(0)
        for i in range(nt):
            ta, tb = arb(i) / nt, arb(i + 1) / nt
            tball = (ta + tb) / 2 + arb(0, 1) * (tb - ta) / 2
            for k in range(nc):
                ca, cb = arb(k) / nc, arb(k + 1) / nc
                cball = (ca + cb) / 2 + arb(0, 1) * (cb - ca) / 2
                Dmax = (1 + (1 - tball) * tau * self.v1_box(tball, cball, arb(1))).upper()
                low = ta * ta * (1 - tb) ** 3 * (1 - cb * cb) / Dmax ** 7
                tot += low.lower() * (tb - ta) * (cb - ca)
        val = (arb(9) / 4 * tau * tot).lower()
        return arb_to_fmpq(val)

    def certify_zero(self, N, K, keep):
        """T0a/T0b: u = 0 is exact at j = 0; the U-scale limit problem is linear with exact solution v1. The engine
        must give S == 0 and REFUSE (never a PASS of the sign checks)."""
        pc0, qcs0 = self.v1_cheb(N, K)
        fn0 = ExactFn(self, pc0, qcs0, fmpq(0))
        S0, _ = fn0.residual_S()
        rec = dict(j='0/1', j_dec='0', degrees=dict(N=N, K=K), S=qstr(S0), S_float=qfloat(S0), p0_U=qstr(fn0.p0),
                   status='REFUSED/TRIVIAL' if not S0 > 0 else 'ERROR: S > 0 for the exact stage-0 solution',
                   why='j = 0: u = 0 exact (delta = 0); U = v1 gives S == 0, a non-strict certificate')
        return rec, (fn0 if keep else None)

    # ---- grid, K1 cells, K2 cover
    def grid_nodes(self, j0=Fraction(1), ratio=Fraction(51, 50), j_end=10, extra=None):
        """j_{k+1} = ratio j_k from j0 to the first node >= j_end, plus extra nodes (exact rationals)."""
        extra = self.EXTRA_NODES if extra is None else extra
        nodes = []
        j = Fraction(j0)
        while True:
            nodes.append(j)
            if j >= j_end:
                break
            j = j * Fraction(ratio)
        for e in extra:
            if Fraction(e) not in nodes:
                nodes.append(Fraction(e))
        return sorted(nodes)

    @staticmethod
    def cell_record(na, nb):
        """K1 cell between certified node records na (j_a) < nb (j_b), arb at 256 bits, outward rounding.
        Uses delta(j)/sqrt(j) nondecreasing (DGC Lemma 3.1 + scaling; Bowen-York specific)."""
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

    def build_result(self, node_recs, sha, j_first=1, j_last=10,
                     grid_line='j_{k+1} = (51/50) j_k from 1 to the first node >= 10, plus j = 2, 4, 6, 8, 10'):
        """nodes + K1 cells + K2 no-gap check, in p2a.py's result format."""
        nodes = sorted(node_recs, key=lambda r: Fraction(*map(int, r['j'].split('/'))))
        clean = [{k: v for k, v in r.items() if not k.startswith('_')} for r in nodes]
        all_cert = all(r['status'] == 'CERTIFIED' for r in clean)
        cells = []
        if all_cert:
            for a, b in zip(clean[:-1], clean[1:]):
                cells.append(self.cell_record(a, b))
        js = [Fraction(*map(int, r['j'].split('/'))) for r in clean]
        no_gaps = bool(all_cert and js and js[0] == j_first and js[-1] >= j_last
                       and all(x < y for x, y in zip(js[:-1], js[1:]))
                       and all(cells[i]['j_hi'] == cells[i + 1]['j_lo'] for i in range(len(cells) - 1))
                       and len(cells) == len(clean) - 1 and all(c['nonempty'] for c in cells))
        return dict(instrument_sha256=sha, object=self.object_line, grid=grid_line,
                    nodes=clean, cells=cells, all_nodes_certified=all_cert, no_gaps=no_gaps,
                    all_t0c_ok=all(r.get('t0c', {}).get('ok', False) for r in clean),
                    max_cell_rel_width=max((c['rel_width'] for c in cells), default=None),
                    max_cell_rel_width_uniform=max((c['rel_width_uniform'] for c in cells), default=None))


BOWEN_YORK = BowenYorkSpinningPuncture()


def certify_bowen_york(j, **kw):
    """Convenience: point certificate for the Bowen-York spinning puncture at exact rational j (mu = 1, J = j/4)."""
    return certify_node(BOWEN_YORK, j, **kw)
