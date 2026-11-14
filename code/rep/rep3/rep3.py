"""
rep3.py - blind replication seat REP3 for BHSLICE: certified Bowen-York spinning puncture at HIGH spin (j = 300, 1000).

Object (DESIGN.md 1.1-1.3):  Delta u = -(9/4) J^2 sin^2(theta) / (r^6 (1 + 1/(2r) + u)^7),  u -> 0 at infinity,
                             mu = 1,  j := 4 J,  delta := M_ADM - 1 = 2 lim r u.

Formulation of THIS seat (derivation in REP3-NOTES.md):
  * radial map y = a/(r + a) with a = 9/20 (deliberately NOT a = 1/2 and NOT a = 1); Chebyshev T_n(2y - 1);
  * angular variable v = cos^2(theta), basis tau_k(v) = T_k(2v - 1) = cos(2 k theta);
  * u = t w, t := j^2,  w = y p(y) + y^2 (1 - y) sum_{k=1..K} q_k(y) tau_k(v);
  * divided operator (r^2 (Delta u - F) = t y^2 (1-y) Gh[w], positive prefactor on 0 < y < 1):
        Gh[w] = (1-y) w_yy - 2 w_y + Ang[w]/(y^2 (1-y)) + CS (1-v) y^2 (1-y)^2 D^-7,
        Ang = 4 v (1-v) w_vv + (2 - 6 v) w_v,  D = 1 - b y + t (1-y) w,  b = 1 - 1/(2a) = -1/9,  CS = s1 / a^4,
        s1 = 9/64 (= (9/4)/16, the source constant over J^2 = t/16);
  * solver: Newton collocation (double) with continuation in j, prolongation to the final degree, iterative refinement
    with arb residuals (closed-form cos(n phi) basis values, exact-rational derivative coefficients: no ball recurrences);
    coefficients rounded to dyadic rationals -> u_N is an explicit function with exact coefficients.
  * certificate (this seat's engine, different from an exact global residual polynomial):
        Q := D_N^7 Lh w_N + CS (1-v) y^2 (1-y)^2 is a polynomial of known degree (8 Ny + 30 in y, 8 K in v).
        On each y-piece B = [y0, y1] x [0, 1], Q is sampled in arb at a tensor Chebyshev-Gauss grid with MORE nodes
        than its degree, so the discrete Chebyshev transform returns its local Chebyshev coefficients exactly (up to
        ball radii); sup_B |Q| <= sum |c_kl| (|T_k| <= 1).  Also inf_B D_N >= max(1 - b y0, c_00(D) - sum' |c_kl(D)|).
        sup_B |Gh[w_N]| <= sup_B |Q| / inf_B D_N^7;  eta := (1/2) max_B (...) (1 + 2^-20), rounded up to a dyadic.
        Positivity: V_- := w_-/y = p + y (1-y) sum q_k tau_k - eta > 0 by the same exact-degree Chebyshev form
        (c_00 - sum' |c| > 0) with adaptive bisection.  Then u_- >= 0, D_- >= 1 - b y >= 1, D_N >= D_-, and since
        D^-7 is decreasing and Lh(y) = -2:  Gh[w_N + eta y] <= Gh[w_N] - 2 eta < 0 < Gh[w_N] + 2 eta <= Gh[w_N - eta y]
        on the CLOSED square (y = 0 is infinity, y = 1 the puncture).  Phi_+- >= 1 + 1/(2r) (kappa = 1); u_+- are
        Lipschitz at 0 and -> 0 at infinity; Lemma C (DESIGN 2.1) gives u_- <= u <= u_+.
        Slope route (r y -> a):  delta in 2 a t [p(0) - eta, p(0) + eta]   (exact rationals).
Run:  python rep3.py              (everything; writes rep3_result.json).  Stage cache: rep3_stages.json (JSON strings,
      never pickle); exit code 3 = "rerun to resume" if a time budget stops it, 0 = finished.
"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
import json, math, sys, time
from fractions import Fraction
import numpy as np
from scipy.linalg import lu_factor, lu_solve
from flint import arb, arb_mat, fmpq, fmpq_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
A_MAP = fmpq(9, 20)
B_MAP = 1 - 1 / (2 * A_MAP)          # -1/9
S1 = fmpq(9, 64)


def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)


def qf(q):
    return float(Fraction(int(q.p), int(q.q)))


def qstr(q):
    return "%s/%s" % (q.p, q.q)


def parse_q(s):
    a, b = s.split('/')
    return fmpq(int(a), int(b))


def arb_upper_q(x):
    """exact rational upper bound of an arb ball."""
    m, e = x.mid().man_exp(); r, er = x.rad().mid().man_exp()
    m, e, r, er = int(m), int(e), int(r), int(er)
    mq = fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))
    rq = fmpq(r * 2 ** er) if er >= 0 else fmpq(r, 2 ** (-er))
    return mq + rq


def arb_lower_q(x):
    return -arb_upper_q(-x)


# =============================================================================================
# 1. double collocation
# =============================================================================================
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
    def __init__(self, Ny, K, a=A_MAP):
        self.Ny, self.K = Ny, K
        self.a = qf(a); self.b = 1 - 1 / (2 * self.a)
        y = gauss01(Ny + 1); v = gauss01(K + 1)
        self.y, self.v = y, v
        T, T1, T2 = cheb_cf(2 * y - 1, Ny)
        Ty, Tyy = 2 * T1, 4 * T2
        yy = y[:, None]; om = 1 - yy
        f = yy * T; f1 = T + yy * Ty; f2 = 2 * Ty + yy * Tyy
        RF = om * f2 - 2 * f1
        h = yy ** 2 - yy ** 3; h1 = 2 * yy - 3 * yy ** 2; h2 = 2 - 6 * yy
        g = h * T; g1 = h1 * T + h * Ty; g2 = h2 * T + 2 * h1 * Ty + h * Tyy
        RG = om * g2 - 2 * g1
        Tv, Tv1, Tv2 = cheb_cf(2 * v - 1, K)
        Tv = Tv[:, 1:]; Tv1 = 2 * Tv1[:, 1:]; Tv2 = 4 * Tv2[:, 1:]
        vv = v[:, None]
        AK = 4 * vv * (1 - vv) * Tv2 + (2 - 6 * vv) * Tv1
        n1 = Ny + 1
        L = np.zeros((n1 * (K + 1), n1 * (K + 1))); U = np.zeros_like(L)
        for jv in range(K + 1):
            rows = np.arange(n1) * (K + 1) + jv
            L[rows, :n1] = RF; U[rows, :n1] = f
            for k in range(1, K + 1):
                L[rows, k * n1:(k + 1) * n1] = RG * Tv[jv, k - 1] + T * AK[jv, k - 1]
                U[rows, k * n1:(k + 1) * n1] = g * Tv[jv, k - 1]
        self.L, self.U = L, U
        self.Y = np.repeat(y, K + 1); self.V = np.tile(v, n1)

    def F(self, w, t, cs):
        D = 1 - self.b * self.Y + t * (1 - self.Y) * (self.U @ w)
        sh = cs * (1 - self.V) * self.Y ** 2 * (1 - self.Y) ** 2
        return self.L @ w + sh * D ** -7, D

    def jac(self, w, t, cs):
        D = 1 - self.b * self.Y + t * (1 - self.Y) * (self.U @ w)
        sh = cs * (1 - self.V) * self.Y ** 2 * (1 - self.Y) ** 2
        return self.L - (7 * t * sh * (1 - self.Y) * D ** -8)[:, None] * self.U

    def newton(self, t, cs, w, iters=40):
        for it in range(iters):
            R, D = self.F(w, t, cs)
            dw = np.linalg.solve(self.jac(w, t, cs), R)
            w = w - dw
            if np.max(np.abs(dw)) <= 1e-16 * max(1e-300, np.max(np.abs(w))):
                break
        R, D = self.F(w, t, cs)
        return w, float(np.max(np.abs(R))), it + 1


def prolong(w, Ny, K, Ny2, K2):
    n1, m1 = Ny + 1, Ny2 + 1
    out = np.zeros(m1 * (K2 + 1))
    for k in range(min(K, K2) + 1):
        nn = min(n1, m1)
        out[k * m1:k * m1 + nn] = w[k * n1:k * n1 + nn]
    return out


def continuation(js_targets, cs, Ny=80, K=8, ratio=1.12, w0=None, j_start=0.5):
    """coarse continuation; returns dict j -> coarse solution."""
    disc = Disc(Ny, K)
    w = np.zeros(disc.L.shape[0]) if w0 is None else w0
    if w0 is None:
        w, res, it = disc.newton(0.0, cs, w)
    out = {}
    x = j_start
    for jt in sorted(js_targets):
        while x < jt:
            w, res, it = disc.newton(x * x, cs, w)
            x *= ratio
        w, res, it = disc.newton(jt * jt, cs, w)
        out[jt] = (w.copy(), res)
        x = jt * ratio
    return out, (Ny, K)


# =============================================================================================
# 2. arb evaluation of w, Lh w (separable), refinement
# =============================================================================================
def cheb_deriv(c):
    """Chebyshev coefficients (in X) of the X-derivative; works for fmpq or arb lists."""
    N = len(c) - 1
    if N == 0:
        return [c[0] * 0]
    d = [c[0] * 0] * (N + 2)
    for k in range(N - 1, -1, -1):
        d[k] = d[k + 2] + 2 * (k + 1) * c[k + 1]
    d[0] = d[0] / 2
    return d[:N]


def ydiff(c):
    return [2 * x for x in cheb_deriv(c)]


def pad(c, n):
    return list(c) + [c[0] * 0] * (n - len(c))


def tau_polys(K):
    """tau_k(v) = T_k(2v-1) and A_k = 4v(1-v) tau'' + (2-6v) tau' as exact fmpq_poly, k = 0..K."""
    x = fmpq_poly([-1, 2]); V = fmpq_poly([0, 1])
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


class Fields:
    """w = y p + y^2(1-y) sum q_k tau_k, with Chebyshev coefficient lists p, q[1..K] (fmpq or arb)."""

    def __init__(self, p, q):
        self.p = list(p); self.q = [list(x) for x in q]
        self.Ny = len(p) - 1; self.K = len(q)
        n1 = self.Ny + 1
        cols = []
        for c in [self.p] + self.q:
            d1 = ydiff(c); d2 = ydiff(d1)
            cols += [pad(c, n1), pad(d1, n1), pad(d2, n1)]
        self.cols = cols            # 3 (K+1) columns

    def radial(self, ys):
        """per y: lists for w-columns [f0, g_1..g_K] and Lh-columns [Rad f0, Rad g_k, q_k] (arb)."""
        n1 = self.Ny + 1
        rows = []
        for y in ys:
            ph = (2 * y - 1).acos()
            rows.append([(ph * n).cos() for n in range(n1)])
        TY = arb_mat(rows)
        C = arb_mat([[arb(col[n]) for col in self.cols] for n in range(n1)])
        P = TY * C                  # len(ys) x 3(K+1)
        W = []; LW = []
        K = self.K
        for i, y in enumerate(ys):
            om = 1 - y
            p, p1, p2 = P[i, 0], P[i, 1], P[i, 2]
            f0 = y * p; f1 = p + y * p1; f2 = 2 * p1 + y * p2
            wr = [f0]; lr = [om * f2 - 2 * f1]
            h = y * y * om; h1 = 2 * y - 3 * y * y; h2 = 2 - 6 * y
            qs = []
            for k in range(1, K + 1):
                q, q1, q2 = P[i, 3 * k], P[i, 3 * k + 1], P[i, 3 * k + 2]
                g = h * q; g1 = h1 * q + h * q1; g2 = h2 * q + 2 * h1 * q1 + h * q2
                wr.append(g); lr.append(om * g2 - 2 * g1); qs.append(q)
            W.append(wr); LW.append(lr + qs)
        return arb_mat(W), arb_mat(LW)

    def angular(self, vs):
        T, A = tau_polys(self.K)
        Aw = [[horner(T[k], v) for v in vs] for k in range(self.K + 1)]
        Aw[0] = [arb(1)] * len(vs)
        Al = [[arb(1)] * len(vs)] + [[horner(T[k], v) for v in vs] for k in range(1, self.K + 1)] + \
             [[horner(A[k], v) for v in vs] for k in range(1, self.K + 1)]
        return arb_mat(Aw), arb_mat(Al)

    def grid(self, ys, vs):
        Wr, Lr = self.radial(ys)
        Aw, Al = self.angular(vs)
        return Wr * Aw, Lr * Al


def split(w, Ny, K):
    n1 = Ny + 1
    return w[:n1], [w[k * n1:(k + 1) * n1] for k in range(1, K + 1)]


def arb_nodes01(N):
    pi = arb.pi()
    return [(1 + (pi * fmpq(2 * i + 1, 2 * N)).cos()) / 2 for i in range(N)]


def refine(disc, wf, t_q, cs_q, prec, maxit=25):
    Ny, K = disc.Ny, disc.K
    ctx.prec = prec
    ys = arb_nodes01(Ny + 1); vs = arb_nodes01(K + 1)
    t = arb(t_q); cs = arb(cs_q); b = arb(B_MAP)
    lu = lu_factor(disc.jac(wf, qf(t_q), qf(cs_q)))
    w = [arb(float(x)) for x in wf]
    hist = []
    target = 2.0 ** (-(prec - 40))
    for it in range(maxit):
        p, q = split(w, Ny, K)
        W, LW = Fields(p, q).grid(ys, vs)
        R = []
        for i in range(Ny + 1):
            y = ys[i]
            sh = cs * y * y * (1 - y) * (1 - y)
            for jv in range(K + 1):
                D = 1 - b * y + t * (1 - y) * W[i, jv]
                R.append(LW[i, jv] + sh * (1 - vs[jv]) / D ** 7)
        Rf = np.array([float(r.mid()) for r in R])
        res = float(np.max(np.abs(Rf)))
        hist.append(res)
        if res < target or (len(hist) > 3 and res > 0.5 * hist[-2]):
            break
        dw = lu_solve(lu, Rf)
        w = [(w[i] - arb(float(dw[i]))).mid() for i in range(len(w))]
    return w, hist


def to_dyadic(x, bits):
    m, e = x.mid().man_exp()
    m = int(m); e = int(e)
    if e >= -bits:
        return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))
    sh = -bits - e
    return fmpq((m + (1 << (sh - 1))) >> sh, 2 ** bits)


# =============================================================================================
# 3. certificate engine: exact-degree Chebyshev forms on pieces
# =============================================================================================
_COS = {}


def dct_mat(N):
    """C[k][i] = (2 - delta_k0)/N cos(k (2i+1) pi / (2N)) as arb_mat (N x N); exact for degree < N."""
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


def nodes_on(lo, hi, N):
    """Chebyshev-Gauss nodes on [lo, hi] (ordered as in dct_mat: x_i = mid + half cos((2i+1)pi/2N))."""
    pi = arb.pi(); mid = arb((lo + hi) / 2); half = arb((hi - lo) / 2)
    return [mid + half * (pi * fmpq(2 * i + 1, 2 * N)).cos() for i in range(N)]


def cheb_coeffs(vals, Ny_, Mv_):
    """vals: arb_mat (Ny_ x Mv_) on the tensor node grid -> local Chebyshev coefficients arb_mat."""
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
    def __init__(self, wq, Ny, K):
        self.wq = wq; self.Ny = Ny; self.K = K
        p, q = split(wq, Ny, K)
        self.F = Fields(p, q)
        self.p0 = sum((wq[n] * (1 if n % 2 == 0 else -1) for n in range(Ny + 1)), fmpq(0))


def V_grid(ans, ylo, yhi, vlo, vhi, shift):
    """V = w/y - shift on the exact-degree grid of the box; returns (vals, N, M)."""
    N = ans.Ny + 3; M = ans.K + 1
    ys = nodes_on(ylo, yhi, N); vs = nodes_on(vlo, vhi, M)
    W, _ = ans.F.grid(ys, vs)
    sh = arb(shift)
    vals = arb_mat(N, M, [W[i, j] / ys[i] - sh for i in range(N) for j in range(M)])
    return vals, N, M


def certify_positive(ans, shift, init_pieces=4, max_boxes=600):
    stack = [(fmpq(i, init_pieces), fmpq(i + 1, init_pieces), fmpq(0), fmpq(1), 0) for i in range(init_pieces)]
    nbox = 0; minlb = None
    while stack:
        ylo, yhi, vlo, vhi, dep = stack.pop()
        nbox += 1
        vals, N, M = V_grid(ans, ylo, yhi, vlo, vhi, shift)
        lb = chebform_lower(vals, N, M)
        if lb > 0:
            lbq = arb_lower_q(lb)
            minlb = lbq if minlb is None or lbq < minlb else minlb
            continue
        # a node value <= 0 means genuinely not positive
        if any(vals[i, j] < 0 for i in range(N) for j in range(M)):
            return False, nbox, None
        if nbox > max_boxes or dep > 14:
            return False, nbox, None
        if (yhi - ylo) >= (vhi - vlo) / 4:
            ym = (ylo + yhi) / 2
            stack += [(ylo, ym, vlo, vhi, dep + 1), (ym, yhi, vlo, vhi, dep + 1)]
        else:
            vm = (vlo + vhi) / 2
            stack += [(ylo, yhi, vlo, vm, dep + 1), (ylo, yhi, vm, vhi, dep + 1)]
    return True, nbox, minlb


def piece_bound(ans, ylo, yhi, t, cs):
    """upper bound of sup |Gh[w_N]| on [ylo,yhi] x [0,1] (rational), plus diagnostics."""
    Ny, K = ans.Ny, ans.K
    if t == 0:
        dQy = 7 * 1 + Ny + 2; dQv = K
    else:
        dQy = 8 * Ny + 30; dQv = 8 * K
    N = dQy + 1; M = dQv + 1
    ys = nodes_on(ylo, yhi, N); vs = nodes_on(fmpq(0), fmpq(1), M)
    W, LW = ans.F.grid(ys, vs)
    ta = arb(t); ca = arb(cs); b = arb(B_MAP)
    flat = []
    for i in range(N):
        y = ys[i]; om = 1 - y
        c0 = 1 - b * y; c1 = ta * om; sh = ca * y * y * om * om
        for jv in range(M):
            D = c0 + c1 * W[i, jv]
            D2 = D * D; D7 = D2 * D2 * D2 * D
            flat.append(D7 * LW[i, jv] + sh * (1 - vs[jv]))
    Qv = arb_mat(N, M, flat)
    SQ = abs_sum(cheb_coeffs(Qv, N, M))
    # lower bound of D_N on the piece
    ND = Ny + 5; MD = K + 1
    ysd = nodes_on(ylo, yhi, ND); vsd = nodes_on(fmpq(0), fmpq(1), MD)
    Wd, _ = ans.F.grid(ysd, vsd)
    Dv = arb_mat(ND, MD, [1 - b * ysd[i] + ta * (1 - ysd[i]) * Wd[i, jv] for i in range(ND) for jv in range(MD)])
    Dlow_cheb = arb_lower_q(chebform_lower(Dv, ND, MD))
    Dlow = max(Dlow_cheb, 1 - B_MAP * ylo)      # second term valid given u_- >= 0 (checked separately)
    SQq = arb_upper_q(SQ)
    bound = SQq / Dlow ** 7
    return bound, dict(ylo=str(ylo), yhi=str(yhi), sup_Q_le="%.3e" % qf(SQq), Dlow="%.6f" % qf(Dlow),
                       Dlow_cheb="%.6f" % qf(Dlow_cheb), bound="%.3e" % qf(bound))


def round_up_dyadic(q, bits=40):
    if q <= 0:
        return fmpq(0)
    e = math.floor(math.log2(qf(q))) - bits
    sc = fmpq(2) ** (-e)
    v = q * sc
    r = fmpq(-((-v.p) // v.q)) / sc
    assert r >= q
    return r


def global_abs_w(ans):
    N = ans.Ny + 4; M = ans.K + 1
    ys = nodes_on(fmpq(0), fmpq(1), N); vs = nodes_on(fmpq(0), fmpq(1), M)
    W, _ = ans.F.grid(ys, vs)
    return arb_upper_q(abs_sum(cheb_coeffs(W, N, M)))


def make_pieces(n_pieces):
    """n_pieces uniform y-pieces; the first one is split geometrically 4 more times toward y = 0 (infinity)."""
    edges = [fmpq(i, n_pieces) for i in range(n_pieces + 1)]
    h = fmpq(1, n_pieces)
    extra = [h / 2 ** m for m in range(4, 0, -1)]
    return [fmpq(0)] + extra + edges[1:]


def certify(wq, Ny, K, t, cs, n_pieces=16, prec=256, t_slab=None, eta_override=None):
    """t_slab: certify uniformly for all t' in [0, t_slab] (requires t == 0 fit)."""
    ctx.prec = prec
    tm = time.time()
    ans = Ansatz(wq, Ny, K)
    pieces = []; sup = fmpq(0)
    E = make_pieces(n_pieces)
    assert E[0] == 0 and E[-1] == 1 and all(E[i] < E[i + 1] for i in range(len(E) - 1))
    for i in range(len(E) - 1):
        bd, info = piece_bound(ans, E[i], E[i + 1], t, cs)
        pieces.append(info)
        sup = max(sup, bd)
    extra = fmpq(0)
    if t_slab is not None:
        Wsup = global_abs_w(ans)
        extra = 7 * t_slab * cs * Wsup / 16      # |S (D_t^-7 - D_0^-7)| <= (cs/16) * 7 t sup|w_N|  (D >= 1)
    sup_tot = sup + extra
    need = sup_tot * (1 + fmpq(1, 2 ** 20)) / 2
    eta = round_up_dyadic(need) if eta_override is None else eta_override
    strict = bool(eta_override is None and eta > 0 and 2 * eta > sup_tot)
    pos_ok, nbox, minlb = certify_positive(ans, eta)
    return dict(sup_Gh_le=sup_tot, eta=eta, p0=ans.p0, pos_ok=bool(pos_ok), pos_boxes=nbox, pos_minlb=minlb,
                strict=strict, certified=bool(strict and pos_ok), pieces=pieces, slab_extra=extra,
                secs=round(time.time() - tm, 1))


# =============================================================================================
# 4. exact-rational sample of the sign conditions (independent of the arb engine)
# =============================================================================================
def cheb_eval_q(c, X):
    b1 = fmpq(0); b2 = fmpq(0)
    for ck in reversed(c[1:]):
        b1, b2 = 2 * X * b1 - b2 + ck, b1
    return X * b1 - b2 + c[0]


class ExactEval:
    def __init__(self, wq, Ny, K):
        p, q = split(wq, Ny, K)
        self.P = [p, ydiff(p), ydiff(ydiff(p))]
        self.Q = [[c, ydiff(c), ydiff(ydiff(c))] for c in q]
        self.T, self.A = tau_polys(K)
        self.K = K

    def at(self, y, v):
        X = 2 * y - 1
        p, p1, p2 = [cheb_eval_q(c, X) for c in self.P]
        om = 1 - y
        w = y * p; f1 = p + y * p1; f2 = 2 * p1 + y * p2
        L = om * f2 - 2 * f1
        h = y * y * om; h1 = 2 * y - 3 * y * y; h2 = 2 - 6 * y
        for k in range(1, self.K + 1):
            q, q1, q2 = [cheb_eval_q(c, X) for c in self.Q[k - 1]]
            tk = self.T[k](v); ak = self.A[k](v)
            w += h * q * tk
            L += (om * (h2 * q + 2 * h1 * q1 + h * q2) - 2 * (h1 * q + h * q1)) * tk + q * ak
        return w, L


def sample_signs(wq, Ny, K, t, cs, eta, n_y=24, n_v=6):
    E = ExactEval(wq, Ny, K)
    mx_plus = None; mn_minus = None; mx_abs = fmpq(0); npos = nneg = 0
    for i in range(n_y + 1):
        y = fmpq(i, n_y)
        for l in range(n_v + 1):
            v = fmpq(l, n_v)
            w, L = E.at(y, v)
            src = cs * (1 - v) * y * y * (1 - y) ** 2
            DN = 1 - B_MAP * y + t * (1 - y) * w
            Gn = L + src / DN ** 7
            Gp = L - 2 * eta + src / (DN + t * (1 - y) * eta * y) ** 7
            Gm = L + 2 * eta + src / (DN - t * (1 - y) * eta * y) ** 7
            mx_plus = Gp if mx_plus is None or Gp > mx_plus else mx_plus
            mn_minus = Gm if mn_minus is None or Gm < mn_minus else mn_minus
            mx_abs = max(mx_abs, abs(Gn)); npos += int(Gn > 0); nneg += int(Gn < 0)
    return dict(max_Gh_plus="%.4e" % qf(mx_plus), min_Gh_minus="%.4e" % qf(mn_minus),
                max_abs_Gh_N="%.4e" % qf(mx_abs), n_GhN_pos=npos, n_GhN_neg=nneg,
                signs_ok=bool(mx_plus < 0 < mn_minus))


# =============================================================================================
# 5. output helpers
# =============================================================================================
def dec_floor(q, nd):
    if q == 0:
        return "0"
    e = math.floor(math.log10(abs(qf(q))))
    sc = nd - 1 - e
    v = q * fmpq(10) ** sc if sc >= 0 else q / fmpq(10) ** (-sc)
    iv = int(v.p // v.q)
    sg = "-" if iv < 0 else ""
    s = str(abs(iv))
    if sc <= 0:
        return sg + s + "0" * (-sc)
    if len(s) <= sc:
        s = "0" * (sc - len(s) + 1) + s
    return sg + s[:-sc] + "." + s[-sc:]


def dec_ceil(q, nd):
    s = dec_floor(-q, nd)
    return s[1:] if s.startswith("-") else ("-" + s if s != "0" else s)


def integral_route(wq, Ny, K, t, cs):
    """double: delta = 2 (cs a) t int y^2 (1-y)^3 <(1-v) D^-7>_c dy  (a^4 cs = s1;  2 s1 t / a^3 = 2 cs a t)."""
    a = qf(A_MAP); b = qf(B_MAP)
    w = np.array([qf(x) for x in wq]); n1 = Ny + 1
    xg, wg = np.polynomial.legendre.leggauss(400)
    yq = (xg + 1) / 2; wy = wg / 2
    cq = xg; wc = wg / 2; vq = cq ** 2
    T, _, _ = cheb_cf(2 * yq - 1, Ny)
    Tv, _, _ = cheb_cf(2 * vq - 1, K)
    p = T @ w[:n1]
    qv = np.stack([T @ w[k * n1:(k + 1) * n1] for k in range(1, K + 1)])   # K x ny
    tot = 0.0
    for iy in range(len(yq)):
        y = yq[iy]
        ww = y * p[iy] + y * y * (1 - y) * (qv[:, iy] @ Tv[:, 1:].T)
        D = 1 - b * y + t * (1 - y) * ww
        tot += wy[iy] * y ** 2 * (1 - y) ** 3 * np.sum(wc * (1 - vq) * D ** -7)
    return 2 * cs * a * t * tot


# =============================================================================================
# 6. drivers
# =============================================================================================
def solve_fine(wc, resc, Ny, K, j_q, cs_q, prec, label):
    disc = Disc(Ny, K)
    w0 = prolong(wc, resc[0], resc[1], Ny, K)
    t_f = qf(j_q) ** 2
    wf, rf, it = disc.newton(t_f, qf(cs_q), w0)
    wa, hist = refine(disc, wf, j_q * j_q, cs_q, prec)
    wq = [to_dyadic(x, prec - 32) for x in wa]
    return wq, rf, hist


def run_point(wq, Ny, K, j_q, cs_q, prec, label, n_pieces=16, extra=None):
    t_q = j_q * j_q
    c = certify(wq, Ny, K, t_q, cs_q, n_pieces=n_pieces, prec=prec)
    a = A_MAP
    lo = 2 * a * t_q * (c['p0'] - c['eta']); hi = 2 * a * t_q * (c['p0'] + c['eta'])
    n1 = Ny + 1
    tail_p = max(abs(qf(x)) for x in wq[n1 - 4:n1]) / max(abs(qf(x)) for x in wq[:n1])
    tail_qK = max(abs(qf(wq[K * n1 + n])) for n in range(n1)) / max(abs(qf(x)) for x in wq[:n1])
    out = dict(label=label, j=str(j_q), Ny=Ny, K=K, prec_bits=prec, n_pieces=n_pieces,
               rel_tail_p_last4="%.2e" % tail_p, rel_max_qK="%.2e" % tail_qK,
               sup_Gh_N_le="%.4e" % qf(c['sup_Gh_le']), eta="%.4e" % qf(c['eta']), eta_exact=qstr(c['eta']),
               p0_exact=qstr(c['p0']), positivity_certified=c['pos_ok'], positivity_boxes=c['pos_boxes'],
               positivity_min_lower_bound=("%.4e" % qf(c['pos_minlb'])) if c['pos_minlb'] is not None else None,
               strict=c['strict'], certified=c['certified'], pieces=c['pieces'],
               delta_lo=dec_floor(lo, 30), delta_hi=dec_ceil(hi, 30), delta_lo_exact=qstr(lo), delta_hi_exact=qstr(hi),
               rel_width="%.3e" % qf((hi - lo) / lo) if lo > 0 else "n/a",
               integral_route_double=repr(integral_route(wq, Ny, K, qf(t_q), qf(cs_q))), secs_cert=c['secs'])
    if extra:
        out.update(extra)
    log("%s Ny=%d K=%d prec=%d: sup|Gh|<=%.3e eta %.3e pos=%s(%d) rel.w %s | delta in [%s, %s] cert=%s (%.1fs)"
        % (label, Ny, K, prec, qf(c['sup_Gh_le']), qf(c['eta']), c['pos_ok'], c['pos_boxes'], out['rel_width'],
           out['delta_lo'][:22], out['delta_hi'][:22], c['certified'], c['secs']))
    return out, (lo, hi), c


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "probe":
    # probe: python rep3.py probe j Ny K prec npieces
    j = int(sys.argv[2]); Ny = int(sys.argv[3]); K = int(sys.argv[4]); prec = int(sys.argv[5]); npc = int(sys.argv[6])
    cs = S1 / A_MAP ** 4
    sol, res = continuation([j], qf(cs))
    log("coarse done", sol[j][1])
    wq, rf, hist = solve_fine(sol[j][0], res, Ny, K, fmpq(j), cs, prec, "probe")
    log("fine newton res %.2e refine %s" % (rf, ["%.1e" % h for h in hist]))
    out, e, c = run_point(wq, Ny, K, fmpq(j), cs, prec, "probe j=%d" % j, n_pieces=npc)
    for pc in out['pieces']:
        print(pc)
    smp = sample_signs(wq, Ny, K, fmpq(j) ** 2, cs, c['eta'], n_y=12, n_v=4)
    log("exact sample:", smp)
    sys.exit(0)


# =============================================================================================
# 7. controls
# =============================================================================================
def identity_controls():
    import random
    out = {}
    ctx.prec = 192
    rnd = random.Random(11)
    # (a) DCT on a piece returns known local Chebyshev coefficients (degree < node count)
    Nn, Mn = 12, 8
    c = [[fmpq(rnd.randint(-50, 50), rnd.randint(1, 9)) for l in range(6)] for k in range(10)]
    pi = arb.pi()
    vals = []
    for i in range(Nn):
        th = pi * fmpq(2 * i + 1, 2 * Nn)
        for m in range(Mn):
            ph = pi * fmpq(2 * m + 1, 2 * Mn)
            vals.append(sum((arb(c[k][l]) * (th * k).cos() * (ph * l).cos() for k in range(10) for l in range(6)), arb(0)))
    C = cheb_coeffs(arb_mat(Nn, Mn, vals), Nn, Mn)
    ok = all(C[k, l].contains(c[k][l] if (k < 10 and l < 6) else 0) for k in range(Nn) for l in range(Mn))
    out['dct_recovers_known_coefficients'] = bool(ok)
    # (b) angular operator: A_k(c^2) = [(1-c^2) g']' with g = tau_k(c^2), and tau_k(c^2) = T_2k(c)
    ok = True
    Cc = fmpq_poly([0, 1]); C2 = Cc * Cc
    T, A = tau_polys(8)
    Tc = [fmpq_poly([1]), Cc]
    for n in range(1, 16):
        Tc.append(2 * Cc * Tc[n] - Tc[n - 1])
    for k in range(1, 9):
        g = T[k](C2)
        ok = ok and (((1 - C2) * g.derivative()).derivative() == A[k](C2)) and (g == Tc[2 * k])
    out['angular_operator_identity'] = bool(ok)
    # (c) radial operator in y = a/(r+a) on r^2 P2: y^2[(1-y)^2 R']' + R Ang(P2) = 0, exact at rational points
    F = Fraction; a = F(9, 20)
    good = True
    for (y, v) in ((F(1, 3), F(1, 5)), (F(5, 7), F(2, 3)), (F(1, 9), F(0))):
        R = a * a * (1 - y) ** 2 / y ** 2
        R1 = -2 * a * a * (1 - y) / y ** 3
        R2 = 2 * a * a * (3 - 2 * y) / y ** 4
        P2 = (3 * v - 1) / 2; AngP2 = (2 - 6 * v) * F(3, 2)
        rad = y * y * ((1 - y) ** 2 * R2 - 2 * (1 - y) * R1)
        good = good and (rad * P2 + R * AngP2 == 0)
        # finite-difference-free consistency of R1, R2 with R (exact rational derivative check by polynomial identity)
        h = F(1, 10 ** 30)
        good = good and abs((a * a * (1 - y - h) ** 2 / (y + h) ** 2 - R) / h - R1) < F(1, 10 ** 20)
    out['radial_operator_harmonic_r2P2'] = bool(good)
    # (d) source transformation: r^2 (9/4) J^2 (1-c^2)/(r^6 Phi^7) = t (s1/a^4) (1-v) y^4 (1-y)^3 / D^7, D = (1-y) Phi
    good = True
    for (y, v, uu, J) in ((F(1, 3), F(1, 5), F(2, 7), F(5, 2)), (F(4, 5), F(0), F(1, 9), F(250)), (F(1, 50), F(1, 2), F(3), F(75))):
        r = a * (1 - y) / y
        Phi = 1 + 1 / (2 * r) + uu
        D = 1 - F(-1, 9) * y + (1 - y) * uu
        t = 16 * J * J
        good = good and (D == (1 - y) * Phi)
        good = good and (r * r * F(9, 4) * J * J * (1 - v) / (r ** 6 * Phi ** 7)
                         == t * (F(9, 64) / a ** 4) * (1 - v) * y ** 4 * (1 - y) ** 3 / D ** 7)
    out['source_transformation'] = bool(good)
    out['B_MAP_is_minus_1_9'] = bool(B_MAP == fmpq(-1, 9))
    # (e) Lh(y) = -2 (w = y <-> p = T_0), through the exact evaluator
    Ny, K = 6, 2
    wq = [fmpq(0)] * ((Ny + 1) * (K + 1)); wq[0] = fmpq(1)
    E = ExactEval(wq, Ny, K)
    out['Lh_of_y_is_minus2'] = bool(all(E.at(fmpq(i, 7), fmpq(l, 3))[1] == -2 for i in range(8) for l in range(4)))
    # (f) arb field evaluator (production path) contains the exact rational evaluation
    Ny, K = 24, 3
    wq = [fmpq(rnd.randint(-2 ** 20, 2 ** 20), 2 ** 24) for _ in range((Ny + 1) * (K + 1))]
    E = ExactEval(wq, Ny, K); An = Ansatz(wq, Ny, K)
    good = True
    for (y, v) in ((fmpq(1, 3), fmpq(2, 7)), (fmpq(0), fmpq(1)), (fmpq(1), fmpq(0)), (fmpq(19, 20), fmpq(1, 2))):
        we, Le = E.at(y, v)
        if y == 0 or y == 1:
            continue   # arb acos(+-1) path not used by the engine (nodes are interior)
        W, LW = An.F.grid([arb(y)], [arb(v)])
        good = good and W[0, 0].contains(we) and LW[0, 0].contains(Le)
    out['arb_evaluator_contains_exact'] = bool(good)
    out['all_ok'] = all(out.values())
    log("identity controls:", out)
    return out


def toy_J0():
    """J = 0: u = 0 is the exact solution (delta = 0); the engine must refuse a strict PASS (eta = 0)."""
    Ny, K, prec = 24, 3, 128
    disc = Disc(Ny, K)
    wf, rf, it = disc.newton(1.0, 0.0, np.zeros(disc.L.shape[0]))
    wa, hist = refine(disc, wf, fmpq(1), fmpq(0), prec)
    wq = [to_dyadic(x, prec - 32) for x in wa]
    c = certify(wq, Ny, K, fmpq(1), fmpq(0), n_pieces=4, prec=prec)
    allzero = all(x == 0 for x in wq)
    out = dict(all_coefficients_zero=bool(allzero), sup_Gh_bound=qstr(c['sup_Gh_le']), eta_exact=qstr(c['eta']),
               strict_certificate=c['strict'], certified=c['certified'], delta="0" if allzero else "nonzero",
               ok=bool(allzero and c['sup_Gh_le'] == 0 and not c['certified']))
    log("toy J=0:", out)
    return out


def v1_gnpp(r, c):
    """GNPP (16)-(17) via DESIGN 1.5, a_design = 1 (mu = 2); for mu = 1 use R = 2 r."""
    A1 = (r ** 4 + 5 * r ** 3 + 10 * r ** 2 + 5 * r + 1) / (40 * (r + 1) ** 5)
    B1 = -r ** 2 / (20 * (r + 1) ** 5)
    return A1 + B1 * (3 * c * c - 1) / 2


def k1_slab(Ny, K, prec, j0_q):
    cs = S1 / A_MAP ** 4
    disc = Disc(Ny, K)
    wf, rf, it = disc.newton(0.0, qf(cs), np.zeros(disc.L.shape[0]))
    wa, hist = refine(disc, wf, fmpq(0), cs, prec)
    wq = [to_dyadic(x, prec - 32) for x in wa]
    t0 = j0_q * j0_q
    c = certify(wq, Ny, K, fmpq(0), cs, n_pieces=8, prec=prec, t_slab=t0)
    lo = 2 * A_MAP * (c['p0'] - c['eta']); hi = 2 * A_MAP * (c['p0'] + c['eta'])
    cont = bool(lo <= fmpq(1, 40) <= hi and c['certified'])
    # GNPP closed-form comparison (double): w_N(y, v) vs v1(2 r, c), r = a (1-y)/y
    E = ExactEval(wq, Ny, K)
    mx = 0.0
    for i in range(1, 20):
        y = fmpq(i, 20)
        for l in range(5):
            v = fmpq(l, 4)
            w, _ = E.at(y, v)
            r = qf(A_MAP) * (1 - qf(y)) / qf(y)
            mx = max(mx, abs(qf(w) - v1_gnpp(2 * r, math.sqrt(qf(v)))))
    out = dict(j0=str(j0_q), Ny=Ny, K=K, prec_bits=prec, refine_history=["%.1e" % h for h in hist],
               sup_Gh_le="%.3e" % qf(c['sup_Gh_le']), slab_extra="%.3e" % qf(c['slab_extra']),
               eta="%.3e" % qf(c['eta']), certified=c['certified'], positivity_certified=c['pos_ok'],
               delta_over_j2_lo=dec_floor(lo, 34), delta_over_j2_hi=dec_ceil(hi, 34),
               lo_exact=qstr(lo), hi_exact=qstr(hi), rel_width="%.3e" % qf((hi - lo) / lo),
               contains_1_over_40=cont, max_abs_diff_vs_GNPP_v1_double="%.2e" % mx)
    log("K1 slab j in (0, %s] Ny=%d K=%d: delta/j^2 in [%s, %s] contains 1/40: %s; |w - v1_GNPP| <= %.1e"
        % (j0_q, Ny, K, out['delta_over_j2_lo'][:30], out['delta_over_j2_hi'][:30], cont, mx))
    return out


def summary(lo, hi, j_q):
    j = j_q
    chi_lo = j / (4 * (1 + hi) ** 2); chi_hi = j / (4 * (1 + lo) ** 2)
    return dict(delta_lo=dec_floor(lo, 30), delta_hi=dec_ceil(hi, 30), delta_lo_exact=qstr(lo), delta_hi_exact=qstr(hi),
                abs_width="%.3e" % qf(hi - lo), rel_width="%.3e" % qf((hi - lo) / lo),
                M_ADM_lo=dec_floor(1 + lo, 30), M_ADM_hi=dec_ceil(1 + hi, 30),
                chi_lo=dec_floor(chi_lo, 30), chi_hi=dec_ceil(chi_hi, 30),
                inside_stage0_upper_j2_over_40=bool(hi <= j * j / 40),
                inside_RN_bound_M2_le_1_plus_3J=bool((1 + hi) ** 2 <= 1 + 3 * (j / 4)),
                chi_below_1=bool(chi_hi < 1))


RES_A = (128, 10, 224)
RES_B = (176, 14, 256)
NPIECES = 16
STAGE_FILE = os.path.join(HERE, "rep3_stages.json")
BUDGET = float(os.environ.get("REP3_BUDGET", "420"))


def sha_self():
    import hashlib
    return hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()


def coarse(cs_q, js=(300, 1000), w0=None):
    sol, res = continuation(list(js), qf(cs_q), w0=w0)
    return sol, res


def target_stage(j, R, cs_q, label):
    Ny, K, prec = R
    sol, res = coarse(S1 / A_MAP ** 4)
    wc = sol[j][0]
    if cs_q != S1 / A_MAP ** 4:
        d = Disc(*res)
        wc, rr, it = d.newton(float(j) ** 2, qf(cs_q), wc)
    wq, rf, hist = solve_fine(wc, res, Ny, K, fmpq(j), cs_q, prec, label)
    out, e, c = run_point(wq, Ny, K, fmpq(j), cs_q, prec, label, n_pieces=NPIECES,
                          extra=dict(fine_newton_residual="%.2e" % rf, refine_history=["%.1e" % h for h in hist]))
    return out, wq, c


def nulls_stage():
    cs = S1 / A_MAP ** 4
    Ny, K, prec = RES_A
    out, wq, c = target_stage(1000, RES_A, cs, "j=1000 (A, for sample/nulls)")
    t = fmpq(1000) ** 2
    res = {}
    smp = sample_signs(wq, Ny, K, t, cs, c['eta'])
    smp['rigorous_bound_sup_Gh'] = "%.4e" % qf(c['sup_Gh_le'])
    res['exact_sample_j1000_A'] = smp
    log("exact rational sample j=1000 (A):", smp)
    # null 1: bumped ansatz (p_0 coefficient + 1e-9 relative)
    wb = list(wq); bump = wq[0] * fmpq(1, 10 ** 9); wb[0] += bump
    smpb = sample_signs(wb, Ny, K, t, cs, c['eta'])
    ob, eb, cb = run_point(wb, Ny, K, fmpq(1000), cs, prec, "null bump 1e-9", n_pieces=NPIECES)
    res['null_bumped_ansatz'] = dict(bump_rel="1e-9", honest_eta="%.4e" % qf(c['eta']),
                                     sample_with_honest_eta_signs_ok=smpb['signs_ok'],
                                     sample_with_honest_eta=smpb, rigorous_eta_bumped=ob['eta'],
                                     rigorous_rel_width_bumped=ob['rel_width'],
                                     eta_inflation_factor="%.3e" % (qf(cb['eta']) / qf(c['eta'])),
                                     ok=bool((not smpb['signs_ok']) and qf(cb['eta']) > 100 * qf(c['eta'])))
    log("null bump:", res['null_bumped_ansatz']['ok'], res['null_bumped_ansatz']['eta_inflation_factor'])
    # null 2: positivity shifted below the sampled minimum must FAIL; honest shift passes
    ctx.prec = prec
    E = ExactEval(wq, Ny, K)
    vmin = None
    for i in range(1, 65):
        y = fmpq(i, 64)
        for l in range(17):
            w, _ = E.at(y, fmpq(l, 16))
            val = w / y
            vmin = val if vmin is None or val < vmin else vmin
    ans = Ansatz(wq, Ny, K)
    okd = certify_positive(ans, vmin + fmpq(1, 10 ** 6))
    oku = certify_positive(ans, c['eta'])
    res['null_positivity'] = dict(sampled_min_V="%.6e" % qf(vmin), shifted_down_passes=bool(okd[0]),
                                  honest_passes=bool(oku[0]), ok=bool((not okd[0]) and oku[0]))
    log("null positivity:", res['null_positivity'])
    return res


def main():
    stages = {}
    me = sha_self()
    if os.path.exists(STAGE_FILE):
        st = json.load(open(STAGE_FILE))
        if st.get('instrument_sha256') == me:
            stages = st['stages']
        else:
            log("stage cache from a different instrument version: ignored")
    order = ['identity', 'toy', 'K1', 'j300_A', 'j300_B', 'j1000_A', 'j1000_B', 'nulls', 'K2_A', 'K2_B']
    cs = S1 / A_MAP ** 4
    cs_k2 = cs * (1 + fmpq(1, 1000))
    for name in order:
        if name in stages:
            continue
        if time.time() - T0 > BUDGET:
            json.dump(dict(instrument_sha256=me, stages=stages), open(STAGE_FILE, 'w'), indent=1, sort_keys=True)
            log("time budget reached before stage %s: saved; exit 3 (rerun to resume)" % name)
            sys.exit(3)
        if name == 'identity':
            r = identity_controls()
        elif name == 'toy':
            r = toy_J0()
        elif name == 'K1':
            j0 = fmpq(1, 2 ** 40)
            r = dict(slab_A=k1_slab(40, 2, 192, j0), slab_B=k1_slab(56, 3, 256, j0))
            r['contains_1_over_40'] = bool(r['slab_A']['contains_1_over_40'] and r['slab_B']['contains_1_over_40'])
        elif name.startswith('j'):
            j = int(name[1:].split('_')[0]); R = RES_A if name.endswith('A') else RES_B
            r = target_stage(j, R, cs, name)[0]
        elif name == 'nulls':
            r = nulls_stage()
        elif name.startswith('K2'):
            R = RES_A if name.endswith('A') else RES_B
            r = target_stage(1000, R, cs_k2, name)[0]
        stages[name] = r
        json.dump(dict(instrument_sha256=me, stages=stages), open(STAGE_FILE, 'w'), indent=1, sort_keys=True)

    res = dict(seat="rep3", instrument_sha256=me,
               object="Bowen-York spinning puncture, mu = 1, j = 4J, delta = M_ADM - 1 = 2 lim r u",
               method=("y = a/(r+a), a = 9/20, Chebyshev in 2y-1 x Chebyshev in 2cos^2(theta)-1; Gauss collocation + "
                       "continuation in j (ratio 1.12) + prolongation + arb refinement; dyadic coefficients; certificate: "
                       "Q = D^7 Gh polynomial sampled in arb on exact-degree Chebyshev grids per y-piece, sup|Q| <= sum|c|, "
                       "sup|Gh| <= max_B sup_B|Q| / inf_B D^7; eta = sup/2 (1+2^-20); u_- >= 0 by Chebyshev-form positivity; "
                       "slope route delta in 2 a j^2 [p(0) -+ eta]"))
    res['identity_controls'] = stages['identity']
    res['toy_J0'] = stages['toy']
    res['K1'] = stages['K1']
    for j in (300, 1000):
        rA, rB = stages['j%d_A' % j], stages['j%d_B' % j]
        eA = (parse_q(rA['delta_lo_exact']), parse_q(rA['delta_hi_exact']))
        eB = (parse_q(rB['delta_lo_exact']), parse_q(rB['delta_hi_exact']))
        lo, hi = max(eA[0], eB[0]), min(eA[1], eB[1])
        s = summary(lo, hi, fmpq(j))
        s['all_certified'] = bool(rA['certified'] and rB['certified'])
        s['A_contains_B'] = bool(eA[0] <= eB[0] and eB[1] <= eA[1])
        s['resolutions_overlap'] = bool(lo <= hi)
        s['rel_width_A'] = rA['rel_width']; s['rel_width_B'] = rB['rel_width']
        s['meets_1e-6'] = bool(qf((hi - lo) / lo) <= 1e-6)
        res['delta_j%d' % j] = s
        res['runs_j%d' % j] = [rA, rB]
        log("FINAL j=%d: delta in [%s, %s] rel width %s certified %s A>B %s" % (j, s['delta_lo'], s['delta_hi'],
                                                                              s['rel_width'], s['all_certified'], s['A_contains_B']))
    res['sample_and_nulls'] = stages['nulls']
    kA, kB = stages['K2_A'], stages['K2_B']
    pA = (parse_q(kA['delta_lo_exact']), parse_q(kA['delta_hi_exact']))
    pB = (parse_q(kB['delta_lo_exact']), parse_q(kB['delta_hi_exact']))
    plo, phi = max(pA[0], pB[0]), min(pA[1], pB[1])
    blo = parse_q(res['delta_j1000']['delta_lo_exact']); bhi = parse_q(res['delta_j1000']['delta_hi_exact'])
    res['K2'] = dict(perturbation="source 9/4 -> 9/4 (1 + 1e-3) at j = 1000", runs=[kA, kB],
                     all_certified=bool(kA['certified'] and kB['certified']),
                     delta_perturbed_lo=dec_floor(plo, 30), delta_perturbed_hi=dec_ceil(phi, 30),
                     shift_lo=dec_floor(plo - bhi, 22), shift_hi=dec_ceil(phi - blo, 22),
                     rel_shift_lo=dec_floor((plo - bhi) / bhi, 16), rel_shift_hi=dec_ceil((phi - blo) / blo, 16))
    log("K2: shift of delta(1000) in [%s, %s], relative [%s, %s]" % (res['K2']['shift_lo'], res['K2']['shift_hi'],
                                                                   res['K2']['rel_shift_lo'], res['K2']['rel_shift_hi']))
    res['expected_false'] = [
        "toy_J0/strict_certificate (J = 0: eta = 0, the engine must refuse a strict PASS)",
        "toy_J0/certified (same: u = 0 is exact, delta = 0, but no strict sub/super pair exists)",
        "sample_and_nulls/null_bumped_ansatz/sample_with_honest_eta_signs_ok (null: bumped ansatz must violate the honest-eta signs)",
        "sample_and_nulls/null_bumped_ansatz/sample_with_honest_eta/signs_ok (same null)",
        "sample_and_nulls/null_positivity/shifted_down_passes (null: V shifted below its minimum must fail)",
    ]
    res['all_flags_ok'] = bool(res['identity_controls']['all_ok'] and res['toy_J0']['ok'] and res['K1']['contains_1_over_40']
                               and res['delta_j300']['all_certified'] and res['delta_j1000']['all_certified']
                               and res['K2']['all_certified'] and res['sample_and_nulls']['null_bumped_ansatz']['ok']
                               and res['sample_and_nulls']['null_positivity']['ok']
                               and res['sample_and_nulls']['exact_sample_j1000_A']['signs_ok'])
    res['total_seconds_this_invocation'] = round(time.time() - T0, 1)
    with open(os.path.join(HERE, "rep3_result.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    log("wrote rep3_result.json; all_flags_ok =", res['all_flags_ok'])


if __name__ == "__main__":
    main()
