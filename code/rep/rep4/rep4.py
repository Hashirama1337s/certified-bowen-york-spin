"""
rep4.py - blind replication seat REP4 for BHSLICE: rigorous enclosures of delta(j) = M_ADM - 1 for the Bowen-York
spinning puncture at j = 3000 and j = 10000.

Object:  Delta_flat u = -(9/4) J^2 sin^2(theta) / (r^6 Phi^7),  Phi = 1 + 1/(2r) + u,  u regular at 0, u -> 0 at infinity,
         mu = 1, j := 4J, delta := M_ADM - 1 = 2 lim r u.

Formulation of THIS seat (derived in REPORT-REP4.md; checked by the exact identity controls below):
  * square-root radial map  Y = 1/(1 + sqrt(r/a)),  a = 1/2   (r = a ((1-Y)/Y)^2; Y = 0 is infinity, Y = 1 the puncture;
    for a = 1/2 the puncture's inversion r -> 1/(4r) becomes Y -> 1-Y);
  * angular variable v = cos^2(theta), basis tau_k(v) = T_k(2v-1) = cos(2k theta);
  * u = Y^2 W,  W = P(Y) + Y(1-Y) sum_{k=1..K} q_k(Y) tau_k(v),  P, q_k Chebyshev in 2Y-1 of degree N;
  * r^2 (Delta u + S Phi^-7) = (1-Y) Y^3 Gh[W],  with
        Gh[W] = -(1/4) [6 W + (6Y-3) W_Y - Y(1-Y) W_YY] + sum_k q_k A_k(v) + CS (1-v) Y^5 (1-Y)^5 / E^7,
        E = 2a (1-Y)^2 (1 + Y^2 W) + Y^2  (= 2a (1-Y)^2 Phi),  CS = 18 a^3 j^2  (times (1+eps) for K2),
        A_k = 4v(1-v) tau_k'' + (2-6v) tau_k'.
    Gh[1] - (source) = -3/2 :  W -> W + eta  is  u -> u + eta Y^2,  and Y^2 is strictly superharmonic.
  * certificate: sub/super pair W_N -+ eta:  Gh[W_N + eta] <= Gh[W_N] - (3/2) eta < 0 < Gh[W_N] + (3/2) eta <= Gh[W_N - eta]
    on the closed square, provided sup |Gh[W_N]| < (3/2) eta and E_- := E_N - 2a(1-Y)^2 Y^2 eta > 0.
    sup |Gh[W_N]| <= max_pieces  sup_B |Q| / inf_B E_N^7,  Q := E_N^7 Lin[W_N] + CS (1-v) Y^5 (1-Y)^5, a polynomial of
    Y-degree <= 8N + 44 and v-degree <= 8K; on each Y-piece Q is sampled in arb on a tensor Chebyshev-Gauss grid with
    MORE nodes than its degree, so the discrete Chebyshev transform returns its local coefficients exactly (balls);
    sup |Q| <= sum |c|.  inf E_N: same exact-degree form (c00 - sum'|c|) with adaptive bisection.
    Comparison (Kato / weak maximum principle, derived in the report) gives u_- <= u <= u_+, so
        delta in [2a (W_N(0) - eta), 2a (W_N(0) + eta)]   (r Y^2 -> a).
Run:  python rep4.py            (all stages; writes result_rep4.json; stage cache rep4_stages.json, exit 3 = rerun)
      python rep4.py probe j N K prec npieces
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
A_MAP = fmpq(1, 2)


def CS_of(j_q, eps=fmpq(0)):
    return 18 * A_MAP ** 3 * j_q * j_q * (1 + eps)


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
    def __init__(self, N, K):
        self.N, self.K = N, K
        a = qf(A_MAP)
        Y = gauss01(N + 1); v = gauss01(K + 1)
        T, T1, T2 = cheb_cf(2 * Y - 1, N)
        TY, TYY = 2 * T1, 4 * T2
        y = Y[:, None]

        def LR(f, f1, f2):
            return -0.25 * (6 * f + (6 * y - 3) * f1 - y * (1 - y) * f2)
        RF = LR(T, TY, TYY)
        h = y - y * y; h1 = 1 - 2 * y; h2 = -2.0
        g = h * T; g1 = h1 * T + h * TY; g2 = h2 * T + 2 * h1 * TY + h * TYY
        RG = LR(g, g1, g2)
        Tv, Tv1, Tv2 = cheb_cf(2 * v - 1, K)
        Tv = Tv[:, 1:]; Tv1 = 2 * Tv1[:, 1:]; Tv2 = 4 * Tv2[:, 1:]
        vv = v[:, None]
        AK = 4 * vv * (1 - vv) * Tv2 + (2 - 6 * vv) * Tv1
        n1 = N + 1
        L = np.zeros((n1 * (K + 1), n1 * (K + 1))); U = np.zeros_like(L)
        for jv in range(K + 1):
            rows = np.arange(n1) * (K + 1) + jv
            L[rows, :n1] = RF; U[rows, :n1] = T
            for k in range(1, K + 1):
                L[rows, k * n1:(k + 1) * n1] = RG * Tv[jv, k - 1] + T * AK[jv, k - 1]
                U[rows, k * n1:(k + 1) * n1] = g * Tv[jv, k - 1]
        self.L, self.U = L, U
        self.Y = np.repeat(Y, K + 1); self.V = np.tile(v, n1)
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
    """coarse continuation in j (CS = CS_scale * j^2); returns coarse solution at j_target."""
    disc = Disc(N, K)
    w = np.zeros(disc.L.shape[0])
    x = min(j_start, j_target)
    while x < j_target:
        w, res, it = disc.newton(CS_scale * x * x, w)
        x *= ratio
    w, res, it = disc.newton(CS_scale * j_target ** 2, w)
    return w, (N, K), res


# =============================================================================================
# 2. arb evaluation of W and Lin[W] (separable), refinement
# =============================================================================================
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
    return [2 * x for x in cheb_deriv(c)]


def pad(c, n):
    return list(c) + [c[0] * 0] * (n - len(c))


def tau_polys(K):
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


def LR_arb(y, f, f1, f2):
    return -(6 * f + (6 * y - 3) * f1 - y * (1 - y) * f2) / 4


class Fields:
    """W = P + Y(1-Y) sum q_k tau_k with Chebyshev coefficient lists P, q[1..K]."""

    def __init__(self, p, q):
        self.p = list(p); self.q = [list(x) for x in q]
        self.N = len(p) - 1; self.K = len(q)
        n1 = self.N + 1
        cols = []
        for c in [self.p] + self.q:
            d1 = ydiff(c); d2 = ydiff(d1)
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
        Wr = []; Lr = []
        K = self.K
        for i, y in enumerate(ys):
            p, p1, p2 = P[i, 0], P[i, 1], P[i, 2]
            wr = [p]; lr = [LR_arb(y, p, p1, p2)]
            h = y - y * y; h1 = 1 - 2 * y
            qs = []
            for k in range(1, K + 1):
                q, q1, q2 = P[i, 3 * k], P[i, 3 * k + 1], P[i, 3 * k + 2]
                g = h * q; g1 = h1 * q + h * q1; g2 = -2 * q + 2 * h1 * q1 + h * q2
                wr.append(g); lr.append(LR_arb(y, g, g1, g2)); qs.append(q)
            Wr.append(wr); Lr.append(lr + qs)
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
    N, K = disc.N, disc.K
    ctx.prec = prec
    ys = arb_nodes01(N + 1); vs = arb_nodes01(K + 1)
    cs = arb(CS_q); a2 = 2 * arb(A_MAP)
    lu = lu_factor(disc.jac(wf, qf(CS_q)))
    w = [arb(float(x)) for x in wf]
    hist = []
    target = 2.0 ** (-(prec - 40))
    for it in range(maxit):
        p, q = split(w, N, K)
        W, LW = Fields(p, q).grid(ys, vs)
        R = []
        for i in range(N + 1):
            y = ys[i]; om = 1 - y
            c2 = a2 * om * om; y2 = y * y
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
    return w, hist


def to_dyadic(x, bits):
    m, e = x.mid().man_exp()
    m = int(m); e = int(e)
    if e >= -bits:
        return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))
    sh = -bits - e
    return fmpq((m + (1 << (sh - 1))) >> sh, 2 ** bits)


# =============================================================================================
# 3. certificate engine
# =============================================================================================
_COS = {}


def dct_mat(N):
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
    pi = arb.pi(); mid = arb((lo + hi) / 2); half = arb((hi - lo) / 2)
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
    def __init__(self, wq, N, K):
        self.wq = wq; self.N = N; self.K = K
        p, q = split(wq, N, K)
        self.F = Fields(p, q)
        self.W0 = sum((wq[n] * (1 if n % 2 == 0 else -1) for n in range(N + 1)), fmpq(0))


def E_grid(ans, ylo, yhi, vlo, vhi):
    N = ans.N + 7; M = ans.K + 1          # E has Y-degree N+6, v-degree K
    ys = nodes_on(ylo, yhi, N); vs = nodes_on(vlo, vhi, M)
    W, _ = ans.F.grid(ys, vs)
    a2 = 2 * arb(A_MAP)
    vals = arb_mat(N, M, [a2 * (1 - ys[i]) ** 2 * (1 + ys[i] ** 2 * W[i, j]) + ys[i] ** 2
                          for i in range(N) for j in range(M)])
    return vals, N, M


def E_lower(ans, ylo, yhi, max_boxes=400, rel=fmpq(9, 10)):
    """rigorous lower bound of E_N on [ylo,yhi] x [0,1]; bisect until each box's form bound is > 0 and
    >= rel * (its smallest node value).  Returns (lower bound or None, boxes)."""
    stack = [(ylo, yhi, fmpq(0), fmpq(1), 0)]
    nbox = 0; lbmin = None
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
            if lb > 0:
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


def piece_bound(ans, ylo, yhi, CS_q):
    N, K = ans.N, ans.K
    if CS_q == 0:
        dQy = N + 2; dQv = K
    else:
        dQy = 8 * N + 44; dQv = 8 * K
    NQ = dQy + 1; MQ = dQv + 1
    ys = nodes_on(ylo, yhi, NQ); vs = nodes_on(fmpq(0), fmpq(1), MQ)
    W, LW = ans.F.grid(ys, vs)
    cs = arb(CS_q); a2 = 2 * arb(A_MAP)
    flat = []
    for i in range(NQ):
        y = ys[i]; om = 1 - y; y2 = y * y
        c2 = a2 * om * om
        s0 = cs * y2 * y2 * y * om ** 5
        for jv in range(MQ):
            E = c2 * (1 + y2 * W[i, jv]) + y2
            E2 = E * E; E7 = E2 * E2 * E2 * E
            flat.append(E7 * LW[i, jv] + s0 * (1 - vs[jv]))
    Qv = arb_mat(NQ, MQ, flat)
    SQ = arb_upper_q(abs_sum(cheb_coeffs(Qv, NQ, MQ)))
    Elow, nb = E_lower(ans, ylo, yhi)
    if Elow is None or Elow <= 0:
        return None, dict(ylo=str(ylo), yhi=str(yhi), sup_Q_le="%.3e" % qf(SQ), Elow="FAIL")
    bound = SQ / Elow ** 7
    return (bound, Elow), dict(ylo=str(ylo), yhi=str(yhi), sup_Q_le="%.3e" % qf(SQ), Elow="%.6f" % qf(Elow),
                               E_boxes=nb, bound="%.3e" % qf(bound))


def round_up_dyadic(q, bits=40):
    if q <= 0:
        return fmpq(0)
    e = math.floor(math.log2(qf(q))) - bits
    sc = fmpq(2) ** (-e)
    v = q * sc
    r = fmpq(-((-v.p) // v.q)) / sc
    assert r >= q
    return r


def make_pieces(n):
    """n uniform Y-pieces, the end pieces split geometrically 3 more times toward Y = 0 and Y = 1."""
    h = fmpq(1, n)
    left = [h / 2 ** m for m in range(3, 0, -1)]
    right = [1 - h / 2 ** m for m in range(1, 4)]
    mid = [fmpq(i, n) for i in range(1, n)]
    return [fmpq(0)] + left + mid + right + [fmpq(1)]


def certify(wq, N, K, CS_q, n_pieces=8, prec=256):
    ctx.prec = prec
    tm = time.time()
    ans = Ansatz(wq, N, K)
    pieces = []; sup = fmpq(0); Emin = None; ok = True
    Ed = make_pieces(n_pieces)
    assert Ed[0] == 0 and Ed[-1] == 1 and all(Ed[i] < Ed[i + 1] for i in range(len(Ed) - 1))
    for i in range(len(Ed) - 1):
        bd, info = piece_bound(ans, Ed[i], Ed[i + 1], CS_q)
        pieces.append(info)
        if bd is None:
            ok = False
            continue
        sup = max(sup, bd[0])
        Emin = bd[1] if Emin is None or bd[1] < Emin else Emin
    need = sup * (1 + fmpq(1, 2 ** 20)) * fmpq(2, 3)
    eta = round_up_dyadic(need)
    strict = bool(ok and eta > 0 and fmpq(3, 2) * eta > sup)
    # E_- = E_N - 2a (1-Y)^2 Y^2 eta >= Emin - (2a/16) eta   ((1-Y) Y <= 1/4)
    Em_low = (Emin - 2 * A_MAP * eta / 16) if Emin is not None else None
    Em_ok = bool(Em_low is not None and Em_low > 0)
    return dict(sup_Gh_le=sup, eta=eta, W0=ans.W0, strict=strict, E_N_min_lower=Emin, E_minus_lower=Em_low,
                E_minus_positive=Em_ok, certified=bool(strict and Em_ok), pieces=pieces, all_pieces_ok=ok,
                secs=round(time.time() - tm, 1))


# =============================================================================================
# 4. exact-rational evaluation (independent of the arb engine)
# =============================================================================================
def cheb_eval_q(c, X):
    b1 = fmpq(0); b2 = fmpq(0)
    for ck in reversed(c[1:]):
        b1, b2 = 2 * X * b1 - b2 + ck, b1
    return X * b1 - b2 + c[0]


class ExactEval:
    def __init__(self, wq, N, K):
        p, q = split(wq, N, K)
        self.P = [p, ydiff(p), ydiff(ydiff(p))]
        self.Q = [[c, ydiff(c), ydiff(ydiff(c))] for c in q]
        self.T, self.A = tau_polys(K)
        self.K = K

    def at(self, y, v):
        """returns (W, Lin[W]) exactly."""
        X = 2 * y - 1
        W, W1, W2 = [cheb_eval_q(c, X) for c in self.P]
        h = y - y * y; h1 = 1 - 2 * y
        angA = fmpq(0)
        for k in range(1, self.K + 1):
            q, q1, q2 = [cheb_eval_q(c, X) for c in self.Q[k - 1]]
            tk = self.T[k](v)
            W += h * q * tk; W1 += (h1 * q + h * q1) * tk; W2 += (-2 * q + 2 * h1 * q1 + h * q2) * tk
            angA += q * self.A[k](v)
        L = -(6 * W + (6 * y - 3) * W1 - y * (1 - y) * W2) / 4 + angA
        return W, L


def sample_signs(wq, N, K, CS_q, eta, n_y=24, n_v=6):
    Ex = ExactEval(wq, N, K)
    a2 = 2 * A_MAP
    mx_plus = None; mn_minus = None; mx_abs = fmpq(0); npos = nneg = 0
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


def integral_route(wq, N, K, CS):
    """double: delta = 2 lim r u = 2 * 2a CS int_0^1 Y^5 (1-Y)^7 <(1-v) E^-7>_c dY  (Gauss's law)."""
    a = qf(A_MAP)
    w = np.array([qf(x) for x in wq]); n1 = N + 1
    xg, wg = np.polynomial.legendre.leggauss(600)
    yq = (xg + 1) / 2; wy = wg / 2
    cq = xg; wc = wg / 2; vq = cq ** 2
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


# =============================================================================================
# 6. drivers
# =============================================================================================
def solve_point(j_q, eps, N, K, prec):
    CS_q = CS_of(j_q, eps)
    scale = qf(CS_q) / qf(j_q) ** 2
    wc, res, rc = continuation(qf(j_q), scale)
    disc = Disc(N, K)
    w0 = prolong(wc, res[0], res[1], N, K)
    wf, rf, it = disc.newton(qf(CS_q), w0)
    wa, hist = refine(disc, wf, CS_q, prec)
    # dyadic grid relative to the coefficient scale (K1 has |W| ~ j^2 ~ 2^-128)
    e_sc = math.ceil(math.log2(max(abs(float(x.mid())) for x in wa)))
    wq = [to_dyadic(x, prec - 32 - e_sc) for x in wa]
    return wq, rf, hist, CS_q


def run_point(j_q, eps, R, label, n_pieces=8, keep=False):
    N, K, prec = R
    wq, rf, hist, CS_q = solve_point(j_q, eps, N, K, prec)
    c = certify(wq, N, K, CS_q, n_pieces=n_pieces, prec=prec)
    a = A_MAP
    lo = 2 * a * (c['W0'] - c['eta']); hi = 2 * a * (c['W0'] + c['eta'])
    n1 = N + 1
    sc = max(abs(qf(x)) for x in wq[:n1])
    tail_p = max(abs(qf(x)) for x in wq[n1 - 4:n1]) / sc
    tail_q = max(abs(qf(wq[k * n1 + n])) for k in range(1, K + 1) for n in range(n1 - 4, n1)) / sc if K else 0.0
    qK = max(abs(qf(wq[K * n1 + n])) for n in range(n1)) / sc if K else 0.0
    out = dict(label=label, j=str(j_q), source_eps=str(eps), N=N, K=K, prec_bits=prec, n_pieces=n_pieces,
               rel_tail_P_last4="%.2e" % tail_p, rel_tail_q_last4="%.2e" % tail_q, rel_max_qK="%.2e" % qK,
               fine_newton_residual="%.2e" % rf, refine_history=["%.1e" % h for h in hist],
               sup_Gh_N_le="%.4e" % qf(c['sup_Gh_le']), eta="%.4e" % qf(c['eta']), eta_exact=qstr(c['eta']),
               W0_exact=qstr(c['W0']), E_N_min_lower=("%.6f" % qf(c['E_N_min_lower'])) if c['E_N_min_lower'] else None,
               E_minus_positive=c['E_minus_positive'], strict=c['strict'], certified=c['certified'],
               pieces=c['pieces'], delta_lo=dec_floor(lo, 34), delta_hi=dec_ceil(hi, 34),
               delta_lo_exact=qstr(lo), delta_hi_exact=qstr(hi),
               rel_width="%.3e" % qf((hi - lo) / lo) if lo > 0 else "n/a",
               integral_route_double=repr(integral_route(wq, N, K, qf(CS_q))), secs_cert=c['secs'])
    log("%s N=%d K=%d prec=%d: sup|Gh|<=%.3e eta %.3e Emin>=%s rel.w %s | delta in [%s, %s] cert=%s (%.1fs) tails %s %s"
        % (label, N, K, prec, qf(c['sup_Gh_le']), qf(c['eta']), out['E_N_min_lower'], out['rel_width'],
           out['delta_lo'][:24], out['delta_hi'][:24], c['certified'], c['secs'], out['rel_tail_P_last4'],
           out['rel_tail_q_last4']))
    if keep:
        return out, (lo, hi), c, wq, CS_q
    return out, (lo, hi), c


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "probe":
    j = fmpq(int(sys.argv[2])); N = int(sys.argv[3]); K = int(sys.argv[4]); prec = int(sys.argv[5])
    npc = int(sys.argv[6])
    out, e, c = run_point(j, fmpq(0), (N, K, prec), "probe j=%s" % j, n_pieces=npc)
    for pc in out['pieces']:
        print(pc)
    print("refine", out['refine_history'], "integral route", out['integral_route_double'])
    sys.exit(0)


# =============================================================================================
# 7. controls
# =============================================================================================
def identity_controls():
    import random
    import sympy as sp
    out = {}
    ctx.prec = 192
    rnd = random.Random(4)
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
    out['dct_recovers_known_coefficients'] = bool(all(C[k, l].contains(c[k][l] if (k < 10 and l < 6) else 0)
                                                      for k in range(Nn) for l in range(Mn)))
    # (b) angular operator: A_k(c^2) = [(1-c^2) g']' with g = tau_k(c^2) = T_2k(c)
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
    # (c) radial operator, symbolic: for u = Y^2 W(Y) with Y = 1/(1+sqrt(r/a)), (r^2 u_r)_r = (1-Y) Y^3 LR[W]
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
        rv = a * sv ** 2
        good = good and sp.simplify((lhs - rhs).subs(r, rv)) == 0
    # and Y^2 is strictly superharmonic: (r^2 (Y^2)_r)_r = -(3/2)(1-Y)Y^3
    g2 = Yr ** 2
    lg = sp.diff(r ** 2 * sp.diff(g2, r), r)
    for sv in (sp.Rational(1, 3), sp.Rational(5, 2)):
        good = good and sp.simplify((lg + sp.Rational(3, 2) * (1 - Yr) * Yr ** 3).subs(r, a * sv ** 2)) == 0
    out['radial_operator_symbolic'] = bool(good)
    # (d) source transformation, exact rationals: r^2 (9/4) J^2 (1-v) / (r^6 Phi^7) = (1-Y)Y^3 CS (1-v) Y^5(1-Y)^5/E^7
    F = Fraction; aF = F(1, 2)
    good = True
    for (Y, v, W, J) in ((F(1, 3), F(1, 5), F(2, 7), F(5, 2)), (F(4, 5), F(0), F(1, 9), F(2500)),
                         (F(1, 50), F(1, 2), F(3), F(750)), (F(9, 10), F(3, 4), F(-1, 5), F(1, 3))):
        rr = aF * ((1 - Y) / Y) ** 2
        u_ = Y * Y * W
        Phi = 1 + 1 / (2 * rr) + u_
        E = 2 * aF * (1 - Y) ** 2 * (1 + Y * Y * W) + Y * Y
        j = 4 * J
        CS = 18 * aF ** 3 * j * j
        good = good and (E == 2 * aF * (1 - Y) ** 2 * Phi)
        good = good and (rr * rr * F(9, 4) * J * J * (1 - v) / (rr ** 6 * Phi ** 7)
                         == (1 - Y) * Y ** 3 * CS * (1 - v) * Y ** 5 * (1 - Y) ** 5 / E ** 7)
        # mass normalisation: r Y^2 = a (1-Y)^2  ->  lim r u = a W(0)
        good = good and (rr * Y * Y == aF * (1 - Y) ** 2)
    out['source_transformation_and_rY2'] = bool(good)
    # (e) Lin[1] = -3/2 through the exact evaluator
    N, K = 6, 2
    wq = [fmpq(0)] * ((N + 1) * (K + 1)); wq[0] = fmpq(1)
    Ex = ExactEval(wq, N, K)
    out['Lin_of_1_is_minus_3_2'] = bool(all(Ex.at(fmpq(i, 7), fmpq(l, 3))[1] == fmpq(-3, 2)
                                            for i in range(8) for l in range(4)))
    # (f) arb field evaluator contains the exact rational evaluation
    N, K = 24, 3
    wq = [fmpq(rnd.randint(-2 ** 20, 2 ** 20), 2 ** 24) for _ in range((N + 1) * (K + 1))]
    Ex = ExactEval(wq, N, K); An = Ansatz(wq, N, K)
    good = True
    for (y, v) in ((fmpq(1, 3), fmpq(2, 7)), (fmpq(19, 20), fmpq(1, 2)), (fmpq(1, 100), fmpq(1))):
        we, Le = Ex.at(y, v)
        W, LW = An.F.grid([arb(y)], [arb(v)])
        good = good and W[0, 0].contains(we) and LW[0, 0].contains(Le)
    out['arb_evaluator_contains_exact'] = bool(good)
    # (g) double operator agrees with the exact evaluator (catches a slip in either)
    d = Disc(N, K)
    wf = np.array([qf(x) for x in wq])
    LWd = d.L @ wf
    mx = 0.0
    for idx in range(0, len(LWd), 7):
        y = d.Y[idx]; v = d.V[idx]
        W, LW = An.F.grid([arb(y)], [arb(v)])
        mx = max(mx, abs(float(LW[0, 0].mid()) - LWd[idx]))
    out['double_operator_matches_arb'] = bool(mx < 1e-9)
    out['all_ok'] = all(out.values())
    log("identity controls:", out)
    return out


def toy_J0():
    """J = 0: u = 0 is the exact solution (delta = 0)."""
    N, K, prec = 24, 3, 128
    disc = Disc(N, K)
    wf, rf, it = disc.newton(0.0, np.zeros(disc.L.shape[0]))
    wa, hist = refine(disc, wf, fmpq(0), prec)
    wq = [to_dyadic(x, prec - 32) for x in wa]
    c = certify(wq, N, K, fmpq(0), n_pieces=2, prec=prec)
    allzero = all(x == 0 for x in wq)
    delta_lo = 2 * A_MAP * (c['W0'] - c['eta']); delta_hi = 2 * A_MAP * (c['W0'] + c['eta'])
    out = dict(all_coefficients_zero=bool(allzero), sup_Gh_bound=qstr(c['sup_Gh_le']), eta_exact=qstr(c['eta']),
               delta_lo=qstr(delta_lo), delta_hi=qstr(delta_hi), strict_certificate=c['strict'],
               u_is_exact_solution=bool(allzero and c['sup_Gh_le'] == 0),
               ok=bool(allzero and c['sup_Gh_le'] == 0 and delta_lo == 0 and delta_hi == 0))
    log("toy J=0:", out)
    return out


def k1_point(R, j_q):
    out, (lo, hi), c = run_point(j_q, fmpq(0), R, "K1 j=%s" % j_q, n_pieces=4)
    t = j_q * j_q
    l2, h2 = lo / t, hi / t
    out['delta_over_j2_lo'] = dec_floor(l2, 34); out['delta_over_j2_hi'] = dec_ceil(h2, 34)
    out['delta_over_j2_lo_exact'] = qstr(l2); out['delta_over_j2_hi_exact'] = qstr(h2)
    out['contains_1_over_40'] = bool(c['certified'] and l2 <= fmpq(1, 40) <= h2)
    log("K1 j=%s: delta/j^2 in [%s, %s] contains 1/40: %s" % (j_q, out['delta_over_j2_lo'], out['delta_over_j2_hi'],
                                                            out['contains_1_over_40']))
    return out


def nulls_stage(R):
    j = fmpq(10000)
    out, e, c, wq, CS_q = run_point(j, fmpq(0), R, "j=10000 (sample/nulls)", keep=True)
    N, K, prec = R
    res = {}
    smp = sample_signs(wq, N, K, CS_q, c['eta'])
    smp['rigorous_bound_sup_Gh'] = "%.4e" % qf(c['sup_Gh_le'])
    res['exact_sample_j10000'] = smp
    log("exact rational sample j=10000:", smp)
    # null 1: bumped ansatz (P_0 coefficient + 1e-9 relative) must violate the honest-eta signs and inflate eta
    wb = list(wq); wb[0] += wq[0] * fmpq(1, 10 ** 9)
    smpb = sample_signs(wb, N, K, CS_q, c['eta'])
    cb = certify(wb, N, K, CS_q, n_pieces=8, prec=prec)
    res['null_bumped_ansatz'] = dict(bump_rel="1e-9", honest_eta="%.4e" % qf(c['eta']),
                                     sample_with_honest_eta_signs_ok=smpb['signs_ok'],
                                     rigorous_eta_bumped="%.4e" % qf(cb['eta']),
                                     eta_inflation_factor="%.3e" % (qf(cb['eta']) / qf(c['eta'])),
                                     ok=bool((not smpb['signs_ok']) and qf(cb['eta']) > 100 * qf(c['eta'])))
    log("null bump:", res['null_bumped_ansatz'])
    # null 2: the solution at the WRONG source constant (K2's 1e-3) certified against the honest source must fail
    wq2, rf2, h2, CS2 = solve_point(j, fmpq(1, 1000), N, K, prec)
    smp2 = sample_signs(wq2, N, K, CS_q, c['eta'])
    res['null_wrong_source_solution'] = dict(sample_with_honest_eta_signs_ok=smp2['signs_ok'],
                                             ok=bool(not smp2['signs_ok']))
    log("null wrong source:", res['null_wrong_source_solution'])
    return res


def summary(lo, hi):
    return dict(delta_lo=dec_floor(lo, 34), delta_hi=dec_ceil(hi, 34), delta_lo_exact=qstr(lo), delta_hi_exact=qstr(hi),
                abs_width="%.3e" % qf(hi - lo), rel_width="%.3e" % qf((hi - lo) / lo),
                M_ADM_lo=dec_floor(1 + lo, 34), M_ADM_hi=dec_ceil(1 + hi, 34))


RES_A = (140, 16, 224)
RES_B = (200, 24, 288)
STAGE_FILE = os.path.join(HERE, "rep4_stages.json")
BUDGET = float(os.environ.get("REP4_BUDGET", "300"))


def sha_self():
    import hashlib
    return hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()


def main():
    stages = {}
    me = sha_self()
    if os.path.exists(STAGE_FILE):
        st = json.load(open(STAGE_FILE))
        if st.get('instrument_sha256') == me:
            stages = st['stages']
        else:
            log("stage cache from a different instrument version: ignored")
    order = ['identity', 'toy', 'K1_A', 'K1_B', 'j3000_A', 'j3000_B', 'j10000_A', 'j10000_B', 'K2_A', 'K2_B', 'nulls']
    for name in order:
        if name in stages:
            continue
        if time.time() - T0 > BUDGET:
            json.dump(dict(instrument_sha256=me, stages=stages), open(STAGE_FILE, 'w'), indent=1, sort_keys=True)
            log("time budget reached before stage %s: saved; exit 3 (rerun to resume)" % name)
            sys.exit(3)
        R = RES_A if name.endswith('A') else RES_B
        if name == 'identity':
            r = identity_controls()
        elif name == 'toy':
            r = toy_J0()
        elif name.startswith('K1'):
            r = k1_point((40, 4, 192) if name.endswith('A') else (56, 6, 256), fmpq(1, 2 ** 64))
        elif name.startswith('j'):
            r = run_point(fmpq(int(name[1:].split('_')[0])), fmpq(0), R, name)[0]
        elif name.startswith('K2'):
            r = run_point(fmpq(10000), fmpq(1, 1000), R, name)[0]
        elif name == 'nulls':
            r = nulls_stage(RES_A)
        stages[name] = r
        json.dump(dict(instrument_sha256=me, stages=stages), open(STAGE_FILE, 'w'), indent=1, sort_keys=True)

    res = dict(seat="rep4", instrument_sha256=me,
               object="Bowen-York spinning puncture, maximal conformally flat slice, mu = 1, j = 4J, "
                      "delta = M_ADM - 1 = 2 lim r u",
               method=("Y = 1/(1+sqrt(2r)) (square-root map, a = 1/2), u = Y^2 W, W = P(Y) + Y(1-Y) sum q_k(Y) "
                       "T_k(2cos^2 theta - 1); Chebyshev in 2Y-1; Gauss collocation + continuation in j + arb "
                       "refinement; dyadic coefficients; certificate: Q = E^7 Gh polynomial sampled in arb on "
                       "exact-degree Chebyshev grids per Y-piece, sup|Q| <= sum|c|, inf E by exact-degree Chebyshev "
                       "form with bisection; eta = (2/3) sup|Gh| (1+2^-20); comparison with u +- eta Y^2 "
                       "(Y^2 strictly superharmonic); delta in 2a [W(0) -+ eta]"))
    res['identity_controls'] = stages['identity']
    res['toy_J0'] = stages['toy']
    k1A, k1B = stages['K1_A'], stages['K1_B']
    l2 = max(parse_q(k1A['delta_over_j2_lo_exact']), parse_q(k1B['delta_over_j2_lo_exact']))
    h2 = min(parse_q(k1A['delta_over_j2_hi_exact']), parse_q(k1B['delta_over_j2_hi_exact']))
    res['K1'] = dict(j="2^-64", runs=[k1A, k1B], delta_over_j2_lo=dec_floor(l2, 34), delta_over_j2_hi=dec_ceil(h2, 34),
                     contains_1_over_40=bool(k1A['contains_1_over_40'] and k1B['contains_1_over_40']
                                             and l2 <= fmpq(1, 40) <= h2))
    for j in (3000, 10000):
        rA, rB = stages['j%d_A' % j], stages['j%d_B' % j]
        eA = (parse_q(rA['delta_lo_exact']), parse_q(rA['delta_hi_exact']))
        eB = (parse_q(rB['delta_lo_exact']), parse_q(rB['delta_hi_exact']))
        lo, hi = max(eA[0], eB[0]), min(eA[1], eB[1])
        s = summary(lo, hi)
        s['all_certified'] = bool(rA['certified'] and rB['certified'])
        s['A_contains_B'] = bool(eA[0] <= eB[0] and eB[1] <= eA[1])
        s['resolutions_overlap'] = bool(lo <= hi)
        s['rel_width_A'] = rA['rel_width']; s['rel_width_B'] = rB['rel_width']
        res['delta_j%d' % j] = s
        res['runs_j%d' % j] = [rA, rB]
        log("FINAL j=%d: delta in [%s, %s] rel width %s certified %s A>B %s" % (j, s['delta_lo'], s['delta_hi'],
                                                                              s['rel_width'], s['all_certified'],
                                                                              s['A_contains_B']))
    kA, kB = stages['K2_A'], stages['K2_B']
    pA = (parse_q(kA['delta_lo_exact']), parse_q(kA['delta_hi_exact']))
    pB = (parse_q(kB['delta_lo_exact']), parse_q(kB['delta_hi_exact']))
    plo, phi = max(pA[0], pB[0]), min(pA[1], pB[1])
    blo = parse_q(res['delta_j10000']['delta_lo_exact']); bhi = parse_q(res['delta_j10000']['delta_hi_exact'])
    res['K2'] = dict(perturbation="source 9/4 -> 9/4 (1 + 1e-3) at j = 10000", runs=[kA, kB],
                     all_certified=bool(kA['certified'] and kB['certified']),
                     delta_perturbed_lo=dec_floor(plo, 34), delta_perturbed_hi=dec_ceil(phi, 34),
                     delta_perturbed_rel_width="%.3e" % qf((phi - plo) / plo),
                     shift_lo=dec_floor(plo - bhi, 22), shift_hi=dec_ceil(phi - blo, 22),
                     rel_shift_lo=dec_floor((plo - bhi) / bhi, 16), rel_shift_hi=dec_ceil((phi - blo) / blo, 16),
                     A_contains_B=bool(pA[0] <= pB[0] and pB[1] <= pA[1]))
    log("K2: delta' in [%s, %s]; shift in [%s, %s]" % (res['K2']['delta_perturbed_lo'], res['K2']['delta_perturbed_hi'],
                                                      res['K2']['shift_lo'], res['K2']['shift_hi']))
    res['sample_and_nulls'] = stages['nulls']
    res['expected_false'] = [
        "sample_and_nulls/null_bumped_ansatz/sample_with_honest_eta_signs_ok",
        "sample_and_nulls/null_wrong_source_solution/sample_with_honest_eta_signs_ok",
        "toy_J0/strict_certificate",
    ]
    res['all_flags_ok'] = bool(res['identity_controls']['all_ok'] and res['toy_J0']['ok'] and res['K1']['contains_1_over_40']
                               and res['delta_j3000']['all_certified'] and res['delta_j10000']['all_certified']
                               and res['delta_j3000']['A_contains_B'] and res['delta_j10000']['A_contains_B']
                               and res['K2']['all_certified'] and res['sample_and_nulls']['null_bumped_ansatz']['ok']
                               and res['sample_and_nulls']['null_wrong_source_solution']['ok']
                               and res['sample_and_nulls']['exact_sample_j10000']['signs_ok'])
    res['total_seconds_last_invocation'] = round(time.time() - T0, 1)
    with open(os.path.join(HERE, "result_rep4.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    log("wrote result_rep4.json; all_flags_ok =", res['all_flags_ok'])


if __name__ == "__main__":
    main()
