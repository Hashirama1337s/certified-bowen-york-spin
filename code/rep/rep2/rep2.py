"""
rep2.py - blind replication seat REP2 for BHSLICE (Bowen-York spinning puncture, bare mass mu = 1, large spin).

Object (DESIGN.md 1.1-1.3):  Delta u = -(9/4) J^2 sin^2(theta) / (r^6 (1 + 1/(2r) + u)^7),  u -> 0 at infinity,
                             delta := M_ADM - 1 = 2 lim r u,  j := 4 J.

Formulation of THIS seat (derived in REP2-NOTES.md):
  * radial y = 1/(1 + r) (scale mu), Chebyshev T_n(2y - 1) on [0, 1];
  * angular variable v = cos^2(theta) in [0, 1]; basis tau_k(v) = T_k(2v - 1) = cos(2 k theta) (Fourier cosine);
  * u = t w, t := j^2,  w = y p(y) + y^2 (1 - y) sum_{k=1..K} q_k(y) tau_k(v);
  * divided operator  Gh[w] = (1-y) w_yy - 2 w_y + Ang[w]/(y^2 (1-y)) + s1 (1-v) y^2 (1-y)^2 D^-7,
        Ang = 4 v (1-v) w_vv + (2 - 6 v) w_v,  D = 1 - y/2 + t (1-y) w,  s1 = 9/64;
    r^2 (Delta u - F) = t y^2 (1-y) Gh[w]   (positive prefactor on 0 < y < 1);
  * solver: Newton collocation on Chebyshev-Gauss (interior) nodes in y x interior equispaced theta nodes, continuation
    in j, then iterative refinement with arb residuals;
  * certificate: exact rationals (python-flint fmpq_poly, Kronecker packing in (y, v)).
        Q := D_N^7 Lh w_N + s1 (1-v) y^2 (1-y)^2   (exact polynomial),  Gh[w_N] = Q / D_N^7,
        S := sum |shifted-Chebyshev coefficients of Q on [0,1]^2|  >= sup |Q|  (exact rational),
        eta := 64 S (1 + 2^-20) rounded up,  w+- := w_N +- eta y;
    positivity fix: certify V_- := w_-/y = p(y) + y(1-y) sum q_k tau_k(v) - eta > 0 on [0,1]^2 by adaptive
    subdivision with exact Taylor (centred-form) lower bounds.  Then u_- >= 0, D_N >= D_- >= 1 - y/2 >= 1/2,
    so sup|Gh[w_N]| <= 128 S < 2 eta, Gh[w+] <= Gh[w_N] - 2 eta < 0 < Gh[w_N] + 2 eta <= Gh[w-],
    Phi_+- >= 1 + 1/(2r) (kappa = 1), and Lemma C gives u_- <= u <= u_+, so delta in 2 t [p(0) - eta, p(0) + eta].

Run:  python rep2.py            (all: controls + targets; writes rep2_result.json)
      python rep2.py probe J Ny K   (double solve only, prints tails)
Never pickles flint objects; all exact numbers are written as "p/q" strings or outward-rounded decimals.
"""
import json, math, os, sys, time
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat, fmpq, fmpz, fmpq_poly, fmpz_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()
S1 = fmpq(9, 64)


def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)


def qf(q):
    return float(Fraction(int(q.p), int(q.q)))


# =============================================================================================
# 1. collocation (double) with Chebyshev-Gauss nodes
# =============================================================================================
def cheb_tab(x, N):
    """T_n, T_n', T_n'' (d/dx) at points x (numpy), n = 0..N."""
    x = np.asarray(x, dtype=float)
    T = np.zeros((N + 1, x.size)); d1 = np.zeros_like(T); d2 = np.zeros_like(T)
    T[0] = 1.0
    if N >= 1:
        T[1] = x; d1[1] = 1.0
    for n in range(1, N):
        T[n + 1] = 2 * x * T[n] - T[n - 1]
        d1[n + 1] = 2 * T[n] + 2 * x * d1[n] - d1[n - 1]
        d2[n + 1] = 4 * d1[n] + 2 * x * d2[n] - d2[n - 1]
    return T, d1, d2


def nodes(Ny, K):
    i = np.arange(Ny + 1)
    y = (1 + np.cos(np.pi * (i + 0.5) / (Ny + 1))) / 2          # Chebyshev-Gauss roots mapped to (0,1)
    jj = np.arange(K + 1)
    th = np.pi * (jj + 0.5) / (2 * (K + 1))                      # interior points of (0, pi/2)
    v = np.cos(th) ** 2
    return y, v


def radial_blocks(y, Ny):
    """matrices (point x n): RF = Lrad(y T_n), FV = y T_n, RG = Lrad(h T_n), GV = h T_n, TY = T_n."""
    T, d1, d2 = cheb_tab(2 * y - 1, Ny)
    T = T.T; Ty1 = 2 * d1.T; Ty2 = 4 * d2.T
    yy = y[:, None]; om = 1 - yy
    f1 = T + yy * Ty1; f2 = 2 * Ty1 + yy * Ty2
    RF = om * f2 - 2 * f1
    FV = yy * T
    h = yy ** 2 - yy ** 3; h1 = 2 * yy - 3 * yy ** 2; h2 = 2 - 6 * yy
    g1 = h1 * T + h * Ty1; g2 = h2 * T + 2 * h1 * Ty1 + h * Ty2
    RG = om * g2 - 2 * g1
    GV = h * T
    return RF, FV, RG, GV, T


def angular_blocks(v, K):
    """TK[j,k] = tau_k(v_j) (k = 1..K), AK[j,k] = (4v(1-v) d2 + (2-6v) d1) tau_k at v_j."""
    T, d1, d2 = cheb_tab(2 * v - 1, K)
    T = T.T[:, 1:]; t1 = 2 * d1.T[:, 1:]; t2 = 4 * d2.T[:, 1:]
    vv = v[:, None]
    A = 4 * vv * (1 - vv) * t2 + (2 - 6 * vv) * t1
    return T, A


class Disc:
    def __init__(self, Ny, K):
        self.Ny, self.K = Ny, K
        y, v = nodes(Ny, K)
        self.y, self.v = y, v
        RF, FV, RG, GV, TY = radial_blocks(y, Ny)
        TK, AK = angular_blocks(v, K)
        n1 = Ny + 1
        npt = n1 * (K + 1)
        nun = n1 * (K + 1)
        L = np.zeros((npt, nun)); U = np.zeros((npt, nun))
        # point index = i*(K+1) + j  (y_i, v_j);  unknown index: p_n -> n, q_{k,n} -> k*n1 + n
        for i in range(n1):
            for j in range(K + 1):
                r = i * (K + 1) + j
                L[r, :n1] = RF[i]; U[r, :n1] = FV[i]
                for k in range(1, K + 1):
                    L[r, k * n1:(k + 1) * n1] = RG[i] * TK[j, k - 1] + TY[i] * AK[j, k - 1]
                    U[r, k * n1:(k + 1) * n1] = GV[i] * TK[j, k - 1]
        self.L, self.U = L, U
        Y = np.repeat(y, K + 1); V = np.tile(v, n1)
        self.Y, self.V = Y, V
        self.sh = (1 - V) * Y ** 2 * (1 - Y) ** 2

    def F(self, w, t, s1):
        D = 1 - self.Y / 2 + t * (1 - self.Y) * (self.U @ w)
        return self.L @ w + s1 * self.sh * D ** -7, D

    def jac(self, w, t, s1):
        D = 1 - self.Y / 2 + t * (1 - self.Y) * (self.U @ w)
        return self.L - (7 * s1 * t * self.sh * (1 - self.Y) * D ** -8)[:, None] * self.U

    def newton(self, t, s1, w0=None, iters=60, tol=1e-14):
        w = np.zeros(self.L.shape[1]) if w0 is None else w0.copy()
        for it in range(iters):
            F, D = self.F(w, t, s1)
            dw = np.linalg.solve(self.jac(w, t, s1), F)
            w -= dw
            if np.max(np.abs(dw)) < tol * max(1.0, np.max(np.abs(w))):
                break
        F, D = self.F(w, t, s1)
        return w, float(np.max(np.abs(F))), it + 1


def solve_continuation(disc, j_target, s1, ratio=1.15):
    """continuation in j from the linear (t = 0) problem; returns double solution at j_target."""
    w, res, _ = disc.newton(0.0, s1)
    if j_target == 0:
        return w, res
    js = []
    jj = min(0.5, j_target)
    while jj < j_target:
        js.append(jj); jj *= ratio
    js.append(j_target)
    for jv in js:
        w, res, it = disc.newton(jv * jv, s1, w)
    return w, res


# =============================================================================================
# 2. arb refinement (separable residual evaluation)
# =============================================================================================
def arb_blocks(Ny, K, prec):
    ctx.prec = prec
    one = arb(1)
    pi = arb.pi()
    ys = [(one + (pi * (2 * i + 1) / (2 * (Ny + 1))).cos()) / 2 for i in range(Ny + 1)]
    vs = [((pi * (2 * j + 1) / (4 * (K + 1))).cos()) ** 2 for j in range(K + 1)]

    def tab(x, N):
        T = [one, x]; d1 = [0 * one, one]; d2 = [0 * one, 0 * one]
        for n in range(1, N):
            T.append(2 * x * T[n] - T[n - 1])
            d1.append(2 * T[n] + 2 * x * d1[n] - d1[n - 1])
            d2.append(4 * d1[n] + 2 * x * d2[n] - d2[n - 1])
        return T[:N + 1], d1[:N + 1], d2[:N + 1]

    RF, FV, RG, GV, TY = [], [], [], [], []
    for y in ys:
        T, a1, a2 = tab(2 * y - 1, Ny)
        T1 = [2 * a for a in a1]; T2 = [4 * a for a in a2]
        om = 1 - y
        RF.append([om * (2 * T1[n] + y * T2[n]) - 2 * (T[n] + y * T1[n]) for n in range(Ny + 1)])
        FV.append([y * T[n] for n in range(Ny + 1)])
        h = y * y - y * y * y; h1 = 2 * y - 3 * y * y; h2 = 2 - 6 * y
        RG.append([om * (h2 * T[n] + 2 * h1 * T1[n] + h * T2[n]) - 2 * (h1 * T[n] + h * T1[n]) for n in range(Ny + 1)])
        GV.append([h * T[n] for n in range(Ny + 1)])
        TY.append(list(T))
    TK, AK = [], []
    for v in vs:
        T, a1, a2 = tab(2 * v - 1, K)
        TK.append(T[1:])
        AK.append([4 * v * (1 - v) * 4 * a2[k] + (2 - 6 * v) * 2 * a1[k] for k in range(1, K + 1)])
    M = lambda rows: arb_mat(rows)
    return dict(RF=M(RF), FV=M(FV), RG=M(RG), GV=M(GV), TY=M(TY), TK=M(TK), AK=M(AK), ys=ys, vs=vs)


def arb_residual(B, w, Ny, K, t, s1):
    n1 = Ny + 1
    P = arb_mat([[w[n]] for n in range(n1)])
    Qm = arb_mat([[w[k * n1 + n] for k in range(1, K + 1)] for n in range(n1)])  # n1 x K
    Lp = B['RF'] * P                     # n1 x 1
    LG = B['RG'] * Qm * B['TK'].transpose() + B['TY'] * Qm * B['AK'].transpose()  # n1 x (K+1)
    Up = B['FV'] * P
    UG = B['GV'] * Qm * B['TK'].transpose()
    out = []
    for i in range(n1):
        y = B['ys'][i]
        for j in range(K + 1):
            v = B['vs'][j]
            wv = Up[i, 0] + UG[i, j]
            D = 1 - y / 2 + t * (1 - y) * wv
            out.append(Lp[i, 0] + LG[i, j] + s1 * (1 - v) * y * y * (1 - y) * (1 - y) / D ** 7)
    return out


def refine(disc, wf, t_q, s1_q, prec, maxit=20):
    Ny, K = disc.Ny, disc.K
    B = arb_blocks(Ny, K, prec)
    t = arb(t_q); s1 = arb(s1_q)
    w = [arb(float(x)) for x in wf]
    tf, s1f = qf(t_q), qf(s1_q)
    hist = []
    target = 2.0 ** (-(prec - 24))
    for it in range(maxit):
        R = arb_residual(B, w, Ny, K, t, s1)
        Rf = np.array([float(r.mid()) for r in R])
        res = float(np.max(np.abs(Rf)))
        hist.append(res)
        if res < target:
            break
        wd = np.array([float(x.mid()) for x in w])
        dw = np.linalg.solve(disc.jac(wd, tf, s1f), Rf)
        w = [(w[i] - arb(float(dw[i]))).mid() for i in range(len(w))]
    return w, hist


def to_dyadic(a, bits):
    """round the exact midpoint of an arb to the nearest multiple of 2^-bits (exact fmpq)."""
    m, e = a.mid().man_exp()
    m = int(m); e = int(e)
    if e >= -bits:
        return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))
    sh = -bits - e
    return fmpq((m + (1 << (sh - 1))) >> sh, 2 ** bits)


# =============================================================================================
# 3. exact certificate
# =============================================================================================
Yv = fmpq_poly([0, 1])
ONE = fmpq_poly([1])


def shifted_cheb(N):
    """T_n(2z - 1) as fmpq_poly, n = 0..N."""
    x = fmpq_poly([-1, 2])
    T = [ONE, x]
    for n in range(1, N):
        T.append(2 * x * T[n] - T[n - 1])
    return T[:N + 1]


def build_polys(wq, Ny, K):
    """returns w-terms, Lh-terms as lists of (poly_y, poly_v), and p(0), plus V-terms (w/y)."""
    n1 = Ny + 1
    Ty = shifted_cheb(Ny)
    Tv = shifted_cheb(K)
    pp = sum((wq[n] * Ty[n] for n in range(n1)), fmpq_poly([0]))
    h = fmpq_poly([0, 0, 1, -1])
    Rad = lambda f: (1 - Yv) * f.derivative().derivative() - 2 * f.derivative()
    f0 = Yv * pp
    w_terms = [(f0, ONE)]
    Lh_terms = [(Rad(f0), ONE)]
    V_terms = [(pp, ONE)]
    for k in range(1, K + 1):
        qk = sum((wq[k * n1 + n] * Ty[n] for n in range(n1)), fmpq_poly([0]))
        tk = Tv[k]
        Ak = 4 * Yv * (1 - Yv) * tk.derivative().derivative() + (2 - 6 * Yv) * tk.derivative()
        g = h * qk
        w_terms.append((g, tk))
        Lh_terms.append((Rad(g), tk))
        Lh_terms.append((qk, Ak))
        V_terms.append((Yv * (1 - Yv) * qk, tk))
    p0 = sum((wq[n] * (1 if n % 2 == 0 else -1) for n in range(n1)), fmpq(0))
    return w_terms, Lh_terms, V_terms, p0


def pack(terms, W):
    """Kronecker packing: y^i v^j -> z^(i W + j)."""
    dy = max(f.degree() for f, g in terms)
    co = [fmpq(0)] * ((dy + 1) * W)
    for f, g in terms:
        gc = g.coeffs()
        assert len(gc) <= W
        for i, fi in enumerate(f.coeffs()):
            if fi == 0:
                continue
            for jv, gj in enumerate(gc):
                if gj != 0:
                    co[i * W + jv] += fi * gj
    return fmpq_poly(co)


_E = {}


def E_int(d):
    """integer matrix E (d+1 x d+1): z^i = 2^-(2d) sum_k E[k][i] T_k(2z-1)   (z^i = cos^{2i}(phi/2) expansion)."""
    if d in _E:
        return _E[d]
    rows = [[0] * (d + 1) for _ in range(d + 1)]
    for i in range(d + 1):
        sc = 1 << (2 * (d - i))
        rows[0][i] = sc * math.comb(2 * i, i)
        for k in range(1, i + 1):
            rows[k][i] = sc * 2 * math.comb(2 * i, i - k)
    M = fmpz_mat(rows)
    _E[d] = M
    return M


def unpack_int(P, W):
    num = P.numer().coeffs()
    den = int(P.denom())
    L = len(num)
    dy = (L - 1) // W
    dv = 0
    flat = [int(c) for c in num] + [0] * ((dy + 1) * W - L)
    for i in range(dy + 1):
        for jv in range(W - 1, -1, -1):
            if flat[i * W + jv] != 0:
                dv = max(dv, jv); break
    M = fmpz_mat(dy + 1, dv + 1, [flat[i * W + jv] for i in range(dy + 1) for jv in range(dv + 1)])
    return M, den, dy, dv


def cheb_abs_sum(P, W):
    """exact S = sum |c_kl|, P(y,v) = sum c_kl T_k(2y-1) T_l(2v-1); sup_{[0,1]^2} |P| <= S."""
    if P == 0:
        return fmpq(0)
    M, den, dy, dv = unpack_int(P, W)
    C = E_int(dy) * M * E_int(dv).transpose()
    tot = 0
    for i in range(C.nrows()):
        for jv in range(C.ncols()):
            x = int(C[i, jv])
            tot += x if x >= 0 else -x
    return fmpq(tot, den * (1 << (2 * dy + 2 * dv)))


def round_up_dyadic(q, bits=64):
    if q <= 0:
        return fmpq(0)
    e = math.floor(math.log2(qf(q))) - bits
    sc = fmpq(2) ** (-e)
    v = q * sc
    r = fmpq(-((-v.p) // v.q)) / sc
    assert r >= q
    return r


# --- positivity of V(y,v) = sum f(y) g(v) - c on [0,1]^2 by centred Taylor forms with subdivision ---
def taylor_box(terms, y0, ry, v0, rv):
    """coefficient table a[a][b] of V(y0 + ry s, v0 + rv t) (monomials s^a t^b), exact; returns lower bound."""
    sy = fmpq_poly([y0, ry]); sv = fmpq_poly([v0, rv])
    const = fmpq(0); rest = fmpq(0)
    acc = {}
    for f, g in terms:
        fs = f(sy).coeffs(); gs = g(sv).coeffs()
        for a, fa in enumerate(fs):
            if fa == 0:
                continue
            for b, gb in enumerate(gs):
                if gb != 0:
                    acc[(a, b)] = acc.get((a, b), fmpq(0)) + fa * gb
    for (a, b), c in acc.items():
        if a == 0 and b == 0:
            const = c
        else:
            rest += abs(c)
    return const - rest, const


def certify_positive(terms, shift, max_boxes=20000):
    """certify sum f g - shift > 0 on [0,1]^2.  Returns (ok, n_boxes, min lower bound, min centre value)."""
    terms = list(terms) + [(fmpq_poly([-shift]), ONE)]
    stack = [(fmpq(1, 2), fmpq(1, 2), fmpq(1, 2), fmpq(1, 2))]
    nbox = 0
    minlb = None; mincen = None
    while stack:
        y0, ry, v0, rv = stack.pop()
        lb, cen = taylor_box(terms, y0, ry, v0, rv)
        nbox += 1
        mincen = cen if mincen is None or cen < mincen else mincen
        if cen <= 0:
            return False, nbox, lb, cen          # a point value <= 0: genuinely not positive
        if lb > 0:
            minlb = lb if minlb is None or lb < minlb else minlb
            continue
        if nbox > max_boxes:
            return False, nbox, lb, cen
        if ry >= rv:
            h = ry / 2
            stack += [(y0 - h, h, v0, rv), (y0 + h, h, v0, rv)]
        else:
            h = rv / 2
            stack += [(y0, ry, v0 - h, h), (y0, ry, v0 + h, h)]
    return True, nbox, minlb, mincen


def certify(wq, Ny, K, t, s1, t_slab=None, eta_override=None):
    """exact certificate. t_slab: uniform over t in [0, t_slab] (wq must be t-independent)."""
    tm = time.time()
    w_terms, Lh_terms, V_terms, p0 = build_polys(wq, Ny, K)
    W = 8 * K + 2
    LhP = pack(Lh_terms, W)
    ShP = pack([(s1 * Yv ** 2 * (1 - Yv) ** 2, fmpq_poly([1, -1]))], W)
    D0P = pack([(fmpq_poly([1, fmpq(-1, 2)]), ONE)], W)
    EP = pack([((1 - Yv) * f, g) for f, g in w_terms], W)
    Sparts = None
    if t_slab is None:
        DP = D0P + t * EP
        Q = DP ** 7 * LhP + ShP
        S = cheb_abs_sum(Q, W)
        Qdeg = Q.degree()
    else:
        S = fmpq(0); Sparts = []
        for m in range(8):
            Qm = (D0P ** 7 * LhP + ShP) if m == 0 else math.comb(7, m) * D0P ** (7 - m) * EP ** m * LhP
            Sm = cheb_abs_sum(Qm, W)
            Sparts.append(Sm)
            S += Sm * t_slab ** m
        Qdeg = None
    tS = time.time() - tm
    if eta_override is not None:
        eta = eta_override
    else:
        eta = round_up_dyadic(64 * S * (1 + fmpq(1, 2 ** 20)))
    strict = eta > 0 and eta >= 64 * S * (1 + fmpq(1, 2 ** 20)) and eta_override is None
    pos_ok, nbox, minlb, mincen = certify_positive(V_terms, eta)
    return dict(S=S, eta=eta, p0=p0, pos_ok=pos_ok, pos_boxes=nbox, pos_minlb=minlb, pos_mincen=mincen,
                strict=bool(strict), certified=bool(strict and pos_ok), Sparts=Sparts,
                secs_S=round(tS, 2), secs=round(time.time() - tm, 2), Qdeg=Qdeg)


# =============================================================================================
# 4. output helpers (outward rounding of exact rationals)
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


def qstr(q):
    return "%s/%s" % (q.p, q.q)


# =============================================================================================
# 5. double cross-check: integral route  delta = 2 (9 J^2/4) int y^2 (1-y)^3 <(1-v)... >  (non-rigorous)
# =============================================================================================
def integral_route(w, Ny, K, t, s1):
    xg, wg = np.polynomial.legendre.leggauss(160)
    yq = (xg + 1) / 2; wy = wg / 2
    cq = xg; wc = wg / 2                     # c in [-1,1], <.>_c = (1/2) int dc
    vq = cq ** 2
    T, _, _ = cheb_tab(2 * yq - 1, Ny)
    Tv, _, _ = cheb_tab(2 * vq - 1, K)
    n1 = Ny + 1
    p = w[:n1] @ T
    tot = 0.0
    for iy in range(len(yq)):
        y = yq[iy]
        wv = y * p[iy] + y * y * (1 - y) * sum((w[k * n1:(k + 1) * n1] @ T[:, iy]) * Tv[k] for k in range(1, K + 1))
        D = 1 - y / 2 + t * (1 - y) * wv
        tot += wy[iy] * y ** 2 * (1 - y) ** 3 * np.sum(wc * (1 - vq) * D ** -7)
    return 2 * s1 * t * tot


# =============================================================================================
# 6. drivers
# =============================================================================================
def run_point(j_q, Ny, K, prec, s1=S1, label="", bump=None, keep=False):
    t_q = j_q * j_q
    disc = Disc(Ny, K)
    wf, resf = solve_continuation(disc, qf(j_q), qf(s1))
    wa, hist = refine(disc, wf, t_q, s1, prec)
    wq = [to_dyadic(a, prec - 24) for a in wa]
    if bump is not None:
        wq[0] += bump
    c = certify(wq, Ny, K, t_q, s1)
    lo = 2 * t_q * (c['p0'] - c['eta'])
    hi = 2 * t_q * (c['p0'] + c['eta'])
    n1 = Ny + 1
    tail_p = max(abs(qf(x)) for x in wq[n1 - 4:n1])
    tail_q = max(abs(qf(wq[K * n1 + n])) for n in range(n1))
    integ = integral_route(np.array([qf(x) for x in wq]), Ny, K, qf(t_q), qf(s1))
    out = dict(label=label, j=str(j_q), J=str(j_q / 4), Ny=Ny, K=K, prec_bits=prec, float_newton_residual=resf,
               refine_history=["%.2e" % h for h in hist], tail_p_last4=tail_p, tail_qK=tail_q,
               S="%.6e" % qf(c['S']), eta="%.6e" % qf(c['eta']), eta_exact=qstr(c['eta']),
               p0_exact=qstr(c['p0']), positivity_certified=c['pos_ok'], positivity_boxes=c['pos_boxes'],
               positivity_min_lower_bound=("%.6e" % qf(c['pos_minlb'])) if c['pos_minlb'] is not None else None,
               strict=c['strict'], certified=c['certified'],
               delta_lo=dec_floor(lo, 34), delta_hi=dec_ceil(hi, 34),
               rel_width="%.3e" % qf((hi - lo) / lo) if lo > 0 else "n/a",
               integral_route_double=repr(integ), secs_S=c['secs_S'], secs_cert=c['secs'], Qdeg=c['Qdeg'])
    log("%s j=%s Ny=%d K=%d prec=%d: tails p %.1e qK %.1e | eta %.3e pos=%s(%d boxes) rel.w %s | delta in [%s, %s] cert=%s (S %.1fs, cert %.1fs)"
        % (label, j_q, Ny, K, prec, tail_p, tail_q, qf(c['eta']), c['pos_ok'], c['pos_boxes'], out['rel_width'],
           out['delta_lo'][:24], out['delta_hi'][:24], c['certified'], c['secs_S'], c['secs']))
    if keep:
        return out, (lo, hi), (wq, c)
    return out, (lo, hi)


def probe(j, Ny, K):
    disc = Disc(Ny, K)
    t1 = time.time()
    w, res = solve_continuation(disc, j, 9 / 64)
    n1 = Ny + 1
    p = w[:n1]
    print("j=%g Ny=%d K=%d newton res %.2e (%.1fs)" % (j, Ny, K, res, time.time() - t1))
    print(" |p_n| last 6:", " ".join("%.1e" % abs(x) for x in p[-6:]))
    for k in range(1, K + 1):
        q = w[k * n1:(k + 1) * n1]
        print(" k=%2d max|q| %.2e  tail %.1e" % (k, np.max(np.abs(q)), np.max(np.abs(q[-4:]))))
    print(" delta(double, slope) = %.15g  integral = %.15g" % (2 * j * j * np.sum(p * (-1.0) ** np.arange(n1)),
                                                          integral_route(w, Ny, K, j * j, 9 / 64)))


# =============================================================================================
# 7. controls
# =============================================================================================
def eval_terms(terms, y, v):
    return sum((f(y) * g(v) for f, g in terms), fmpq(0))


def sample_signs(wq, Ny, K, t, s1, eta, n_y=40, n_v=10):
    """exact rational evaluation of Gh[w+], Gh[w-], Gh[w_N] on a rational grid incl. y = 0, 1 and v = 0, 1."""
    w_terms, Lh_terms, V_terms, p0 = build_polys(wq, Ny, K)
    mx_plus = None; mn_minus = None; mx_abs = fmpq(0); npos = 0; nneg = 0
    for i in range(n_y + 1):
        y = fmpq(i, n_y)
        for l in range(n_v + 1):
            v = fmpq(l, n_v)
            Lw = eval_terms(Lh_terms, y, v)
            w = eval_terms(w_terms, y, v)
            src = s1 * (1 - v) * y * y * (1 - y) * (1 - y)
            DN = 1 - y / 2 + t * (1 - y) * w
            Gn = Lw + src / DN ** 7
            Gp = Lw - 2 * eta + src / (DN + t * (1 - y) * eta * y) ** 7
            Gm = Lw + 2 * eta + src / (DN - t * (1 - y) * eta * y) ** 7
            mx_plus = Gp if mx_plus is None or Gp > mx_plus else mx_plus
            mn_minus = Gm if mn_minus is None or Gm < mn_minus else mn_minus
            mx_abs = max(mx_abs, abs(Gn))
            npos += int(Gn > 0); nneg += int(Gn < 0)
    return dict(max_Gh_plus=qf(mx_plus), min_Gh_minus=qf(mn_minus), max_abs_Gh_N=qf(mx_abs),
                n_GhN_pos=int(npos), n_GhN_neg=int(nneg), signs_ok=bool(mx_plus < 0 < mn_minus))


def identity_controls():
    out = {}
    # (a) shifted-Chebyshev abs-sum on a polynomial built from known Chebyshev coefficients
    import random
    rnd = random.Random(7)
    Ty = shifted_cheb(9); Tv = shifted_cheb(5)
    terms = []; tot = fmpq(0)
    for k in range(10):
        for l in range(6):
            c = fmpq(rnd.randint(-50, 50), rnd.randint(1, 9))
            tot += abs(c)
            terms.append((c * Ty[k], Tv[l]))
    out['cheb_abs_sum_exact'] = bool(cheb_abs_sum(pack(terms, 7), 7) == tot)
    # (b) angular operator in v equals [(1-c^2) g_c]_c for g = tau_k(c^2), and tau_k(c^2) = T_2k(c)
    ok = True
    C = fmpq_poly([0, 1]); C2 = C * C
    for k in range(1, 9):
        tk = shifted_cheb(k)[k]
        g = tk(C2)
        lhs = ((1 - C2) * g.derivative()).derivative()
        Ak = 4 * Yv * (1 - Yv) * tk.derivative().derivative() + (2 - 6 * Yv) * tk.derivative()
        ok = ok and (lhs == Ak(C2))
        Tc = [ONE, C]
        for n in range(1, 2 * k):
            Tc.append(2 * C * Tc[n] - Tc[n - 1])
        ok = ok and (g == Tc[2 * k])
    out['angular_operator_identity'] = bool(ok)
    # (c) radial + angular operator on the harmonic function r^2 P2(c): y^2[(1-y)^2 u_y]_y + Ang u = 0 (exact)
    F = Fraction
    good = True
    for (y, v) in ((F(1, 3), F(1, 5)), (F(5, 7), F(2, 3))):
        a = 1 / y - 1                     # r
        R = a * a; R1 = 2 * a * (-1 / y ** 2); R2 = 2 / y ** 4 + 2 * a * 2 / y ** 3
        Sv = (3 * v - 1) / 2; AngS = (2 - 6 * v) * F(3, 2)
        rad = y * y * ((1 - y) ** 2 * R2 - 2 * (1 - y) * R1)
        good = good and (rad * Sv + R * AngS == 0)
    out['radial_operator_harmonic_r2P2'] = bool(good)
    # (d) source transformation: r^2 * (9/4) J^2 (1-c^2) / (r^6 Phi^7) == (9/4) J^2 (1-v) y^4 (1-y)^3 / D^7
    good = True
    for (y, v, uu, J) in ((F(1, 3), F(1, 5), F(2, 7), F(5, 2)), (F(4, 5), F(0), F(1, 9), F(3, 2))):
        r = (1 - y) / y
        Phi = 1 + 1 / (2 * r) + uu
        D = 1 - y / 2 + (1 - y) * uu
        good = good and (r * r * F(9, 4) * J * J * (1 - v) / (r ** 6 * Phi ** 7)
                         == F(9, 4) * J * J * (1 - v) * y ** 4 * (1 - y) ** 3 / D ** 7)
        good = good and (D == (1 - y) * Phi)
    out['source_transformation'] = bool(good)
    # (e) Lh(y) = -2
    out['Lh_of_y_is_minus2'] = bool((1 - Yv) * Yv.derivative().derivative() - 2 * Yv.derivative() == -2)
    out['all_ok'] = all(out.values())
    log("identity controls:", out)
    return out


def k1_slab(j0_q, Ny, K, prec):
    disc = Disc(Ny, K)
    wf, resf, _ = disc.newton(0.0, qf(S1))
    wa, hist = refine(disc, wf, fmpq(0), S1, prec)
    wq = [to_dyadic(a, prec - 24) for a in wa]
    t0 = j0_q * j0_q
    c = certify(wq, Ny, K, None, S1, t_slab=t0)
    lo = 2 * (c['p0'] - c['eta']); hi = 2 * (c['p0'] + c['eta'])
    cont = bool(lo <= fmpq(1, 40) <= hi and c['certified'])
    out = dict(j0=str(j0_q), Ny=Ny, K=K, prec_bits=prec, eta="%.3e" % qf(c['eta']),
               S_parts=["%.3e" % qf(s) for s in c['Sparts']], certified=c['certified'],
               positivity_certified=c['pos_ok'], delta_over_j2_lo=dec_floor(lo, 30), delta_over_j2_hi=dec_ceil(hi, 30),
               rel_width="%.3e" % qf((hi - lo) / lo), contains_1_over_40=cont)
    log("K1 slab j in (0, %s] Ny=%d K=%d: delta/j^2 in [%s, %s], contains 1/40: %s"
        % (j0_q, Ny, K, out['delta_over_j2_lo'], out['delta_over_j2_hi'], cont))
    return out


def toy_J0():
    """J = 0 (source coefficient 0): the engine must return u = 0 exactly and refuse a strict PASS (eta = 0)."""
    Ny, K = 24, 4
    disc = Disc(Ny, K)
    wf, resf = solve_continuation(disc, 1.0, 0.0)
    wa, hist = refine(disc, wf, fmpq(1), fmpq(0), 128)
    wq = [to_dyadic(a, 104) for a in wa]
    c = certify(wq, Ny, K, fmpq(1), fmpq(0))
    allzero = all(x == 0 for x in wq)
    out = dict(all_coefficients_zero=bool(allzero), S_exact=qstr(c['S']), eta_exact=qstr(c['eta']),
               strict_certificate=c['strict'], certified=c['certified'],
               u_is_exact_zero=bool(allzero and c['S'] == 0), delta="0" if allzero else "nonzero",
               ok=bool(allzero and c['S'] == 0 and not c['certified']))
    log("toy J=0:", out)
    return out


def intersect(encs):
    lo = max(e[0] for e in encs); hi = min(e[1] for e in encs)
    return lo, hi


def summary(lo, hi, j_q):
    j = j_q
    chi_lo = j / (4 * (1 + hi) ** 2); chi_hi = j / (4 * (1 + lo) ** 2)
    return dict(delta_lo=dec_floor(lo, 30), delta_hi=dec_ceil(hi, 30),
                delta_lo_exact=qstr(lo), delta_hi_exact=qstr(hi),
                abs_width="%.3e" % qf(hi - lo), rel_width="%.3e" % qf((hi - lo) / lo),
                M_ADM_lo=dec_floor(1 + lo, 30), M_ADM_hi=dec_ceil(1 + hi, 30),
                chi_lo=dec_floor(chi_lo, 30), chi_hi=dec_ceil(chi_hi, 30),
                inside_stage0_upper_j2_over_40=bool(hi <= j * j / 40),
                inside_RN_bound_M2_le_1_plus_3J=bool((1 + hi) ** 2 <= 1 + 3 * (j / 4)),
                chi_below_1=bool(chi_hi < 1))


def parse_q(s):
    a, b = s.split('/')
    return fmpq(int(a), int(b))


# =============================================================================================
# 8. main
# =============================================================================================
RES_A = (64, 12, 224)
RES_B = (88, 16, 288)


def main():
    res = dict(seat="rep2", object="Bowen-York spinning puncture, mu = 1, j = 4J, delta = M_ADM - 1",
               method=("y = 1/(1+r) Chebyshev x Fourier-cosine (T_k(2cos^2-1)) in v = cos^2 theta; Gauss-node collocation "
                       "+ arb refinement; exact-rational residual Q, shifted-Chebyshev abs-sum; u_- >= 0 by centred "
                       "Taylor-form subdivision (D >= 1/2, kappa = 1); slope route delta in 2 j^2 [p(0) -+ eta]"))
    res['identity_controls'] = identity_controls()
    res['toy_J0'] = toy_J0()

    # K1
    res['K1'] = dict(slab_A=k1_slab(fmpq(1, 2 ** 24), 32, 4, 192), slab_B=k1_slab(fmpq(1, 2 ** 24), 44, 6, 224))
    res['K1']['contains_1_over_40'] = bool(res['K1']['slab_A']['contains_1_over_40']
                                           and res['K1']['slab_B']['contains_1_over_40'])

    # cross-check against REP1's certified j = 1 enclosure (REP1-NOTES.md, Results table, J = 0.25 row)
    rep1_lo = fmpq(238361051533259143533505592110, 10 ** 31)
    rep1_hi = fmpq(238361051533259143533844732458, 10 ** 31)
    o, e1 = run_point(fmpq(1), 48, 8, 224, label="j=1 (REP1 overlap)")
    res['crosscheck_j1_vs_REP1'] = dict(run=o, overlaps_REP1=bool(max(e1[0], rep1_lo) <= min(e1[1], rep1_hi)))
    log("j=1 overlaps REP1:", res['crosscheck_j1_vs_REP1']['overlaps_REP1'])

    # targets
    for jv in (6, 10):
        j_q = fmpq(jv)
        runs = []; encs = []
        for (Ny, K, prec) in (RES_A, RES_B):
            o, e = run_point(j_q, Ny, K, prec, label="j=%d" % jv)
            runs.append(o); encs.append(e)
        lo, hi = intersect(encs)
        s = summary(lo, hi, j_q)
        s['all_certified'] = all(r['certified'] for r in runs)
        s['coarse_contains_fine'] = bool(encs[0][0] <= encs[1][0] and encs[1][1] <= encs[0][1])
        s['resolutions_overlap'] = bool(lo <= hi)
        res['delta_j%d' % jv] = s
        res['runs_j%d' % jv] = runs
        log("FINAL j=%d: delta in [%s, %s] rel width %s certified %s" % (jv, s['delta_lo'], s['delta_hi'],
                                                                         s['rel_width'], s['all_certified']))

    # falsification sample (exact rational evaluation) at j = 10, resolution A
    o, e, (wq, c) = run_point(fmpq(10), *RES_A, label="j=10 sample", keep=True)
    smp = sample_signs(wq, RES_A[0], RES_A[1], fmpq(100), S1, c['eta'])
    smp['bound_128S'] = 128 * qf(c['S']); smp['two_eta'] = 2 * qf(c['eta'])
    res['sample_signs_j10'] = smp
    log("exact sample j=10:", smp)

    # nulls at j = 10 (resolution A)
    nulls = []
    for b in (fmpq(1, 10 ** 6), fmpq(1, 10 ** 9)):
        on, en = run_point(fmpq(10), *RES_A, label="null bump %s" % b, bump=b)
        nulls.append(dict(bump=str(b), eta=on['eta'], rel_width=on['rel_width'], certified=on['certified'],
                          delta_lo=on['delta_lo'], delta_hi=on['delta_hi']))
    res['null_bumped_ansatz'] = nulls
    # positivity null: shift V down past its sampled minimum -> must FAIL; normal shift passes
    w_terms, Lh_terms, V_terms, p0 = build_polys(wq, RES_A[0], RES_A[1])
    vmin = min(eval_terms(V_terms, fmpq(i, 64), fmpq(l, 16)) for i in range(65) for l in range(17))
    okd = certify_positive(V_terms, vmin + fmpq(1, 10 ** 6))
    oku = certify_positive(V_terms, c['eta'])
    res['null_positivity'] = dict(sampled_min_V="%.6e" % qf(vmin), shifted_down_passes=okd[0],
                                  normal_passes=oku[0], normal_boxes=oku[1], ok=bool((not okd[0]) and oku[0]))
    log("positivity null:", res['null_positivity'])

    # K2: source 9/4 -> 9/4 (1 + 1e-3) at j = 10
    s1p = S1 * (1 + fmpq(1, 1000))
    k2runs = []; k2e = []
    for (Ny, K, prec) in (RES_A, RES_B):
        o, e = run_point(fmpq(10), Ny, K, prec, s1=s1p, label="K2 j=10")
        k2runs.append(o); k2e.append(e)
    plo, phi = intersect(k2e)
    blo = parse_q(res['delta_j10']['delta_lo_exact']); bhi = parse_q(res['delta_j10']['delta_hi_exact'])
    res['K2'] = dict(runs=k2runs, delta_perturbed_lo=dec_floor(plo, 30), delta_perturbed_hi=dec_ceil(phi, 30),
                     all_certified=all(r['certified'] for r in k2runs),
                     shift_lo=dec_floor(plo - bhi, 22), shift_hi=dec_ceil(phi - blo, 22),
                     rel_shift_lo=dec_floor((plo - bhi) / bhi, 16), rel_shift_hi=dec_ceil((phi - blo) / blo, 16))
    log("K2: shift of delta(10) in [%s, %s], relative [%s, %s]" % (res['K2']['shift_lo'], res['K2']['shift_hi'],
                                                                 res['K2']['rel_shift_lo'], res['K2']['rel_shift_hi']))
    res['total_seconds'] = round(time.time() - T0, 1)
    with open(os.path.join(HERE, "rep2_result.json"), "w") as f:
        json.dump(res, f, indent=1, sort_keys=True)
    log("wrote rep2_result.json")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "probe":
        probe(float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
        sys.exit(0)
    main()
