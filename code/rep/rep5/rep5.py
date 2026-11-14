"""REP5 (BHSLICE blind replication seat): certified slope s(j) = dM_ADM/dmu (J fixed, mu = 1) and delta(j) = M_ADM - 1 of
the Bowen-York spinning puncture, by an extension of the s-coordinate engine idea to the linear tangent problem.

Object (brief): Delta u = -S Phi^-7, Phi = 1 + mu/(2r) + u, S = 9 J^2 sin^2(theta)/(4 r^6), u -> 0, M = mu + 2 lim r u,
j = 4J/mu^2, mu = 1, J = j/4. Tangent: w = du/dmu, Delta w = 7 S Phi^-8 (1/(2r) + w), slope = 1 + 2 lim r w.

Coordinate (derived in RESULT notes; same as the engine): s = a/(r + a), a = 1/2, c = cos(theta), t = c^2.
  r = (1-s)/(2s), 1/(2r) = s/(1-s), Phi = D/(1-s) with D = 1 + (1-s) u.
  Delta f(r) P_l(c) = (4 s^4/(1-s)) [ (1-s) f'' - 2 f' - l(l+1) f/(s^2 (1-s)) ] P_l(c)        (f' = df/ds)
  S Phi^-7 = 144 J^2 (1-c^2) s^6 (1-s) D^-7,  S Phi^-8 = 144 J^2 (1-c^2) s^6 (1-s)^2 D^-8.
  Point problem:   Lh u + tau src D^-7 = 0,               tau src = (9/4) j^2 (1-c^2) s^2 (1-s)^2.
  Tangent problem: Lh w = lam tau src (s + (1-s) w) D^-8,  lam = 7 (lam = 6 is the KN control).
  2 lim r f = p_f(0) for f = s p_f(s) + O(s^2) (a = 1/2): delta = p_u(0), slope = 1 + p_w(0).
Ansatz (engine's): f = s p(s) + sum_{k=1..K} s^2 (1-s) q_k(s) P_{2k}(c), p, q_k Chebyshev in X = 2s - 1 (degree N).

Certificates (all inequalities on the closed square s in [0,1], t in [0,1]):
  point:   Q = D_N^7 Lh u_N + tau src is a polynomial; S_u >= sum |Chebyshev coefficients of Q| (arb, outward) >= sup|Q|.
           e_u > S_u/2; u_pm = u_N pm e_u s; Bernstein (exact rationals, engine) u_-/s > 0  =>  D_- >= 1, D_N >= 1.
           Lh(s) = -2 and src >= 0, D^-7 decreasing => u_+ super, u_- sub => u_- <= u <= u_+ (comparison).
           delta in [p_u(0) - e_u, p_u(0) + e_u].
  tangent: G[psi] = D^8 Lh psi - lam tau src Z_psi, Z_psi = s + (1-s) psi (so Delta psi - V z_psi = 4 s^4 G/((1-s) D^8)).
           UPPER bound on the true slope needs Delta psi <= V z_psi for the TRUE V = 7 S Phi^-8. If z_psi >= 0 this
           follows from Delta psi <= V_lo z_psi with V_lo = 7 S Phi_+^-8 <= V (Phi_+ from the point SUPER-solution).
           So: fit psi_N against D_+, S_G >= sum|cheb(G_+[psi_N])|, e > S_G/2, psi = psi_N + e s:
             G_+[psi] = G_+[psi_N] - 2 e D_+^8 - lam tau src e s (1-s) <= S_G - 2e < 0 (D_+ >= 1),
             Bernstein Z_{psi_N}/s > 0 (then Z_psi >= Z_{psi_N} > 0)  =>  (Delta - V)(psi - w) <= 0  =>  psi >= w
             =>  slope <= 1 + p_psi(0) + e.
           LOWER bound: phi = phi_N - e' s fitted against D_- (Phi_- from the point SUB-solution), G_-[phi] >= 0 and
             Z_phi >= 0 => Delta phi >= V_hi z_phi >= V z_phi => phi <= w => slope >= 1 + p_phi(0) - e'.
           BRIEF OBJECT (super-solution with V = 7 S Phi_-^-8): psi = phi_N + e' s, G_-[psi] <= 0. It bounds from above
             the slope of the comparison problem with potential V_hi >= V, whose slope is <= the true slope; it is NOT
             an upper bound on the true slope (reported separately, labelled).
No pickle anywhere; exact rationals are stored as "p/q" strings; arb only inside functions (prec restored).
"""
import os
import sys
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import json
import math
import time
import hashlib
from fractions import Fraction

sys.path.insert(0, '<LAB>/POCKETLAB')
from gravity_engine import lichnerowicz_axisym as E   # used for: bernstein_positive, cheb_s, leg_c, round_up_dyadic, dec_*
import numpy as np
import scipy.linalg as sla
from flint import arb, arb_mat, fmpq, fmpq_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
PREC = 256
RB = 180                      # fitted coefficients rounded to multiples of 2^-RB
WARM_N, WARM_K = 120, 14      # fit-only continuation degree for large N


# ============================================================ exact helpers
def Fq(x):
    x = Fraction(x)
    return fmpq(int(x.numerator), int(x.denominator))


def qs(q):
    return "%d/%d" % (int(q.p), int(q.q))


def qp(s):
    a, b = s.split('/')
    return fmpq(int(a), int(b))


def qf(q):
    return float(Fraction(int(q.p), int(q.q)))


def dyad(x, bits=RB):
    """float or fmpq -> nearest multiple of 2^-bits (fmpq)"""
    f = Fraction(int(x.p), int(x.q)) if isinstance(x, fmpq) else Fraction(x)
    v = f * (1 << bits)
    n = math.floor(v + Fraction(1, 2))
    return fmpq(n, 1 << bits)


def arb_up_fmpq(a):
    """exact fmpq >= every point of the arb ball a"""
    m, e = a.upper().mid().man_exp()
    m, e = int(m), int(e)
    return fmpq(m * 2 ** e) if e >= 0 else fmpq(m, 2 ** (-e))


class Prec:
    def __init__(self, p):
        self.p = p

    def __enter__(self):
        self.old = ctx.prec
        ctx.prec = self.p

    def __exit__(self, *a):
        ctx.prec = self.old


# ============================================================ Chebyshev coefficient algebra (X = 2s - 1), generic field
def c_add(a, b):
    n = max(len(a), len(b))
    z = a[0] * 0 if a else b[0] * 0
    return [(a[i] if i < len(a) else z) + (b[i] if i < len(b) else z) for i in range(n)]


def c_scale(a, k):
    return [k * x for x in a]


def c_mulX(a):
    z = a[0] * 0
    out = [z] * (len(a) + 1)
    for k, x in enumerate(a):
        if k == 0:
            out[1] = out[1] + x
        else:
            out[k + 1] = out[k + 1] + x / 2
            out[k - 1] = out[k - 1] + x / 2
    return out


def c_muls(a):          # s = (1 + X)/2
    return [x / 2 for x in c_add(a, c_mulX(a))]


def c_muloms(a):        # 1 - s = (1 - X)/2
    return [x / 2 for x in c_add(a, [-y for y in c_mulX(a)])]


def c_derX(a):
    n = len(a) - 1
    z = a[0] * 0
    if n == 0:
        return [z]
    b = [z] * (n + 2)
    for k in range(n, 0, -1):
        b[k - 1] = b[k + 1] + 2 * k * a[k]
    b[0] = b[0] / 2
    return b[:n]


def c_ders(a):
    return [2 * x for x in c_derX(a)]


def c_Rop(f):           # (1-s) f'' - 2 f'
    f1 = c_ders(f)
    f2 = c_ders(f1)
    return c_add(c_muloms(f2), c_scale(f1, -2))


def c_at0(a):           # value at s = 0 (X = -1)
    return sum((x if k % 2 == 0 else -x) for k, x in enumerate(a))


def radial_pieces(pc, qcs):
    """value and Lh coefficient lists for each Legendre mode k = 0..K (mode k <-> P_{2k}(c))."""
    vals, lhs = [], []
    v0 = c_muls(pc)
    vals.append(v0)
    lhs.append(c_Rop(v0))
    for k, q in enumerate(qcs, 1):
        l = 2 * k
        g = c_muls(c_muls(c_muloms(q)))
        vals.append(g)
        lhs.append(c_add(c_Rop(g), c_scale(q, -l * (l + 1))))
    return vals, lhs


def deg_of(a):
    d = 0
    for k, x in enumerate(a):
        if x != 0:
            d = k
    return d


# ============================================================ float operators for the fit (never part of a proof)
_FOPS = {}


def float_ops(N, K):
    key = (N, K)
    if key in _FOPS:
        return _FOPS[key]
    s = np.array([(1 + math.cos(math.pi * i / N)) / 2 for i in range(N + 1)])
    s[0], s[N] = 1.0, 0.0
    cm = np.array([math.cos(math.pi * m / (2 * K)) for m in range(K + 1)]) if K > 0 else np.array([1.0])
    if K > 0:
        cm[0], cm[K] = 1.0, 0.0
    Lmax = N + 5
    CosF = np.array([[math.cos(math.pi * ((k * i) % (2 * N)) / N) for k in range(Lmax)] for i in range(N + 1)])
    PL = np.zeros((len(cm), K + 1))
    for m, c in enumerate(cm):
        P = [1.0, c]
        for l in range(1, 2 * K):
            P.append(((2 * l + 1) * c * P[l] - l * P[l - 1]) / (l + 1))
        for k in range(K + 1):
            PL[m, k] = P[2 * k]
    def vals_of(lst):
        M = np.zeros((Lmax, len(lst)))
        for j, a in enumerate(lst):
            M[:len(a), j] = a
        return CosF @ M
    Lblocks, Ublocks = [], []
    for k in range(K + 1):
        vl, ll = [], []
        for n in range(N + 1):
            e = [0.0] * (N + 1)
            e[n] = 1.0
            if k == 0:
                v = c_muls(e)
                lh = c_Rop(v)
            else:
                l = 2 * k
                v = c_muls(c_muls(c_muloms(e)))
                lh = c_add(c_Rop(v), c_scale(e, -l * (l + 1)))
            vl.append(v)
            ll.append(lh)
        Vk, Lk = vals_of(vl), vals_of(ll)
        Lblocks.append(np.kron(Lk, PL[:, k:k + 1]))
        Ublocks.append(np.kron(Vk, PL[:, k:k + 1]))
    L = np.hstack(Lblocks)
    U = np.hstack(Ublocks)
    sv = np.repeat(s, K + 1)
    cv = np.tile(cm, N + 1)
    om = 1 - sv
    src = 9.0 / 4.0 * (1 - cv * cv) * sv * sv * om * om
    out = dict(L=L, U=U, s=sv, om=om, src=src, N=N, K=K)
    _FOPS.clear()
    _FOPS[key] = out
    return out


def newton_point(ops, tau, w, itmax=60):
    L, U, om, src = ops['L'], ops['U'], ops['om'], ops['src']
    for it in range(itmax):
        D = 1 + om * (U @ w)
        F = L @ w + tau * src * D ** -7
        J = L - (7 * tau * src * om * D ** -8)[:, None] * U
        dw = np.linalg.solve(J, F)
        # damped step: keep D > 0
        lam = 1.0
        while lam > 1e-4:
            wn = w - lam * dw
            if np.all(1 + om * (U @ wn) > 0):
                break
            lam /= 2
        w = wn
        if np.max(np.abs(lam * dw)) < 1e-14 * max(1.0, np.max(np.abs(w))):
            break
    D = 1 + om * (U @ w)
    return w, float(np.max(np.abs(L @ w + tau * src * D ** -7)))


# ============================================================ arb grid machinery (proof path and refinement residuals)
_COS = {}


def cos_table(n):
    """arb cos(pi m / n), m = 0..2n-1 (balls)"""
    if (n, ctx.prec) not in _COS:
        tab = [arb.cos_pi_fmpq(fmpq(m, n)) for m in range(2 * n)]
        _COS.clear()
        _COS[(n, ctx.prec)] = tab
    return _COS[(n, ctx.prec)]


def legendre_even_vals(n_t, K):
    """arb matrix (K+1) x (n_t+1): P_{2k}(c_m), c_m = cos(pi m/(2 n_t)) (so t_m = c_m^2 are Lobatto points of [0,1])."""
    M = arb_mat(K + 1, n_t + 1)
    for m in range(n_t + 1):
        c = arb.cos_pi_fmpq(fmpq(m, 2 * n_t)) if n_t > 0 else arb(1)
        if n_t > 0 and m == n_t:
            c = arb(0)
        if m == 0:
            c = arb(1)
        P = [arb(1), c]
        for l in range(1, 2 * K):
            P.append(((2 * l + 1) * c * P[l] - l * P[l - 1]) / (l + 1))
        for k in range(K + 1):
            M[k, m] = P[2 * k]
    return M


def grid_values(lists, n_s, PLA):
    """lists: K+1 exact (fmpq) Chebyshev coefficient lists (one per Legendre mode). Returns arb_mat (n_s+1) x (n_t+1)
    of sum_k f_k(s_i) P_{2k}(c_m) at s_i = (1 + cos(pi i/n_s))/2."""
    tab = cos_table(n_s)
    Lmax = max(len(a) for a in lists)
    K1 = len(lists)
    C = arb_mat(Lmax, K1)
    for k, a in enumerate(lists):
        for i, x in enumerate(a):
            if x != 0:
                C[i, k] = arb(x)
    Cs = arb_mat(n_s + 1, Lmax)
    for i in range(n_s + 1):
        for k in range(Lmax):
            Cs[i, k] = tab[(k * i) % (2 * n_s)]
    return (Cs * C) * PLA


def s_nodes(n_s):
    tab = cos_table(n_s)
    out = [(1 + tab[i]) / 2 for i in range(n_s + 1)]
    out[0] = arb(1)
    out[n_s] = arb(0)
    return out


def t_nodes(n_t):
    if n_t == 0:
        return [arb(1)]
    out = [(1 + arb.cos_pi_fmpq(fmpq(m, n_t))) / 2 for m in range(n_t + 1)]
    out[0] = arb(1)
    out[n_t] = arb(0)
    return out


_DCT = {}


def dct_mat(n):
    """arb (n+1)x(n+1): coefficients a_k of sum a_k T_k from values at x_i = cos(pi i/n) (exact for degree <= n)."""
    if (n, ctx.prec) in _DCT:
        return _DCT[(n, ctx.prec)]
    if n == 0:
        M = arb_mat(1, 1)
        M[0, 0] = arb(1)
        return M
    tab = cos_table(n)
    M = arb_mat(n + 1, n + 1)
    two_n = arb(2) / n
    for k in range(n + 1):
        vk = arb(1) / 2 if k in (0, n) else arb(1)
        for i in range(n + 1):
            wi = 0.5 if i in (0, n) else 1
            M[k, i] = two_n * vk * wi * tab[(k * i) % (2 * n)]
    if len(_DCT) > 2:
        _DCT.clear()
    _DCT[(n, ctx.prec)] = M
    return M


def cheb_abs_sum(F, n_s, n_t):
    """upper bound (exact fmpq) of sum |Chebyshev coefficients| of the polynomial with values F on the Lobatto grid."""
    A = dct_mat(n_s) * F * dct_mat(n_t).transpose()
    tot = arb(0)
    for k in range(n_s + 1):
        for m in range(n_t + 1):
            tot += abs(A[k, m])
    return arb_up_fmpq(tot), A


def mat_absmax_mid(F):
    best = 0.0
    for i in range(F.nrows()):
        for j in range(F.ncols()):
            v = abs(float(F[i, j].mid()))
            if v > best:
                best = v
    return best


# ============================================================ problem pieces on a grid
def point_residual_grid(pc, qcs, tau, n_s, n_t, form='Q'):
    """arb grid of Q = D^7 Lh u + tau src (form 'Q') or of R = Lh u + tau src D^-7 (form 'R')."""
    K = len(qcs)
    vals, lhs = radial_pieces(pc, qcs)
    PLA = legendre_even_vals(n_t, K)
    Ug = grid_values(vals, n_s, PLA)
    Lg = grid_values(lhs, n_s, PLA)
    sn, tn = s_nodes(n_s), t_nodes(n_t)
    taua = arb(tau)
    F = arb_mat(n_s + 1, n_t + 1)
    for i in range(n_s + 1):
        si = sn[i]
        oi = 1 - si
        srs = arb(9) / 4 * si * si * oi * oi
        for m in range(n_t + 1):
            D = 1 + oi * Ug[i, m]
            ts = taua * srs * (1 - tn[m])
            D2 = D * D
            D4 = D2 * D2
            if form == 'Q':
                F[i, m] = D4 * D2 * D * Lg[i, m] + ts
            else:
                F[i, m] = Lg[i, m] + ts / (D4 * D2 * D)
    return F


def tangent_residual_grid(upc, uqcs, pc, qcs, tau, lam, n_s, n_t, form='G'):
    """G = D^8 Lh psi - lam tau src (s + (1-s) psi) ('G') or R = G / D^8 ('R'); D from u given by (upc, uqcs)."""
    Ku, K = len(uqcs), len(qcs)
    KK = max(Ku, K)
    uv, _ = radial_pieces(upc, uqcs)
    pv, pl = radial_pieces(pc, qcs)
    z = [fmpq(0)]
    uv = uv + [z] * (KK - Ku)
    pv = pv + [z] * (KK - K)
    pl = pl + [z] * (KK - K)
    PLA = legendre_even_vals(n_t, KK)
    Ug = grid_values(uv, n_s, PLA)
    Pg = grid_values(pv, n_s, PLA)
    Lg = grid_values(pl, n_s, PLA)
    sn, tn = s_nodes(n_s), t_nodes(n_t)
    taua = arb(tau) * lam
    F = arb_mat(n_s + 1, n_t + 1)
    for i in range(n_s + 1):
        si = sn[i]
        oi = 1 - si
        srs = arb(9) / 4 * si * si * oi * oi
        for m in range(n_t + 1):
            D = 1 + oi * Ug[i, m]
            D2 = D * D
            D4 = D2 * D2
            D8 = D4 * D4
            rhs = taua * srs * (1 - tn[m]) * (si + oi * Pg[i, m])
            F[i, m] = D8 * Lg[i, m] - rhs if form == 'G' else Lg[i, m] - rhs / D8
    return F


def degs_point(pc, qcs):
    Np, Nq = len(pc) - 1, (len(qcs[0]) - 1 if qcs else 0)
    du = max(Np + 1, Nq + 3 if qcs else 0)
    dlu = max(Np, Nq + 2 if qcs else 0)
    K = len(qcs)
    return 7 * (du + 1) + dlu, 8 * K


def degs_tangent(upc, uqcs, pc, qcs):
    Np, Nq = len(upc) - 1, (len(uqcs[0]) - 1 if uqcs else 0)
    du = max(Np + 1, Nq + 3 if uqcs else 0)
    Mp, Mq = len(pc) - 1, (len(qcs[0]) - 1 if qcs else 0)
    dpsi = max(Mp + 1, Mq + 3 if qcs else 0)
    dlp = max(Mp, Mq + 2 if qcs else 0)
    ds = max(8 * (du + 1) + dlp, 4 + 1 + dpsi)
    dt = max(8 * len(uqcs) + len(qcs), 1 + len(qcs))
    return ds, dt


# ============================================================ fits
def split(vec, N, K):
    return list(vec[:N + 1]), [list(vec[(N + 1) * k:(N + 1) * (k + 1)]) for k in range(1, K + 1)]


def join(pc, qcs):
    out = list(pc)
    for q in qcs:
        out += list(q)
    return out


def fit_point(j, N, K, refine=True, log=print):
    """u-scale coefficients (exact dyadic fmpq) of the point problem at j (fit only)."""
    jf = float(j)
    sched, x = [], min(1.0, jf)
    while x < jf:
        sched.append(x)
        x *= 2 ** 0.5
    sched.append(jf)
    if N > WARM_N:
        # continuation at a lower degree (fit only), then zero-padded warm start at (N, K) and Newton at the final j
        N0, K0 = WARM_N, min(K, WARM_K)
        ops0 = float_ops(N0, K0)
        w0 = np.zeros((N0 + 1) * (K0 + 1))
        for jj in sched:
            w0, _ = newton_point(ops0, jj * jj, w0)
        w = np.zeros((N + 1) * (K + 1))
        for k in range(K0 + 1):
            w[(N + 1) * k:(N + 1) * k + N0 + 1] = w0[(N0 + 1) * k:(N0 + 1) * (k + 1)]
        ops = float_ops(N, K)
        w, resf = newton_point(ops, jf * jf, w)
    else:
        ops = float_ops(N, K)
        w = np.zeros((N + 1) * (K + 1))
        for jj in sched:
            w, resf = newton_point(ops, jj * jj, w)
    tau = Fq(j) ** 2
    coef = [dyad(v) for v in w]
    hist = [resf]
    if refine:
        D = 1 + ops['om'] * (ops['U'] @ w)
        J = ops['L'] - (7 * jf * jf * ops['src'] * ops['om'] * D ** -8)[:, None] * ops['U']
        lu = sla.lu_factor(J)
        with Prec(PREC):
            for it in range(16):
                pc, qcs = split(coef, N, K)
                F = point_residual_grid(pc, qcs, tau, N, K, form='R')
                Ff = np.array([float(F[i, m].mid()) for i in range(N + 1) for m in range(K + 1)])
                r = float(np.max(np.abs(Ff)))
                hist.append(r)
                if r < 1e-60 or (len(hist) > 3 and r > 0.5 * hist[-2]):
                    break
                d = sla.lu_solve(lu, Ff)
                coef = [dyad(c - Fq(float(v))) for c, v in zip(coef, d)]
    return coef, hist


def fit_tangent(j, N, K, upc, uqcs, lam=7, refine=True):
    """tangent coefficients against the D built from (upc, uqcs) (fit only)."""
    ops = float_ops(N, K)
    Ku = len(uqcs)
    # D at collocation points (float) from u coefficients (u may have other degrees): evaluate via arb grid at N, K grid
    with Prec(PREC):
        KK = max(K, Ku)
        uv, _ = radial_pieces(upc, uqcs)
        uv = uv + [[fmpq(0)]] * (KK - Ku)
        PLA = legendre_even_vals(K, KK)
        Ug = grid_values(uv, N, PLA)
        uf = np.array([float(Ug[i, m].mid()) for i in range(N + 1) for m in range(K + 1)])
    om, src, s = ops['om'], ops['src'], ops['s']
    jf = float(j)
    D = 1 + om * uf
    V = lam * jf * jf * src * D ** -8
    J = ops['L'] - (V * om)[:, None] * ops['U']
    rhs = V * s
    lu = sla.lu_factor(J)
    w = sla.lu_solve(lu, rhs)
    coef = [dyad(v) for v in w]
    tau = Fq(j) ** 2
    hist = [float(np.max(np.abs(J @ w - rhs)))]
    if refine:
        with Prec(PREC):
            for it in range(16):
                pc, qcs = split(coef, N, K)
                F = tangent_residual_grid(upc, uqcs, pc, qcs, tau, lam, N, K, form='R')
                Ff = np.array([float(F[i, m].mid()) for i in range(N + 1) for m in range(K + 1)])
                r = float(np.max(np.abs(Ff)))
                hist.append(r)
                if r < 1e-60 or (len(hist) > 3 and r > 0.5 * hist[-2]):
                    break
                d = sla.lu_solve(lu, Ff)
                coef = [dyad(c - Fq(float(v))) for c, v in zip(coef, d)]
    return coef, hist


# ============================================================ Bernstein positivity (exact; engine routine)
def cheb_to_poly(a):
    T = E.cheb_s(len(a) - 1)
    p = fmpq_poly([0])
    for k, x in enumerate(a):
        if x != 0:
            p += x * T[k]
    return p


def positivity_terms(pc, qcs, const=None, oms_factor=False, shift=fmpq(0)):
    """terms of  [const +] (1-s)^{oms} * (p(s) - shift + s(1-s) sum q_k P_{2k}(c))  as [(poly s, poly c)]"""
    K = len(qcs)
    Pl = E.leg_c(2 * K)
    fac = E.OMS if oms_factor else E.ONE
    p = cheb_to_poly(pc) - shift
    terms = [(fac * p + (const if const is not None else 0), E.ONE)]
    for k, q in enumerate(qcs, 1):
        qq = cheb_to_poly(q)
        if qq != 0:
            terms.append((fac * E.SP * E.OMS * qq, Pl[2 * k]))
    return terms


# ============================================================ certificates
def certify_point(j, coef, N, K, log=print):
    t0 = time.time()
    tau = Fq(j) ** 2
    pc, qcs = split(coef, N, K)
    n_s, n_t = degs_point(pc, qcs)
    with Prec(PREC):
        F = point_residual_grid(pc, qcs, tau, n_s, n_t, form='Q')
        S, _ = cheb_abs_sum(F, n_s, n_t)
    t1 = time.time()
    e = E.round_up_dyadic(S * (1 + fmpq(1, 2 ** 20)) / 2) if S > 0 else fmpq(0)
    pos = E.bernstein_positive(positivity_terms(pc, qcs, shift=e))
    t2 = time.time()
    p0 = c_at0(pc)
    ok = bool(pos['status'] == 'PASS' and (S > 0) and 2 * e > S)
    rec = dict(j=qs(Fq(j)), N=N, K=K, grid=[n_s, n_t], S=qs(S), S_float=qf(S), e=qs(e), e_float=qf(e),
               p0=qs(p0), positivity=pos, certified=ok,
               delta_lo=qs(p0 - e), delta_hi=qs(p0 + e),
               delta_lo_dec=E.dec_floor(p0 - e, 30), delta_hi_dec=E.dec_ceil(p0 + e, 30),
               time_S=t1 - t0, time_bern=t2 - t1)
    log("  point j=%s N=%d K=%d grid %s: S=%.3e e=%.3e pos=%s delta in [%s, %s] (%.1fs + %.1fs)" % (
        rec['j'], N, K, rec['grid'], qf(S), qf(e), pos['status'], rec['delta_lo_dec'][:24], rec['delta_hi_dec'][:24],
        t1 - t0, t2 - t1))
    return rec


def u_shift(coef, N, K, e):
    pc, qcs = split(coef, N, K)
    pc = list(pc)
    pc[0] = pc[0] + e
    return pc, qcs


def certify_tangent(j, ucoef, Nu, Ku, eu, sign, tcoef, N, K, lam=7, log=print):
    """sign +1: D_+ (u_N + e_u s), super-solution psi_N + e s => upper bound of slope.
    sign -1: D_- (u_N - e_u s), sub-solution phi_N - e s => lower bound; and the brief object phi_N + e s (super w.r.t.
    V_hi) whose bound 1 + p(0) + e is reported as 'brief_object_upper' (not a bound on the true slope)."""
    t0 = time.time()
    tau = Fq(j) ** 2
    upc, uqcs = u_shift(ucoef, Nu, Ku, sign * eu)
    pc, qcs = split(tcoef, N, K)
    n_s, n_t = degs_tangent(upc, uqcs, pc, qcs)
    with Prec(PREC):
        F = tangent_residual_grid(upc, uqcs, pc, qcs, tau, lam, n_s, n_t, form='G')
        S, _ = cheb_abs_sum(F, n_s, n_t)
    t1 = time.time()
    e = E.round_up_dyadic(S * (1 + fmpq(1, 2 ** 20)) / 2) if S > 0 else fmpq(0)
    # Z/s positivity: 1 + (1-s)(p - shift + s(1-s) sum q P)
    shift = fmpq(0) if sign > 0 else e
    pos = E.bernstein_positive(positivity_terms(pc, qcs, const=E.ONE, oms_factor=True, shift=shift))
    t2 = time.time()
    p0 = c_at0(pc)
    strict = bool(S > 0 and 2 * e > S)
    rec = dict(j=qs(Fq(j)), lam=lam, sign=sign, N=N, K=K, grid=[n_s, n_t], S=qs(S), S_float=qf(S), e=qs(e),
               e_float=qf(e), p0=qs(p0), positivity_Z=pos, strict=strict, time_S=t1 - t0, time_bern=t2 - t1)
    if sign > 0:
        ub = 1 + p0 + e
        rec.update(slope_upper=qs(ub), slope_upper_dec=E.dec_ceil(ub, 30),
                   certified=bool(strict and pos['status'] == 'PASS'))
        log("  tangent+ j=%s lam=%d N=%d K=%d grid %s: S=%.3e e=%.3e posZ=%s slope <= %s (%.1fs + %.1fs)" % (
            rec['j'], lam, N, K, rec['grid'], qf(S), qf(e), pos['status'], rec['slope_upper_dec'][:24], t1 - t0, t2 - t1))
    else:
        lb = 1 + p0 - e
        bo = 1 + p0 + e
        rec.update(slope_lower=qs(lb), slope_lower_dec=E.dec_floor(lb, 30),
                   brief_object_upper=qs(bo), brief_object_upper_dec=E.dec_ceil(bo, 30),
                   certified=bool(strict and pos['status'] == 'PASS'), brief_object_certified=strict)
        log("  tangent- j=%s lam=%d N=%d K=%d grid %s: S=%.3e e=%.3e posZ=%s slope >= %s ; brief-object <= %s (%.1fs + %.1fs)" % (
            rec['j'], lam, N, K, rec['grid'], qf(S), qf(e), pos['status'], rec['slope_lower_dec'][:24],
            rec['brief_object_upper_dec'][:24], t1 - t0, t2 - t1))
    return rec


# ============================================================ float probe (never a proof)
def probe(j, N, K, Nt=None, Kt=None, refine=False):
    t0 = time.time()
    coef, hist = fit_point(j, N, K, refine=refine)
    pc, qcs = split(coef, N, K)
    p0 = qf(c_at0(pc))
    tails = dict(p_tail=max(abs(qf(x)) for x in pc[-4:]), q_tail=max(abs(qf(x)) for q in qcs for x in q[-4:]),
                 qK=max(abs(qf(x)) for x in qcs[-1]) if qcs else 0)
    # off-grid residual R at a finer grid (arb, but only a diagnostic)
    with Prec(128):
        F = point_residual_grid(pc, qcs, Fq(j) ** 2, 2 * N + 7, 2 * K + 3, form='R')
        rmax = mat_absmax_mid(F)
    Nt = Nt or N
    Kt = Kt or K
    tc, th = fit_tangent(j, Nt, Kt, pc, qcs, refine=refine)
    tp, tq = split(tc, Nt, Kt)
    slope = 1 + qf(c_at0(tp))
    with Prec(128):
        Ft = tangent_residual_grid(pc, qcs, tp, tq, Fq(j) ** 2, 7, 2 * Nt + 7, 2 * Kt + 3, form='R')
        trmax = mat_absmax_mid(Ft)
    ttails = dict(p_tail=max(abs(qf(x)) for x in tp[-4:]), q_tail=max(abs(qf(x)) for q in tq for x in q[-4:]))
    return dict(j=str(j), N=N, K=K, delta=p0, fit_hist=hist[:3] + hist[-1:], tails=tails, offgrid_R=rmax,
                slope=slope, tangent_offgrid_R=trmax, ttails=ttails, time=time.time() - t0)


def sha256_self():
    with open(os.path.abspath(__file__), 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()




# ============================================================ orchestration with JSON checkpoints
CKDIR = os.path.join(HERE, 'ckpt')


def ck_load(tag):
    p = os.path.join(CKDIR, tag + '.json')
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return {}


def ck_save(tag, d):
    os.makedirs(CKDIR, exist_ok=True)
    p = os.path.join(CKDIR, tag + '.json')
    with open(p + '.tmp', 'w') as f:
        json.dump(d, f)
    os.replace(p + '.tmp', p)


def coef_str(c):
    return [qs(x) for x in c]


def coef_parse(l):
    return [qp(x) for x in l]


def run_j(tag, j, Nu, Ku, Nt, Kt, kn=False, stop_after=None, log=print):
    """all stages for one j; resumable (each finished stage is checkpointed)."""
    d = ck_load(tag)
    d.setdefault('j', qs(Fq(j)))
    d.setdefault('degrees', dict(Nu=Nu, Ku=Ku, Nt=Nt, Kt=Kt))
    assert d['degrees'] == dict(Nu=Nu, Ku=Ku, Nt=Nt, Kt=Kt), 'checkpoint degrees differ'
    d['instrument_sha256_last'] = sha256_self()
    stages = ['ufit', 'ucert', 'tfit+', 'tcert+', 'tfit-', 'tcert-'] + (['knfit', 'kncert'] if kn else [])
    for st in stages:
        if st in d:
            continue
        t0 = time.time()
        if st == 'ufit':
            c, h = fit_point(j, Nu, Ku)
            d[st] = dict(coef=coef_str(c), hist=h)
        elif st == 'ucert':
            d[st] = certify_point(j, coef_parse(d['ufit']['coef']), Nu, Ku, log=log)
        elif st in ('tfit+', 'tfit-', 'knfit'):
            sign = -1 if st == 'tfit-' else 1
            lam = 6 if st == 'knfit' else 7
            uc = coef_parse(d['ufit']['coef'])
            eu = qp(d['ucert']['e'])
            upc, uqcs = u_shift(uc, Nu, Ku, sign * eu)
            c, h = fit_tangent(j, Nt, Kt, upc, uqcs, lam=lam)
            d[st] = dict(coef=coef_str(c), hist=h, lam=lam, sign=sign)
        elif st in ('tcert+', 'tcert-', 'kncert'):
            fst = {'tcert+': 'tfit+', 'tcert-': 'tfit-', 'kncert': 'knfit'}[st]
            sign, lam = d[fst]['sign'], d[fst]['lam']
            d[st] = certify_tangent(j, coef_parse(d['ufit']['coef']), Nu, Ku, qp(d['ucert']['e']), sign,
                                    coef_parse(d[fst]['coef']), Nt, Kt, lam=lam, log=log)
        d.setdefault('stage_time', {})[st] = time.time() - t0
        ck_save(tag, d)
        log("  [%s] stage %s done in %.1fs" % (tag, st, time.time() - t0))
        if stop_after == st:
            return d
    return d


def summarize(d):
    out = dict(j=d['j'], degrees=d['degrees'])
    u = d.get('ucert', {})
    out['point_certified'] = u.get('certified')
    out['delta_lo'], out['delta_hi'] = u.get('delta_lo'), u.get('delta_hi')
    out['delta_lo_dec'], out['delta_hi_dec'] = u.get('delta_lo_dec'), u.get('delta_hi_dec')
    tp, tm = d.get('tcert+', {}), d.get('tcert-', {})
    out['slope_upper'] = tp.get('slope_upper')
    out['slope_upper_dec'] = tp.get('slope_upper_dec')
    out['slope_upper_certified'] = tp.get('certified')
    out['slope_lower'] = tm.get('slope_lower')
    out['slope_lower_dec'] = tm.get('slope_lower_dec')
    out['slope_lower_certified'] = tm.get('certified')
    out['brief_object_upper_dec'] = tm.get('brief_object_upper_dec')
    if 'kncert' in d:
        out['KN_slope_upper_dec'] = d['kncert'].get('slope_upper_dec')
        out['KN_certified'] = d['kncert'].get('certified')
    out['e_u'] = u.get('e_float')
    out['e_t+'] = tp.get('e_float')
    out['e_t-'] = tm.get('e_float')
    out['positivity'] = dict(u=u.get('positivity', {}).get('status'), Zp=tp.get('positivity_Z', {}).get('status'),
                             Zm=tm.get('positivity_Z', {}).get('status'))
    out['stage_time'] = d.get('stage_time')
    return out


# ============================================================ controls
def control_dct():
    """C0: the DCT recovers known Chebyshev coefficients of a 2-D polynomial (exact identity, degree <= grid);
    null: a polynomial of degree above the grid is NOT recovered (aliasing shows), so the control can fail."""
    import random
    rnd = random.Random(5)
    ns, nt = 13, 5
    A = [[Fraction(rnd.randint(-99, 99), rnd.randint(1, 50)) for m in range(nt + 1)] for k in range(ns + 1)]
    with Prec(PREC):
        sn, tn = s_nodes(ns), t_nodes(nt)
        F = arb_mat(ns + 1, nt + 1)
        for i in range(ns + 1):
            X = 2 * sn[i] - 1
            Ts = [arb(1), X]
            for k in range(1, ns):
                Ts.append(2 * X * Ts[k] - Ts[k - 1])
            for m in range(nt + 1):
                Y = 2 * tn[m] - 1
                Tt = [arb(1), Y]
                for k in range(1, nt):
                    Tt.append(2 * Y * Tt[k] - Tt[k - 1])
                F[i, m] = sum((arb(Fq(A[k][mm])) * Ts[k] * Tt[mm] for k in range(ns + 1) for mm in range(nt + 1)), arb(0))
        S, C = cheb_abs_sum(F, ns, nt)
        ok = all(C[k, m].contains(arb(Fq(A[k][m]))) and C[k, m].rad() < 1e-60
                 for k in range(ns + 1) for m in range(nt + 1))
        Sx = sum(abs(A[k][m]) for k in range(ns + 1) for m in range(nt + 1))
        ok = ok and (S >= Fq(Sx)) and qf(S - Fq(Sx)) < 1e-50
        F2 = arb_mat(ns + 1, nt + 1)
        for i in range(ns + 1):
            X = 2 * sn[i] - 1
            T = [arb(1), X]
            for k in range(1, ns + 2):
                T.append(2 * X * T[k] - T[k - 1])
            for m in range(nt + 1):
                F2[i, m] = T[ns + 2]
        S2, C2 = cheb_abs_sum(F2, ns, nt)
        alias_detected = (not C2[ns - 2, 0].contains(arb(0)))
    return dict(name='C0 DCT exact recovery', ok=bool(ok), S=qf(S), S_exact=float(Sx),
                aliasing_null_detected=bool(alias_detected))


def control_operator():
    """C1: my Chebyshev-coefficient divided operator vs the engine's monomial Rop, Lh(s) = -2, and the Laplacian
    identity Delta(a/(r+a)) = -2 a^2/(r (r+a)^3) = (4 s^4/(1-s)) * (-2) at a = 1/2 (rational r)."""
    import random
    rnd = random.Random(7)
    ok = True
    for trial in range(5):
        a = [Fq(Fraction(rnd.randint(-50, 50), rnd.randint(1, 9))) for _ in range(9)]
        mine = cheb_to_poly(c_Rop(a))
        eng = E.SCompactBasis.Rop(cheb_to_poly(a))
        ok = ok and (mine == eng)
        g = c_muls(c_muls(c_muloms(a)))
        ok = ok and (cheb_to_poly(g) == E.HP * cheb_to_poly(a))
    ok2 = cheb_to_poly(c_Rop([fmpq(1, 2), fmpq(1, 2)])) == fmpq_poly([-2])
    ok3 = True
    for r in (Fraction(1, 3), Fraction(2), Fraction(7, 5)):
        a_ = Fraction(1, 2)
        s_ = a_ / (r + a_)
        ok3 = ok3 and (-2 * a_ ** 2 / (r * (r + a_) ** 3) == 4 * s_ ** 4 / (1 - s_) * (-2))
    a = [Fq(Fraction(rnd.randint(-50, 50), 7)) for _ in range(12)]
    ok4 = cheb_to_poly(c_ders(a)) == cheb_to_poly(a).derivative()
    return dict(name='C1 operator algebra', ok=bool(ok and ok2 and ok3 and ok4), Rop_vs_engine=bool(ok),
                Lh_s_eq_m2=bool(ok2), laplacian_identity=bool(ok3), derivative=bool(ok4))


def control_engine_j2(log=print):
    """C2: positive control against the engine's certified point enclosure at j = 2 (different residual-bound method:
    engine = exact rational Chebyshev of Q; mine = arb DCT on a Lobatto grid)."""
    old = ctx.prec
    try:
        rec_e, _ = E.certify_bowen_york(2, log=lambda *a: None)
    finally:
        ctx.prec = old
    coef, h = fit_point(2, 40, 8)
    mine = certify_point(2, coef, 40, 8, log=log)
    elo, ehi = qp(rec_e['delta_lo']), qp(rec_e['delta_hi'])
    mlo, mhi = qp(mine['delta_lo']), qp(mine['delta_hi'])
    overlap = bool(max(elo, mlo) <= min(ehi, mhi))
    return dict(name='C2 engine j=2 overlap', ok=bool(overlap and mine['certified'] and rec_e['status'] == 'CERTIFIED'),
                engine=[rec_e['delta_lo_dec'], rec_e['delta_hi_dec']], mine=[mine['delta_lo_dec'], mine['delta_hi_dec']],
                engine_eta_float=rec_e['eta_float'], mine_e=mine['e_float'])


def toy_J0(log=print):
    """T: J = 0. Point problem: u = 0 exactly. Tangent: V = 0, w = 0, slope exactly 1. Same code path."""
    N, K = 12, 2
    zero = [fmpq(0)] * ((N + 1) * (K + 1))
    rp = certify_point(0, zero, N, K, log=log)
    c, h = fit_tangent(0, N, K, *split(zero, N, K))
    allzero = all(x == 0 for x in c)
    rt = certify_tangent(0, zero, N, K, fmpq(0), +1, c, N, K, log=log)
    rm = certify_tangent(0, zero, N, K, fmpq(0), -1, c, N, K, log=log)
    ok = (rp['S'] == '0/1' and allzero and rt['S'] == '0/1' and rt['slope_upper'] == '1/1'
          and rm['slope_lower'] == '1/1')
    return dict(name='T toy J=0', ok=bool(ok), point_S=rp['S'], tangent_fit_all_zero=bool(allzero),
                slope_upper=rt['slope_upper'], slope_lower=rm['slope_lower'],
                note='S == 0 exactly: an exact, non-strict certificate; psi = 0 solves Delta psi = V z_psi with V = 0')


def series_coeffs():
    """exact small-j coefficients from series.py (derived in this seat; re-run there)."""
    return [Fraction(1, 40), Fraction(-133, 105600), Fraction(8923, 84864000)]


def slope_series(j, nterms=3):
    """s(j) = 1 + delta - 2 j delta'(j) = 1 + sum_k (1 - 4k) c_k j^(2k)  (scaling M = mu G(4J/mu^2))."""
    cs = series_coeffs()
    jj = Fraction(j)
    return 1 + sum((1 - 4 * (k + 1)) * cs[k] * jj ** (2 * k + 2) for k in range(nterms))


def delta_series(j, nterms=3):
    cs = series_coeffs()
    jj = Fraction(j)
    return sum(cs[k] * jj ** (2 * k + 2) for k in range(nterms))


# ============================================================ secant bounds from convexity (independent rigorous route)
def secant_bounds(tag, j, h, N, K, log=print):
    """M(mu) at fixed J is convex in mu (Phi_t = (1-t) Phi_0 + t Phi_1 is a super-solution for mu_t because x^-7 is
    convex, so u_{mu_t} <= (1-t) u_0 + t u_1 and lim r u is convex; mu enters M linearly). Hence
      [M(1) - M(1-h)]/h <= s(j) <= [M(1+h) - M(1)]/h,   M(mu) = mu (1 + delta(j/mu^2)) at J = j/4.
    Uses only certified point enclosures (no tangent problem)."""
    j, h = Fraction(j), Fraction(h)
    d = ck_load(tag)
    out = {}
    for name, jj in (('j', j), ('jp', j / (1 + h) ** 2), ('jm', j / (1 - h) ** 2)):
        if name in d:
            out[name] = d[name]
            continue
        c, hist = fit_point(jj, N, K)
        rec = certify_point(jj, c, N, K, log=log)
        rec.pop('positivity', None)
        rec['positivity_status'] = 'PASS' if rec['certified'] else 'FAIL'
        d[name] = rec
        out[name] = rec
        ck_save(tag, d)
    hq = Fq(h)
    d0lo, d0hi = qp(out['j']['delta_lo']), qp(out['j']['delta_hi'])
    dplo, dphi = qp(out['jp']['delta_lo']), qp(out['jp']['delta_hi'])
    dmlo, dmhi = qp(out['jm']['delta_lo']), qp(out['jm']['delta_hi'])
    upper = ((1 + hq) * (1 + dphi) - (1 + d0lo)) / hq
    lower = ((1 + d0lo) - (1 - hq) * (1 + dmhi)) / hq
    res = dict(j=qs(Fq(j)), h=qs(hq), N=N, K=K, all_certified=bool(all(out[k]['certified'] for k in out)),
               secant_upper=qs(upper), secant_lower=qs(lower),
               secant_upper_dec=E.dec_ceil(upper, 25), secant_lower_dec=E.dec_floor(lower, 25),
               nodes={k: dict(j=out[k]['j'], delta_lo_dec=out[k]['delta_lo_dec'], delta_hi_dec=out[k]['delta_hi_dec'],
                              e_float=out[k]['e_float']) for k in out})
    d['result'] = res
    ck_save(tag, d)
    log("  secant j=%s h=%s: %s <= s <= %s" % (res['j'], res['h'], res['secant_lower_dec'], res['secant_upper_dec']))
    return res


def k1t(js=('1/4', '1/2', '1'), N=40, K=8, log=print):
    """K1T: certified slope enclosures at small j vs the series s = 1 + sum (1-4k) c_k j^(2k) derived in series.py."""
    rows = []
    for js_ in js:
        j = Fraction(js_)
        d = run_j('k1t_' + js_.replace('/', '_'), j, N, K, N, K, log=log)
        lo, hi = qp(d['tcert-']['slope_lower']), qp(d['tcert+']['slope_upper'])
        s2, s3 = Fq(slope_series(j, 2)), Fq(slope_series(j, 3))
        dl, dh = qp(d['ucert']['delta_lo']), qp(d['ucert']['delta_hi'])
        rows.append(dict(j=js_, slope_lo_dec=E.dec_floor(lo, 25), slope_hi_dec=E.dec_ceil(hi, 25),
                         series2=E.dec_floor(s2, 25), series3=E.dec_floor(s3, 25),
                         resid2_over_j6=[qf((lo - s2) / Fq(j) ** 6), qf((hi - s2) / Fq(j) ** 6)],
                         resid3_over_j8=[qf((lo - s3) / Fq(j) ** 8), qf((hi - s3) / Fq(j) ** 8)],
                         delta_lo_dec=d['ucert']['delta_lo_dec'], delta_hi_dec=d['ucert']['delta_hi_dec'],
                         delta_resid3_over_j8=[qf((dl - Fq(delta_series(j))) / Fq(j) ** 8),
                                               qf((dh - Fq(delta_series(j))) / Fq(j) ** 8)],
                         certified=bool(d['tcert-']['certified'] and d['tcert+']['certified'] and d['ucert']['certified'])))
        log(json.dumps(rows[-1]))
    return dict(name='K1T small-j series', rows=rows, minus_11_c3=-11 * float(series_coeffs()[2]),
                note='resid2/j^6 must tend to -11 c3 as j -> 0; resid3/j^8 must stay bounded (next term -15 c4)')


# ============================================================ final collector
def _peak(logname):
    p = os.path.join(HERE, logname)
    if not os.path.exists(p):
        return None
    pk = None
    for line in open(p, encoding='utf-8', errors='replace'):
        if line.startswith('memcap:') and 'peak' in line:
            pk = line.strip()
    return pk


def collect():
    sha = sha256_self()
    res = dict(instrument='<LAB>/BHSLICE/rep5/rep5.py', instrument_sha256=sha,
               object='Bowen-York spinning puncture, maximal conformally flat slice, mu = 1, J = j/4; '
                      'slope s(j) = dM_ADM/dmu at fixed J; delta = M_ADM - 1',
               formulation='s = a/(r+a), a = 1/2; Chebyshev(2s-1) x even Legendre(cos theta); residual bound = sum of '
                           '|Chebyshev coefficients| of the exact polynomial residual, enclosed by an arb DCT on a '
                           'Lobatto grid of size >= its degree; positivity by exact Bernstein (engine routine)',
               notes=dict(
                   slope_upper='rigorous: super-solution psi_N + e s of the tangent problem w.r.t. V_lo = 7 S Phi_+^-8 '
                               '(Phi_+ from the certified point SUPER-solution) plus Z = s + (1-s) psi >= 0',
                   slope_lower='rigorous: sub-solution phi_N - e s w.r.t. V_hi = 7 S Phi_-^-8 plus Z >= 0',
                   brief_object_upper='the brief literal object: super-solution w.r.t. V_hi = 7 S Phi_-^-8. It bounds '
                                      'from above the slope of the problem with potential V_hi, which is <= the true '
                                      'slope; it is NOT an upper bound on s(j)'))
    for key, main, partner in (('j64', 'F_j64_main', 'F_j64_partner'), ('j4096', 'F_j4096_main', 'F_j4096_partner')):
        dm, dp = ck_load(main), ck_load(partner)
        if not dm:
            continue
        sm = summarize(dm)
        sp_ = summarize(dp) if dp else None
        ent = dict(j=sm['j'], degrees_main=sm['degrees'],
                   delta_lo=sm['delta_lo_dec'], delta_hi=sm['delta_hi_dec'],
                   delta_lo_exact=sm['delta_lo'], delta_hi_exact=sm['delta_hi'],
                   slope_upper=sm['slope_upper_dec'], slope_lower=sm['slope_lower_dec'],
                   slope_upper_exact=sm['slope_upper'], slope_lower_exact=sm['slope_lower'],
                   brief_object_upper=sm['brief_object_upper_dec'],
                   certified=bool(sm['point_certified'] and sm['slope_upper_certified'] and sm['slope_lower_certified']),
                   positivity=sm['positivity'], e=dict(u=sm['e_u'], tplus=sm['e_t+'], tminus=sm['e_t-']),
                   stage_time_s=sm['stage_time'], instrument_sha256_at_run=dm.get('instrument_sha256_last'))
        if 'KN_slope_upper_dec' in sm:
            ent['KN_lambda6_slope_upper'] = sm['KN_slope_upper_dec']
            ent['KN_certified'] = sm['KN_certified']
        if sp_:
            ent['partner'] = dict(degrees=sp_['degrees'], delta_lo=sp_['delta_lo_dec'], delta_hi=sp_['delta_hi_dec'],
                                  slope_upper=sp_['slope_upper_dec'], slope_lower=sp_['slope_lower_dec'],
                                  KN_lambda6_slope_upper=sp_.get('KN_slope_upper_dec'),
                                  instrument_sha256_at_run=dp.get('instrument_sha256_last'))
            mlo, mhi = qp(sm['delta_lo']), qp(sm['delta_hi'])
            plo, phi = qp(sp_['delta_lo']), qp(sp_['delta_hi'])
            slo, shi = qp(sm['slope_lower']), qp(sm['slope_upper'])
            tlo, thi = qp(sp_['slope_lower']), qp(sp_['slope_upper'])
            ent['convergence'] = dict(delta_intervals_overlap=bool(max(mlo, plo) <= min(mhi, phi)),
                                      main_delta_inside_partner=bool(plo <= mlo and mhi <= phi),
                                      slope_intervals_overlap=bool(max(slo, tlo) <= min(shi, thi)),
                                      main_slope_inside_partner=bool(tlo <= slo and shi <= thi))
        ent['peak_memory'] = dict(main=_peak('run_' + main + '.log'), partner=_peak('run_' + partner + '.log'))
        res[key] = ent
    for key, tag in (('secant_j64', 'S_j64'), ('secant_j4096', 'S_j4096')):
        d = ck_load(tag)
        if 'result' in d:
            r = d['result']
            r['peak_memory'] = _peak('run_' + tag + '.log')
            ent = res.get(key[7:])
            if ent:
                r['tangent_enclosure_inside_secant'] = bool(qp(r['secant_lower']) <= qp(ent['slope_lower_exact'])
                                                            and qp(ent['slope_upper_exact']) <= qp(r['secant_upper']))
            res[key] = r
    for fn, key in (('controls_rep5.json', 'controls'), ('k1t_rep5.json', 'K1T')):
        p = os.path.join(HERE, fn)
        if os.path.exists(p):
            res[key] = json.load(open(p))
    res['peak_memory_controls'] = _peak('run_controls.log')
    res['peak_memory_k1t'] = _peak('run_k1t.log')
    if 'j64' in res and 'KN_lambda6_slope_upper' in res['j64']:
        res['KN'] = dict(j='64', lambda_=6, slope_upper=res['j64']['KN_lambda6_slope_upper'],
                         note='exponent 7 -> 6 in the tangent operator only (point problem unchanged); rigorous upper '
                              'bound of the modified problem slope; must differ from the true slope enclosure')
    with open(os.path.join(HERE, 'result_rep5.json'), 'w') as f:
        json.dump(res, f, indent=1)
    return res


def main():
    mode = sys.argv[1]
    if mode == 'probe':
        for NK in sys.argv[3:]:
            N, K = map(int, NK.split(','))
            print(json.dumps(probe(Fraction(sys.argv[2]), N, K)), flush=True)
    elif mode == 'controls':
        out = dict(instrument_sha256=sha256_self())
        for f in (control_dct, control_operator, toy_J0, control_engine_j2):
            t0 = time.time()
            r = f()
            r['time'] = time.time() - t0
            print(json.dumps(r), flush=True)
            out[r['name']] = r
        with open(os.path.join(HERE, 'controls_rep5.json'), 'w') as f:
            json.dump(out, f, indent=1)
    elif mode == 'secant':
        r = secant_bounds(sys.argv[2], Fraction(sys.argv[3]), Fraction(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))
        print(json.dumps(r, indent=1), flush=True)
    elif mode == 'k1t':
        r = k1t()
        with open(os.path.join(HERE, 'k1t_rep5.json'), 'w') as f:
            json.dump(r, f, indent=1)
        print(json.dumps(r, indent=1), flush=True)
    elif mode == 'collect':
        r = collect()
        print(json.dumps({k: v for k, v in r.items() if k in ('j64', 'j4096', 'KN', 'secant_j64', 'secant_j4096')},
                         indent=1), flush=True)
    elif mode == 'run':
        tag, j = sys.argv[2], Fraction(sys.argv[3])
        Nu, Ku, Nt, Kt = map(int, sys.argv[4:8])
        kn = 'kn' in sys.argv[8:]
        stop = None
        for a in sys.argv[8:]:
            if a.startswith('stop='):
                stop = a[5:]
        d = run_j(tag, j, Nu, Ku, Nt, Kt, kn=kn, stop_after=stop)
        print(json.dumps(summarize(d), indent=1), flush=True)


if __name__ == '__main__':
    main()
