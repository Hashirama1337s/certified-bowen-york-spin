"""
rep1.py  -  blind replication seat REP1 for BHSLICE (Bowen-York spinning puncture, bare mass mu = 1).

Object: Delta u = -(9/4) J^2 sin^2(theta) / (r^6 (1 + 1/(2r) + u)^7) on R^3, u -> 0 at infinity,
        M_ADM = 1 + 2 lim_{r->oo} r u,  delta := M_ADM - 1,  j := 4 J.

Discretisation of THIS seat (chosen to differ from a Chebyshev-in-s x Legendre-in-cos(theta) fit):
  * radial compactification y = mu/(r + mu) = 1/(1 + r)   (scale mu, not mu/2), Chebyshev T_n(2y - 1);
  * angular basis cos(2 k theta) = T_{2k}(cos theta)  (Fourier cosine in theta, not Legendre);
  * unknown w with u = j^2 w = tau w, tau := j^2:
        w = y p(y) + y^2 (1 - y) sum_{k=1..K} q_k(y) cos(2 k theta);
  * Newton collocation (double), then iterative refinement with residuals in arb (high precision);
  * certificate: EXACT rational polynomial arithmetic (python-flint fmpq_poly / fmpz_mat), no rounding at all
    until one final coefficient-sum inequality  sup |sum C_kl T_k T_l| <= sum |C_kl|.

Equation in (y, c = cos theta)  (derived in REP1-NOTES.md):
    r^2 (Delta u - F) / (tau y^2 (1-y)) = Gh[w] := Lh w + s1 (1 - c^2) y^2 (1 - y)^2 D^-7,
    Lh w = (1 - y) w_yy - 2 w_y + ((1 - c^2) w_c)_c / (y^2 (1 - y)),
    D    = 1 - y/2 + tau (1 - y) w   (= (1 - y) Phi),      s1 = 9/64 (source coefficient 9/4 times J^2/j^2 = 1/16).
Sub/super pair: w+- = w_N +- eta y.  Lh(y) = -2 and D+ >= D_N >= 1/2, so
    Gh[w+] <= Gh[w_N] - 2 eta,   Gh[w-] >= Gh[w_N] + 2 eta   (needs D- > 0),
and Gh[w_N] = Q / D_N^7 with Q := D_N^7 Lh w_N + s1 (1-c^2) y^2 (1-y)^2 an exact polynomial.
With D >= 1/2:  sup|Gh[w_N]| <= 128 sup|Q| <= 128 S,  S = sum of |Chebyshev coefficients of Q| (exact rational).
eta := 64 S (1 + 2^-20) makes w+ a strict super-solution and w- a strict sub-solution on the closed square.
Comparison lemma => u- <= u <= u+  => lim r u in tau [p(0) - eta, p(0) + eta]  => delta in 2 tau [p(0) -+ eta].

Run:  python rep1.py   (prints and writes rep1_result.json next to this file).  Never pickles flint objects.
"""
import json, math, os, sys, time
from fractions import Fraction
import numpy as np
from flint import arb, arb_mat, fmpq, fmpz, fmpq_poly, fmpz_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
T0 = time.time()


def log(*a):
    print("[%7.1fs]" % (time.time() - T0), *a, flush=True)


# ----------------------------------------------------------------------------------------------
# Chebyshev values and derivatives by recurrence (works for float and arb)
# ----------------------------------------------------------------------------------------------
def cheb_all(x, N, one):
    T = [one, x]
    d1 = [0 * one, one]
    d2 = [0 * one, 0 * one]
    for n in range(1, N):
        T.append(2 * x * T[n] - T[n - 1])
        d1.append(2 * T[n] + 2 * x * d1[n] - d1[n - 1])
        d2.append(4 * d1[n] + 2 * x * d2[n] - d2[n - 1])
    return T[:N + 1], d1[:N + 1], d2[:N + 1]


def build_ops(Ny, K, mode):
    """Collocation operators. mode 'f' (float) or 'a' (arb). Returns L, U (lists of rows), point data."""
    if mode == 'f':
        one = 1.0
        pi = math.pi
        cosf = math.cos
    else:
        one = arb(1)
        pi = arb.pi()
        cosf = lambda t: t.cos()
    ys = [(one + cosf(pi * i / Ny)) / 2 for i in range(Ny + 1)]
    ths = [pi * j / (2 * K) for j in range(K + 1)]
    cs = [cosf(t) for t in ths]
    # exact endpoints
    ys[0] = one * 1
    ys[Ny] = 0 * one
    cs[0] = one * 1
    cs[K] = 0 * one
    ang = []
    for c in cs:
        Tc, dTc, _ = cheb_all(c, 2 * K, one)
        T2 = [Tc[2 * k] for k in range(K + 1)]
        A = [-c * dTc[2 * k] - 4 * k * k * Tc[2 * k] for k in range(K + 1)]
        ang.append((T2, A))
    rows_L, rows_U, pts = [], [], []
    for y in ys:
        x = 2 * y - one
        T, d1, d2 = cheb_all(x, Ny, one)
        Ty = T
        Ty1 = [2 * v for v in d1]
        Ty2 = [4 * v for v in d2]
        om = one - y
        # p-basis f = y T
        f1 = [Ty[n] + y * Ty1[n] for n in range(Ny + 1)]
        f2 = [2 * Ty1[n] + y * Ty2[n] for n in range(Ny + 1)]
        Rf = [om * f2[n] - 2 * f1[n] for n in range(Ny + 1)]
        fv = [y * Ty[n] for n in range(Ny + 1)]
        # q-basis g = h T, h = y^2 - y^3
        h = y * y - y * y * y
        h1 = 2 * y - 3 * y * y
        h2 = 2 * one - 6 * y
        g1 = [h1 * Ty[n] + h * Ty1[n] for n in range(Ny + 1)]
        g2 = [h2 * Ty[n] + 2 * h1 * Ty1[n] + h * Ty2[n] for n in range(Ny + 1)]
        Rg = [om * g2[n] - 2 * g1[n] for n in range(Ny + 1)]
        gv = [h * Ty[n] for n in range(Ny + 1)]
        for jc, c in enumerate(cs):
            T2, A = ang[jc]
            rl = list(Rf)
            ru = list(fv)
            for k in range(1, K + 1):
                rl += [Rg[n] * T2[k] + Ty[n] * A[k] for n in range(Ny + 1)]
                ru += [gv[n] * T2[k] for n in range(Ny + 1)]
            rows_L.append(rl)
            rows_U.append(ru)
            sh = (one - c * c) * y * y * om * om
            pts.append((y, c, sh))
    return rows_L, rows_U, pts


def solve_float(Ny, K, tau, s1, w0=None, iters=40):
    L, U, pts = build_ops(Ny, K, 'f')
    L = np.array(L)
    U = np.array(U)
    y = np.array([p[0] for p in pts])
    sh = np.array([p[2] for p in pts])
    n = L.shape[1]
    w = np.zeros(n) if w0 is None else w0.copy()
    for it in range(iters):
        D = 1 - y / 2 + tau * (1 - y) * (U @ w)
        F = L @ w + s1 * sh * D ** -7
        J = L - (7 * s1 * tau * sh * (1 - y) * D ** -8)[:, None] * U
        dw = np.linalg.solve(J, F)
        w = w - dw
        if np.max(np.abs(dw)) < 1e-15:
            break
    D = 1 - y / 2 + tau * (1 - y) * (U @ w)
    F = L @ w + s1 * sh * D ** -7
    return w, float(np.max(np.abs(F))), (L, U, y, sh)


def refine_arb(Ny, K, tau_q, s1_q, wf, fl, prec, maxit=12):
    """Iterative refinement: residual in arb at prec bits, correction from the double Jacobian."""
    ctx.prec = prec
    L, U, y, sh = fl
    Lr, Ur, pts = build_ops(Ny, K, 'a')
    La = arb_mat(Lr)
    Ua = arb_mat(Ur)
    n = len(wf)
    tau = arb(tau_q)
    s1 = arb(s1_q)
    w = arb_mat([[arb(float(v))] for v in wf])
    tf = fq_float(tau_q)
    s1f = fq_float(s1_q)
    hist = []
    for it in range(maxit):
        Lw = La * w
        Uw = Ua * w
        F = []
        for i, (yy, cc, ss) in enumerate(pts):
            D = 1 - yy / 2 + tau * (1 - yy) * Uw[i, 0]
            F.append(Lw[i, 0] + s1 * ss / D ** 7)
        Ff = np.array([float(f.mid()) for f in F])
        res = float(np.max(np.abs(Ff)))
        hist.append(res)
        if res < 2.0 ** (-(prec - 40)):
            break
        wd = np.array([float(w[i, 0].mid()) for i in range(n)])
        D = 1 - y / 2 + tf * (1 - y) * (U @ wd)
        J = L - (7 * s1f * tf * sh * (1 - y) * D ** -8)[:, None] * U
        dw = np.linalg.solve(J, Ff)
        w = arb_mat([[(w[i, 0] - arb(float(dw[i]))).mid()] for i in range(n)])
    return [w[i, 0].mid() for i in range(n)], hist


def to_fmpq(a, bits):
    """exact dyadic rounding of an exact arb midpoint to a multiple of 2^-bits."""
    m, e = a.man_exp()
    m = int(m)
    e = int(e)
    if e >= -bits:
        return fmpq(m * 2 ** (e + bits), 2 ** bits) if e + bits >= 0 else fmpq(m, 2 ** (-e))
    sh = -bits - e
    mm = (m + (1 << (sh - 1))) >> sh  # round to nearest
    return fmpq(mm, 2 ** bits)


# ----------------------------------------------------------------------------------------------
# Exact polynomial certificate
# ----------------------------------------------------------------------------------------------
Yp = fmpq_poly([0, 1])
OMY = fmpq_poly([1, -1])
Cp = fmpq_poly([0, 1])


def cheb_polys(var, N):
    T = [fmpq_poly([1]), var]
    for n in range(1, N):
        T.append(2 * var * T[n] - T[n - 1])
    return T[:N + 1]


def Rop(f):
    return OMY * f.derivative().derivative() - 2 * f.derivative()


def pack(terms, W):
    dy = max(f.degree() for f, g in terms if f.degree() >= 0)
    co = [fmpq(0)] * ((dy + 1) * W)
    for f, g in terms:
        gc = g.coeffs()
        assert len(gc) <= W
        for i, fi in enumerate(f.coeffs()):
            if fi == 0:
                continue
            base = i * W
            for j, gj in enumerate(gc):
                if gj != 0:
                    co[base + j] += fi * gj
    return fmpq_poly(co)


_AX = {}


def Ax_int(d):
    """integer matrix A with x^m = 2^-d sum_k A[k,m] T_k(x), k,m = 0..d."""
    if d in _AX:
        return _AX[d]
    A = [[0] * (d + 1) for _ in range(d + 1)]
    for m in range(d + 1):
        sc = 1 << (d - m)
        for l in range(m + 1):
            A[abs(m - 2 * l)][m] += sc * math.comb(m, l)
    _AX[d] = A
    return A


def cheb_abs_sum(P, W):
    """Exact rational S = sum |Chebyshev coefficients| of the packed bivariate P(y,c) on [0,1]x[-1,1]."""
    num = P.numer().coeffs()
    den = int(P.denom())
    L = len(num)
    dy = (L - 1) // W
    rows = dy + 1
    flat = [0] * (rows * W)
    for idx, v in enumerate(num):
        flat[idx] = int(v)
    # drop all-zero top columns in c
    dc = 0
    for i in range(rows):
        for j in range(W - 1, -1, -1):
            if flat[i * W + j] != 0:
                dc = max(dc, j)
                break
    M = fmpz_mat(rows, dc + 1, [flat[i * W + j] for i in range(rows) for j in range(dc + 1)])
    Ax = Ax_int(dy)
    # B[m][i] = 2^(dy-i) C(i,m)
    B = fmpz_mat(dy + 1, dy + 1, [((1 << (dy - i)) * math.comb(i, m) if m <= i else 0)
                                   for m in range(dy + 1) for i in range(dy + 1)])
    Ay = fmpz_mat(dy + 1, dy + 1, [Ax[k][m] for k in range(dy + 1) for m in range(dy + 1)]) * B
    Ac = Ax_int(dc)
    AcT = fmpz_mat(dc + 1, dc + 1, [Ac[l][j] for j in range(dc + 1) for l in range(dc + 1)])
    C = Ay * M * AcT
    tot = 0
    for i in range(C.nrows()):
        for j in range(C.ncols()):
            v = int(C[i, j])
            tot += -v if v < 0 else v
    return fmpq(tot, den * (1 << (2 * dy + dc)))


def build_exact(wq, Ny, K):
    p = wq[:Ny + 1]
    qs = [wq[(Ny + 1) * k:(Ny + 1) * (k + 1)] for k in range(1, K + 1)]
    Ty = cheb_polys(fmpq_poly([-1, 2]), Ny)
    Tc = cheb_polys(Cp, 2 * K)
    one = fmpq_poly([1])
    pp = sum((p[n] * Ty[n] for n in range(Ny + 1)), fmpq_poly([0]))
    f0 = Yp * pp
    h = fmpq_poly([0, 0, 1, -1])
    Lh_terms = [(Rop(f0), one)]
    w_terms = [(f0, one)]
    omw_terms = [(OMY * f0, one)]
    for k in range(1, K + 1):
        qk = sum((qs[k - 1][n] * Ty[n] for n in range(Ny + 1)), fmpq_poly([0]))
        gk = h * qk
        T2 = Tc[2 * k]
        Ak = -Cp * T2.derivative() - 4 * k * k * T2
        Lh_terms += [(Rop(gk), T2), (qk, Ak)]
        w_terms.append((gk, T2))
        omw_terms.append((OMY * gk, T2))
    p_at_0 = sum((p[n] * (1 if n % 2 == 0 else -1) for n in range(Ny + 1)), fmpq(0))  # T_n(-1)
    return Lh_terms, w_terms, omw_terms, p_at_0


def certify(wq, Ny, K, tau, s1, tau_slab=None):
    """Returns dict with exact rationals: eta, p(0), and checks. tau_slab: certify uniformly for tau in [0,tau_slab]
    (then wq must be the tau-independent ansatz and tau is ignored)."""
    t = time.time()
    Lh_terms, w_terms, omw_terms, p0 = build_exact(wq, Ny, K)
    W = 16 * K + 3
    one = fmpq_poly([1])
    Dbase = fmpq_poly([1, fmpq(-1, 2)])
    Sh = (s1 * Yp ** 2 * OMY ** 2, fmpq_poly([1, 0, -1]))
    LhP = pack(Lh_terms, W)
    ShP = pack([Sh], W)
    # sup |w| bound (for D >= 1/2):  D = 1/2 + (1-y)(1/2 + tau w)
    Sw = cheb_abs_sum(pack(w_terms, 2 * K + 1), 2 * K + 1)
    if tau_slab is None:
        DP = pack([(Dbase, one)] + [(tau * f, g) for f, g in omw_terms], W)
        Q = DP ** 7 * LhP + ShP
        S = cheb_abs_sum(Q, W)
        Sparts = None
        taumax = tau
    else:
        D0P = pack([(Dbase, one)], W)
        omwP = pack(omw_terms, W)
        S = fmpq(0)
        Sparts = []
        for m in range(8):
            if m == 0:
                Qm = D0P ** 7 * LhP + ShP
            else:
                Qm = math.comb(7, m) * D0P ** (7 - m) * omwP ** m * LhP
            Sm = cheb_abs_sum(Qm, W)
            Sparts.append(Sm)
            S += Sm * tau_slab ** m
        taumax = tau_slab
    eta_exact = 64 * S * (1 + fmpq(1, 2 ** 20))
    # round eta up to a dyadic with 64 significant bits
    if eta_exact > 0:
        e = math.floor(math.log2(fq_float(eta_exact))) - 64
        sc = fmpq(2) ** (-e) if e < 0 else fmpq(1, 2 ** e)
        v = eta_exact * sc
        eta = fmpq(-((-v.p) // v.q)) / sc
    else:
        eta = fmpq(0)
    assert eta >= eta_exact
    # D >= 1/2 for w_N and w- (w- = w_N - eta y): need 1/2 - taumax (Sw + eta) > 0
    D_ok = (fmpq(1, 2) - taumax * (Sw + eta)) > 0
    return dict(S=S, eta=eta, p0=p0, Sw=Sw, D_ok=bool(D_ok), Sparts=Sparts, secs=time.time() - t,
                degQ=None)


# ----------------------------------------------------------------------------------------------
# decimal output with outward rounding (exact rationals)
# ----------------------------------------------------------------------------------------------
def dec_floor(q, nd):
    s = 1 if q >= 0 else -1
    a = abs(q)
    if a == 0:
        return "0"
    e = math.floor(math.log10(fq_float(a)))
    scale = nd - 1 - e
    if scale >= 0:
        v = q * fmpq(10) ** scale
    else:
        v = q / fmpq(10) ** (-scale)
    iv = v.p // v.q  # floor
    return fmt_int(int(iv), scale)


def dec_ceil(q, nd):
    return neg(dec_floor(-q, nd))


def neg(s):
    return s[1:] if s.startswith("-") else ("-" + s if s != "0" else s)


def fmt_int(iv, scale):
    sgn = "-" if iv < 0 else ""
    iv = abs(iv)
    s = str(iv)
    if scale <= 0:
        return sgn + s + "0" * (-scale)
    if len(s) <= scale:
        s = "0" * (scale - len(s) + 1) + s
    return sgn + s[:-scale] + "." + s[-scale:]


def fq_float(q):
    return float(Fraction(int(q.p), int(q.q)))


# ----------------------------------------------------------------------------------------------
# non-rigorous cross-check: integral route (double, Gauss-Legendre)
# ----------------------------------------------------------------------------------------------
def integral_route(w, Ny, K, tau, s1):
    xg, wg = np.polynomial.legendre.leggauss(120)
    yq = (xg + 1) / 2
    wy = wg / 2
    cq, wc = xg, wg / 2
    p = w[:Ny + 1]
    tot = 0.0
    for iy, yv in enumerate(yq):
        xx = 2 * yv - 1
        Tn = np.polynomial.chebyshev.chebvander(np.array([xx]), Ny)[0]
        val_p = Tn @ p
        qv = [Tn @ w[(Ny + 1) * k:(Ny + 1) * (k + 1)] for k in range(1, K + 1)]
        th = np.arccos(cq)
        wv = yv * val_p + yv ** 2 * (1 - yv) * sum(qv[k - 1] * np.cos(2 * k * th) for k in range(1, K + 1))
        D = 1 - yv / 2 + tau * (1 - yv) * wv
        inner = np.sum(wc * (1 - cq ** 2) * D ** -7)
        tot += wy[iy] * yv ** 2 * (1 - yv) ** 3 * inner
    # delta = 2 lim r u = 2 * (9 J^2/4) * tot, 9J^2/4 = s1*tau
    return 2 * s1 * tau * tot


# ----------------------------------------------------------------------------------------------
def run_point(j_q, Ny, K, prec, s1, label):
    tau = j_q * j_q
    tf = fq_float(tau)
    s1f = fq_float(s1)
    wf, resf, fl = solve_float(Ny, K, tf, s1f)
    wa, hist = refine_arb(Ny, K, tau, s1, wf, fl, prec)
    wq = [to_fmpq(a, prec - 16) for a in wa]
    c = certify(wq, Ny, K, tau, s1)
    p0, eta = c['p0'], c['eta']
    lo = 2 * tau * (p0 - eta)
    hi = 2 * tau * (p0 + eta)
    ok = c['D_ok'] and eta > 0
    rel = fq_float((hi - lo) / lo)
    integ = integral_route(np.array([fq_float(q) for q in wq]), Ny, K, tf, s1f)
    out = dict(label=label, j=str(j_q), J=str(j_q / 4), Ny=Ny, K=K, prec_bits=prec,
               float_collocation_residual=resf, arb_refine_residual_history=hist,
               S_cheb_abs_sum_Q=fq_float(c['S']), eta=fq_float(eta), D_ge_half_certified=c['D_ok'],
               certified=bool(ok), delta_lo=dec_floor(lo, 30), delta_hi=dec_ceil(hi, 30),
               rel_width=rel, integral_route_double=repr(integ), cert_seconds=round(c['secs'], 2))
    log("%s j=%s Ny=%d K=%d prec=%d: eta=%.3e rel.width=%.3e  delta in [%s, %s]  certified=%s (%.1fs)"
        % (label, j_q, Ny, K, prec, fq_float(eta), rel, out['delta_lo'][:22], out['delta_hi'][:22], ok, c['secs']))
    return out, (lo, hi)


def run_slab(j0_q, Ny, K, prec, s1):
    """uniform certificate for all j in (0, j0]: w = solution of the tau = 0 (linear) problem."""
    tau0 = j0_q * j0_q
    wf, resf, fl = solve_float(Ny, K, 0.0, fq_float(s1))
    wa, hist = refine_arb(Ny, K, fmpq(0), s1, wf, fl, prec)
    wq = [to_fmpq(a, prec - 16) for a in wa]
    c = certify(wq, Ny, K, None, s1, tau_slab=tau0)
    p0, eta = c['p0'], c['eta']
    lo = 2 * (p0 - eta)
    hi = 2 * (p0 + eta)
    ok = c['D_ok'] and eta > 0
    contains = bool(lo <= fmpq(1, 40) <= hi) and ok
    out = dict(j0=str(j0_q), Ny=Ny, K=K, prec_bits=prec, eta=fq_float(eta),
               S_parts=[fq_float(s) for s in c['Sparts']], certified=bool(ok),
               delta_over_j2_lo=dec_floor(lo, 30), delta_over_j2_hi=dec_ceil(hi, 30),
               rel_width=fq_float((hi - lo) / lo), contains_1_over_40=contains,
               arb_refine_residual_history=hist)
    log("K1 slab j in (0,%s] Ny=%d K=%d: delta/j^2 in [%s, %s]  contains 1/40: %s"
        % (j0_q, Ny, K, out['delta_over_j2_lo'][:24], out['delta_over_j2_hi'][:24], contains))
    return out, (lo, hi)


def null_test(j_q, Ny, K, prec, s1, bump):
    """perturbed ansatz (p_0 coefficient + bump): the certifier must respond with eta >~ bump-size."""
    tau = j_q * j_q
    wf, resf, fl = solve_float(Ny, K, fq_float(tau), fq_float(s1))
    wa, hist = refine_arb(Ny, K, tau, s1, wf, fl, prec)
    wq = [to_fmpq(a, prec - 16) for a in wa]
    wq[0] += bump
    c = certify(wq, Ny, K, tau, s1)
    return dict(bump=str(bump), eta=fq_float(c['eta']), S=fq_float(c['S']))


def main():
    s1 = fmpq(9, 64)
    res = dict(seat="rep1", object="Bowen-York spinning puncture, mu = 1, delta = M_ADM - 1, j = 4J",
               method="y = 1/(1+r) Chebyshev x cos(2k theta) collocation; arb refinement; exact fmpq certificate")
    # controls: analytic j->0 coefficient from the integral route (double quadrature, positive control of formula)
    res['control_integral_tau0_double'] = repr(integral_route(np.zeros((12 + 1) * 3), 12, 2, 1.0, 9 / 64))
    log("control: integral route at tau=0, delta/j^2 (double) =", res['control_integral_tau0_double'])

    # main targets, two resolutions + two precisions
    runs = []
    encl = {}
    for (Jlabel, j_q) in (("J=0.0625", fmpq(1, 4)), ("J=0.25", fmpq(1))):
        for (Ny, K, prec) in ((40, 8, 224), (56, 12, 320)):
            o, (lo, hi) = run_point(j_q, Ny, K, prec, s1, Jlabel)
            runs.append(o)
            encl.setdefault(Jlabel, []).append((lo, hi, o))
    res['runs'] = runs
    final = {}
    for Jl, lst in encl.items():
        lo = max(x[0] for x in lst)
        hi = min(x[1] for x in lst)
        overlap = lo <= hi
        final[Jl] = dict(delta_lo=dec_floor(lo, 30), delta_hi=dec_ceil(hi, 30),
                         M_ADM_lo=dec_floor(1 + lo, 30), M_ADM_hi=dec_ceil(1 + hi, 30),
                         abs_width=fq_float(hi - lo), rel_width=fq_float((hi - lo) / lo),
                         all_certified=all(x[2]['certified'] for x in lst), resolutions_overlap=bool(overlap))
        log("FINAL %s: delta in [%s, %s], rel width %.2e" % (Jl, final[Jl]['delta_lo'], final[Jl]['delta_hi'],
                                                           final[Jl]['rel_width']))
    res['delta_J0.0625'] = final['J=0.0625']
    res['delta_J0.25'] = final['J=0.25']

    # K1: uniform slab j in (0, j0] and small fixed j
    k1 = {}
    sl, _ = run_slab(fmpq(1, 2 ** 16), 28, 6, 192, s1)
    k1['slab'] = sl
    sl2, _ = run_slab(fmpq(1, 2 ** 16), 40, 8, 256, s1)
    k1['slab_resolution2'] = sl2
    small = []
    for j_q in (fmpq(1, 64), fmpq(1, 256)):
        o, (lo, hi) = run_point(j_q, 36, 6, 224, s1, "small-j")
        o['delta_over_j2_lo'] = dec_floor(lo / (j_q * j_q), 30)
        o['delta_over_j2_hi'] = dec_ceil(hi / (j_q * j_q), 30)
        small.append(o)
    k1['small_j_points'] = small
    k1['contains_1_over_40'] = bool(sl['contains_1_over_40'] and sl2['contains_1_over_40'])
    res['K1'] = k1

    # K2: source 9/4 -> 9/4 (1 + 1e-3) at J = 0.25
    s1p = s1 * (1 + fmpq(1, 1000))
    k2runs = []
    pl = []
    for (Ny, K, prec) in ((40, 8, 224), (56, 12, 320)):
        o, (lo, hi) = run_point(fmpq(1), Ny, K, prec, s1p, "K2 perturbed")
        k2runs.append(o)
        pl.append((lo, hi))
    plo = max(x[0] for x in pl)
    phi = min(x[1] for x in pl)
    blo = max(x[0] for x in encl['J=0.25'])
    bhi = min(x[1] for x in encl['J=0.25'])
    res['K2'] = dict(perturbed_runs=k2runs, delta_perturbed_lo=dec_floor(plo, 30), delta_perturbed_hi=dec_ceil(phi, 30),
                     shift_lo=dec_floor(plo - bhi, 20), shift_hi=dec_ceil(phi - blo, 20),
                     rel_shift_lo=dec_floor((plo - bhi) / bhi, 12), rel_shift_hi=dec_ceil((phi - blo) / blo, 12))
    log("K2: delta(0.25) moves by [%s, %s] (relative [%s, %s])" % (res['K2']['shift_lo'], res['K2']['shift_hi'],
                                                                  res['K2']['rel_shift_lo'], res['K2']['rel_shift_hi']))

    # null: perturbed ansatz must be penalised by the certifier
    res['null_bumped_ansatz'] = [null_test(fmpq(1), 40, 8, 224, s1, fmpq(1, 10 ** 6)),
                                 null_test(fmpq(1), 40, 8, 224, s1, fmpq(1, 10 ** 9))]
    log("null (bumped ansatz):", res['null_bumped_ansatz'])
    res['total_seconds'] = round(time.time() - T0, 1)
    with open(os.path.join(HERE, "rep1_result.json"), "w") as f:
        json.dump(res, f, indent=1)
    log("wrote rep1_result.json")


if __name__ == "__main__":
    main()
