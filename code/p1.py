"""BHSLICE P1 certifier -- computer-certified (arb ball arithmetic) sub/super-solution sandwich for the Hamiltonian
constraint of a single spinning Bowen-York puncture, bare mass mu = 1, J in [0, 1/4] (j = 4J in [0, 1]).

Problem (DESIGN.md 1.4): s = a/(r+a), c = cos(theta), a = mu/2, j = J/a^2 = 4J. Write u = y U with y = j^2. Then
   G[u] = y s^2 (1-s) Gh[U],
   Gh[U] := sum_l [ (1-s) U_l'' - 2 U_l' - l(l+1) U_l/(s^2 (1-s)) ] P_l(c) + Sh(U),
   Sh(U) := (9/4) (1-c^2) s^2 (1-s)^2 (1 + (1-s) y U)^-7,
with U = s p(s) P_0 + sum_{l=2,4,..} s^2 (1-s) q_l(s) P_l(c) (p, q_l polynomials; every term above is a polynomial on the
closed square). Super-solution <=> Gh[U+] < 0, sub-solution <=> Gh[U-] > 0 (for y > 0; the factor y s^2 (1-s) is >= 0).
Sandwich: U+- = U_N +- eta0 s, U_N a Chebyshev fit (in X = 2s-1 and in the slab variable eta of y) of double-precision
collocation solutions. Comparison lemma: DESIGN 2.1 (Lemma C). u- >= 0 is certified as V := U-/s > 0.

Certified bound of an analytic F on a box B = [sa,sb] x [-1,1] x [ya,yb] (CLW seat-E method, adapted):
 (1) node values of F at tensor Chebyshev (2nd kind) points, in arb (exact nodes enclosed);
 (2) Chebyshev coefficients by the exact DCT matrix (arb); polynomial parts of F are reproduced exactly (node counts
     exceed their degrees); the analytic part Sh (or the quadrature integrand) adds
     |F - p| <= Lam * sum_v 4 M rho_v^-n_v/(rho_v - 1),  Lam = prod_v (1 + Lambda_{n_v}),  Lambda_n <= (2/pi) log(n+1) + 1
     (Trefethen ATAP Thm 8.2 per variable + Lebesgue constants, Thm 15.2, as cited in the CLW engine), with M a rigorous
     bound of the analytic part on the complex polyellipse (|T_n(x)| <= rho_x^n, |P_l(c)| <= rho_c^l on E_rho_c,
     |1-c^2| <= R_c^2, |(1+w)^-7| <= (1-|w|)^-7);
 (3) sup F <= c_0 + sum_{k != 0} |c_k| + E, inf F >= c_0 - sum |c_k| - E.
Failure of (3) triggers a split of the box (work queue); any node value with a certified wrong sign is reported as a
certified VIOLATION (exact point evaluation, no interpolation involved).

Usage (every run under memcap):
  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 p1.py --probe            (timing/feasibility probe)
  ... p1.py --controls controls_p1.json                                              (controls, < 10 min)
  ... p1.py --main main_p1.json                                                      (refuses without GO-P1; resumable)
No pickle anywhere; arb values are serialised as outward-rounded float pairs plus display strings.
"""
import sys, os, json, time, math, hashlib, argparse
import numpy as np
from numpy.polynomial import legendre as npleg, chebyshev as npch
from flint import arb, arb_mat, ctx

ctx.prec = 128
HERE = os.path.dirname(os.path.abspath(__file__))

# ------------------------------------------------------------------ configuration (justified in BUILD-P1.md)
NS_SOLVE, NL_SOLVE = 48, 8          # double collocation solve: s-degree, Legendre modes (l <= 16)
NB = 28                             # kept degree of the reduced s-polynomials p, q_l
LMAX = 12                           # kept Legendre degree
J_EDGES = [0.0, 0.1] + [round(0.1 + 0.05 * k, 10) for k in range(1, 19)]   # 19 slabs (DESIGN sec 3)
Y_EDGES = [j * j for j in J_EDGES]  # float (dyadic) slab edges in y = j^2; consecutive slabs share the same float
ETA_CAP = 1.0e-9                    # need-based cap: slope-route relative width 2 eta0/p(0) ~ 80 eta0 <= 1e-7
ETA_FLOOR = 2.0e-13
MAX_DEPTH = 6                       # box splits per slab (work queue bound)
REPORT_J = ['0.0025', '0.025', '0.0625', '0.125', '0.25']
TIME_BUDGET = 520.0
PI = arb.pi()


# ------------------------------------------------------------------ double-precision collocation solver (not a proof)
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    cc = np.ones(N + 1); cc[0] = cc[-1] = 2; cc *= (-1) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(cc, 1 / cc) / (X - X.T + np.eye(N + 1))
    D -= np.diag(D.sum(1))
    return x, D


def wfun(l, s):
    if l == 0:
        return s, np.ones_like(s), np.zeros_like(s)
    return s ** 2 * (1 - s), 2 * s - 3 * s ** 2, 2 - 6 * s


def setup(Ns=NS_SOLVE, NL=NL_SOLVE):
    x, D = cheb(Ns)
    s = (1 - x) / 2
    Ds = -2 * D; Ds2 = Ds @ Ds
    Lv = 2 * np.arange(NL + 1)
    cq, wq = npleg.leggauss(4 * NL + 24)
    P = np.array([npleg.legval(cq, np.eye(l + 1)[l]) for l in Lv])
    proj = (2 * Lv[:, None] + 1) / 2 * P * wq[None, :]
    A, W = [], []
    for l in Lv:
        w, w1, w2 = wfun(l, s)
        Al = np.diag((1 - s) * w2) + np.diag(2 * (1 - s) * w1) @ Ds + np.diag((1 - s) * w) @ Ds2 - 2 * np.diag(w1) - 2 * np.diag(w) @ Ds
        if l > 0:
            Al -= l * (l + 1) * np.eye(Ns + 1)
        A.append(Al); W.append(w)
    return dict(s=s, Lv=Lv, cq=cq, P=P, proj=proj, A=A, W=np.array(W), Tm=npch.chebvander(2 * s - 1, Ns), Ns=Ns, NL=NL)


def solveU(y, S, Q0=None, tol=1e-14, itmax=25):
    s, P, proj, W = S['s'], S['P'], S['proj'], S['W']
    nl, ns = S['NL'] + 1, S['Ns'] + 1
    Q = np.zeros((nl, ns)) if Q0 is None else Q0.copy()
    pref = 2.25 * (1 - S['cq'] ** 2)[None, :] * (s ** 2 * (1 - s) ** 2)[:, None]
    for it in range(itmax):
        U = (W * Q).T @ P
        base = 1 + (1 - s)[:, None] * y * U
        Sh = pref * base ** -7
        dSh = -7 * (1 - s)[:, None] * y * pref * base ** -8
        R = np.array([S['A'][i] @ Q[i] + Sh @ proj[i] for i in range(nl)])
        J = np.zeros((nl * ns, nl * ns))
        for i in range(nl):
            J[i * ns:(i + 1) * ns, i * ns:(i + 1) * ns] += S['A'][i]
            for m in range(nl):
                J[i * ns:(i + 1) * ns, m * ns:(m + 1) * ns] += np.diag(((dSh * P[m][None, :]) @ proj[i]) * W[m])
        dQ = np.linalg.solve(J, -R.ravel()).reshape(nl, ns)
        Q += dQ
        if np.max(np.abs(dQ)) < tol:
            break
    return Q


def fit_slab(ya, yb, ny, S=None):
    """Chebyshev fit of U on the slab y in [ya, yb]: a[l] of shape (NB+1, ny+1), U_l = w_l(s) sum a[l][n,g] T_n(2s-1) T_g(eta)."""
    S = S or setup()
    eta = np.cos(np.pi * np.arange(ny + 1) / ny)
    ym, yh = (ya + yb) / 2, (yb - ya) / 2
    B = []
    Q = None
    for e in eta[::-1]:
        Q = solveU(ym + yh * e, S, Q)
        B.append(np.linalg.solve(S['Tm'], Q.T).T)      # (nl, Ns+1)
    B = np.array(B[::-1])                               # (ny+1, nl, Ns+1) at eta nodes
    Vy = npch.chebvander(eta, ny)
    C = np.linalg.solve(Vy, B.reshape(ny + 1, -1)).reshape(B.shape)   # (ny+1 [g], nl, Ns+1)
    a = {}
    for i, l in enumerate(S['Lv']):
        if l > LMAX:
            continue
        a[int(l)] = np.ascontiguousarray(C[:, i, :NB + 1].T)            # (NB+1, ny+1)
    return dict(ya=ya, yb=yb, ny=ny, a=a)


# ------------------------------------------------------------------ float evaluation (only to choose eta0; never a certificate)
def float_Gh(model, s, c, y):
    ya, yb = model['ya'], model['yb']
    eta = 0.0 if yb == ya else (y - (ya + yb) / 2) / ((yb - ya) / 2)
    Tg = npch.chebvander(np.array([eta]), model['ny'])[0]
    X = 2 * s - 1
    U = np.zeros((len(s), len(c))); Lh = np.zeros_like(U)
    for l, a in model['a'].items():
        b = a @ Tg
        q = npch.chebval(X, b); q1 = 2 * npch.chebval(X, npch.chebder(b)); q2 = 4 * npch.chebval(X, npch.chebder(b, 2))
        w, w1, w2 = wfun(l, s)
        Ul = w * q
        L = (1 - s) * (w2 * q + 2 * w1 * q1 + w * q2) - 2 * (w1 * q + w * q1) - (l * (l + 1) * q if l > 0 else 0)
        Pl = npleg.legval(c, np.eye(l + 1)[l])
        U += np.outer(Ul, Pl); Lh += np.outer(L, Pl)
    Sh = 2.25 * (1 - c ** 2)[None, :] * (s ** 2 * (1 - s) ** 2)[:, None] * (1 + (1 - s)[:, None] * y * U) ** -7
    return Lh + Sh


def float_resid(model, ns=241, nc=41, nyg=None):
    s = np.linspace(0, 1, ns); c = np.linspace(-1, 1, nc)
    ya, yb = model['ya'], model['yb']
    ys = np.linspace(ya, yb, nyg or (2 * model['ny'] + 3))
    return max(float(np.max(np.abs(float_Gh(model, s, c, y)))) for y in ys)


# ------------------------------------------------------------------ arb helpers
def A(x):
    return x if isinstance(x, arb) else arb(x)


def is_pt(ybox):
    """a y-'box' is a point when both ends are the same object (arb point) or equal floats."""
    return ybox[0] is ybox[1] or (not isinstance(ybox[0], arb) and not isinstance(ybox[1], arb) and ybox[0] == ybox[1])


def nodes(n):
    return [arb(0)] if n == 0 else [(arb(i) / n).cos_pi() for i in range(n + 1)]


_DCT = {}


def dct_mat(n):
    if n not in _DCT:
        if n == 0:
            _DCT[n] = arb_mat([[1]])
        else:
            M = []
            for jj in range(n + 1):
                row = []
                for i in range(n + 1):
                    v = (arb(i * jj % (2 * n)) / n).cos_pi() * 2 / n
                    if i in (0, n): v = v / 2
                    if jj in (0, n): v = v / 2
                    row.append(v)
                M.append(row)
            _DCT[n] = arb_mat(M)
    return _DCT[n]


def to_mat(T2):
    r = T2.shape[0]; cols = int(np.prod(T2.shape[1:])) if T2.ndim > 1 else 1
    return arb_mat(r, cols, list(T2.reshape(r, cols).ravel()))


def from_mat(B, shape):
    return np.array([B[i, j] for i in range(B.nrows()) for j in range(B.ncols())], dtype=object).reshape(shape)


def mode_mul(Mat, T, axis):
    T2 = np.moveaxis(T, axis, 0)
    sh = T2.shape
    B = Mat * to_mat(T2)
    return np.moveaxis(from_mat(B, (Mat.nrows(),) + sh[1:]), 0, axis)


def cheb_T012(Xs, N):
    E0, E1, E2 = [], [], []
    for X in Xs:
        T = [arb(1), X]; D1 = [arb(0), arb(1)]; D2 = [arb(0), arb(0)]
        for n in range(1, N):
            T.append(2 * X * T[n] - T[n - 1])
            D1.append(2 * T[n] + 2 * X * D1[n] - D1[n - 1])
            D2.append(4 * D1[n] + 2 * X * D2[n] - D2[n - 1])
        E0.append(T[:N + 1]); E1.append([2 * d for d in D1[:N + 1]]); E2.append([4 * d for d in D2[:N + 1]])
    return arb_mat(E0), arb_mat(E1), arb_mat(E2)


def cheb_T(xs, N):
    out = []
    for x in xs:
        T = [arb(1), x]
        for n in range(1, N):
            T.append(2 * x * T[n] - T[n - 1])
        out.append(T[:N + 1])
    return arb_mat(out)


def leg_P(cs, L):
    out = []
    for c in cs:
        P = [arb(1), c]
        for l in range(1, L):
            P.append(((2 * l + 1) * c * P[l] - l * P[l - 1]) / (l + 1))
        out.append(P[:L + 1])
    return out


def lam(n):
    return 1 + (2 / math.pi) * math.log(n + 1) + 1 if n > 0 else 1.0


def finite_all(arr):
    return all(x.is_finite() for x in arr.flat)


def outward(x):
    lo = float(x.lower()); hi = float(x.upper())
    lo = np.nextafter(np.nextafter(lo, -np.inf), -np.inf); hi = np.nextafter(np.nextafter(hi, np.inf), np.inf)
    return [float(lo), float(hi)]


def fromw(p):
    lo, hi = arb(p[0]), arb(p[1])
    return (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)   # ball containing [lo, hi]


def dstr(x, d=20):
    return x.str(d, radius=True)


# ------------------------------------------------------------------ model in arb form
class Model:
    def __init__(self, m, shift0=None):
        """m: fitted dict; shift0: list (len ny+1) of floats added to a[0][0, g] (U += s * shift(eta))."""
        self.ya, self.yb, self.ny = m['ya'], m['yb'], m['ny']
        self.ls = sorted(m['a'].keys())
        self.af = {l: np.array(m['a'][l], dtype=float).copy() for l in self.ls}
        if shift0 is not None:
            for g, v in enumerate(shift0):
                self.af[0][0, g] += v
        self.Am = {l: arb_mat([[arb(float(v)) for v in row] for row in self.af[l]]) for l in self.ls}
        self.Ym = (A(self.ya) + A(self.yb)) / 2; self.Yh = (A(self.yb) - A(self.ya)) / 2
        self.finite = all(np.all(np.isfinite(self.af[l])) for l in self.ls)

    def eta_of(self, y):
        return arb(0) if self.yb == self.ya else (A(y) - self.Ym) / self.Yh

    def p0(self, y):
        """p(s=0; y) = U_s(0; y) exactly (arb)."""
        Tg = cheb_T([self.eta_of(y)], self.ny)
        a = self.af[0]
        tot = arb(0)
        for n in range(a.shape[0]):
            for g in range(a.shape[1]):
                tot += arb(float(a[n, g])) * (1 if n % 2 == 0 else -1) * Tg[0, g]
        return tot

    def p0_coefs_eta(self):
        """float coefficients (in T_g(eta)) of p(0; y): used only to build the N1 control function (exact as floats)."""
        a = self.af[0]
        return [float(sum(a[n, g] * (1 if n % 2 == 0 else -1) for n in range(a.shape[0]))) for g in range(a.shape[1])]


# ------------------------------------------------------------------ rigorous polyellipse bounds
def rho_glob(m, h, rho, K=96):
    """Upper bound of the Bernstein-ellipse parameter of x = m + h z over z in the closed ellipse E_rho (m, h real arb)."""
    if h == 0:
        return arb(1)
    rho = A(rho)
    R = (rho + 1 / rho) / 2; r = (rho - 1 / rho) / 2
    Dmax = 1.0
    for k in range(K):
        th = PI * arb((k + 0.5) / K, 0.5 / K + 1e-17)
        a = m + h * R * th.cos(); b = h * r * th.sin()
        D = (((a - 1) ** 2 + b ** 2).sqrt() + ((a + 1) ** 2 + b ** 2).sqrt()) / 2
        Dmax = max(Dmax, float(np.nextafter(float(D.upper()), np.inf)))
    Du = arb(Dmax)
    return (Du + (Du * Du - 1).sqrt()).upper()


def poly_sum(af, rx, ry, rc_pow):
    """sum_{n,g} |a[n,g]| rx^n ry^g (arb), times rc_pow."""
    tot = arb(0)
    nN, nG = af.shape
    rxp = [rx ** n for n in range(nN)]; ryp = [ry ** g for g in range(nG)]
    absa = np.abs(af)
    for n in range(nN):
        rowsum = arb(0)
        for g in range(nG):
            if absa[n, g] != 0:
                rowsum += arb(float(absa[n, g])) * ryp[g]
        tot += rowsum * rxp[n]
    return tot * rc_pow


def ell_M(model, kind, sbox, ybox, rs, rc, ry):
    """Rigorous bound M of the analytic part of `kind` on the polyellipse (s in E_rs mapped to sbox, c in E_rc,
    y in E_ry mapped to ybox, or y the point ybox[0]). Returns arb or None if 1 + (1-s) y U may vanish."""
    sa, sb = sbox; ya, yb = ybox
    sm = (A(sa) + A(sb)) / 2; sh = (A(sb) - A(sa)) / 2
    Rs = (A(rs) + 1 / A(rs)) / 2; Rc = (A(rc) + 1 / A(rc)) / 2; Ry = (A(ry) + 1 / A(ry)) / 2
    rx = rho_glob(2 * sm - 1, 2 * sh, rs)
    if is_pt(ybox):
        rye = arb(1); ymax = abs(A(ya)).upper()
    else:
        ym = (A(ya) + A(yb)) / 2; yh = (A(yb) - A(ya)) / 2
        rye = rho_glob(model.eta_of(ym), yh / model.Yh, ry)
        ymax = (abs(ym) + yh * Ry).upper()
    smax = (abs(sm) + sh * Rs).upper(); omax = (abs(1 - sm) + sh * Rs).upper()
    Ub = smax * poly_sum(model.af[0], rx, rye, arb(1))
    for l in model.ls:
        if l > 0:
            Ub += smax * smax * omax * poly_sum(model.af[l], rx, rye, A(rc) ** l)
    Ub = Ub.upper()
    W = (omax * ymax * Ub).upper()
    if not (W < 1):
        return None
    if kind == 'Gh':
        M = arb(9) / 4 * Rc * Rc * smax ** 2 * omax ** 2 * (1 - W) ** -7
    elif kind == 'I':
        M = Rc * Rc * smax ** 2 * omax ** 3 * (1 - W) ** -7
    elif kind == 'I4':
        M = Rc * Rc * smax ** 2 * omax ** 3 * 7 * omax * Ub * (1 - W) ** -8
    else:
        M = arb(0)
    return M.upper()


def interp_err(M, ns, nc, ny, rs, rc, ry, ybox_is_point):
    if M == 0:
        return arb(0)
    Lam = lam(ns) * lam(nc) * (1.0 if ybox_is_point else lam(ny))
    e = 4 * M * A(rs) ** (-ns) / (A(rs) - 1) + 4 * M * A(rc) ** (-nc) / (A(rc) - 1)
    if not ybox_is_point:
        e += 4 * M * A(ry) ** (-ny) / (A(ry) - 1)
    return (arb(Lam) * arb(1 + 1e-12) * e).upper()


RS_C = [1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 16.0, 24.0]
RC_C = [1.5, 2.0, 3.0, 4.0, 6.0, 8.0]
RY_C = [2.0, 4.0, 8.0, 16.0, 32.0, 64.0, 128.0]


def plan(model, kind, sbox, ybox, tolE):
    """choose (rho, n) per variable: fewest nodes with E <= tolE; polynomial-degree minima for exact reproduction."""
    point = is_pt(ybox)
    if kind in ('Gh', 'V'):
        mins = (NB + 2, LMAX, 0 if point else model.ny)
    else:
        mins = (4, 4, 0 if point else model.ny)
    if kind == 'V':
        return dict(ns=mins[0], nc=mins[1], ny=mins[2], rs=2.0, rc=2.0, ry=2.0, M=arb(0), E=arb(0))
    best = None
    Mcache = {}
    for rs in RS_C:
        for rc in RC_C:
            for ry in (RY_C if not point else [2.0]):
                key = (rs, rc, ry)
                M = ell_M(model, kind, sbox, ybox, rs, rc, ry)
                if M is None:
                    continue
                Mf = float(M)
                Lam = 40.0 if point else 250.0          # generous (1+Lambda)^k for n <= ~200
                part = tolE / (3 * Lam)
                need = lambda r: max(0, math.ceil(math.log(4 * Mf / ((r - 1) * part)) / math.log(r)))
                ns = max(mins[0], need(rs)); nc = max(mins[1], need(rc)); ny = 0 if point else max(mins[2], need(ry))
                if nc % 2: nc += 1
                cost = (ns + 1) * (nc + 1) * (ny + 1)
                if ns > 400 or nc > 200 or ny > 60:
                    continue
                if best is None or cost < best[0]:
                    best = (cost, dict(ns=ns, nc=nc, ny=ny, rs=rs, rc=rc, ry=ry, M=M))
    if best is None:
        return None
    p = best[1]
    for _ in range(20):
        E = interp_err(p['M'], p['ns'], p['nc'], p['ny'], p['rs'], p['rc'], p['ry'], point)
        if E < tolE:
            break
        p['ns'] += 4; p['nc'] += 2
        if not point: p['ny'] += 1
    p['E'] = E
    return p


# ------------------------------------------------------------------ node values of F on a box
def box_values(model, kind, sbox, ybox, ns, nc, ny):
    sa, sb = sbox; ya, yb = ybox
    point = is_pt(ybox)
    sm = (A(sa) + A(sb)) / 2; sh = (A(sb) - A(sa)) / 2
    us = nodes(ns); ss = [sm + sh * u for u in us]
    Xs = [2 * s - 1 for s in ss]
    if point:
        yv = [A(ya)]; ny = 0
    else:
        ym = (A(ya) + A(yb)) / 2; yh = (A(yb) - A(ya)) / 2
        yv = [ym + yh * v for v in nodes(ny)]
    etas = [model.eta_of(y) for y in yv]
    cs = nodes(nc)
    E0, E1, E2 = cheb_T012(Xs, NB)
    Ty = cheb_T(etas, model.ny).transpose()          # (ny_model+1, ny+1)
    Pc = leg_P(cs, LMAX)                              # [k][l]
    nS, nY, nC = ns + 1, ny + 1, nc + 1
    L = len(model.ls)
    Ucols = np.empty((nS, nY, L), dtype=object); Lcols = np.empty((nS, nY, L), dtype=object); Vcols = np.empty((nS, nY, L), dtype=object)
    for li, l in enumerate(model.ls):
        Qv = model.Am[l] * Ty
        q = E0 * Qv; q1 = E1 * Qv; q2 = E2 * Qv
        for i in range(nS):
            s = ss[i]
            if l == 0:
                w, w1, w2 = s, arb(1), arb(0)
            else:
                w, w1, w2 = s * s * (1 - s), 2 * s - 3 * s * s, 2 - 6 * s
            for m in range(nY):
                qq, qa, qb = q[i, m], q1[i, m], q2[i, m]
                Ucols[i, m, li] = w * qq
                Lv = (1 - s) * (w2 * qq + 2 * w1 * qa + w * qb) - 2 * (w1 * qq + w * qa)
                if l > 0:
                    Lv -= l * (l + 1) * qq
                Lcols[i, m, li] = Lv
                Vcols[i, m, li] = qq if l == 0 else s * (1 - s) * qq
    Pm = arb_mat([[Pc[k][l] for k in range(nC)] for l in model.ls])     # (L, nC)
    def contract(cols):
        return from_mat(to_mat(cols.reshape(nS * nY, L)) * Pm, (nS, nY, nC))
    if kind == 'V':
        return contract(Vcols)
    U = contract(Ucols)
    F = np.empty((nS, nY, nC), dtype=object)
    if kind == 'Gh':
        Lh = contract(Lcols)
        for i in range(nS):
            s = ss[i]; f1 = arb(9) / 4 * s * s * (1 - s) ** 2; om = 1 - s
            for m in range(nY):
                y = yv[m]
                for k in range(nC):
                    F[i, m, k] = Lh[i, m, k] + f1 * (1 - cs[k] ** 2) * (1 + om * y * U[i, m, k]) ** -7
    elif kind == 'I':
        for i in range(nS):
            t = ss[i]; f1 = t * t * (1 - t) ** 3; om = 1 - t
            for m in range(nY):
                y = yv[m]
                for k in range(nC):
                    F[i, m, k] = f1 * (1 - cs[k] ** 2) * (1 + om * y * U[i, m, k]) ** -7
    elif kind == 'I4':
        for i in range(nS):
            t = ss[i]; f1 = t * t * (1 - t) ** 3; om = 1 - t
            for m in range(nY):
                y = yv[m]
                for k in range(nC):
                    if y == 0:
                        g = -7 * om * U[i, m, k]
                    else:
                        g = ((1 + om * y * U[i, m, k]) ** -7 - 1) / y
                    F[i, m, k] = f1 * (1 - cs[k] ** 2) * g
    return F


def coefs(F, ns, ny, nc):
    C = mode_mul(dct_mat(ns), F, 0)
    if ny > 0:
        C = mode_mul(dct_mat(ny), C, 1)
    return mode_mul(dct_mat(nc), C, 2)


# ------------------------------------------------------------------ sign certification with a work queue
def check_box(model, kind, sign, sbox, ybox, tolE):
    """sign=-1: certify sup F < 0; sign=+1: certify inf F > 0. Returns dict(status=PASS/FAIL/VIOLATION/INVALID, ...)."""
    if not model.finite:
        return dict(status='INVALID', why='non-finite coefficients')
    p = plan(model, kind, sbox, ybox, tolE)
    if p is None:
        return dict(status='INVALID', why='no admissible polyellipse (1+(1-s)yU may vanish) or node cap')
    ns, nc, ny = p['ns'], p['nc'], p['ny']
    F = box_values(model, kind, sbox, ybox, ns, nc, ny)
    if not finite_all(F):
        return dict(status='INVALID', why='non-finite node value')
    # certified violation at an exact node?
    for x in F.flat:
        if (sign < 0 and x > 0) or (sign > 0 and x < 0):
            return dict(status='VIOLATION', node_value=dstr(x, 8), ns=ns, nc=nc, ny=ny)
    C = coefs(F, ns, ny, nc)
    if not finite_all(C):
        return dict(status='INVALID', why='non-finite coefficient')
    c0 = C[0, 0, 0]
    tot = arb(0)
    for x in C.flat:
        tot += abs(x)
    rest = tot - abs(c0)
    E = p['E']
    if sign < 0:
        bnd = c0 + rest + E
        ok = bnd < 0
        val = bnd.upper()
    else:
        bnd = c0 - rest - E
        ok = bnd > 0
        val = bnd.lower()
    return dict(status='PASS' if ok else 'FAIL', bound=float(val), c0=float(c0.mid()), rest=float(rest.upper()),
                E=float(E), ns=ns, nc=nc, ny=ny, rs=p['rs'], rc=p['rc'], ry=p['ry'])


def certify(model, kind, sign, ybox, tolE, s_init=None, max_depth=MAX_DEPTH, log=None):
    """work queue over s-subintervals (then y-halves). Every popped box either passes or is split in two halves;
    the queue is bounded by depth, so it shrinks to empty or stops with FAIL."""
    queue = [(b, ybox, 0) for b in (s_init or [(0.0, 1.0)])]
    passed = []
    nb = 0
    while queue:
        sbox, yb_, depth = queue.pop(0)
        r = check_box(model, kind, sign, sbox, yb_, tolE)
        nb += 1
        if log: log(f"    box s={sbox} y=[{yb_[0]:.6g},{yb_[1]:.6g}] {kind}{'<0' if sign < 0 else '>0'}: {r['status']} "
                    + (f"bound={r.get('bound'):.3e} c0={r.get('c0'):.3e} rest={r.get('rest'):.2e} E={r.get('E'):.1e} n=({r.get('ns')},{r.get('nc')},{r.get('ny')})" if 'bound' in r else str({k: v for k, v in r.items() if k != 'status'})))
        if r['status'] == 'PASS':
            passed.append((sbox, yb_, r)); continue
        if r['status'] in ('VIOLATION', 'INVALID') or depth >= max_depth:
            return dict(status=r['status'] if r['status'] != 'FAIL' else 'FAIL', at=dict(s=sbox, y=yb_), detail=r, boxes=nb)
        sa, sb = sbox
        if depth % 3 == 2 and yb_[0] != yb_[1]:
            ymid = (yb_[0] + yb_[1]) / 2
            queue += [(sbox, (yb_[0], ymid), depth + 1), (sbox, (ymid, yb_[1]), depth + 1)]
        else:
            smid = (sa + sb) / 2
            queue += [((sa, smid), yb_, depth + 1), ((smid, sb), yb_, depth + 1)]
    return dict(status='PASS', boxes=nb, pieces=len(passed),
                worst=max(r['bound'] for _, _, r in passed) if sign < 0 else min(r['bound'] for _, _, r in passed))


# ------------------------------------------------------------------ certified quadrature (integral route)
def quad(model, kind, ybox, tolE, tsplit=4):
    """certified enclosure of Q := (1/2) int_0^1 int_{-1}^1 F dc dt for kind 'I' (F = t^2(1-t)^3(1-c^2)(1+(1-t)yU)^-7)
    or 'I4'; y a point (is_pt(ybox)) or a box (then the enclosure holds for every y in it)."""
    if not model.finite:
        return None
    tot = arb(0)
    edges = [k / tsplit for k in range(tsplit + 1)]
    for a_, b_ in zip(edges[:-1], edges[1:]):
        p = plan(model, kind, (a_, b_), ybox, tolE)
        if p is None:
            return None
        ns, nc, ny = p['ns'], p['nc'], p['ny']
        F = box_values(model, kind, (a_, b_), ybox, ns, nc, ny)
        C = coefs(F, ns, ny, nc)
        if not finite_all(C):
            return None
        wt = lambda k: arb(0) if k % 2 else arb(2) / (1 - k * k)
        th = (A(b_) - A(a_)) / 2
        # integrate in t (axis 0) and c (axis 2): leaves a polynomial in eta (axis 1)
        d = []
        for m in range(ny + 1):
            acc = arb(0)
            for k in range(0, ns + 1, 2):
                wk = wt(k)
                for kc in range(0, nc + 1, 2):
                    acc += C[k, m, kc] * wk * wt(kc)
            d.append(acc * th)
        val = d[0]
        if ny > 0:
            r = arb(0)
            for x in d[1:]:
                r += abs(x)
            val = val + arb(0, 1) * r.upper()
        val = val + arb(0, 1) * (p['E'] * 2 * th * 2).upper()
        tot += val
    return tot / 2


# ------------------------------------------------------------------ one slab
def slab_job(k, log, eta_override=None, ny=None):
    ya, yb = Y_EDGES[k], Y_EDGES[k + 1]
    ny = ny or (4 if k == 0 else 6)
    t0 = time.time()
    fm = fit_slab(ya, yb, ny)
    R0 = float_resid(fm)
    eta0 = eta_override or max(ETA_FLOOR, 8 * R0)
    log(f"  slab {k}: j in [{J_EDGES[k]}, {J_EDGES[k+1]}], y in [{ya:.6g}, {yb:.6g}], ny={ny}, float max|Gh(U_N)| = {R0:.3e}, eta0 = {eta0:.3e}  (fit {time.time()-t0:.1f}s)")
    res = None
    while eta0 <= ETA_CAP:
        mp = Model(fm, [eta0] + [0.0] * ny); mm = Model(fm, [-eta0] + [0.0] * ny)
        tolE = 0.05 * eta0
        sup = certify(mp, 'Gh', -1, (ya, yb), tolE, log=log)
        sub = certify(mm, 'Gh', +1, (ya, yb), tolE, log=log) if sup['status'] == 'PASS' else dict(status='SKIPPED')
        pos = certify(mm, 'V', +1, (ya, yb), 0.0, log=log) if sub['status'] == 'PASS' else dict(status='SKIPPED')
        res = dict(sup=sup, sub=sub, pos=pos)
        if all(r['status'] == 'PASS' for r in res.values()):
            break
        if any(r['status'] in ('VIOLATION', 'INVALID') for r in res.values()) and eta0 * 4 > ETA_CAP:
            break
        eta0 *= 4
        log(f"  slab {k}: retry with eta0 = {eta0:.3e}")
    ok = res is not None and all(r['status'] == 'PASS' for r in res.values())
    return dict(k=k, j=[J_EDGES[k], J_EDGES[k + 1]], y=[ya, yb], ny=ny, eta0=eta0, R0_float=R0, ok=ok, res=res,
                model={str(l): fm['a'][l].tolist() for l in fm['a']}, seconds=time.time() - t0)


def model_from_rec(rec):
    return dict(ya=rec['y'][0], yb=rec['y'][1], ny=rec['ny'], a={int(l): np.array(v) for l, v in rec['model'].items()})


def enclose(rec, y, log=None, tolq=1e-16):
    """certified enclosure of delta at the exact point y (arb) inside the slab of rec (only if the slab passed)."""
    fm = model_from_rec(rec); eta0 = rec['eta0']; ny = rec['ny']
    mp = Model(fm, [eta0] + [0.0] * ny); mm = Model(fm, [-eta0] + [0.0] * ny)
    yv = y if isinstance(y, arb) else arb(y)
    if not (arb(rec['y'][0]) <= yv and yv <= arb(rec['y'][1])):
        raise ValueError('point not certified inside this slab')
    p0 = mp.p0(yv) - eta0                      # p(0; y) of U_N (mp carries +eta0)
    slope = yv * (p0 + arb(0, 1) * eta0)       # delta = y U_s(0) in y [p0 - eta0, p0 + eta0]
    out = dict(slope=slope)
    Ip = quad(mp, 'I', (yv, yv), tolq); Im = quad(mm, 'I', (yv, yv), tolq)
    if Ip is not None and Im is not None:
        lo = (arb(9) / 4 * yv * Ip).lower(); hi = (arb(9) / 4 * yv * Im).upper()   # I decreasing in U
        out['integral'] = (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)
    d = slope
    if 'integral' in out:
        lo = max(slope.lower(), out['integral'].lower()); hi = min(slope.upper(), out['integral'].upper())
        lo = arb(lo); hi = arb(hi)
        if not (lo <= hi):
            out['EMPTY'] = True
        d = (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)
    out['delta'] = d
    return out


def find_slab(recs, y):
    for r in recs:
        if r.get('ok') and arb(r['y'][0]) <= y and y <= arb(r['y'][1]):
            return r
    return None


def report_points(recs):
    """the five J values and every slab endpoint; y must be an exact float for the integral route, so the J values
    0.0025 ... 0.25 are used at the float y nearest to (4J)^2 when (4J)^2 is not a float, and that is stated."""
    pts = []
    for Js in REPORT_J:
        Jx = arb(Js); yx = (4 * Jx) ** 2
        yf = float(yx.mid())
        pts.append(dict(label=f'J={Js}', J=Jx, y=yx, yf=yf))
    for k, yv in enumerate(Y_EDGES):
        yx = arb(yv)
        pts.append(dict(label=f'edge{k} j={J_EDGES[k]}', J=yx.sqrt() / 4 if yv > 0 else arb(0), y=yx, yf=yv))
    return pts


def eval_points(recs, log):
    out = []
    for p in report_points(recs):
        if p['yf'] == 0:
            continue
        cands = [r for r in recs if r.get('ok') and arb(r['y'][0]) <= p['y'] and p['y'] <= arb(r['y'][1])]
        best = None
        for r in cands:
            e = enclose(r, p['y'], log)           # exact point (arb ball of the decimal / dyadic y)
            if best is None:
                best = e
            else:
                lo = max(best['delta'].lower(), e['delta'].lower()); hi = min(best['delta'].upper(), e['delta'].upper())
                best['delta'] = (arb(lo) + arb(hi)) / 2 + arb(0, 1) * ((arb(hi) - arb(lo)) / 2)
        if best is None:
            continue
        d = best['delta']
        Jx = p['J']
        chi = Jx / (1 + d) ** 2
        rel = float(((d.upper() - d.lower()) / d.mid()))
        out.append(dict(label=p['label'], y=p['yf'], delta=outward(d), delta_str=dstr(d, 16), rel_width=rel,
                        slope=outward(best['slope']), integral=outward(best['integral']) if 'integral' in best else None,
                        slope_rel=float((best['slope'].upper() - best['slope'].lower()) / best['slope'].mid()),
                        integral_rel=float((best['integral'].upper() - best['integral'].lower()) / best['integral'].mid()) if 'integral' in best else None,
                        chi=outward(chi), chi_str=dstr(chi, 16), M_ADM=outward(1 + d), empty=best.get('EMPTY', False)))
        log(f"  {p['label']:<22} delta = {dstr(d, 16)}  rel.width {rel:.2e}  chi = {dstr(chi, 14)}")
    return out


def first_slab_limits(rec, log):
    fm = model_from_rec(rec); eta0 = rec['eta0']; ny = rec['ny']
    mp = Model(fm, [eta0] + [0.0] * ny); mm = Model(fm, [-eta0] + [0.0] * ny)
    pz = mp.p0(arb(0)) - eta0
    lim2 = pz + arb(0, 1) * eta0                       # delta/j^2 as j -> 0 (slope route; integral route gives 1/40 trivially)
    # uniform over the first slab: delta/j^2 in [min_y p0(y) - eta0, max_y p0(y) + eta0]
    a0 = fm['a'][0]
    cfa = [sum((arb(float(a0[n, g])) * (1 if n % 2 == 0 else -1) for n in range(a0.shape[0])), arb(0)) for g in range(a0.shape[1])]
    unif = cfa[0] + arb(0, 1) * (sum((abs(c) for c in cfa[1:]), arb(0)) + eta0).upper()
    I4p = quad(mp, 'I4', (0.0, 0.0), 1e-17); I4m = quad(mm, 'I4', (0.0, 0.0), 1e-17)
    d4 = None
    if I4p is not None and I4m is not None:
        lo = (arb(9) / 4 * I4p).lower(); hi = (arb(9) / 4 * I4m).upper()
        d4 = (lo + hi) / 2 + arb(0, 1) * ((hi - lo) / 2)
    out = dict(delta_over_j2_at_0=outward(lim2), delta_over_j2_at_0_str=dstr(lim2, 16), contains_1_40=bool(lim2.contains(arb(1) / 40)),
               delta_over_j2_uniform_slab=outward(unif), delta4_limit=outward(d4) if d4 is not None else None,
               delta4_limit_str=dstr(d4, 12) if d4 is not None else None,
               delta4_contains_m133_105600=bool(d4.contains(arb(-133) / 105600)) if d4 is not None else None)
    log(f"  first slab: delta/j^2 -> {dstr(lim2, 16)} (contains 1/40: {out['contains_1_40']}); uniform on slab: {dstr(unif, 10)}")
    if d4 is not None:
        log(f"  first slab: (delta - j^2/40)/j^4 -> {dstr(d4, 12)}  (contains -133/105600 = {-133/105600:.12e}: {out['delta4_contains_m133_105600']})")
    return out


def sha_self():
    return hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()


# ------------------------------------------------------------------ stage 0 (control N3)
def stage0_model():
    """v1 = A1 + B1 P2 in our basis, exact rationals as floats where exact (1/40 is not a float: we keep it exact
    through the polynomial identity in sympy and certify positivity of v1/s with a float-coefficient ball model)."""
    import sympy as sp
    s, c = sp.symbols('s c')
    p = ((1 - s) ** 4 + 5 * (1 - s) ** 3 * s + 10 * (1 - s) ** 2 * s ** 2 + 5 * (1 - s) * s ** 3 + s ** 4) / 40
    q2 = -(1 - s) * s / 20
    U0 = s * p; U2 = s ** 2 * (1 - s) * q2
    Lh0 = (1 - s) * sp.diff(U0, s, 2) - 2 * sp.diff(U0, s)
    Lh2 = (1 - s) * sp.diff(U2, s, 2) - 2 * sp.diff(U2, s) - 6 * q2
    P2 = (3 * c ** 2 - 1) / 2
    ident = sp.expand(Lh0 + Lh2 * P2 + sp.Rational(9, 4) * (1 - c ** 2) * s ** 2 * (1 - s) ** 2)
    # sanity: v1 matches GNPP A1, B1 with r = (1-s)/s
    r = sp.symbols('r', positive=True)
    A1 = (r ** 4 + 5 * r ** 3 + 10 * r ** 2 + 5 * r + 1) / (40 * (r + 1) ** 5); B1 = -r ** 2 / (20 * (r + 1) ** 5)
    chkA = sp.simplify(A1.subs(r, (1 - s) / s) - U0); chkB = sp.simplify(B1.subs(r, (1 - s) / s) - U2)
    # Chebyshev coefficients (exact rationals) of p and q2 in X = 2s - 1, then as arb-exact? use fmpq -> arb balls
    X = sp.symbols('X')
    def chebco(expr, deg):
        e = sp.expand(expr.subs(s, (X + 1) / 2))
        poly = sp.Poly(e, X)
        mono = [poly.coeff_monomial(X ** k) for k in range(deg + 1)]
        # monomial -> Chebyshev exact: X^k = 2^(1-k) sum binom...
        ch = [sp.Rational(0)] * (deg + 1)
        for k, mk in enumerate(mono):
            if mk == 0: continue
            for jx in range(k + 1):
                if (k - jx) % 2: continue
                coef = sp.binomial(k, (k - jx) // 2) / sp.Integer(2) ** (k - 1 if k > 0 else 0)
                if jx == 0 and k > 0: coef = coef / 2
                ch[jx] += mk * coef
        return ch
    return dict(ident_zero=(ident == 0), A1_match=(chkA == 0), B1_match=(chkB == 0), p_cheb=chebco(p, 4), q2_cheb=chebco(q2, 2))


def stage0_checks(log):
    import sympy as sp
    st = stage0_model()
    log(f"  N3 exact identity Lh[v1] = -(9/4)(1-c^2)s^2(1-s)^2 (sympy, exact): {st['ident_zero']}; v1 == GNPP A1,B1: {st['A1_match'] and st['B1_match']}")
    # v1/s = p + s(1-s) q2 P2 > 0 on the closed square: certify with the 'V' machinery. Coefficients are rationals; we
    # enclose them by two float models (rounded down / up) and certify both lower-bounds: V is linear in the coefficients
    # and we use the coefficientwise ball instead: build a Model whose coefficients are float approximations, then add the
    # rigorous rounding slack |rational - float| * (max |basis| = 1 per term) to the bound.
    a0 = np.zeros((NB + 1, 1)); a2 = np.zeros((NB + 1, 1))
    slack = 0.0
    for n, v in enumerate(st['p_cheb']):
        a0[n, 0] = float(v); slack += abs(float(sp.Abs(v - sp.Rational(a0[n, 0])))) * 1.0000001
    for n, v in enumerate(st['q2_cheb']):
        a2[n, 0] = float(v); slack += abs(float(sp.Abs(v - sp.Rational(a2[n, 0])))) * 1.0000001
    fm = dict(ya=0.0, yb=1.0, ny=0, a={0: a0, 2: a2})
    m = Model(fm)
    m.ny = 0
    r = certify(m, 'V', +1, (0.0, 0.0), 0.0, log=log)
    # V = p + s(1-s) q2 P2: |V_float - V_exact| <= slack (|T_n| <= 1, |s(1-s) P2| <= 1/4)
    okpos = r['status'] == 'PASS' and r['worst'] > slack
    log(f"  N3 v1/s > 0 on [0,1]x[-1,1]: {r['status']} (min bound {r.get('worst')}, rounding slack {slack:.2e}) -> {okpos}")
    # Stage-0 delta bracket at the report points: [I(v1) y (9/4)..., y/40]
    br = []
    for Js in REPORT_J:
        yx = (4 * arb(Js)) ** 2; yf = float(yx.mid())
        Iv = quad(m, 'I', (yx, yx), 1e-15)
        # rounding slack of v1 coefficients moves U by <= slack*s... integrand monotone: widen by derivative bound
        # |dI/dU| <= 7 y * (1/2) int t^2 (1-t)^4 (1-c^2) dt dc <= 7 y / 105 * (4/3) / 2; tiny, included:
        widen = 7 * yf * slack
        lo = (arb(9) / 4 * yx * Iv).lower() - widen
        br.append(dict(J=Js, lo=float(lo), hi=float((yx / 40).upper())))
    return dict(identity=bool(st['ident_zero']), gnpp_match=bool(st['A1_match'] and st['B1_match']), v1_pos=okpos,
                v1_pos_detail=r, brackets=br, ok=bool(st['ident_zero'] and okpos))


# ------------------------------------------------------------------ modes
def run_probe(log):
    t = time.time()
    rec = slab_job(18, log)
    log(f"probe slab 18: ok={rec['ok']} seconds={rec['seconds']:.1f}")
    t2 = time.time()
    pts = enclose(rec, arb(1.0), log)
    log(f"probe enclosure at y=1: delta = {dstr(pts['delta'], 16)} slope {dstr(pts['slope'], 12)} integral {dstr(pts.get('integral', arb(0)), 16)}  ({time.time()-t2:.1f}s)")
    rec0 = slab_job(0, log)
    log(f"probe slab 0: ok={rec0['ok']} seconds={rec0['seconds']:.1f}")
    first_slab_limits(rec0, log)
    log(f"probe total {time.time()-t:.1f}s")


def run_controls(out, log):
    T0 = time.time()
    res = dict(sha256_p1=sha_self(), started=time.strftime('%Y-%m-%d %H:%M:%S'))
    # NaN / inf guards
    fm = fit_slab(Y_EDGES[0], Y_EDGES[1], 4)
    g = {}
    for bad in (float('nan'), float('inf'), 1e300):
        fb = dict(ya=fm['ya'], yb=fm['yb'], ny=fm['ny'], a={l: v.copy() for l, v in fm['a'].items()})
        fb['a'][2][3, 1] = bad
        r = check_box(Model(fb), 'Gh', -1, (0.0, 1.0), (fm['ya'], fm['yb']), 1e-12)
        g[str(bad)] = r['status']
    res['nan_inf_guard'] = dict(statuses=g, ok=all(v == 'INVALID' for v in g.values()))
    log(f"guard: {g}")
    # slabs needed for the 5 J values: 0 (j=0.01, 0.1), 4 or 5 (j=0.25), 9 (j=0.5), 18 (j=1)
    recs = []
    for k in (0, 4, 9, 18):
        rec = slab_job(k, log)
        recs.append(rec)
        log(f"slab {k}: ok={rec['ok']} eta0={rec['eta0']:.3e} ({rec['seconds']:.1f}s)")
    res['slabs'] = [{kk: v for kk, v in r.items() if kk != 'model'} for r in recs]
    # N1: sub-solution with slope coefficient raised by 1e-3 relative; N2: u+ - 1e-6 j^2 s
    n1 = {}
    for rec in (recs[0], recs[-1]):
        fmk = model_from_rec(rec); eta0 = rec['eta0']; ny = rec['ny']
        cf = Model(fmk).p0_coefs_eta()
        shift = [-eta0 + 1e-3 * (cf[0] - eta0)] + [1e-3 * c for c in cf[1:]]
        r = certify(Model(fmk, shift), 'Gh', +1, (rec['y'][0], rec['y'][1]), 0.05 * eta0, log=log)
        n1[rec['k']] = r['status']
    res['N1'] = dict(statuses=n1, ok=all(v == 'VIOLATION' for v in n1.values()))
    rec = recs[0]; fmk = model_from_rec(rec); eta0 = rec['eta0']
    r = certify(Model(fmk, [eta0 - 1e-6] + [0.0] * rec['ny']), 'Gh', -1, (rec['y'][0], rec['y'][1]), 0.05 * eta0, log=log)
    res['N2'] = dict(status=r['status'], ok=r['status'] == 'VIOLATION')
    log(f"N1: {n1}  N2: {r['status']}")
    # N3 stage 0
    res['N3'] = stage0_checks(log)
    # known answer: our own double solve at the 5 J values (full double, Ns=48/NL=8 collocation) and probe0/solve_out.json
    S = setup()
    ref = {}
    try:
        so = json.load(open(os.path.join(HERE, 'probe0', 'solve_out.json')))['64_10']
        ref_probe0 = {round(rw['j'], 10): rw['delta'] for rw in so}
    except Exception:
        ref_probe0 = {}
    pts = eval_points(recs, log)
    ka = []
    for p in pts:
        if not p['label'].startswith('J='):
            continue
        Js = p['label'][2:]; j = 4 * float(Js); y = p['y']
        Q = solveU(y, S); b = np.linalg.solve(S['Tm'], Q.T).T
        d_slope = y * npch.chebval(-1.0, b[0])
        # double integral route on the same double solution (its sensitivity to errors in U carries a factor ~ y)
        tq, wt_ = npleg.leggauss(160); tq = (tq + 1) / 2; wt_ = wt_ / 2
        cq, wc = npleg.leggauss(48)
        Ut = np.zeros((len(tq), len(cq)))
        for i, l in enumerate(S['Lv']):
            w, _, _ = wfun(l, tq)
            Ut += np.outer(w * npch.chebval(2 * tq - 1, b[i]), npleg.legval(cq, np.eye(l + 1)[l]))
        # split form y/40 + correction: the plain Gauss sum of the full integrand carries ~1e-14 relative roundoff
        corr = (tq ** 2 * (1 - tq) ** 3)[:, None] * (1 - cq ** 2)[None, :] * ((1 + (1 - tq)[:, None] * y * Ut) ** -7 - 1)
        d_int = y / 40 + 2.25 * y * float(wt_ @ corr @ (wc / 2))
        d_p0 = ref_probe0.get(round(j, 10))
        lo, hi = p['delta']
        tol = 2e-15 * abs(d_int)                       # double floor of the 2-D quadrature sum (stated in BUILD-P1.md)
        st0 = next(b_ for b_ in res['N3']['brackets'] if b_['J'] == Js)
        ka.append(dict(J=Js, enclosure=p['delta'], rel_width=p['rel_width'],
                       double_integral_route=d_int, inside_raw=lo <= d_int <= hi, inside_within_double_floor=lo - tol <= d_int <= hi + tol,
                       double_slope_route=d_slope, slope_route_enclosure=p['slope'], inside_slope_enclosure=p['slope'][0] <= d_slope <= p['slope'][1],
                       double_probe0_slope=d_p0, probe0_inside_slope_enclosure=(d_p0 is not None and p['slope'][0] <= d_p0 <= p['slope'][1]),
                       probe0_dist_rel=(None if d_p0 is None else (0.0 if lo <= d_p0 <= hi else min(abs(d_p0 - lo), abs(d_p0 - hi)) / d_p0)),
                       inside_stage0=st0['lo'] <= lo and hi <= st0['hi'], stage0=[st0['lo'], st0['hi']], chi=p['chi']))
        kk = ka[-1]
        kk['ok'] = bool(kk['inside_within_double_floor'] and kk['inside_slope_enclosure'] and kk['probe0_inside_slope_enclosure'] and kk['inside_stage0'])
        log(f"known-answer J={Js}: enclosure [{lo:.16e}, {hi:.16e}] rel {p['rel_width']:.2e}; double integral route {d_int:.16e} "
            f"inside raw={kk['inside_raw']} within 2e-15 floor={kk['inside_within_double_floor']}; double slope routes (here {d_slope:.16e}, probe0 {d_p0}) "
            f"inside slope enclosure={kk['inside_slope_enclosure']},{kk['probe0_inside_slope_enclosure']}; probe0 rel dist to final enclosure {kk['probe0_dist_rel']:.1e}; inside stage0 {kk['inside_stage0']}")
    res['known_answer'] = ka
    res['points'] = pts
    res['first_slab'] = first_slab_limits(recs[0], log)
    res['seconds'] = time.time() - T0
    res['all_ok'] = bool(res['nan_inf_guard']['ok'] and res['N1']['ok'] and res['N2']['ok'] and res['N3']['ok']
                         and all(r['ok'] for r in recs) and all(x['ok'] for x in ka) and len(ka) == 5 and res['first_slab']['contains_1_40'] and res['first_slab']['delta4_contains_m133_105600'])
    json.dump(res, open(out, 'w'), indent=1, default=str)
    log(f"controls done in {res['seconds']:.1f}s; all_ok={res['all_ok']}; sha256(p1.py)={res['sha256_p1']}")


def run_main(out, log):
    if not os.path.exists(os.path.join(HERE, 'GO-P1')):
        log("REFUSED: <LAB>/BHSLICE/GO-P1 does not exist (the sealed bar must be armed first).")
        sys.exit(2)
    T0 = time.time()
    state = json.load(open(out)) if os.path.exists(out) else dict(sha256_p1=sha_self(), slabs={}, invocations=0)
    if state.get('sha256_p1') != sha_self():
        log(f"WARNING: p1.py changed since the checkpoint was started (was {state.get('sha256_p1')}); refusing to mix.")
        sys.exit(3)
    state['invocations'] = state.get('invocations', 0) + 1
    worst = 30.0
    for k in range(len(Y_EDGES) - 1):
        if str(k) in state['slabs']:
            continue
        if time.time() - T0 + 1.5 * worst > TIME_BUDGET - 20:
            log(f"time budget reached; {len(state['slabs'])}/19 slabs done; rerun to resume")
            json.dump(state, open(out, 'w'), indent=1, default=str)
            return
        rec = slab_job(k, log)
        worst = max(worst, rec['seconds'])
        state['slabs'][str(k)] = rec
        json.dump(state, open(out, 'w'), indent=1, default=str)
    recs = [state['slabs'][str(k)] for k in range(len(Y_EDGES) - 1)]
    if 'points' not in state:
        if time.time() - T0 > TIME_BUDGET - 200:
            json.dump(state, open(out, 'w'), indent=1, default=str); log("slabs done; rerun for enclosures"); return
        state['points'] = eval_points(recs, log)
        state['first_slab'] = first_slab_limits(recs[0], log)
    state['all_slabs_ok'] = all(r['ok'] for r in recs)
    state['finished'] = True
    json.dump(state, open(out, 'w'), indent=1, default=str)
    log(f"main finished: all slabs ok = {state['all_slabs_ok']}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--probe', action='store_true')
    ap.add_argument('--controls')
    ap.add_argument('--main')
    args = ap.parse_args()
    logf = None
    def log(msg):
        print(msg, flush=True)
        if logf: logf.write(msg + '\n'); logf.flush()
    if args.controls:
        logf = open(os.path.splitext(args.controls)[0] + '.log', 'w')
        run_controls(args.controls, log)
    elif args.main:
        logf = open(os.path.splitext(args.main)[0] + '.log', 'a')
        run_main(args.main, log)
    elif args.probe:
        run_probe(log)
