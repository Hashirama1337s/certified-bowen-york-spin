"""BHSLICE B3T instrument: certified spin of the EXTREME (mu = 0) Bowen-York data, computed directly on the throat.

BAR-B3T.md, seal 10a65e940a. Design: scout/B3-EXTREME.md (route B). Certificate machinery: the REP6 polyellipse /
Cauchy patch certificate (rep6/rep6.py, logic re-implemented here; rep6/ is not imported and not edited).

Object (J = 1, scale-free; chi(0) = J / M_0^2):
  extreme puncture  Delta u = -S (1 + u)^-p,  S = coef * ang(c) * r^-(p+5)/2,
  Bowen-York: coef = 9/4, p = 7, ang = 1 - c^2 (S = (9/4) J^2 sin^2 theta r^-6);
  control RN (extreme Reissner-Nordstrom, q = 1): coef = 1/4, p = 3, ang = 1 (exact u = sqrt(1 + 1/r) - 1, M = 1).
Throat coordinates: sigma = log r, x = xi = tanh(sigma/8), so r^(1/4) = (1 + x)/(1 - x), s := sqrt r = ((1+x)/(1-x))^2,
  c = cos theta.  (NB scout/B3-EXTREME.md sec. 3 prints "s = sqrt r = ((1+xi)/(1-xi))^4"; that expression is r, not
  sqrt r; the probe code b3probe.py uses r = ((1+x)/(1-x))^4 correctly. Logged in BUILD-B3T.md.)
Ansatz (B3-EXTREME sec. 3, rational prefactor):  w = r^-1/2 G / (1 + sqrt r),  G = g0 + (1 - x)^2 H(x, c),
  H = sum_{n<=N, k<=K} h_nk T_n(x) P_2k(c).
Normalised residual (derived in BUILD-B3T.md sec. 2; checked against brute-force mpmath Laplacians by control C0):
  Delta w + S (1 + w)^-p = r^-5/2 (1 + s)^-1 E_G,   E_G = Q (1 + 3 s) / (4 (1 + s)^2)   (positive factors), and
  Q = [A2 G_xx - A1 P1 + 4 (1+x^2)^2 L_S H + coef ang 2^(p+3) (1+x^2)^(p+3) (1-x)^(2p-4) / D^p] / (1 + x + x^2) - G,
  A2 = (1+x^2)^2 (1+x)^2 / 16,  A1 = (1+x^2)(1+x)(x(1+x^2) + 2(1+x)^2) / 8,  P1 = (1-x) H_x - 2 H,
  G_xx = (1-x)^2 H_xx - 4 (1-x) H_x + 2 H,  D = 2 (1+x)^2 (1+x^2) + G (1-x)^4,  L_S = (1-c^2) d_cc - 2 c d_c.
Shift: G -> G + eta (H fixed) changes Q by exactly -eta in the linear part, and the source only decreases (D grows,
  (1-x)^4 >= 0). So w_+ := r^-1/2 (G_N + eta)/(1 + sqrt r) is a super-solution (Q <= 0 everywhere) if eta >= sup Q_N
  on the closed square [-1,1] x [0,1] (G even in c).
Class B3-5: C^2 off 0 (structural: x smooth in r on (0, oo), G polynomial in c^2 = z^2/r^2), w -> 0 at infinity
  (structural), 1 + w > 0 and w -> +oo at 0 from the certified inf (G_N + eta) > 0.
Mass bounds (B3-5 (ii); J = 1):
  M_0 <= M_up := 2 lim r w_+ = 2 (g0 + eta)  (exact: G_N(1, c) = g0);
  M_0 >= M_lo := (1/(2 pi)) int S (1 + w_+)^-p dx = 2 int_0^1 dc int_{-1}^{1} dx f(x, c),
  f = coef ang 2^(p+3) (1+x)(1+x^2)^p (1-x)^(2p-3) / D_+^p,  D_+ = D with G -> G_N + eta.
Rigorous quadrature (BUILD-B3T.md sec. 4): on each patch the tensor Chebyshev-Lobatto interpolant p of f is
  integrated exactly (int_{-1}^1 T_j = 2/(1 - j^2), j even), and |int (f - p)| <= area * ||f - p||_inf with the
  polyellipse/Cauchy bound ||f - p|| <= 8 M (s1 s2 - t1 t2); the patch contribution is bounded below by
  hw_x hw_c (int p - 4 err).
chi(0) in [1/M_up^2, 1/M_lo^2], exact Fractions, decimal strings rounded outward.

Usage (every run: py -3.11 <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 b3t.py ..., from <LAB>/BHSLICE):
  b3t.py --controls controls_b3t.json      (resumable; exit 0 = all parts done, 3 = rerun the same command)
  b3t.py --main result_b3t.json            (refuses without GO-B3T and a ledger GO row for seal 10a65e940a naming
                                            this file's sha256; resumable; exit 0 finished, 3 = rerun, other = error)
Exit codes: 0 finished, 3 rerun to resume, 1 controls failed, 2 no GO file, 4 checkpoint refused, 6 self-sha not
pinned, 7 certificate failed.
"""
import sys
sys.dont_write_bytecode = True
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # single-thread LAPACK => deterministic
    os.environ[_v] = '1'
import json, time, math, hashlib, argparse, shutil
from fractions import Fraction
import numpy as np
from numpy.polynomial import chebyshev as Ch
from numpy.polynomial import legendre as Lg
from flint import arb, acb, arb_poly, acb_poly, arb_mat, acb_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-B3T')
SEAL_ID = '10a65e940a'
INSTRUMENT_NAME = 'b3t.py'
SCRATCH = os.path.join(HERE, 'seats', 'build-b3t')
TIME_BUDGET = 520.0
EXIT_FINISHED, EXIT_RESUME, EXIT_CONTROLS_FAILED, EXIT_NO_GO, EXIT_CKPT_REFUSED, EXIT_NOT_PINNED, EXIT_CERT_FAILED = \
    0, 3, 1, 2, 4, 6, 7

PREC = 256
ctx.prec = PREC
NSAMP = 64                     # Lobatto degree per patch (n1 = n2)
M1, M2 = 24, 48                # arcs: upper half x-ellipse (conjugate symmetry), full c-ellipse
RHO_COMBOS = [(2.0, 2.0), (1.6, 1.6), (3.0, 3.0), (2.0, 1.4), (1.4, 2.0), (1.3, 1.3), (1.15, 1.15)]
MIN_PATCH = Fraction(1, 4096)
QTOL = Fraction(1, 10 ** 14)   # per-patch bound on ||f - p||_inf accepted for the quadrature
NX_INIT = 16                   # initial patches in x over [-1, 1] (c in [0, 1] whole)

# main configuration (fixed before any certified number; chosen from float probes, BUILD-B3T.md sec. 5)
MAIN_CFG = dict(model='BY', fac='1', N=112, K=12, eta_factor='4')
P2E_WINDOW = (Fraction('0.9282004538'), Fraction('0.9282004687'))   # BAR-B3T ceiling (P2e certified window, rounded)
T4_PROBE = Fraction('0.9282004587652')

MODELS = {'BY': dict(coef=Fraction(9, 4), p=7, ang='1-c2'),
          'RN': dict(coef=Fraction(1, 4), p=3, ang='1')}


# ================================================================== small helpers
def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def sha_self():
    return sha_file(os.path.abspath(__file__))


def fstr(q):
    q = Fraction(q)
    return f'{q.numerator}/{q.denominator}'


def fhex(x):
    return float(x).hex()


def unhex(s):
    return float.fromhex(s)


def atomic_json(path, obj):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=True)
    for k in range(50):                    # Windows: a reader (indexer, viewer) can hold the target briefly
        try:
            os.replace(tmp, path)
            return
        except PermissionError:
            time.sleep(0.2)
    os.replace(tmp, path)


def arb_q(q):
    q = Fraction(q)
    return arb(q.numerator) / arb(q.denominator)


def arb_point_frac(p):
    m, e = p.mid().man_exp()
    return Fraction(int(m)) * (Fraction(2) ** int(e))


def ball_upper(a):
    return arb_point_frac(a.upper())


def ball_lower(a):
    return arb_point_frac(a.lower())


def dyadic_up(q, bits=64):
    q = Fraction(q)
    n = -((-q.numerator * 2 ** bits) // q.denominator)
    return Fraction(n, 2 ** bits)


def dec_round(q, digits, up):
    """decimal string with `digits` significant digits, rounded up (up=True) or down; asserted outward."""
    q = Fraction(q)
    if q == 0:
        return '0'
    neg = q < 0
    a = -q if neg else q
    e = int(math.floor(math.log10(a.numerator) - math.log10(a.denominator)))
    while Fraction(10) ** e > a:
        e -= 1
    while Fraction(10) ** (e + 1) <= a:
        e += 1
    scale = Fraction(10) ** (digits - 1 - e)
    v = a * scale
    mag_up = (up and not neg) or ((not up) and neg)
    n = v.numerator // v.denominator
    if mag_up and n * v.denominator != v.numerator:
        n += 1
    exp10 = digits - 1 - e
    if exp10 > 0:
        txt = str(n).rjust(exp10 + 1, '0')
        txt = txt[:-exp10] + '.' + txt[-exp10:]
    else:
        txt = str(n) + '0' * (-exp10)
    out = ('-' if neg else '') + txt
    assert (Fraction(out) >= q) if up else (Fraction(out) <= q)
    return out


def peak_memory():
    """peak working set and peak commit of THIS process, GB (Windows psapi); None elsewhere."""
    try:
        import ctypes
        import ctypes.wintypes as wt

        class PMC(ctypes.Structure):
            _fields_ = [('cb', wt.DWORD), ('PageFaultCount', wt.DWORD), ('PeakWorkingSetSize', ctypes.c_size_t),
                        ('WorkingSetSize', ctypes.c_size_t), ('QuotaPeakPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaPagedPoolUsage', ctypes.c_size_t), ('QuotaPeakNonPagedPoolUsage', ctypes.c_size_t),
                        ('QuotaNonPagedPoolUsage', ctypes.c_size_t), ('PagefileUsage', ctypes.c_size_t),
                        ('PeakPagefileUsage', ctypes.c_size_t)]
        k = ctypes.WinDLL('kernel32')
        k.GetCurrentProcess.restype = wt.HANDLE
        k.K32GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(PMC), wt.DWORD]
        k.K32GetProcessMemoryInfo.restype = wt.BOOL
        c = PMC()
        c.cb = ctypes.sizeof(c)
        if not k.K32GetProcessMemoryInfo(k.GetCurrentProcess(), ctypes.byref(c), c.cb):
            return None
        return dict(peak_working_set_gb=round(c.PeakWorkingSetSize / 1024 ** 3, 3),
                    peak_commit_gb=round(c.PeakPagefileUsage / 1024 ** 3, 3))
    except Exception:
        return None


class Refusal(Exception):
    pass


# ================================================================== the model (generic scalar code: numpy / arb / acb)
def model_params(model, fac):
    m = MODELS[model]
    return dict(model=model, coef=m['coef'] * Fraction(fac), p=m['p'], ang=m['ang'], fac=str(Fraction(fac)))


def ang_val(mp, c):
    return (1 - c * c) if mp['ang'] == '1-c2' else (c * 0 + 1)


def Q_point(x, c, g0, H, Hx, Hxx, LH, coefc, mp, eta=0):
    """normalised residual Q of G = g0 + eta + (1-x)^2 H (eta enters only through -G and D); returns (Q, G, D, tri).
    coefc = coef * 2^(p+3) as a scalar of the working type."""
    p = mp['p']
    w1 = 1 - x
    x2p = 1 + x * x
    xp1 = 1 + x
    tri = 1 + x + x * x
    w2 = w1 * w1
    G = g0 + eta + w2 * H
    Gxx = w2 * Hxx - 4 * w1 * Hx + 2 * H
    P1 = w1 * Hx - 2 * H
    A2 = x2p * x2p * xp1 * xp1 / 16
    A1 = x2p * xp1 * (x * x2p + 2 * xp1 * xp1) / 8
    D = 2 * xp1 * xp1 * x2p + G * w2 * w2
    src = coefc * ang_val(mp, c) * x2p ** (p + 3) * w1 ** (2 * p - 4) / D ** p
    Q = (A2 * Gxx - A1 * P1 + 4 * x2p * x2p * LH + src) / tri - G
    return Q, G, D, tri


def f_point(x, c, Gp, coefc, mp):
    """mass integrand f = coef ang 2^(p+3) (1+x)(1+x^2)^p (1-x)^(2p-3) / D_+^p for G_+ = Gp; returns (f, D_+)."""
    p = mp['p']
    w1 = 1 - x
    x2p = 1 + x * x
    xp1 = 1 + x
    D = 2 * xp1 * xp1 * x2p + Gp * (w1 * w1) * (w1 * w1)
    return coefc * ang_val(mp, c) * xp1 * x2p ** p * w1 ** (2 * p - 3) / D ** p, D


# ================================================================== float solve (weighted = Q-normalised collocation)
def tbasis(x, N):
    V = np.zeros((len(x), N + 1)); D1 = np.zeros_like(V); D2 = np.zeros_like(V)
    for n in range(N + 1):
        cf = np.zeros(N + 1); cf[n] = 1
        V[:, n] = Ch.chebval(x, cf)
        D1[:, n] = Ch.chebval(x, Ch.chebder(cf, 1))
        D2[:, n] = Ch.chebval(x, Ch.chebder(cf, 2))
    return V, D1, D2


def pbasis(c, K):
    V = np.zeros((len(c), K + 1)); LV = np.zeros_like(V)
    for k in range(K + 1):
        cf = np.zeros(2 * k + 1); cf[2 * k] = 1
        V[:, k] = Lg.legval(c, cf)
        LV[:, k] = -2 * k * (2 * k + 1) * V[:, k]
    return V, LV


class FloatModel:
    def __init__(self, mp, N, K):
        self.mp, self.N, self.K = mp, N, K
        self.coefc = float(mp['coef'] * 2 ** (mp['p'] + 3))

    def mats(self, x, c):
        N, K = self.N, self.K
        Tv, T1, T2 = tbasis(x, N)
        Pv, PL = pbasis(c, K)
        X = np.kron(x, np.ones(len(c))); C = np.kron(np.ones(len(x)), c)
        return dict(Hv=np.kron(Tv, Pv), Hx=np.kron(T1, Pv), Hxx=np.kron(T2, Pv), LH=np.kron(Tv, PL), x=X, c=C)

    def resid(self, M, z, eta=0.0):
        g0, h = z[0], z[1:]
        x, c = M['x'], M['c']
        H = M['Hv'] @ h; Hx = M['Hx'] @ h; Hxx = M['Hxx'] @ h; LH = M['LH'] @ h
        Q, G, D, tri = Q_point(x, c, g0, H, Hx, Hxx, LH, self.coefc, self.mp, eta)
        return Q, G, D, tri

    def jac(self, M, z):
        Q, G, D, tri = self.resid(M, z)
        x, c = M['x'], M['c']
        p = self.mp['p']
        w1 = 1 - x; w2 = w1 * w1; x2p = 1 + x * x; xp1 = 1 + x
        A2 = x2p * x2p * xp1 * xp1 / 16
        A1 = x2p * xp1 * (x * x2p + 2 * xp1 * xp1) / 8
        src = self.coefc * ang_val(self.mp, c) * x2p ** (p + 3) * w1 ** (2 * p - 4) / D ** p
        dsrc = -p * src * w2 * w2 / D                         # d src / d G
        Gxx_b = w2[:, None] * M['Hxx'] - (4 * w1)[:, None] * M['Hx'] + 2 * M['Hv']
        P1_b = w1[:, None] * M['Hx'] - 2 * M['Hv']
        lin = (A2[:, None] * Gxx_b - A1[:, None] * P1_b + (4 * x2p * x2p)[:, None] * M['LH']) / tri[:, None]
        Jh = lin + ((dsrc / tri - 1) * w2)[:, None] * M['Hv']
        Jg = dsrc / tri - 1
        return Q, G, np.hstack([Jg[:, None], Jh])


def colloc_grid(N, K):
    nx = N + 12
    xg = np.cos(np.pi * (np.arange(nx) + 0.5) / nx)
    xg = np.concatenate([[-1.0], np.sort(xg)])
    cg, _ = Lg.leggauss(2 * (K + 4))
    cg = np.sort(cg[cg > 0])
    return xg, cg


def float_solve(mp, N, K, z0=None, log=print):
    fm = FloatModel(mp, N, K)
    xg, cg = colloc_grid(N, K)
    M = fm.mats(xg, cg)
    nunk = 1 + (N + 1) * (K + 1)
    if z0 is None:
        z = np.zeros(nunk)
        if mp['model'] == 'RN':
            z[0] = 0.5; z[1] = 0.5 / 4          # G(-1) ~ 1 = g0 + 4 H(-1)
        else:
            z[0] = 0.52; z[1] = 0.73 / 4        # G(-1) ~ 1.25 (B3-EXTREME probe start)
        z[0] *= float(Fraction(mp['fac'])) ** 0.25
    else:
        z = z0.copy()
    nr = None
    for it in range(80):
        Q, G, J = fm.jac(M, z)
        nr = float(np.max(np.abs(Q)))
        dz, *_ = np.linalg.lstsq(J, -Q, rcond=None)
        lam = 1.0
        while lam > 1e-4:
            zn = z + lam * dz
            Qn, Gn, _, _ = fm.resid(M, zn)
            if np.all(np.isfinite(Qn)) and np.min(Gn) > 0 and np.max(np.abs(Qn)) < max(nr, 1e-300) * (1 - 0.1 * lam) + 1e-14:
                break
            lam /= 2
        z = zn
        if np.max(np.abs(dz)) < 1e-15 or nr < 1e-15:
            break
    return z, dict(iterations=it + 1, colloc_max_abs_Q=nr)


def float_check(mp, N, K, z, eta=0.0, nx=2001, nc=41):
    """sampled (NOT rigorous) sup/inf of Q_N and min G_N on a fine grid of the closed square."""
    fm = FloatModel(mp, N, K)
    xall = -np.cos(np.pi * np.linspace(0, 1, nx))
    c = np.linspace(0, 1, nc)
    sup, inf, gmin, arg = -np.inf, np.inf, np.inf, None
    for j0 in range(0, nx, 250):                       # chunks in x: bounded memory
        M = fm.mats(xall[j0:j0 + 250], c)
        Q, G, D, _ = fm.resid(M, z, eta)
        i = int(np.argmax(Q))
        if Q[i] > sup:
            sup, arg = float(Q[i]), [float(M['x'][i]), float(M['c'][i])]
        inf = min(inf, float(np.min(Q))); gmin = min(gmin, float(np.min(G)))
    return dict(supQ=sup, infQ=inf, Gmin=gmin, argmax=arg)


def float_mass_lo(mp, N, K, z, eta, nq=800, ncq=40):
    """double-precision Gauss-Legendre value of M_lo (probe only)."""
    fm = FloatModel(mp, N, K)
    xq, wq = Lg.leggauss(nq)
    cq, wc = Lg.leggauss(2 * ncq); sel = cq > 0; cq, wc = cq[sel], wc[sel]
    M = fm.mats(xq, cq)
    _, G, _, _ = fm.resid(M, z)
    f, _ = f_point(M['x'], M['c'], G + eta, fm.coefc, mp)
    return float(2 * np.sum(np.kron(wq, wc) * f))


# ================================================================== arb patch machinery (REP6 logic, re-implemented)
def cheb_local(a, b, nmax):
    """T_n(m + hw X) as arb_poly in the local variable X, n = 0..nmax."""
    m = (arb_q(a) + arb_q(b)) / 2; hw = (arb_q(b) - arb_q(a)) / 2
    X = arb_poly([m, hw])
    T = [arb_poly([1]), X]
    for n in range(1, nmax):
        T.append(2 * X * T[n] - T[n - 1])
    return T[:nmax + 1], hw, m


def leg_local(a, b, kmax2):
    m = (arb_q(a) + arb_q(b)) / 2; hw = (arb_q(b) - arb_q(a)) / 2
    X = arb_poly([m, hw])
    P = [arb_poly([1]), X]
    for n in range(1, kmax2):
        P.append(((2 * n + 1) * X * P[n] - n * P[n - 1]) * (1 / arb(n + 1)))
    return P[:kmax2 + 1], hw, m


class Cand:
    """certified candidate G = g0 + (1-x)^2 sum h[n,k] T_n(x) P_2k(c); h exact dyadic floats."""

    def __init__(self, d):
        self.d = d
        self.mp = model_params(d['model'], Fraction(d['fac']))
        self.N, self.K = d['N'], d['K']
        self.g0 = unhex(d['g0'])
        self.h = np.array([unhex(s) for s in d['h']]).reshape(self.N + 1, self.K + 1)
        self.eta = Fraction(d['eta'])
        self.coefc = arb_q(self.mp['coef'] * 2 ** (self.mp['p'] + 3))

    def local(self, xa, xb, ca, cb):
        T, hwx, mx = cheb_local(xa, xb, max(self.N, 1))
        P, hwc, mc = leg_local(ca, cb, max(2 * self.K, 1))
        hk = []
        for k in range(self.K + 1):
            pk = arb_poly([0])
            for n in range(self.N + 1):
                cf = self.h[n, k]
                if cf != 0.0:
                    pk += arb(float(cf)) * T[n]
            hk.append(pk)
        hk1 = [q.derivative() * (1 / hwx) for q in hk]
        hk2 = [q.derivative() * (1 / hwx) for q in hk1]
        Pk = [P[2 * k] for k in range(self.K + 1)]
        lam = [-(2 * k) * (2 * k + 1) for k in range(self.K + 1)]
        return dict(hk=hk, hk1=hk1, hk2=hk2, Pk=Pk, lam=lam, hwx=hwx, hwc=hwc, mx=mx, mc=mc)


def fields(cd, loc, xs, cs, cplx):
    MAT = acb_mat if cplx else arb_mat
    POLY = (lambda q: acb_poly(q)) if cplx else (lambda q: q)
    hk = [POLY(q) for q in loc['hk']]; hk1 = [POLY(q) for q in loc['hk1']]; hk2 = [POLY(q) for q in loc['hk2']]
    Pk = [POLY(q) for q in loc['Pk']]
    K1 = cd.K + 1
    Hx0 = MAT([[hk[k](x) for k in range(K1)] for x in xs])
    Hx1 = MAT([[hk1[k](x) for k in range(K1)] for x in xs])
    Hx2 = MAT([[hk2[k](x) for k in range(K1)] for x in xs])
    PyT = MAT([[Pk[k](y) for y in cs] for k in range(K1)])
    PLT = MAT([[loc['lam'][k] * Pk[k](y) for y in cs] for k in range(K1)])
    xi = [loc['mx'] + loc['hwx'] * x for x in xs]
    cc = [loc['mc'] + loc['hwc'] * y for y in cs]
    return Hx0 * PyT, Hx1 * PyT, Hx2 * PyT, Hx0 * PLT, xi, cc


def patch_grid(cd, xa, xb, ca, cb, xs, cs, cplx):
    """Q_N, G_N, f (with G_+ = G_N + eta) on the grid; plus D_N, D_+, tri for the analyticity checks."""
    loc = cd.local(xa, xb, ca, cb)
    H, H1, H2, LH, xi, cc = fields(cd, loc, xs, cs, cplx)
    g0 = arb(cd.g0)
    eta = arb_q(cd.eta)
    out = dict(Q=[], G=[], f=[], D=[], Dp=[], tri=[])
    for i, x in enumerate(xi):
        rq = []; rg = []; rf = []; rd = []; rdp = []
        for l, c in enumerate(cc):
            Q, G, D, tri = Q_point(x, c, g0, H[i, l], H1[i, l], H2[i, l], LH[i, l], cd.coefc, cd.mp)
            f, Dp = f_point(x, c, G + eta, cd.coefc, cd.mp)
            rq.append(Q); rg.append(G); rf.append(f); rd.append(D); rdp.append(Dp)
        out['Q'].append(rq); out['G'].append(rg); out['f'].append(rf); out['D'].append(rd); out['Dp'].append(rdp)
        out['tri'].append(1 + x + x * x)
    return out


_DCT = {}


def dct_matrix(n):
    if n in _DCT:
        return _DCT[n]
    rows = []
    for jj in range(n + 1):
        row = []
        for i in range(n + 1):
            v = 2 * arb.cos_pi(arb(jj * i) / n) / n
            if i in (0, n):
                v = v / 2
            if jj in (0, n):
                v = v / 2
            row.append(v)
        rows.append(row)
    _DCT[n] = arb_mat(rows)
    return _DCT[n]


def lobatto(n):
    return [arb.cos_pi(arb(i) / n) for i in range(n + 1)]


def cheb_coeffs(F, n):
    Cm = dct_matrix(n)
    return Cm * arb_mat(F) * Cm.transpose()


def coef_bounds(Co, n):
    c00 = Co[0, 0]
    tot = arb(0)
    for a in range(n + 1):
        for b in range(n + 1):
            if a or b:
                tot += abs(Co[a, b])
    return c00, tot


def cheb_integral(Co, n):
    """int_{-1}^1 int_{-1}^1 sum c_jk T_j(x) T_k(y) dx dy, exactly in ball arithmetic."""
    I = [arb(2) / (1 - j * j) if j % 2 == 0 else arb(0) for j in range(n + 1)]
    tot = arb(0)
    for a in range(0, n + 1, 2):
        for b in range(0, n + 1, 2):
            tot += Co[a, b] * I[a] * I[b]
    return tot


def ellipse_balls(rho, m, full):
    """complex boxes covering the Bernstein ellipse boundary E_rho (upper half if not full)."""
    rho = arb(rho)
    a = (rho + 1 / rho) / 2; b = (rho - 1 / rho) / 2
    span = 2 if full else 1
    dphi = arb(span) / m
    r = (arb.pi() * dphi / 2 * a).upper() * (1 + arb(2) ** -40)
    out = []
    for i in range(m):
        ph = dphi * (i + arb(1) / 2)
        cx = a * arb.cos_pi(ph); cy = b * arb.sin_pi(ph)
        out.append(acb(arb(cx.mid(), (r + cx.rad()).upper()), arb(cy.mid(), (r + cy.rad()).upper())))
    return out


def tail_factor(rho1, rho2, n1, n2):
    r1 = arb(rho1); r2 = arb(rho2)
    s1 = r1 / (r1 - 1); s2 = r2 / (r2 - 1)
    t1 = (1 - r1 ** (-(n1 + 1))) * s1; t2 = (1 - r2 ** (-(n2 + 1))) * s2
    return 8 * (s1 * s2 - t1 * t2)


def certify_patch(cd, xa, xb, ca, cb, n=NSAMP):
    """one patch: rigorous sup Q_N, inf G_N, and a rigorous lower bound of int_patch f. Returns a record or None
    (None = no polyellipse passed the analyticity checks)."""
    nodes = lobatto(n)
    g = patch_grid(cd, xa, xb, ca, cb, nodes, nodes, cplx=False)
    CQ = cheb_coeffs(g['Q'], n); CG = cheb_coeffs(g['G'], n); CF = cheb_coeffs(g['f'], n)
    q00, qtot = coef_bounds(CQ, n)
    g00, gtot = coef_bounds(CG, n)
    pint = cheb_integral(CF, n)
    hwx = (arb_q(xb) - arb_q(xa)) / 2; hwc = (arb_q(cb) - arb_q(ca)) / 2
    best = None
    for (r1, r2) in RHO_COMBOS:
        bx = ellipse_balls(r1, M1, full=False)
        by = ellipse_balls(r2, M2, full=True)
        gc = patch_grid(cd, xa, xb, ca, cb, bx, by, cplx=True)
        ok = all(t.real > 0 for t in gc['tri'])
        MQ = arb(0); MG = arb(0); MF = arb(0)
        for i in range(len(bx)):
            if not ok:
                break
            for l in range(len(by)):
                vq = abs(gc['Q'][i][l]); vg = abs(gc['G'][i][l]); vf = abs(gc['f'][i][l])
                if not (gc['D'][i][l].real > 0) or not (gc['Dp'][i][l].real > 0) or not (
                        vq.is_finite() and vg.is_finite() and vf.is_finite()):
                    ok = False
                    break
                MQ = arb(max(MQ.upper(), vq.upper())); MG = arb(max(MG.upper(), vg.upper()))
                MF = arb(max(MF.upper(), vf.upper()))
        if not ok:
            continue
        tf = tail_factor(r1, r2, n, n)
        eQ = (MQ * tf).upper(); eG = (MG * tf).upper(); eF = (MF * tf).upper()
        sup = (q00 + qtot + eQ).upper()
        ginf = (g00 - gtot - eG).lower()
        qlo = (hwx * hwc * (pint - 4 * eF)).lower()
        qhi = (hwx * hwc * (pint + 4 * eF)).upper()
        if not (sup.is_finite() and ginf.is_finite() and qlo.is_finite()):
            continue
        rec = dict(rho=[r1, r2], MQ=float(MQ.upper()), errQ=float(eQ), errF=float(eF),
                   sup=ball_upper(sup), ginf=ball_lower(ginf), qlo=ball_lower(qlo), qhi=ball_upper(qhi),
                   ptot=float(qtot.upper()))
        if best is None or rec['errQ'] < best['errQ']:
            best = rec
        if rec['errQ'] < 1e-2 * float(qtot.upper()) and rec['errF'] < float(QTOL) / 10:
            break
    return best


# ================================================================== adaptive certification with JSON checkpoint
def patches_initial(nx=NX_INIT):
    pts = [Fraction(-1) + Fraction(2, nx) * i for i in range(nx + 1)]
    return [[fstr(pts[i]), fstr(pts[i + 1]), '0/1', '1/1'] for i in range(nx)]


def patch_good(rec, eta):
    return (rec is not None and rec['sup'] <= eta and rec['ginf'] > 0 and Fraction(rec['errF']) <= QTOL)


def cert_step(st, cd, budget, log, save=None):
    """advance the patch queue in st (dict with todo/done/failed) until empty or budget used (checked before each
    patch; a patch takes seconds). Returns True if the queue is empty. save() is called after every patch.
    The state holds only exact strings (p/q) and no timings, so a resumed run writes the same bytes."""
    t0 = time.time()
    processed = 0
    while st['todo']:
        if processed and time.time() - t0 > budget:      # at least one patch per chunk (progress guaranteed)
            return False
        processed += 1
        p = st['todo'].pop(0)
        xa, xb, ca, cb = [Fraction(s) for s in p]
        rec = certify_patch(cd, xa, xb, ca, cb)
        if patch_good(rec, cd.eta):
            st['done'].append(dict(patch=p, rho=rec['rho'], sup=fstr(rec['sup']), ginf=fstr(rec['ginf']),
                                   qlo=fstr(rec['qlo']), qhi=fstr(rec['qhi']), errQ=fhex(rec['errQ']),
                                   errF=fhex(rec['errF'])))
            log(f"ok    {p} rho {rec['rho']} sup {float(rec['sup']):.3e} ginf {float(rec['ginf']):.5f} "
                f"errQ {rec['errQ']:.1e} errF {rec['errF']:.1e}")
        else:
            w = xb - xa
            if w < MIN_PATCH:
                st['failed'].append(dict(patch=p, rec=None if rec is None else {k: str(v) for k, v in rec.items()}))
                log(f"FAIL  {p} {rec}")
                continue
            if w >= (cb - ca) / 8:
                m = (xa + xb) / 2
                st['todo'][0:0] = [[fstr(xa), fstr(m), p[2], p[3]], [fstr(m), fstr(xb), p[2], p[3]]]
            else:
                m = (ca + cb) / 2
                st['todo'][0:0] = [[p[0], p[1], fstr(ca), fstr(m)], [p[0], p[1], fstr(m), fstr(cb)]]
            log(f"split {p} " + ('no analytic ellipse' if rec is None else
                                 f"sup {float(rec['sup']):.3e} (eta {float(cd.eta):.3e}) ginf {float(rec['ginf']):.4f} "
                                 f"errF {rec['errF']:.1e}"))
        if save:
            save()
    return True


def cover_check(done):
    """exact check that the accepted patches tile [-1,1] x [0,1] (area + pairwise disjoint interiors)."""
    P = [[Fraction(s) for s in d['patch']] for d in done]
    for a, b, c, d in P:
        assert -1 <= a < b <= 1 and 0 <= c < d <= 1, 'patch outside the square'
    srt = sorted(range(len(P)), key=lambda i: P[i][0])
    for ii in range(len(srt)):
        a, b, c, d = P[srt[ii]]
        for kk in range(ii + 1, len(srt)):
            e, f, g, h = P[srt[kk]]
            if e >= b:
                break
            assert min(b, f) <= max(a, e) or min(d, h) <= max(c, g), 'overlap'
    area = sum((b - a) * (d - c) for a, b, c, d in P)
    assert area == 2, f'area {area}'
    return len(P)


def assemble(cand_d, st):
    """from a finished certificate state: exact M_up, M_lo, chi window, all checks asserted."""
    assert not st['todo'] and not st['failed'], 'certification incomplete or failed'
    npatch = cover_check(st['done'])
    g0 = Fraction(unhex(cand_d['g0'])); eta = Fraction(cand_d['eta'])
    sup = max(Fraction(d['sup']) for d in st['done'])
    ginf = min(Fraction(d['ginf']) for d in st['done'])
    assert sup <= eta, 'super-solution inequality not certified'
    assert ginf > 0 and eta >= 0, 'positivity not certified'
    qlo = sum(Fraction(d['qlo']) for d in st['done'])
    qhi = sum(Fraction(d['qhi']) for d in st['done'])
    M_up = 2 * (g0 + eta)
    M_lo = 2 * qlo
    assert M_lo > 0
    chi_lo = 1 / M_up ** 2; chi_hi = 1 / M_lo ** 2
    errF_max = max(unhex(d['errF']) for d in st['done'])
    return dict(patches=npatch, g0=g0, eta=eta, supQ=sup, ginf=ginf, M_up=M_up, M_lo=M_lo, M_int_hi=2 * qhi,
                chi_lo=chi_lo, chi_hi=chi_hi, width=chi_hi - chi_lo, errF_max=errF_max)


def summary(a, digits=14):
    return dict(patches=a['patches'], g0=fstr(a['g0']), eta=fstr(a['eta']), eta_dec=dec_round(a['eta'], 6, True),
                supQ_upper=dec_round(a['supQ'], 6, True) if a['supQ'] != 0 else '0',
                inf_G_N_lower=dec_round(a['ginf'], 10, False),
                M_up=fstr(a['M_up']), M_lo=fstr(a['M_lo']), M_up_dec=dec_round(a['M_up'], digits, True),
                M_lo_dec=dec_round(a['M_lo'], digits, False), M_int_upper_dec=dec_round(a['M_int_hi'], digits, True),
                chi_lo=fstr(a['chi_lo']), chi_hi=fstr(a['chi_hi']), chi_lo_dec=dec_round(a['chi_lo'], digits, False),
                chi_hi_dec=dec_round(a['chi_hi'], digits, True), width_dec=dec_round(a['width'], 4, True),
                quad_err_per_area_max=a['errF_max'])


# ================================================================== candidate build (float solve + a priori eta)
def build_candidate(model, fac, N, K, eta_factor, log):
    mp = model_params(model, Fraction(fac))
    t0 = time.time()
    z, info = float_solve(mp, N, K, log=log)
    ck = float_check(mp, N, K, z)
    supQ = max(ck['supQ'], 0.0)
    eta = dyadic_up(Fraction(eta_factor) * Fraction(supQ) + Fraction(1, 2 ** 80), bits=80)
    d = dict(model=model, fac=fstr(Fraction(fac)), N=N, K=K, g0=fhex(z[0]), h=[fhex(v) for v in z[1:]],
             eta=fstr(eta), eta_rule=f'eta = dyadic_up({eta_factor} * sampled sup Q_N (2001 x 41 grid) + 2^-80)',
             float_info=dict(info, **ck, solve_seconds=round(time.time() - t0, 1), M_float=2 * float(z[0])))
    log(f"candidate {model} fac {fac} N {N} K {K}: g0 {z[0]:.15f} sampled supQ {ck['supQ']:.3e} infQ {ck['infQ']:.3e} "
        f"Gmin {ck['Gmin']:.5f} eta {float(eta):.3e} ({time.time()-t0:.1f} s)")
    return d


def cand_float_z(d):
    return np.concatenate([[unhex(d['g0'])], [unhex(s) for s in d['h']]])


# ================================================================== GO gate and self-sha check
def ledger_path():
    from pathlib import Path
    return str(Path(os.environ.get('LAB_DATA') or (Path.home() / '.pocketlab')) / 'ledger.jsonl')


def go_pinned(sha, ledger=None):
    """True iff the lab ledger has a 'go' row for SEAL_ID naming an instrument b3t.py with exactly this sha256."""
    p = ledger or ledger_path()
    if not os.path.exists(p):
        return False, f'ledger {p} not found'
    for line in open(p, encoding='utf-8'):
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if r.get('kind') == 'go' and r.get('seal_id') == SEAL_ID:
            for i in r.get('instruments', []):
                if str(i.get('path', '')).lower().replace('\\', '/').endswith('/' + INSTRUMENT_NAME) and i.get('sha256') == sha:
                    return True, f"GO row {r.get('at')}"
    return False, f'no GO row for seal {SEAL_ID} names {INSTRUMENT_NAME} with sha256 {sha}'


def main_gate(go_file=GO_FILE, ledger=None):
    if not os.path.exists(go_file):
        return EXIT_NO_GO, f'REFUSED: {go_file} does not exist (the lead arms it with lab.py go after sealing).'
    ok, why = go_pinned(sha_self(), ledger)
    if not ok:
        return EXIT_NOT_PINNED, f'REFUSED: self-sha check failed: {why}'
    return None, 'GO: ' + why


# ================================================================== a full certified run (main and control cases)
def run_case(state_path, cfg, budget, log, sha=None):
    """resumable: candidate (float solve, once) then patch certification. Returns (exit code, assembled or None).
    The state file holds the candidate and the patch queue as exact strings."""
    sha = sha or sha_self()
    t0 = time.time()
    if os.path.exists(state_path):
        st = json.load(open(state_path))
        if st.get('instrument_sha256') != sha or st.get('cfg') != cfg:
            log(f'REFUSED: checkpoint {state_path} was written by {st.get("instrument_sha256")} / cfg {st.get("cfg")}')
            return EXIT_CKPT_REFUSED, None
    else:
        st = dict(instrument_sha256=sha, cfg=cfg, cand=None, todo=patches_initial(), done=[], failed=[],
                  prec=PREC, nsamp=NSAMP)
    if st['cand'] is None:
        st['cand'] = build_candidate(cfg['model'], cfg['fac'], cfg['N'], cfg['K'], cfg['eta_factor'], log)
        st['cand']['float_info'].pop('solve_seconds', None)       # no timings in the state (byte-identical resume)
        atomic_json(state_path, st)
    cd = Cand(st['cand'])
    done = cert_step(st, cd, budget - (time.time() - t0), log, save=lambda: atomic_json(state_path, st))
    atomic_json(state_path, st)
    if not done:
        log(f'chunk budget used: {len(st["done"])} patches done, {len(st["todo"])} queued -> exit 3 (rerun)')
        return EXIT_RESUME, None
    if st['failed']:
        log(f'CERTIFICATE FAILED on {len(st["failed"])} minimal patches')
        return EXIT_CERT_FAILED, None
    return EXIT_FINISHED, assemble(st['cand'], st)


# ================================================================== BAR-B3T scoring
def score(a):
    lo, hi = a['chi_lo'], a['chi_hi']
    T1 = True                                           # assemble() asserted every certificate condition
    T2 = a['width'] <= Fraction(1, 10 ** 9)
    T3 = P2E_WINDOW[0] <= lo and hi <= P2E_WINDOW[1]
    dist = Fraction(0) if lo <= T4_PROBE <= hi else min(abs(T4_PROBE - lo), abs(T4_PROBE - hi))
    T4 = dist <= Fraction(1, 10 ** 9)
    return dict(T1_certificate=T1, T2_width_le_1e9=bool(T2), T2_width=dec_round(a['width'], 4, True),
                T3_inside_P2e=bool(T3), T4_within_1e9_of_probe=bool(T4),
                T4_distance='0' if dist == 0 else dec_round(dist, 4, True))


def run_main(out_path, log):
    rc, why = main_gate()
    log(why)
    if rc is not None:
        return rc
    rc, a = run_case(out_path + '.ckpt', MAIN_CFG, TIME_BUDGET, log)
    if rc != EXIT_FINISHED:
        return rc
    st = json.load(open(out_path + '.ckpt'))
    res = dict(seal=SEAL_ID, instrument_sha256=sha_self(), cfg=MAIN_CFG, prec=PREC, nsamp=NSAMP,
               object='extreme (mu = 0) Bowen-York, J = 1; chi(0) = J / M_0^2',
               window=summary(a), bar=score(a), float_info=st['cand']['float_info'],
               depends_on='B3-5 + B3-2 (scout/B3-EXTREME.md, referee R pending); E1, E3 (LEMMA-REP4); arb ball '
                          'arithmetic; polyellipse/Cauchy patch bound (REP6 logic, unrefereed)')
    atomic_json(out_path, res)
    log(f"RESULT chi(0) in [{res['window']['chi_lo_dec']}, {res['window']['chi_hi_dec']}] width "
        f"{res['window']['width_dec']}; bar {res['bar']}")
    return EXIT_FINISHED


# ================================================================== controls
CTRL_CASES = {
    'RN': dict(model='RN', fac='1', N=64, K=0, eta_factor='2'),
    'NULL': dict(model='BY', fac='1001/1000', N=MAIN_CFG['N'], K=MAIN_CFG['K'], eta_factor=MAIN_CFG['eta_factor']),
}
PARTS = ['C0', 'RN', 'NULL', 'MARGIN', 'COVER', 'NANINF', 'RESUME', 'GATE', 'MEM']


def ctrl_C0():
    """brute-force mpmath Laplacian of w vs the Q formula (both models, 5 points each) and the mass-integrand identity."""
    import mpmath as mpm
    mpm.mp.dps = 40
    out = {}
    worst_all = mpm.mpf(0)
    for model in ('BY', 'RN'):
        mpp = model_params(model, 1)
        rng = np.random.default_rng(1)
        N, K = 5, 2
        h = rng.normal(size=(N + 1, K + 1)) * 0.1
        g0 = mpm.mpf('0.6')
        p = mpp['p']; coef = mpm.mpf(mpp['coef'].numerator) / mpp['coef'].denominator

        def Hf(x, c):
            return sum(h[n, k] * mpm.chebyt(n, x) * mpm.legendre(2 * k, c) for n in range(N + 1) for k in range(K + 1))

        def Gf(x, c):
            return g0 + (1 - x) ** 2 * Hf(x, c)

        def wf(r, th):
            return r ** mpm.mpf(-0.5) * Gf(mpm.tanh(mpm.log(r) / 8), mpm.cos(th)) / (1 + mpm.sqrt(r))

        worst = mpm.mpf(0)
        for (rs, ts) in [('0.003', '0.7'), ('0.4', '1.1'), ('2.5', '0.3'), ('70', '1.4'), ('1', '0.9')]:
            r = mpm.mpf(rs); th = mpm.mpf(ts)
            lap = mpm.diff(lambda rr: wf(rr, th), r, 2) + 2 / r * mpm.diff(lambda rr: wf(rr, th), r, 1) + (
                mpm.diff(lambda tt: wf(r, tt), th, 2)
                + mpm.cos(th) / mpm.sin(th) * mpm.diff(lambda tt: wf(r, tt), th, 1)) / r ** 2
            ang = mpm.sin(th) ** 2 if model == 'BY' else 1
            lhs = lap + coef * ang * r ** (-(p + 5) / mpm.mpf(2)) * (1 + wf(r, th)) ** -p
            x = mpm.tanh(mpm.log(r) / 8); c = mpm.cos(th); s = mpm.sqrt(r)
            Hx = mpm.diff(lambda xx: Hf(xx, c), x, 1); Hxx = mpm.diff(lambda xx: Hf(xx, c), x, 2)
            LH = (1 - c * c) * mpm.diff(lambda cc: Hf(x, cc), c, 2) - 2 * c * mpm.diff(lambda cc: Hf(x, cc), c, 1)
            Q, _, _, _ = Q_point(x, c, g0, Hf(x, c), Hx, Hxx, LH, coef * 2 ** (p + 3), mpp)
            rhs = r ** mpm.mpf(-2.5) / (1 + s) * Q * (1 + 3 * s) / (4 * (1 + s) ** 2)
            worst = max(worst, abs(lhs - rhs) / abs(lhs))
        r, th = mpm.mpf('0.37'), mpm.mpf('0.8')
        x = mpm.tanh(mpm.log(r) / 8); c = mpm.cos(th)
        ang = mpm.sin(th) ** 2 if model == 'BY' else 1
        dens = coef * ang * r ** (-(p + 5) / mpm.mpf(2)) * (1 + wf(r, th)) ** -p * r ** 3 * 8 / (1 - x * x)
        fv, _ = f_point(x, c, Gf(x, c), coef * 2 ** (p + 3), mpp)
        rel_f = abs(dens - fv) / fv
        out[model] = dict(residual_worst_rel=mpm.nstr(worst, 3), integrand_rel=mpm.nstr(rel_f, 3))
        worst_all = max(worst_all, worst, rel_f)
    out['ok'] = bool(worst_all < mpm.mpf(10) ** -30)
    return out


def case_path(name, scratch):
    return os.path.join(scratch, f'ctrl_{name}_state.json')


def point_Q(cd, xq, cq, eta=None):
    """arb value of Q at an exact point (xq, cq) for candidate cd (local polys on a tiny patch starting at xq)."""
    loc = cd.local(xq, xq + Fraction(1, 2 ** 30), Fraction(0), Fraction(1))
    H, H1, H2, LH, xi, cc = fields(cd, loc, [arb(-1)], [arb(-1) + 2 * arb_q(cq)], cplx=False)
    Q, _, _, _ = Q_point(xi[0], cc[0], arb(cd.g0), H[0, 0], H1[0, 0], H2[0, 0], LH[0, 0], cd.coefc, cd.mp,
                         0 if eta is None else arb_q(eta))
    return Q


def ctrl_margin(scratch, log):
    """margin null on the certified NULL candidate: eta' = half the arb residual at the float argmax (below the
    certified sup). (a) the patch certificate must reject; (b) a rigorous point witness Q[G_N + eta'] > 0 must exist;
    (c) the assembly with eta' must refuse."""
    st = json.load(open(case_path('NULL', scratch)))
    cd = Cand(st['cand'])
    xs, cs = st['cand']['float_info']['argmax']
    xq, cq = Fraction(xs), Fraction(cs)
    qn_lo = ball_lower(point_Q(cd, xq, cq))
    assert qn_lo > 0, 'float argmax residual not positive in arb'
    for div in (2, 4, 8, 16, 32, 64):                  # the source also drops with eta: halve until a witness holds
        eta_p = dyadic_up(qn_lo / div)
        witness_lo = ball_lower(point_Q(cd, xq, cq, eta_p))
        if witness_lo > 0:
            break
    cert_sup = max(Fraction(d['sup']) for d in st['done'])
    pt = next(d for d in st['done'] if Fraction(d['patch'][0]) <= xq <= Fraction(d['patch'][1])
              and Fraction(d['patch'][2]) <= cq <= Fraction(d['patch'][3]))
    d2 = dict(st['cand']); d2['eta'] = fstr(eta_p)
    cd2 = Cand(d2)
    rec = certify_patch(cd2, *[Fraction(s) for s in pt['patch']])
    accepted = patch_good(rec, cd2.eta)
    try:
        assemble(d2, st)
        assembled = True
    except AssertionError:
        assembled = False
    res = dict(eta_null=fstr(eta_p), eta_null_dec=dec_round(eta_p, 6, True), certified_sup=dec_round(cert_sup, 6, True),
               eta_null_below_certified_sup=bool(eta_p < cert_sup), witness_point=[xs, cs], witness_patch=pt['patch'],
               witness_Q_shifted_lower=dec_round(witness_lo, 6, False), witness_proves_not_super=bool(witness_lo > 0),
               margin_patch_accepted=bool(accepted), margin_assembly_succeeded=bool(assembled))
    res['margin_null_FAILED_as_required'] = bool(res['eta_null_below_certified_sup'] and res['witness_proves_not_super']
                                                 and not accepted and not assembled)
    log(f"MARGIN: eta' {res['eta_null_dec']} < certified sup {res['certified_sup']}; witness Q >= "
        f"{res['witness_Q_shifted_lower']}; patch accepted {accepted}; assembly {assembled}")
    return res


def ctrl_cover(scratch):
    out = {}
    ok = True
    for name in ('RN', 'NULL'):
        st = json.load(open(case_path(name, scratch)))
        out[f'{name}_patches'] = cover_check(st['done'])
        for tag, lst in (('missing_patch_cover_ok', st['done'][1:]), ('duplicate_patch_cover_ok', st['done'] + st['done'][:1])):
            try:
                cover_check(lst)
                r = True
            except AssertionError:
                r = False
            out[f'{name}_{tag}'] = r
            ok = ok and not r
    out['ok'] = ok
    return out


def ctrl_naninf(scratch):
    st = json.load(open(case_path('RN', scratch)))
    acc = {}
    for tag, bad in (('nan', float('nan')), ('inf', float('inf'))):
        d = dict(st['cand']); h = list(d['h']); h[3] = fhex(bad); d['h'] = h
        rec = certify_patch(Cand(d), Fraction(-1, 2), Fraction(-3, 8), Fraction(0), Fraction(1))
        acc[tag] = patch_good(rec, Fraction(d['eta']))
    def bad(v):
        if isinstance(v, float):
            return not math.isfinite(v)
        if isinstance(v, str):
            return v.strip().lower().lstrip('+-') in ('nan', 'inf', 'infinity')
        if isinstance(v, dict):
            return any(bad(x) for x in v.values())
        if isinstance(v, list):
            return any(bad(x) for x in v)
        return False

    clean = not any(bad(json.load(open(case_path(n, scratch)), parse_constant=lambda s: float(s)))
                    for n in ('RN', 'NULL'))
    return dict(nan_patch_accepted=bool(acc['nan']), inf_patch_accepted=bool(acc['inf']), states_free_of_nan_inf=clean,
                ok=bool(not acc['nan'] and not acc['inf'] and clean))


def ctrl_resume(scratch, log):
    a_path = os.path.join(scratch, 'resume_A.json'); b_path = os.path.join(scratch, 'resume_B.json')
    for p in (a_path, b_path):
        if os.path.exists(p):
            os.remove(p)
    cfg = CTRL_CASES['RN']

    def quiet(s):
        pass

    rcA, aA = run_case(a_path, cfg, 1e9, quiet)
    chunks = 0
    while True:
        rcB, aB = run_case(b_path, cfg, 0.0, quiet)     # budget 0: exactly one patch per chunk, then exit 3
        chunks += 1
        if rcB != EXIT_RESUME or chunks > 500:
            break
    same_state = open(a_path, 'rb').read() == open(b_path, 'rb').read()
    same_result = (json.dumps(summary(aA)) == json.dumps(summary(aB))) if (aA and aB) else False
    c_path = os.path.join(scratch, 'resume_C.json')
    st = json.load(open(b_path)); st['instrument_sha256'] = '0' * 64; st['todo'] = patches_initial()
    atomic_json(c_path, st)
    rcC, _ = run_case(c_path, cfg, 1e9, quiet)
    res = dict(rc_A=rcA, rc_B_final=rcB, chunks_B=chunks, state_bytes_identical=bool(same_state),
               result_identical=bool(same_result), foreign_checkpoint_rc=rcC,
               foreign_checkpoint_refused=bool(rcC == EXIT_CKPT_REFUSED))
    res['ok'] = bool(rcA == 0 and rcB == 0 and chunks >= 2 and same_state and same_result and rcC == EXIT_CKPT_REFUSED)
    log(f'RESUME: {res}')
    return res


def ctrl_gate(scratch):
    fake_go = os.path.join(scratch, 'GO-fake-absent')
    if os.path.exists(fake_go):
        os.remove(fake_go)
    rc_nogo, msg_nogo = main_gate(go_file=fake_go)
    present = os.path.join(scratch, 'GO-fake-present')
    with open(present, 'w') as f:
        f.write('scratch GO file for the self-sha control (not the real GO-B3T)\n')
    me = sha_self()
    rows = {'wrong_sha': dict(kind='go', seal_id=SEAL_ID, instruments=[dict(path='<lab>/bhslice/b3t.py', sha256='0' * 64)]),
            'other_seal': dict(kind='go', seal_id='560847f6a4', instruments=[dict(path='<lab>/bhslice/b3t.py', sha256=me)]),
            'other_file': dict(kind='go', seal_id=SEAL_ID, instruments=[dict(path='<lab>/bhslice/p2e.py', sha256=me)]),
            'right': dict(kind='go', seal_id=SEAL_ID, at='controls', instruments=[dict(path='<lab>/bhslice/b3t.py', sha256=me)])}
    led = {}
    for k, row in rows.items():
        lp = os.path.join(scratch, f'ledger_{k}.jsonl')
        with open(lp, 'w') as f:
            f.write(json.dumps(row) + '\n')
        led[k] = main_gate(go_file=present, ledger=lp)[0]
    # --main without GO-B3T through the real entry point (only if the real GO file is absent): no output written
    out_fake = os.path.join(scratch, 'main_refusal_out.json')
    if os.path.exists(out_fake):
        os.remove(out_fake)
    real_go = os.path.exists(GO_FILE)
    rc_main = run_main(out_fake, lambda s: None) if not real_go else None
    res = dict(no_go_rc=rc_nogo, no_go_msg=msg_nogo, run_main_rc_without_GO=rc_main,
               no_go_output_exists=bool(os.path.exists(out_fake) or os.path.exists(out_fake + '.ckpt')),
               wrong_sha_rc=led['wrong_sha'], other_seal_rc=led['other_seal'], other_file_rc=led['other_file'],
               right_sha_passes=led['right'] is None, real_GO_file_exists=real_go)
    res['ok'] = bool(rc_nogo == EXIT_NO_GO and led['wrong_sha'] == EXIT_NOT_PINNED and led['other_seal'] == EXIT_NOT_PINNED
                     and led['other_file'] == EXIT_NOT_PINNED and led['right'] is None and not res['no_go_output_exists']
                     and (real_go or rc_main == EXIT_NO_GO))
    return res


EXPECTED_FALSE = [
    'parts.NULL.overlaps_P2e_window (the 1e-3 source null must move the window off P2e)',
    'parts.NULL.inside_P2e_window (same: the null window must not sit in the P2e window)',
    'parts.MARGIN.margin_patch_accepted (eta below the certified sup: the patch certificate must reject)',
    'parts.MARGIN.margin_assembly_succeeded (eta below the certified sup: the assembly must refuse)',
    'parts.COVER.RN_missing_patch_cover_ok (one patch dropped: the exact cover check must fail)',
    'parts.COVER.RN_duplicate_patch_cover_ok (one patch duplicated: the overlap check must fail)',
    'parts.COVER.NULL_missing_patch_cover_ok (one patch dropped: the exact cover check must fail)',
    'parts.COVER.NULL_duplicate_patch_cover_ok (one patch duplicated: the overlap check must fail)',
    'parts.NANINF.nan_patch_accepted (mock variant: NaN coefficient injected; the patch must be rejected)',
    'parts.NANINF.inf_patch_accepted (mock variant: inf coefficient injected; the patch must be rejected)',
    'parts.RESUME.foreign_checkpoint_refused is True but its rc is 4, not 0 (mock checkpoint with a zero sha)',
    'parts.GATE.no_go_output_exists (no output without GO-B3T)',
    'parts.GATE.real_GO_file_exists (GO-B3T is armed by the lead after the controls; absent now)',
    'parts.GATE.wrong_sha_rc / other_seal_rc / other_file_rc are 6, not 0 (mock scratch ledgers that must not pin)',
    'parts.GATE.right_sha_passes (mock: scratch ledger + scratch GO file; the real ledger is never written)',
]


def run_controls(out_path, log):
    t0 = time.time()
    sha = sha_self()
    os.makedirs(SCRATCH, exist_ok=True)
    R = None
    if os.path.exists(out_path):
        R = json.load(open(out_path))
        if R.get('instrument_sha256') != sha:
            log(f'controls file {out_path} is from instrument {R.get("instrument_sha256")}: starting fresh')
            R = None
    if R is None:
        for n in ('RN', 'NULL'):
            if os.path.exists(case_path(n, SCRATCH)):
                os.remove(case_path(n, SCRATCH))
        R = dict(instrument_sha256=sha, seal=SEAL_ID, parts={}, oks={}, expected_false=EXPECTED_FALSE,
                 cases=CTRL_CASES, main_cfg=MAIN_CFG, all_ok=False, resume_tested=False, chunk_peaks=[])
    P, oks = R['parts'], R['oks']

    def left():
        return TIME_BUDGET - (time.time() - t0)

    for part in PARTS:
        if part in P:
            continue
        if left() <= 0:
            break
        log(f'--- part {part}')
        if part == 'C0':
            P[part] = ctrl_C0(); oks[part] = P[part]['ok']
        elif part in ('RN', 'NULL'):
            rc, a = run_case(case_path(part, SCRATCH), CTRL_CASES[part], left(), log, sha=sha)
            if rc == EXIT_RESUME:
                break
            if rc != EXIT_FINISHED:
                P[part] = dict(rc=rc, ok=False); oks[part] = False
            else:
                s = summary(a)
                if part == 'RN':
                    s['window_contains_1'] = bool(a['M_lo'] <= 1 <= a['M_up'] and a['chi_lo'] <= 1 <= a['chi_hi'])
                    oks[part] = s['window_contains_1']
                else:
                    lo, hi = a['chi_lo'], a['chi_hi']
                    s['overlaps_P2e_window'] = bool(not (hi < P2E_WINDOW[0] or lo > P2E_WINDOW[1]))
                    s['inside_P2e_window'] = bool(P2E_WINDOW[0] <= lo and hi <= P2E_WINDOW[1])
                    s['expected_shift_note'] = 'chi scales as fac^-1/2 at fixed J: about -4.6e-4 for fac = 1.001'
                    oks[part] = not s['overlaps_P2e_window']
                s['float_info'] = json.load(open(case_path(part, SCRATCH)))['cand']['float_info']
                P[part] = s
        elif part == 'MARGIN':
            P[part] = ctrl_margin(SCRATCH, log); oks[part] = P[part]['margin_null_FAILED_as_required']
        elif part == 'COVER':
            P[part] = ctrl_cover(SCRATCH); oks[part] = P[part]['ok']
        elif part == 'NANINF':
            P[part] = ctrl_naninf(SCRATCH); oks[part] = P[part]['ok']
        elif part == 'RESUME':
            P[part] = ctrl_resume(SCRATCH, log); oks[part] = P[part]['ok']
        elif part == 'GATE':
            P[part] = ctrl_gate(SCRATCH); oks[part] = P[part]['ok']
        elif part == 'MEM':
            m = peak_memory()
            peaks = [x for x in R['chunk_peaks'] if x] + ([m] if m else [])
            P[part] = dict(peak_this_chunk=m, earlier_chunk_peaks=R['chunk_peaks'], cap_gb=8)
            oks[part] = bool(m is not None and all(x['peak_commit_gb'] < 8 and x['peak_working_set_gb'] < 8 for x in peaks))
        log(f'part {part}: ok = {oks.get(part)}')
        R['parts_done'] = [p for p in PARTS if p in P]
        atomic_json(out_path, R)
    R['parts_done'] = [p for p in PARTS if p in P]
    R['resume_tested'] = bool(oks.get('RESUME', False))
    R['all_ok'] = bool(set(PARTS) <= set(P.keys()) and all(oks.values()))
    R['chunk_peaks'].append(peak_memory())
    atomic_json(out_path, R)
    log(f"controls parts done {R['parts_done']}; all_ok = {R['all_ok']}; oks = {oks}")
    if not set(PARTS) <= set(P.keys()):
        return EXIT_RESUME
    return EXIT_FINISHED if R['all_ok'] else EXIT_CONTROLS_FAILED


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--controls')
    ap.add_argument('--main')
    args = ap.parse_args()
    logp = (args.controls or args.main or 'b3t') + '.log'

    def log(s):
        line = f"{time.strftime('%H:%M:%S')} {s}"
        print(line, flush=True)
        with open(logp, 'a', encoding='utf-8') as f:
            f.write(line + '\n')

    if args.controls:
        log(f'--- controls invocation sha256 {sha_self()}')
        sys.exit(run_controls(args.controls, log))
    if args.main:
        log(f'--- main invocation sha256 {sha_self()}')
        sys.exit(run_main(args.main, log))
    ap.print_help()
    sys.exit(EXIT_CONTROLS_FAILED)


if __name__ == '__main__':
    main()
