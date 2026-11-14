"""REP6 blind replication seat (BHSLICE): rigorous sandwich for the Bowen-York spinning puncture at j = 2^21
in a log-tanh compactification, with an analytic-interpolation positivity certificate.

Formulation (derived in RESULT.md, checked by ctrl_formula.py):
  mu = 1, sigma = log(2r), xi = tanh(sigma/8), c = cos(theta), A = 2(1+6xi^2+xi^4), w = 1-xi^2,
  u = 2G/(2r+1),  G = g0 + w^4 H(xi,c),  delta = M_ADM - 1 = 2 g0.
  Rhat = r^2 (Delta u + S Phi^-7) / (e^{-sigma/2} sech^3(sigma/2)) (formula in rep6num.py docstring).
  Super: G_N + eta_p gives Rhat <= Rhat_N - eta_p/2;  sub: G_N - eta_m gives Rhat >= Rhat_N + eta_m/2 (if E_- > 0).
  Tangent: psi = 2Z/(2r+1), Z = z0 + z1 (1+t)/2 + w^4 HZ, t = tanh(sigma/2) = 8 xi (1+xi^2)/A;
  That = r^2 (Delta psi - V1 (1/(2r) + psi)) / (e^{-sigma/2} sech^3); adding kappa to z0 lowers That by >= kappa/2.

Positivity certificate per patch [xa,xb] x [ca,cb] (local variables x, y in [-1,1]):
  f = Rhat (or That, or E) is analytic on a neighbourhood of the closed product of Bernstein ellipses E_rho1 x E_rho2
  (certified Re E > 0, Re A > 0 on the distinguished boundary + minimum principle for the pluriharmonic Re E),
  |f| <= M there (ball arithmetic on a cover of the distinguished boundary + maximum modulus, conj symmetry);
  Chebyshev coefficients |a_jk| <= 4 M rho1^-j rho2^-k (Cauchy on the Joukowski torus);
  tensor Chebyshev-Lobatto interpolant p of degree (n1, n2) has ||f - p|| <= 2 sum_{outside} |a_jk|
  = 8 M (s1 s2 - t1 t2) (aliasing: I T_j = T_j' exactly at Lobatto nodes);  sup p <= c00 + sum' |c_jk|.
Usage: rep6.py solve | certify <case> | final        (exit 0 = finished, 3 = rerun to resume)
"""
import sys, os, json, time
from fractions import Fraction
import numpy as np
from flint import arb, acb, arb_poly, acb_poly, arb_mat, acb_mat, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rep6num as RN

PREC = int(os.environ.get('REP6_PREC', '256'))
ctx.prec = PREC
NSAMP = int(os.environ.get('REP6_NSAMP', '64'))    # n1 = n2
M1 = 24     # arcs on the upper half xi-ellipse
M2 = 48     # arcs on the full c-ellipse
CHUNK_SECONDS = 480


# ---------------------------------------------------------------------------------------------- utilities
def fhex(x):
    return float(x).hex()


def unhex(s):
    return float.fromhex(s)


def arb_exact(x):
    """arb from a Python float or Fraction (dyadic floats are exact)."""
    if isinstance(x, Fraction):
        return arb(x.numerator) / arb(x.denominator)
    return arb(float(x))


def frac_of_float(x):
    return Fraction(float(x))


def _ball_upper(a):
    # a.upper() is an arb point (exact dyadic) >= the ball; convert exactly through man_exp
    u = a.upper()
    return _arb_point_to_frac(u, up=True)


def _ball_lower(a):
    return _arb_point_to_frac(a.lower(), up=False)


def _arb_point_to_frac(p, up):
    m, e = p.mid().man_exp()
    m = int(m); e = int(e)
    return Fraction(m) * (Fraction(2) ** e)


def dec_round(q, digits, up):
    """Decimal string with `digits` significant digits, rounded up (up=True) or down."""
    q = Fraction(q)
    if q == 0:
        return '0'
    import math
    neg = q < 0
    a = -q if neg else q
    e = int(math.floor(math.log10(a.numerator) - math.log10(a.denominator)))
    scale = Fraction(10) ** (digits - 1 - e)
    v = a * scale
    # direction for magnitude
    mag_up = (up and not neg) or ((not up) and neg)
    n = v.numerator // v.denominator
    if mag_up and n * v.denominator != v.numerator:
        n += 1
    s = Fraction(n) / scale
    if neg:
        s = -s
    # format exactly
    exp10 = digits - 1 - e
    if exp10 > 0:
        txt = str(abs(n)).rjust(exp10 + 1, '0')
        txt = txt[:-exp10] + '.' + txt[-exp10:]
    else:
        txt = str(abs(n)) + '0' * (-exp10)
    out = ('-' if neg else '') + txt
    assert (Fraction(out) >= q) if up else (Fraction(out) <= q)
    return out


# ---------------------------------------------------------------------------------------------- local polynomials
def cheb_local(a, b, nmax):
    """T_n(m + hw x) as arb_poly in x, n = 0..nmax."""
    m = (arb_exact(a) + arb_exact(b)) / 2; hw = (arb_exact(b) - arb_exact(a)) / 2
    X = arb_poly([m, hw])
    T = [arb_poly([1]), X]
    for n in range(1, nmax):
        T.append(2 * X * T[n] - T[n - 1])
    return T[:nmax + 1], hw


def leg_local(a, b, kmax2):
    m = (arb_exact(a) + arb_exact(b)) / 2; hw = (arb_exact(b) - arb_exact(a)) / 2
    X = arb_poly([m, hw])
    P = [arb_poly([1]), X]
    for n in range(1, kmax2):
        P.append(((2 * n + 1) * X * P[n] - n * P[n - 1]) * (1 / arb(n + 1)))
    return P[:kmax2 + 1], hw


class Sol:
    """A certified-candidate solution: G = g0 + w^4 sum_{n,k} h[n,k] T_{bn}(xi) P_{2k}(c) (b = 2 even / 1 all)."""

    def __init__(self, g0, h, even=True):
        self.g0 = float(g0); self.h = np.asarray(h, float); self.even = even
        self.N, self.K1 = self.h.shape

    def local(self, xa, xb, ca, cb):
        deg = 2 * (self.N - 1) if self.even else self.N - 1
        T, hwx = cheb_local(xa, xb, max(deg, 1))
        P, hwc = leg_local(ca, cb, max(2 * (self.K1 - 1), 1))
        hk = []
        for k in range(self.K1):
            p = arb_poly([0])
            for n in range(self.N):
                cf = self.h[n, k]
                if cf != 0.0:
                    p += arb(float(cf)) * T[2 * n if self.even else n]
            hk.append(p)
        hk1 = [p.derivative() * (1 / hwx) for p in hk]
        hk2 = [p.derivative() * (1 / hwx) for p in hk1]
        hk3 = [p.derivative() * (1 / hwx) for p in hk2]
        Pk = [P[2 * k] for k in range(self.K1)]
        lam = [-(2 * k) * (2 * k + 1) for k in range(self.K1)]
        return dict(hk=hk, hk1=hk1, hk2=hk2, hk3=hk3, Pk=Pk, lam=lam, hwx=hwx, hwc=hwc,
                    mx=(arb_exact(xa) + arb_exact(xb)) / 2, mc=(arb_exact(ca) + arb_exact(cb)) / 2)


def fields(sol, loc, xs, cs, cplx, extra=False):
    """H, H_xi, H_xixi, Lap_S H on the grid xs (local x) times cs (local y); returns lists of lists + xi, c values.
    extra=True also returns H_xixixi and d_xi Lap_S H."""
    MAT = acb_mat if cplx else arb_mat
    POLY = (lambda p: acb_poly(p)) if cplx else (lambda p: p)
    hk = [POLY(p) for p in loc['hk']]; hk1 = [POLY(p) for p in loc['hk1']]; hk2 = [POLY(p) for p in loc['hk2']]
    Pk = [POLY(p) for p in loc['Pk']]
    K1 = sol.K1
    Hx = MAT([[hk[k](x) for k in range(K1)] for x in xs])
    Hx1 = MAT([[hk1[k](x) for k in range(K1)] for x in xs])
    Hx2 = MAT([[hk2[k](x) for k in range(K1)] for x in xs])
    PyT = MAT([[Pk[k](y) for y in cs] for k in range(K1)])
    PLT = MAT([[loc['lam'][k] * Pk[k](y) for y in cs] for k in range(K1)])
    H = Hx * PyT; H1 = Hx1 * PyT; H2 = Hx2 * PyT; LH = Hx * PLT
    xi = [loc['mx'] + loc['hwx'] * x for x in xs]
    cc = [loc['mc'] + loc['hwc'] * y for y in cs]
    if not extra:
        return H, H1, H2, LH, xi, cc
    hk3 = [POLY(p) for p in loc['hk3']]
    Hx3 = MAT([[hk3[k](x) for k in range(K1)] for x in xs])
    return H, H1, H2, LH, xi, cc, Hx3 * PyT, Hx1 * PLT


def lin_terms(x, A, w, H, H1, H2, LH):
    return (A * A / 256) * (w * w * H2 - 18 * x * w * H1 + (64 * x * x - 8 * w) * H) \
        - (x * (1 + x * x) * A / 4) * (w * H1 - 8 * x * H) + (A * A / 4) * LH


# ---------------------------------------------------------------------------------------------- function families
class MainF:
    """Rhat for G_N (and E_N)."""

    def __init__(self, sol, CS):
        self.sol = sol; self.CS = CS

    def grid(self, xa, xb, ca, cb, xs, cs, cplx):
        sol = self.sol
        loc = sol.local(xa, xb, ca, cb)
        H, H1, H2, LH, xi, cc = fields(sol, loc, xs, cs, cplx)
        g0 = arb(sol.g0)
        R = []; Eg = []; ReA = []
        for i, x in enumerate(xi):
            A = 2 * (1 + 6 * x * x + x ** 4); w = 1 - x * x; w4 = w ** 4
            A2 = A * A; A10 = A2 ** 5; w8 = w4 * w4
            rowR = []; rowE = []
            for l, c in enumerate(cc):
                h = H[i, l]
                G = g0 + w4 * h
                E = A2 + 2 * w4 * G
                lin = lin_terms(x, A, w, h, H1[i, l], H2[i, l], LH[i, l]) - G / 2
                src = self.CS * (1 - c * c) * A10 * w8 / E ** 7
                rowR.append(lin + src); rowE.append(E)
            R.append(rowR); Eg.append(rowE); ReA.append(A)
        return {'R': R, 'E': Eg, 'A': ReA}


class Dual:
    """First-order forward-mode dual number over flint / numpy scalars (value a, xi-derivative b)."""
    __slots__ = ('a', 'b')

    def __init__(self, a, b=0):
        self.a = a; self.b = b

    @staticmethod
    def _c(o):
        return o if isinstance(o, Dual) else Dual(o, 0)

    def __add__(self, o):
        o = Dual._c(o); return Dual(self.a + o.a, self.b + o.b)
    __radd__ = __add__

    def __sub__(self, o):
        o = Dual._c(o); return Dual(self.a - o.a, self.b - o.b)

    def __rsub__(self, o):
        o = Dual._c(o); return Dual(o.a - self.a, o.b - self.b)

    def __mul__(self, o):
        o = Dual._c(o); return Dual(self.a * o.a, self.a * o.b + self.b * o.a)
    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, Dual):
            return Dual(self.a / o.a, (self.b * o.a - self.a * o.b) / (o.a * o.a))
        return Dual(self.a / o, self.b / o)

    def __neg__(self):
        return Dual(-self.a, -self.b)

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 1
        r = self
        for _ in range(n - 1):
            r = r * self
        return r


def lin_N_and_D(x, g0, H, H1, H2, H3, LH, LH1):
    """Lin_N = linear part of Rhat (incl. -G/2) and its xi-derivative."""
    X = Dual(x, 1)
    A = 2 * (1 + 6 * X * X + X ** 4); w = 1 - X * X
    val = lin_terms(X, A, w, Dual(H, H1), Dual(H1, H2), Dual(H2, H3), Dual(LH, LH1)) \
        - (g0 + (w ** 4) * Dual(H, H1)) / 2
    return val.a, val.b


def tan_point(x, c, g0, g0m, zeta0, c63a, c63b, hg, hg1, hg2, hg3, lhg, lhg1, h, h1, h2, lh, dscale=1):
    """Normalised tangent residual That at one point (generic scalar type). Returns (That, E_-, A)."""
    A = 2 * (1 + 6 * x * x + x ** 4); w = 1 - x * x; w4 = w ** 4
    A2 = A * A; A10 = A2 ** 5; w8 = w4 * w4; w12 = w8 * w4
    t = 8 * x * (1 + x * x) / A
    om4 = (1 - x) ** 4
    G = g0 + w4 * hg
    DG = w4 * (w * hg1 / 8 - x * hg)
    ZD = dscale * ((1 + t) * G / 2 - DG)
    LN, dLN = lin_N_and_D(x, g0, hg, hg1, hg2, hg3, lhg, lhg1)
    linD = dscale * ((1 + 3 * t) / 2 * LN - (w / 8) * dLN)
    Ze = zeta0 + w4 * h
    line = lin_terms(x, A, w, h, h1, h2, lh) - Ze / 2
    Em = A2 + 2 * w4 * (g0m + w4 * hg)
    s = 1 - c * c
    E8 = Em ** 8
    Z = ZD + Ze
    pot = c63a * s * A10 * w12 * Z / E8
    frc = c63b * s * A10 * A * w8 * om4 / E8
    return linD + line - pot - frc, Em, A


class TanF:
    """That for psi_+ = 2 Z/(2r+1), Z = Z_D + Z_e,
       Z_D = (1+t) G_N/2 - D G_N  (so 2 Z_D/(2r+1) = -r d_r u_N),  Z_e = zeta0 + w^4 H_e (H_e even in xi),
       against the certified sub-solution u_- = 2 (G_N - eta_m)/(2r+1).
       lin(Z_D) = (1+3t)/2 Lin_N - D Lin_N  (L has constant coefficients in sigma)."""

    def __init__(self, sol, eta_m, esol, zeta0, j2, dscale=1):
        self.sol = sol; self.eta_m = eta_m; self.esol = esol; self.zeta0 = zeta0; self.j2 = j2; self.dscale = dscale

    def grid(self, xa, xb, ca, cb, xs, cs, cplx):
        sol, es = self.sol, self.esol
        loc = sol.local(xa, xb, ca, cb)
        Hg, Hg1, Hg2, LHg, xi, cc, Hg3, LHg1 = fields(sol, loc, xs, cs, cplx, extra=True)
        loce = es.local(xa, xb, ca, cb)
        H, H1, H2, LH, _, _ = fields(es, loce, xs, cs, cplx)
        g0 = arb(sol.g0)
        g0m = g0 - arb_exact(self.eta_m)
        z0 = arb(self.zeta0)
        c63a = arb(63) / 16 * self.j2; c63b = arb(63) / 32 * self.j2
        R = []; Eg = []; ReA = []
        for i, x in enumerate(xi):
            rowR = []; rowE = []
            for l, c in enumerate(cc):
                v, Em, A = tan_point(x, c, g0, g0m, z0, c63a, c63b, Hg[i, l], Hg1[i, l], Hg2[i, l], Hg3[i, l],
                                     LHg[i, l], LHg1[i, l], H[i, l], H1[i, l], H2[i, l], LH[i, l], self.dscale)
                rowR.append(v); rowE.append(Em)
            R.append(rowR); Eg.append(rowE); ReA.append(A)
        return {'R': R, 'E': Eg, 'A': ReA}


# ---------------------------------------------------------------------------------------------- certificate
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


def coef_bounds(F, n):
    """F: (n+1)x(n+1) list of arb samples at Lobatto x Lobatto. Returns (c00, sum'|c|) as arb."""
    Cm = dct_matrix(n)
    Fm = arb_mat(F)
    Co = Cm * Fm * Cm.transpose()
    c00 = Co[0, 0]
    tot = arb(0)
    for a in range(n + 1):
        for b in range(n + 1):
            if a == 0 and b == 0:
                continue
            tot += abs(Co[a, b])
    return c00, tot


def ellipse_balls(rho, m, full):
    """Complex boxes covering the Bernstein ellipse boundary (upper half if not full)."""
    rho = arb(rho)
    a = (rho + 1 / rho) / 2; b = (rho - 1 / rho) / 2
    span = 2 if full else 1
    dphi = arb(span) / m        # in units of pi
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


RHO_COMBOS = [(2.0, 2.0), (1.6, 1.6), (3.0, 3.0), (2.0, 1.4), (1.4, 2.0), (1.3, 1.3), (1.15, 1.15)]


def certify_patch(fam, xa, xb, ca, cb, n=NSAMP):
    nodes = lobatto(n)
    g = fam.grid(xa, xb, ca, cb, nodes, nodes, cplx=False)
    c00, tot = coef_bounds(g['R'], n)
    cE00, totE = coef_bounds(g['E'], n)
    best = None
    for (r1, r2) in RHO_COMBOS:
        bx = ellipse_balls(r1, M1, full=False)
        by = ellipse_balls(r2, M2, full=True)
        gc = fam.grid(xa, xb, ca, cb, bx, by, cplx=True)
        ok = all(a.real > 0 for a in gc['A'])
        M = arb(0); ME = arb(0)
        for i in range(len(bx)):
            if not ok:
                break
            for l in range(len(by)):
                E = gc['E'][i][l]
                v = abs(gc['R'][i][l])
                if not (E.real > 0) or not v.is_finite():
                    ok = False; break
                M = arb(max(M.upper(), v.upper()))
                ME = arb(max(ME.upper(), abs(E).upper()))
        if not ok:
            continue
        tf = tail_factor(r1, r2, n, n)
        err = (M * tf).upper()
        errE = (ME * tf).upper()
        sup = (c00 + tot + err).upper(); inf = (c00 - tot - err).lower()
        Einf = (cE00 - totE - errE).lower()
        rec = dict(rho=[r1, r2], M=float(M.upper()), err=float(err), sup=sup, inf=inf, Einf=Einf,
                   ptot=float(tot.upper()))
        if best is None or float(err) < best['err']:
            best = rec
        if float(err) < 1e-2 * float(tot.upper()):
            break
    return best


# ---------------------------------------------------------------------------------------------- driver
def sol_path(name):
    return os.path.join(HERE, f'sol_{name}.json')


def save_sol(name, d):
    with open(sol_path(name), 'w') as f:
        json.dump(d, f, indent=1)


def load_sol(name):
    with open(sol_path(name)) as f:
        return json.load(f)


def cmd_solve(which):
    t0 = time.time()
    J21 = 2.0 ** 21
    cases = {
        'main': dict(j=J21, N=84, K=11, fac=1.0),
        'mainB': dict(j=J21, N=72, K=11, fac=1.0),
        'KN': dict(j=J21, N=72, K=11, fac=1.001),
        'K1': dict(j=2.0 ** -20, N=48, K=6, fac=1.0),
        'K1b': dict(j=2.0 ** -6, N=48, K=6, fac=1.0),
    }
    for name in which:
        cs = cases[name]
        g0, h, grid = RN.solve(cs['j'], cs['N'], cs['K'], CSfac=cs['fac'], j0=min(cs['j'], 2.0 ** -4))
        d = dict(j=fhex(cs['j']), fac_num=str(Fraction(cs['fac']).limit_denominator(10 ** 6)), N=cs['N'], K=cs['K'],
                 g0=fhex(g0), h=[fhex(v) for v in h], float_delta=2 * g0)
        save_sol(name, d)
        print(name, 'delta(float)', 2 * g0, f'{time.time()-t0:.1f}s', flush=True)


def CS_of(d):
    j = arb(unhex(d['j']))
    fac = Fraction(d['fac_num'])
    return arb(9) / 32 * j * j * arb(fac.numerator) / arb(fac.denominator)


def main_sol(d):
    N, K = d['N'], d['K']
    h = np.array([unhex(s) for s in d['h']]).reshape(N, K + 1)
    return Sol(unhex(d['g0']), h, even=True)


def patches_initial(even, nx):
    lo = 0 if even else -1
    pts = [Fraction(lo) + Fraction(1 - lo, nx) * i for i in range(nx + 1)]
    return [[str(pts[i]), str(pts[i + 1]), '0', '1'] for i in range(nx)]


def run_cert(state_file, fam_builder, even, nx_init, target, time_budget=CHUNK_SECONDS, onesided=False):
    """Adaptive patch certification with JSON checkpoint. Returns (done, state)."""
    if os.path.exists(state_file):
        with open(state_file) as f:
            st = json.load(f)
    else:
        st = dict(todo=patches_initial(even, nx_init), done=[], failed=[], prec=PREC, nsamp=NSAMP)
    fam = fam_builder()
    t0 = time.time()
    while st['todo']:
        if time.time() - t0 > time_budget:
            break
        p = st['todo'].pop(0)
        xa, xb, ca, cb = [Fraction(s) for s in p]
        t1 = time.time()
        rec = certify_patch(fam, xa, xb, ca, cb)
        good = rec is not None and float(rec['sup']) <= target and (onesided or float(rec['inf']) >= -target)
        if good:
            st['done'].append(dict(patch=p, rho=rec['rho'], M=rec['M'], err=rec['err'],
                                   sup=str(_ball_upper(rec['sup'] if isinstance(rec['sup'], arb) else arb(rec['sup']))),
                                   inf=str(_ball_lower(rec['inf'] if isinstance(rec['inf'], arb) else arb(rec['inf']))),
                                   Einf=str(_ball_lower(rec['Einf'] if isinstance(rec['Einf'], arb) else arb(rec['Einf']))),
                                   sec=round(time.time() - t1, 2)))
            print('ok  ', p, 'rho', rec['rho'], 'sup %.3e inf %.3e err %.2e M %.2e Einf %.4g  %.1fs' % (
                float(rec['sup']), float(rec['inf']), rec['err'], rec['M'], float(rec['Einf']), time.time() - t1), flush=True)
        else:
            w = xb - xa
            if w < Fraction(1, 4096):
                st['failed'].append(dict(patch=p, rec=None if rec is None else {k: str(v) for k, v in rec.items()}))
                print('FAIL', p, rec, flush=True)
                continue
            if w >= (cb - ca) / 8:
                m = (xa + xb) / 2
                st['todo'][0:0] = [[str(xa), str(m), p[2], p[3]], [str(m), str(xb), p[2], p[3]]]
            else:
                m = (ca + cb) / 2
                st['todo'][0:0] = [[p[0], p[1], str(ca), str(m)], [p[0], p[1], str(m), str(cb)]]
            print('split', p, None if rec is None else 'sup %.3e inf %.3e err %.2e' % (
                float(rec['sup']), float(rec['inf']), rec['err']), flush=True)
        with open(state_file, 'w') as f:
            json.dump(st, f, indent=0)
    with open(state_file, 'w') as f:
        json.dump(st, f, indent=0)
    return (not st['todo']), st


def cert_main(name, target):
    d = load_sol(name)
    sf = os.path.join(HERE, f'cert_{name}_p{PREC}_n{NSAMP}.json')
    done, st = run_cert(sf, lambda: MainF(main_sol(d), CS_of(d)), True, 16, target)
    return done


def cmd_certify(case, target):
    if case.startswith('tan'):
        return cert_tan(case, target)
    return cert_main(case, target)


# ---------------------------------------------------------------------------------------------- tangent
def solve_tangent(mainname, eta_m, Ne=60, Ke=12, even=True):
    """Float least-squares solve for (zeta0, H_e) in That = 0 against u_- = G_N - eta_m (collocation on [-1,1])."""
    from numpy.polynomial import chebyshev as Ch
    from numpy.polynomial import legendre as Lg
    d = load_sol(mainname)
    N, K = d['N'], d['K']
    g0 = unhex(d['g0']); h = np.array([unhex(s) for s in d['h']]).reshape(N, K + 1)
    j = unhex(d['j']); fac = float(Fraction(d['fac_num'])); j2 = j * j * fac
    nx = 2 * Ne + 20
    xi = np.cos(np.pi * (np.arange(nx) + 0.5) / nx)
    xi = np.concatenate([[1.0], xi, [-1.0]])
    cn, _ = Lg.leggauss(2 * (Ke + 4)); cc = np.concatenate([np.sort(cn[cn > 0]), [1.0]])
    XI, CC = np.meshgrid(xi, cc, indexing='ij'); x = XI.ravel(); c = CC.ravel()
    # G_N fields
    P = np.array([Lg.legval(c, np.eye(2 * K + 1)[2 * k]) for k in range(K + 1)]).T      # (pts, K+1)
    lam = np.array([-(2 * k) * (2 * k + 1) for k in range(K + 1)])
    der = []
    for m in range(4):
        cols = []
        for k in range(K + 1):
            cf = np.zeros(2 * N - 1); cf[0::2] = h[:, k]
            cols.append(Ch.chebval(x, Ch.chebder(cf, m) if m else cf))
        der.append(np.array(cols).T)
    Hg = [np.sum(der[m] * P, axis=1) for m in range(4)]
    LHg = np.sum(der[0] * P * lam, axis=1); LHg1 = np.sum(der[1] * P * lam, axis=1)
    c63a = 63 / 16 * j2; c63b = 63 / 32 * j2
    z = np.zeros_like(x)
    b, Em, A = tan_point(x, c, g0, g0 - eta_m, 0.0, c63a, c63b, Hg[0], Hg[1], Hg[2], Hg[3], LHg, LHg1, z, z, z, z)
    w = 1 - x ** 2
    potZ = c63a * (1 - c ** 2) * A ** 10 * w ** 12 / Em ** 8
    gz = RN.Grid(xi, cc, Ne, Ke, even=even)
    Mx = np.empty((x.size, 1 + gz.B.shape[1]))
    Mx[:, 0] = -0.5 - potZ
    Mx[:, 1:] = gz.Lin - (potZ * w ** 4)[:, None] * gz.B
    sol = np.linalg.lstsq(Mx, -b, rcond=None)[0]
    res = Mx @ sol + b
    return sol[0], sol[1:].reshape(Ne, Ke + 1), np.max(np.abs(res))


def cert_tan(case, target):
    info = json.load(open(os.path.join(HERE, f'tan_{case}.json')))
    d = load_sol(info['main'])
    sol = main_sol(d)
    es = Sol(0.0, np.array([unhex(s) for s in info['he']]).reshape(info['Ne'], info['Ke'] + 1), even=info['even'])
    j = arb(unhex(d['j'])); fac = Fraction(d['fac_num'])
    j2 = j * j * arb(fac.numerator) / arb(fac.denominator)
    eta_m = Fraction(info['eta_m'])
    sf = os.path.join(HERE, f'cert_{case}_p{PREC}_n{NSAMP}.json')
    done, st = run_cert(sf, lambda: TanF(sol, eta_m, es, unhex(info['zeta0']), j2, info['dscale']), False, 32, target, onesided=True)
    return done


def dyadic_up(q, bits=64):
    q = Fraction(q)
    n = -((-q.numerator * 2 ** bits) // q.denominator)
    return Fraction(n, 2 ** bits)


def cert_state(case):
    sf = os.path.join(HERE, f'cert_{case}_p{PREC}_n{NSAMP}.json')
    st = json.load(open(sf))
    assert not st['todo'] and not st['failed'], f'{case}: certification incomplete'
    return st


def cover_check(st, even):
    """Exact check that the accepted patches tile [lo,1] x [0,1] (area + pairwise disjoint interiors)."""
    lo = Fraction(0) if even else Fraction(-1)
    P = [[Fraction(s) for s in p['patch']] for p in st['done']]
    area = sum((b - a) * (d - c) for a, b, c, d in P)
    for a, b, c, d in P:
        assert lo <= a < b <= 1 and 0 <= c < d <= 1
    for i in range(len(P)):
        for k in range(i + 1, len(P)):
            a, b, c, d = P[i]; e, f, g, h = P[k]
            assert min(b, f) <= max(a, e) or min(d, h) <= max(c, g), 'overlap'
    assert area == (1 - lo), f'area {area}'
    return len(P)


def sandwich(case):
    st = cert_state(case)
    npatch = cover_check(st, True)
    d = load_sol(case)
    g0 = Fraction(unhex(d['g0']))
    sup = max(Fraction(p['sup']) for p in st['done'])
    inf = min(Fraction(p['inf']) for p in st['done'])
    Einf = min(Fraction(p['Einf']) for p in st['done'])
    if sup == 0 and inf == 0:                       # exact zero residual (toy): any positive eta works
        eta_p = eta_m = Fraction(1, 2 ** 60)
    else:
        eta_p = dyadic_up(max(2 * sup, Fraction(0)) * (1 + Fraction(1, 2 ** 20)) + Fraction(1, 2 ** 200))
        eta_m = dyadic_up(max(-2 * inf, Fraction(0)) * (1 + Fraction(1, 2 ** 20)) + Fraction(1, 2 ** 200))
    assert eta_p > 2 * sup and eta_m > -2 * inf      # strict super / sub conditions
    Eminus_lower = Einf - 2 * eta_m                  # E_- = E_N - 2 w^4 eta_m >= E_N - 2 eta_m (0 <= w <= 1)
    assert Eminus_lower > 0
    return dict(case=case, N=d['N'], K=d['K'], j=d['j'], fac=d['fac_num'], g0=g0, eta_p=eta_p, eta_m=eta_m,
                sup=sup, inf=inf, Einf=Einf, Eminus_lower=Eminus_lower, patches=npatch,
                delta_lo=2 * (g0 - eta_m), delta_hi=2 * (g0 + eta_p))


def cmd_tangent(case, mainname, Ne, Ke, even, dscale=1):
    sw = sandwich(mainname)
    if dscale == 0:
        zeta0, he, res = 0.0, np.zeros((Ne, Ke + 1)), 0.0
    else:
        zeta0, he, res = solve_tangent(mainname, float(sw['eta_m']), Ne, Ke, even=even)
    info = dict(main=mainname, eta_m=str(sw['eta_m']), Ne=Ne, Ke=Ke, even=even, dscale=dscale,
                zeta0=fhex(zeta0), he=[fhex(v) for v in np.asarray(he).ravel()], float_colloc_res=float(res))
    with open(os.path.join(HERE, f'tan_{case}.json'), 'w') as f:
        json.dump(info, f)
    print(case, 'zeta0', zeta0, 'float colloc res', res)


def tangent_result(case):
    info = json.load(open(os.path.join(HERE, f'tan_{case}.json')))
    st = cert_state(case)
    npatch = cover_check(st, False)
    sw = sandwich(info['main'])
    supT = max(Fraction(p['sup']) for p in st['done'])
    Eminus = min(Fraction(p['Einf']) for p in st['done'])
    assert Eminus > 0
    kappa = Fraction(0) if supT <= 0 else dyadic_up(2 * supT * (1 + Fraction(1, 2 ** 20)) + Fraction(1, 2 ** 200))
    zeta0 = Fraction(unhex(info['zeta0']))
    Zinf = (sw['g0'] if info['dscale'] else 0) + zeta0 + kappa     # Z_+(xi = 1) = lim r psi_+
    return dict(case=case, supT=supT, kappa=kappa, zeta0=zeta0, twoL_hi=1 + 2 * Zinf, patches=npatch,
                Eminus_lower_tan=Eminus)


def make_toy():
    d = dict(j=fhex(0.0), fac_num='1', N=4, K=1, g0=fhex(0.0), h=[fhex(0.0)] * 8, float_delta=0.0)
    save_sol('toy', d)


def fstr(q):
    q = Fraction(q)
    return f'{q.numerator}/{q.denominator}' if q.denominator != 1 else str(q.numerator)


def cmd_final():
    out = dict(seat='REP6', seal='0405395ed8', object='Bowen-York spinning puncture, mu1 = 1, j = 2^21 (J = 2^19)',
               method='log-tanh compactification xi = tanh(log(2r)/8) (a quarter-power map, not the sqrt map), '
                      'inversion-symmetric ansatz u = 2G/(2r+1), G = g0 + (1-xi^2)^4 H(xi, cos theta) with H even '
                      'Chebyshev x even Legendre; positivity by tensor Chebyshev-Lobatto interpolation of the exact '
                      'residual with an analytic-continuation (Bernstein polyellipse, Cauchy) error bound in ball arithmetic.',
               prec_bits=PREC, nsamp=NSAMP)
    sw = sandwich('main')
    out['a_delta_enclosure'] = dict(
        delta_lo=dec_round(sw['delta_lo'], 13, up=False), delta_hi=dec_round(sw['delta_hi'], 13, up=True),
        delta_lo_exact=fstr(sw['delta_lo']), delta_hi_exact=fstr(sw['delta_hi']),
        M_ADM_lo=dec_round(1 + sw['delta_lo'], 13, up=False), M_ADM_hi=dec_round(1 + sw['delta_hi'], 13, up=True),
        g0_candidate=fstr(sw['g0']), eta_plus=fstr(sw['eta_p']), eta_minus=fstr(sw['eta_m']),
        Rhat_sup_upper=dec_round(sw['sup'], 6, True), Rhat_inf_lower=dec_round(sw['inf'], 6, False),
        E_N_inf_lower=dec_round(sw['Einf'], 8, False), E_minus_inf_lower=dec_round(sw['Eminus_lower'], 8, False),
        patches=sw['patches'], N=sw['N'], K=sw['K'], strict_super=True, strict_sub=True, E_minus_positive=True)
    tr = tangent_result('tanmain')
    out['b_tangent'] = dict(
        twoL_hi_exact=fstr(tr['twoL_hi']), twoL_hi=dec_round(tr['twoL_hi'], 12, up=True),
        That_sup_upper=dec_round(tr['supT'], 6, True), kappa=fstr(tr['kappa']), zeta0=fstr(tr['zeta0']),
        Phi1_positive_lower_E=dec_round(tr['Eminus_lower_tan'], 8, False), patches=tr['patches'],
        TL_H_certified=True,
        class_K='structural: psi_+ = 2Z/(2r+1), Z = 2r/(2r+1) G_N - (1-xi^2)^4((1-xi^2)H_xi/8 - xi H) + zeta0 + kappa '
                '+ (1-xi^2)^4 H_e(xi,c); (1-xi^2)^4 = 512 r/(1+(2r)^(1/4))^8, so every non-constant piece is r times a '
                'smooth function of ((2r)^(1/4), cos theta): continuous at 0 with value 2 (zeta0 + kappa), |grad| = O(1), '
                '|Lap| = O(1/r), decays like 1/r at infinity (K1-K3 hold).')
    if 1 + sw['delta_lo'] - max(tr['twoL_hi'], 0) <= 0:
        raise SystemExit('denominator not positive')
    j = Fraction(2) ** 21
    chi = (j / 4) / (1 + sw['delta_lo'] - max(tr['twoL_hi'], Fraction(0))) ** 2
    out['c_chi_tail'] = dict(chi_tail=dec_round(chi, 15, up=True), chi_tail_exact=fstr(chi))
    # controls
    ctr = {}
    for name in ('K1', 'K1b'):
        s = sandwich(name)
        jj = Fraction(unhex(load_sol(name)['j']))
        lo = s['delta_lo'] / jj ** 2; hi = s['delta_hi'] / jj ** 2
        ctr[name] = dict(j=fstr(jj), delta_over_j2_lo=dec_round(lo, 15, False), delta_over_j2_hi=dec_round(hi, 15, True),
                         contains_1_over_40=bool(lo <= Fraction(1, 40) <= hi))
    stkt = json.load(open(os.path.join(HERE, f'cert_tanKT_p{PREC}_n{NSAMP}.json')))
    ctr['KT'] = dict(twoL='1', proof='psi = 0 (class K trivially): That(0) = -(63/32) j^2 (1-c^2) A^11 w^8 (1-xi)^4 / E_-^8 <= 0 '
                     'since every factor is >= 0 and E_- > 0 is certified in (a); i.e. Delta 0 = 0 <= V1/(2r). '
                     '2L = 1 + 2 lim r psi = 1.',
                     machinery_partial=dict(patches_done=len(stkt['done']), patches_left=len(stkt['todo']),
                                            max_sup_bound_on_done=dec_round(max(Fraction(p['sup']) for p in stkt['done']), 6, True),
                                            note='interpolation run stopped by the seat: the bound c00+sum|c| is loose where That(0) '
                                                 'touches 0 at c = 1, so patches kept splitting; not needed for the analytic proof'))
    s = sandwich('KN')
    ctr['KN'] = dict(source_factor='1001/1000', delta_lo=dec_round(s['delta_lo'], 12, False),
                     delta_hi=dec_round(s['delta_hi'], 12, True))
    s = sandwich('toy'); tt = tangent_result('tantoy')
    ctr['Toy'] = dict(delta_lo=fstr(s['delta_lo']), delta_hi=fstr(s['delta_hi']), contains_0=bool(s['delta_lo'] <= 0 <= s['delta_hi']),
                      twoL_hi=fstr(tt['twoL_hi']))
    out['d_controls'] = ctr
    conv = {}
    s = sandwich('mainB')
    conv['mainB_N72_K11'] = dict(delta_lo=dec_round(s['delta_lo'], 13, False), delta_hi=dec_round(s['delta_hi'], 13, True))
    tb = tangent_result('tanB')
    conv['tanB_mainB_Ne100'] = dict(twoL_hi=dec_round(tb['twoL_hi'], 12, True))
    out['convergence'] = conv
    with open(os.path.join(HERE, 'result_rep6.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'solve':
        cmd_solve(sys.argv[2:])
    elif cmd == 'toy':
        make_toy()
    elif cmd == 'certify':
        done = cmd_certify(sys.argv[2], float(sys.argv[3]))
        print('DONE' if done else 'RESUME')
        sys.exit(0 if done else 3)
    elif cmd == 'tangent':
        cmd_tangent(sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), sys.argv[6] == 'even',
                    int(sys.argv[7]) if len(sys.argv) > 7 else 1)
    elif cmd == 'sandwich':
        r = sandwich(sys.argv[2])
        print({k: (float(v) if isinstance(v, Fraction) else v) for k, v in r.items()})
    elif cmd == 'final':
        cmd_final()
    elif cmd == 'tanres':
        r = tangent_result(sys.argv[2])
        print({k: (float(v) if isinstance(v, Fraction) else v) for k, v in r.items()})
