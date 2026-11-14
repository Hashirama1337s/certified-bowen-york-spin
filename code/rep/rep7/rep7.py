"""REP7 blind replication: rigorous enclosure of the ADM mass M_0 of the extreme (mu = 0)
Bowen-York spinning data with J = 1, plus RN and KN controls.

Usage (always under memcap, from <LAB>/BHSLICE/rep7):
  py -3.11 rep7.py run BY|RN|KN      (exit 0 = case finished, 3 = rerun to resume, other = error)
  py -3.11 rep7.py assemble          (writes result_rep7.json)

Method summary (details in RESULT.md / final report):
  * numerical solution: rep7_solve.py, formulation psi^2 = 1 + F/r, F polynomial on a far element in
    x = 1/r and Lobatto elements in t = -log r; eta = cos 2 theta.
  * the approximate solution u_a = sqrt(1 + F_a/r) - 1 is the exact interpolant of the float nodal values
    (continuous across elements by construction); far element corrected so F_a(x=0, eta) is exactly constant.
  * residual R = Delta u_a + S (1+u_a)^-p bounded rigorously (python-flint ball arithmetic) on cells by
    first-kind Chebyshev interpolation + Bernstein-polyellipse remainder (max modulus on the distinguished
    boundary), interface kinks bounded as surface measures, throat tail t > T+ bounded analytically.
  * rho(r) >= sup_theta |R| (as a radial measure); h = Newton potential of rho; v = u_a - h is a sub-solution,
    w = u_a + h a super-solution; comparison gives  M_a - 2 delta <= M_0 <= M_a + 2 delta,
    delta = int_0^inf rho(s) s^2 ds.
"""
import sys, os, json, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from flint import arb, acb, arb_mat, acb_mat, ctx
import scipy.sparse as sp
import scipy.sparse.linalg as spla

# ===================== numerical solver (float64) =====================
# Formulation: psi = sqrt(1 + F/r), g = e^{-t} + F, t = -log r, eta = cos 2theta;
# E(F) = e^{-t}(F_tt+F_t+LF) + F(F_tt+LF) - (F_t^2 + 4(1-eta^2)F_eta^2)/2 - F^2/2 + 2K g^{(3-p)/2} = 0



def bary_weights(x):
    n = len(x)
    w = np.ones(n)
    for j in range(n):
        for k in range(n):
            if k != j:
                w[j] *= (x[j] - x[k])
    return 1.0 / w


def diff_matrix(x):
    n = len(x)
    w = bary_weights(x)
    D = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                D[i, j] = (w[j] / w[i]) / (x[i] - x[j])
        D[i, i] = -np.sum(D[i, :])
    return D


def lobatto(N):
    # increasing order
    return np.cos(np.pi * np.arange(N, -1, -1) / N)


def radau_right(N):
    # Chebyshev-Gauss-Radau nodes including s=+1, excluding s=-1; increasing order
    s = np.cos(2 * np.pi * np.arange(N + 1) / (2 * N + 1))
    return np.sort(s)


class Problem:
    def __init__(self, p=7, kamp=9.0 / 4.0, Nf=20, x1=0.5, tb=None, Nt=20, Ne=20):
        self.p = p
        self.kamp = kamp
        self.q = (3 - p) / 2.0
        self.a = (p - 1) / 2.0
        self.Nf, self.x1, self.Nt, self.Ne = Nf, x1, Nt, Ne
        if tb is None:
            tb = [np.log(x1), 0.0, 0.75, 1.5, 2.5, 3.5, 5.0, 7.0, 9.5, 13.0, 17.0, 21.0, 26.0]
        assert abs(tb[0] - np.log(x1)) < 1e-14
        self.tb = list(tb)
        self.eta = lobatto(Ne)
        self.De = diff_matrix(self.eta)
        self.nE = Ne + 1
        # far element
        self.sf = radau_right(Nf)
        self.Dsf = diff_matrix(self.sf)
        self.xf = x1 * (1 + self.sf) / 2
        # t elements
        self.st = lobatto(Nt)
        self.Dst = diff_matrix(self.st)
        self.elems = []
        off = 0
        self.elems.append(dict(kind='x', off=0, n=Nf + 1))
        off = Nf
        for k in range(len(tb) - 1):
            c = (tb[k] + tb[k + 1]) / 2
            h = (tb[k + 1] - tb[k]) / 2
            self.elems.append(dict(kind='t', off=off, n=Nt + 1, c=c, h=h, t=c + h * self.st))
            off += Nt
        self.nR = off + 1
        self.n = self.nR * self.nE

    def K(self):
        if self.p == 7:
            return self.kamp * (1 - self.eta) / 2
        return self.kamp * np.ones_like(self.eta)

    def elem_eval(self, el, Fe, want_jac=True):
        """residual of element at all its nodes, plus local Jacobian (dense)."""
        nE = self.nE
        eta = self.eta
        K = self.K()[None, :]
        Kv = np.tile(self.K(), el['n'])
        I_e = np.eye(nE)
        De = self.De
        Fe_eta = Fe @ De.T
        Fe_ee = Fe_eta @ De.T
        LF = 4 * (1 - eta ** 2) * Fe_ee - 2 * (1 + 3 * eta) * Fe_eta
        ev = np.tile(eta, el['n'])
        if el['kind'] == 't':
            Ds = self.Dst / el['h']
            Ft = Ds @ Fe
            Ftt = Ds @ Ft
            et = np.exp(-el['t'])[:, None]
            g = et + Fe
            res = (et * (Ftt + Ft + LF) + Fe * (Ftt + LF) - 0.5 * (Ft ** 2 + 4 * (1 - eta ** 2) * Fe_eta ** 2)
                   - 0.5 * Fe ** 2 + 2 * K * g ** self.q)
            if not want_jac:
                return res, None
            I_r = np.eye(el['n'])
            Dt = np.kron(Ds, I_e)
            Dtt = Dt @ Dt
            Dq = np.kron(I_r, De)
            Dqq = Dq @ Dq
            Lop = (4 * (1 - ev ** 2))[:, None] * Dqq - (2 * (1 + 3 * ev))[:, None] * Dq
            etv = np.repeat(np.exp(-el['t']), nE)
            f = Fe.ravel(); ft = Ft.ravel(); ftt = Ftt.ravel(); lf = LF.ravel(); fq = Fe_eta.ravel()
            gv = etv + f
            J = (etv[:, None] * (Dtt + Dt + Lop) + np.diag(ftt + lf) + f[:, None] * (Dtt + Lop)
                 - ft[:, None] * Dt - (4 * (1 - ev ** 2) * fq)[:, None] * Dq - np.diag(f)
                 + np.diag(2 * Kv * self.q * gv ** (self.q - 1)))
            return res, J
        else:
            x = self.xf[:, None]
            Ds = self.Dsf * (2.0 / self.x1)
            Fx = Ds @ Fe
            Fxx = Ds @ Fx
            Ft = x * Fx
            Ftt = x * x * Fxx + x * Fx
            res = ((Ftt + Ft + LF) + x * (Fe * (Ftt + LF) - 0.5 * (Ft ** 2 + 4 * (1 - eta ** 2) * Fe_eta ** 2)
                                         - 0.5 * Fe ** 2) + 2 * K * x ** self.a * (1 + x * Fe) ** self.q)
            if not want_jac:
                return res, None
            I_r = np.eye(el['n'])
            xv = np.repeat(self.xf, nE)
            Dx = np.kron(Ds, I_e)
            Dxx = Dx @ Dx
            Dt = xv[:, None] * Dx
            Dtt = (xv ** 2)[:, None] * Dxx + xv[:, None] * Dx
            Dq = np.kron(I_r, De)
            Dqq = Dq @ Dq
            Lop = (4 * (1 - ev ** 2))[:, None] * Dqq - (2 * (1 + 3 * ev))[:, None] * Dq
            f = Fe.ravel(); ft = Ft.ravel(); ftt = Ftt.ravel(); lf = LF.ravel(); fq = Fe_eta.ravel()
            inner = (np.diag(ftt + lf) + f[:, None] * (Dtt + Lop) - ft[:, None] * Dt
                     - (4 * (1 - ev ** 2) * fq)[:, None] * Dq - np.diag(f))
            J = (Dtt + Dt + Lop) + xv[:, None] * inner + np.diag(
                2 * Kv * self.q * xv ** (self.a + 1) * (1 + xv * f) ** (self.q - 1))
            return res, J

    def dt_rows(self, el, Fe, idx):
        """t-derivative at local radial node idx (all eta), and its linear operator."""
        nE = self.nE
        if el['kind'] == 't':
            row = self.Dst[idx] / el['h']
        else:
            row = self.Dsf[idx] * (2.0 / self.x1) * self.xf[idx]
        val = row @ Fe
        return val, row

    def residual_jac(self, F, want_jac=True):
        nE, nR = self.nE, self.nR
        F2 = F.reshape(nR, nE)
        res = np.zeros(self.n)
        rows, cols, vals = [], [], []
        for ie, el in enumerate(self.elems):
            off, n = el['off'], el['n']
            Fe = F2[off:off + n]
            r, J = self.elem_eval(el, Fe, want_jac)
            if el['kind'] == 'x':
                keep = np.arange(0, n - 1)
            else:
                keep = np.arange(1, n - 1)
            for i in keep:
                gi = off + i
                res[gi * nE:(gi + 1) * nE] = r[i]
            if want_jac:
                ridx = ((off + keep)[:, None] * nE + np.arange(nE)[None, :]).ravel()
                lrows = (keep[:, None] * nE + np.arange(nE)[None, :]).ravel()
                blk = J[lrows, :]
                cidx = off * nE + np.arange(n * nE)
                rows.append(np.repeat(ridx, len(cidx)))
                cols.append(np.tile(cidx, len(ridx)))
                vals.append(blk.ravel())
        # interfaces: derivative continuity
        for ie in range(len(self.elems) - 1):
            L, R = self.elems[ie], self.elems[ie + 1]
            gi = R['off']
            vL, rowL = self.dt_rows(L, F2[L['off']:L['off'] + L['n']], L['n'] - 1)
            vR, rowR = self.dt_rows(R, F2[R['off']:R['off'] + R['n']], 0)
            res[gi * nE:(gi + 1) * nE] = vL - vR
            if want_jac:
                k = np.arange(nE)
                for j in range(L['n']):
                    rows.append(gi * nE + k); cols.append((L['off'] + j) * nE + k); vals.append(np.full(nE, rowL[j]))
                for j in range(R['n']):
                    rows.append(gi * nE + k); cols.append((R['off'] + j) * nE + k); vals.append(np.full(nE, -rowR[j]))
        # throat end Neumann
        L = self.elems[-1]
        gi = L['off'] + L['n'] - 1
        vL, rowL = self.dt_rows(L, F2[L['off']:L['off'] + L['n']], L['n'] - 1)
        res[gi * nE:(gi + 1) * nE] = vL
        if not want_jac:
            return res, None
        k = np.arange(nE)
        for j in range(L['n']):
            rows.append(gi * nE + k); cols.append((L['off'] + j) * nE + k); vals.append(np.full(nE, rowL[j]))
        J = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.n, self.n))
        return res, J

    def solve(self, F0=None, maxit=40, verbose=True):
        F = np.ones(self.n) if F0 is None else F0.copy()
        best = (np.inf, F.copy())
        stall = 0
        for it in range(maxit):
            res, J = self.residual_jac(F)
            nr = np.max(np.abs(res))
            if verbose:
                print(f"  newton {it}: |res|_inf = {nr:.3e}", flush=True)
            stall = 0 if nr < 0.5 * best[0] else stall + 1
            if nr < best[0]:
                best = (nr, F.copy())
            if best[0] < 1e-9 and stall >= 2:
                break
            dF = spla.spsolve(J.tocsc(), -res)
            lam = 1.0
            if nr > 1e-6:
                while lam > 1e-4:
                    Fn = F + lam * dF
                    if np.all(Fn > 0):
                        rn, _ = self.residual_jac(Fn, want_jac=False)
                        if np.max(np.abs(rn)) < nr:
                            break
                    lam /= 2
            F = F + lam * dF
        self.F = best[1]
        self.best_res = best[0]
        return best[1]

    def mass_far(self, F):
        """F at x = 0 for each eta via barycentric interpolation on the far element."""
        F2 = F.reshape(self.nR, self.nE)
        Fe = F2[0:self.Nf + 1]
        s = self.sf
        w = bary_weights(s)
        c = w / (-1.0 - s)
        return (c @ Fe) / np.sum(c)

# ===================== rigorous verification =====================

PREC = 128
CASES = {
    'BY': dict(p=7, kamp_num=9, kamp_den=4),
    'RN': dict(p=3, kamp_num=1, kamp_den=4),
    'KN': dict(p=7, kamp_num=9 * 1001, kamp_den=4 * 1000),
    # convergence control: same object as BY, independent discretisation / verification parameters
    'BY2': dict(p=7, kamp_num=9, kamp_den=4,
                disc=dict(Nf=16, x1=0.4, Nt=18, Ne=18,
                          tb=['log', 0.3, 1.0, 2.0, 3.0, 4.5, 6.5, 9.0, 12.0, 16.0, 20.0, 25.0, 30.0]),
                ver=dict(ns=5, nq=3, n=28, rho1=3.5, rho2=3.5, arcs1=40, arcs2=80, nthroat=20, n1d=28,
                         rho1d=3.5, arcs1d=80)),
}
# BYP: can-fail control of the certificate: the converged BY solution is deliberately spoiled by
# +1e-7 (1 + eta/2) at every node (shifts M_a by ~1e-7); the window must still contain the BY windows.
CASES['BYP'] = dict(p=7, kamp_num=9, kamp_den=4)
PERTURB = {'BYP': 1e-7}
# RN: start Newton away from the exact solution F = 1 so the solver path is exercised
INIT = {'RN': (1.3, 0.2)}
# discretisation (fixed before the target run; chosen by the convergence probes)
DISC = dict(Nf=20, x1=0.5, Nt=20, Ne=20,
            tb=['log', 0.0, 0.75, 1.5, 2.5, 3.5, 5.0, 7.0, 9.5, 13.0, 17.0, 21.0, 26.0])
VER = dict(ns=4, nq=4, n=24, rho1=4.0, rho2=4.0, arcs1=32, arcs2=64, nthroat=16, n1d=24, rho1d=4.0, arcs1d=64)
TIME_BUDGET = 480.0


def ckpt_path(case):
    return os.path.join(HERE, f'ckpt_{case}.json')


def load_ckpt(case):
    p = ckpt_path(case)
    if os.path.exists(p):
        with open(p) as f:
            return json.load(f)
    return {}


def save_ckpt(case, d):
    p = ckpt_path(case)
    with open(p + '.tmp', 'w') as f:
        json.dump(d, f, indent=1)
    os.replace(p + '.tmp', p)


# ----------------------------------------------------------------------------------------------
# exact Chebyshev helpers (ball arithmetic)
# ----------------------------------------------------------------------------------------------

def a(x):
    return arb(x)


def lobatto_theta(N):
    # increasing nodes s_j = cos(pi (N-j)/N)  -> theta fraction (N-j)/N
    return [(N - j, N) for j in range(N + 1)]


def radau_theta(N):
    # increasing nodes s_i = cos(2 pi (N-i)/(2N+1))
    return [(2 * (N - i), 2 * N + 1) for i in range(N + 1)]


def vander_from_theta(th, nb):
    """V[j,k] = T_k(cos(pi p/q)) = cos(pi k p / q), exact rational arguments."""
    V = arb_mat(len(th), nb)
    for j, (p, q) in enumerate(th):
        for k in range(nb):
            V[j, k] = (arb(k * p) / q).cos_pi()
    return V


def identity(n):
    I = arb_mat(n, n)
    for i in range(n):
        I[i, i] = 1
    return I


def cheb_T_real(x, nb):
    """T_0..T_{nb-1} at arb x (tiny radius) by recurrence."""
    T = [arb(1), x]
    for k in range(2, nb):
        T.append(2 * x * T[-1] - T[-2])
    return T[:nb]


def deriv_matrix(nb):
    """coefficient differentiation: b = D a for p = sum a_k T_k, p' = sum b_k T_k."""
    D = arb_mat(nb, nb)
    for k in range(nb):
        ck = 2 if k == 0 else 1
        for j in range(k + 1, nb):
            if (j - k) % 2 == 1:
                D[k, j] = arb(2 * j) / ck
    return D


def first_kind_basis(n, nb):
    """B[i,k] = T_k(cos((2i+1) pi / (2n)))"""
    B = arb_mat(n, nb)
    for i in range(n):
        for k in range(nb):
            B[i, k] = (arb(k * (2 * i + 1)) / (2 * n)).cos_pi()
    return B


def dct_first_kind(n):
    """coefficients c = Cn f of the first-kind interpolant (degree n-1)."""
    C = arb_mat(n, n)
    for k in range(n):
        for i in range(n):
            v = (arb(k * (2 * i + 1)) / (2 * n)).cos_pi() * 2 / n
            if k == 0:
                v = v / 2
            C[k, i] = v
    return C


def local_reexp_matrix(nb, m, h):
    """A with  T_k(m + h sigma) = sum_k' A[k',k] T_k'(sigma)  (exact polynomial identity, enclosed)."""
    th = lobatto_theta(nb - 1)
    Vl = vander_from_theta(th, nb)
    W = arb_mat(nb, nb)
    for j, (p, q) in enumerate(th):
        sig = (arb(p) / q).cos_pi()
        T = cheb_T_real(m + h * sig, nb)
        for k in range(nb):
            W[j, k] = T[k]
    return Vl.solve(W)


def mat_abs_sum(Mx):
    tot = arb(0)
    for i in range(Mx.nrows()):
        for j in range(Mx.ncols()):
            tot += abs(Mx[i, j])
    return arb(tot.upper())


def arc_balls(rho, nar, half):
    """balls (acb) covering |w| = rho; half=True covers only arg in [0, pi]."""
    span = np.pi if half else 2 * np.pi
    out = []
    dphi = span / nar
    R = rho * dphi / 2 * (1 + 1e-9) + 1e-14
    for k in range(nar):
        ph = (k + 0.5) * dphi
        c = complex(rho * np.cos(ph), rho * np.sin(ph))
        out.append(acb(arb(c.real, R), arb(c.imag, R)))
    return out


def cheb_T_arc(W, nb):
    """enclosures of T_k(J(w)), J(w) = (w + 1/w)/2, for all w in the ball W."""
    Wi = 1 / W
    out = [acb(1)]
    pw, pwi = acb(1), acb(1)
    for k in range(1, nb):
        pw = pw * W
        pwi = pwi * Wi
        out.append((pw + pwi) / 2)
    return out, (W + Wi) / 2


# ----------------------------------------------------------------------------------------------
# residual functions (generic over arb / acb)
# ----------------------------------------------------------------------------------------------

def phi_t(case, F, Ft, Ftt, Fq, Fqq, t, eta):
    """Phi = E / (2 g^{3/2}),  R = r^{-5/2} Phi.   returns (Phi, g)."""
    p = CASES[case]['p']
    kamp = arb(CASES[case]['kamp_num']) / CASES[case]['kamp_den']
    LF = 4 * (1 - eta * eta) * Fqq - 2 * (1 + 3 * eta) * Fq
    et = (-t).exp()
    g = et + F
    if p == 7:
        K = kamp * (1 - eta) / 2
        src = 2 * K / (g * g)
    else:
        src = 2 * kamp
    E = et * (Ftt + Ft + LF) + F * (Ftt + LF) - (Ft * Ft + 4 * (1 - eta * eta) * Fq * Fq) / 2 - F * F / 2 + src
    return E / (2 * g * g.sqrt()), g


def phi_x(case, F, Fx, Fxx, Fq, Fqq, x, eta):
    """Q = Ehat / (2 x (1+xF)^{3/2}),  R = x^4 Q = r^{-4} Q.  returns (Q, 1+xF)."""
    p = CASES[case]['p']
    kamp = arb(CASES[case]['kamp_num']) / CASES[case]['kamp_den']
    LF = 4 * (1 - eta * eta) * Fqq - 2 * (1 + 3 * eta) * Fq
    Ft = x * Fx
    Ftt = x * x * Fxx + x * Fx
    G = 1 + x * F
    if p == 7:
        K = kamp * (1 - eta) / 2
        src = 2 * K * x * x * x / (G * G)
    else:
        src = 2 * kamp * x
    Eh = (Ftt + Ft + LF) + x * (F * (Ftt + LF) - (Ft * Ft + 4 * (1 - eta * eta) * Fq * Fq) / 2 - F * F / 2) + src
    return Eh / (2 * x * G * G.sqrt()), G


# ----------------------------------------------------------------------------------------------
# stage 1: numerical solve
# ----------------------------------------------------------------------------------------------

def make_problem(case, disc):
    c = CASES[case]
    tb = [np.log(disc['x1']) if v == 'log' else v for v in disc['tb']]
    return Problem(p=c['p'], kamp=c['kamp_num'] / c['kamp_den'], Nf=disc['Nf'], x1=disc['x1'], tb=tb,
                   Nt=disc['Nt'], Ne=disc['Ne'])


def stage_solve(case, ck):
    t0 = time.time()
    P = make_problem(case, DISC)
    F0 = None
    if case in INIT:
        c0, c1 = INIT[case]
        F0 = np.tile(c0 + c1 * P.eta, P.nR)
    F = P.solve(F0=F0, verbose=True)
    if case in PERTURB:
        F = F + PERTURB[case] * np.tile(1 + 0.5 * P.eta, P.nR)
    m = P.mass_far(F)
    ck['solve'] = dict(F_hex=[float(v).hex() for v in F], M_float=[float(m.min()), float(m.max())],
                       time=time.time() - t0)
    return ck


# ----------------------------------------------------------------------------------------------
# stage 2: rigorous coefficients and cell bounds
# ----------------------------------------------------------------------------------------------

class Rig:
    def __init__(self, case, ck):
        ctx.prec = PREC
        self.case = case
        self.P = make_problem(case, DISC)
        P = self.P
        F = np.array([float.fromhex(h) for h in ck['solve']['F_hex']])
        F2 = F.reshape(P.nR, P.nE)
        nE = P.nE
        self.nbq = nE
        Vq = vander_from_theta(lobatto_theta(P.Ne), nE)
        Vqinv = Vq.solve(identity(nE))
        self.Dq = deriv_matrix(nE)
        self.coef = []
        for el in P.elems:
            n = el['n']
            th = radau_theta(P.Nf) if el['kind'] == 'x' else lobatto_theta(P.Nt)
            Vs = vander_from_theta(th, n)
            Fv = arb_mat(n, nE)
            for i in range(n):
                for k in range(nE):
                    Fv[i, k] = arb(float(F2[el['off'] + i, k]))   # exact dyadic
            C = Vs.solve(Fv) * Vqinv.transpose()
            self.coef.append(C)
        # far-element correction: make F(x=0, eta) exactly constant (q(s) = ((1-s)/2)^2)
        C = self.coef[0]
        n = C.nrows()
        f0 = []
        for l in range(nE):
            v = arb(0)
            for k in range(n):
                v += C[k, l] * (1 if k % 2 == 0 else -1)
            f0.append(v)
        qc = [arb(3) / 8, arb(-1) / 2, arb(1) / 8]
        for l in range(1, nE):
            for k in range(3):
                C[k, l] = C[k, l] - f0[l] * qc[k]
        self.coef[0] = C
        self.Ma = f0[0]
        self.f0_correction = max(abs(float(f0[l].mid())) for l in range(1, nE))
        self.Ds = {}
        self.cache = {}
        # exact element maps t = c + h s (interfaces shared exactly; first boundary = log(x1) exactly)
        tbx = [arb(P.x1).log()] + [arb(float(v)) for v in P.tb[1:]]
        self.ec, self.eh = [None], [None]
        for k in range(len(tbx) - 1):
            self.ec.append((tbx[k] + tbx[k + 1]) / 2)
            self.eh.append((tbx[k + 1] - tbx[k]) / 2)
        self.tbx = tbx

    def dmat(self, nb):
        if nb not in self.Ds:
            self.Ds[nb] = deriv_matrix(nb)
        return self.Ds[nb]

    def cell(self, ie, isub, iq):
        """bound sup |Phi| (t-element) or sup |Q| (far element) on one cell."""
        P = self.P
        el = P.elems[ie]
        C = self.coef[ie]
        ns, nq, n = VER['ns'], VER['nq'], VER['n']
        nb, nbq = C.nrows(), C.ncols()
        sa = arb(-1) + arb(2 * isub) / ns
        hs = arb(1) / ns
        ms = sa + hs
        qa = arb(-1) + arb(2 * iq) / nq
        hq = arb(1) / nq
        mq = qa + hq
        As = local_reexp_matrix(nb, ms, hs)
        Aq = local_reexp_matrix(nbq, mq, hq)
        Cl = As * C * Aq.transpose()
        Ds, Dq = self.dmat(nb), self.dmat(nbq)
        Cs = Ds * Cl
        Css = Ds * Cs
        Cq = Cl * Dq.transpose()
        Cqq = Cq * Dq.transpose()
        if el['kind'] == 't':
            he = self.eh[ie]
            ce = self.ec[ie]
            fs = 1 / (hs * he)
        else:
            x1 = arb(P.x1)
            fs = (2 / x1) / hs
        fq = 1 / hq
        mats = [Cl, Cs * fs, Css * (fs * fs), Cq * fq, Cqq * (fq * fq)]
        # --- real first-kind nodes
        if ('fk', n, nb) not in self.cache:
            self.cache[('fk', n, nb)] = first_kind_basis(n, nb)
        if ('fk', n, nbq) not in self.cache:
            self.cache[('fk', n, nbq)] = first_kind_basis(n, nbq)
        if ('dct', n) not in self.cache:
            self.cache[('dct', n)] = dct_first_kind(n)
        Bs = self.cache[('fk', n, nb)]
        Bq = self.cache[('fk', n, nbq)]
        vals = [Bs * Mx * Bq.transpose() for Mx in mats]
        Phi = arb_mat(n, n)
        for i in range(n):
            sig = (arb(2 * i + 1) / (2 * n)).cos_pi()
            sloc = ms + hs * sig
            for j in range(n):
                tau = (arb(2 * j + 1) / (2 * n)).cos_pi()
                eta = mq + hq * tau
                args = [vals[m][i, j] for m in range(5)]
                if el['kind'] == 't':
                    ph, _ = phi_t(self.case, *args, ce + he * sloc, eta)
                else:
                    ph, _ = phi_x(self.case, *args, x1 * (1 + sloc) / 2, eta)
                Phi[i, j] = ph
        Cn = self.cache[('dct', n)]
        cc = Cn * Phi * Cn.transpose()
        S_nodes = mat_abs_sum(cc)
        maxnode = max(abs(float(Phi[i, j].mid())) for i in range(n) for j in range(n))
        # --- complex distinguished boundary
        rho1, rho2 = VER['rho1'], VER['rho2']
        W1 = arc_balls(rho1, VER['arcs1'], True)
        W2 = arc_balls(rho2, VER['arcs2'], False)
        A1 = acb_mat(len(W1), nb)
        sig1 = []
        for a_, W in enumerate(W1):
            T, sg = cheb_T_arc(W, nb)
            sig1.append(sg)
            for k in range(nb):
                A1[a_, k] = T[k]
        A2 = acb_mat(len(W2), nbq)
        tau2 = []
        for b_, W in enumerate(W2):
            T, sg = cheb_T_arc(W, nbq)
            tau2.append(sg)
            for k in range(nbq):
                A2[b_, k] = T[k]
        cvals = [A1 * acb_mat(Mx) * A2.transpose() for Mx in mats]
        Mmax = arb(0)
        glo = None
        ok = True
        for a_ in range(len(W1)):
            sloc = ms + hs * sig1[a_]
            if el['kind'] == 't':
                tt = ce + he * sloc
            else:
                xx = x1 * (1 + sloc) / 2
            for b_ in range(len(W2)):
                eta = mq + hq * tau2[b_]
                args = [cvals[m][a_, b_] for m in range(5)]
                if el['kind'] == 't':
                    ph, g = phi_t(self.case, *args, tt, eta)
                else:
                    ph, g = phi_x(self.case, *args, xx, eta)
                rl = g.real.lower()
                if not (rl > 0):
                    ok = False
                glo = rl if glo is None else min(glo, rl)
                if ph.is_finite():
                    Mmax = max(Mmax, abs(ph).upper())
                else:
                    ok = False
        r1, r2 = arb(rho1), arb(rho2)
        S1, S2 = r1 / (r1 - 1), r2 / (r2 - 1)
        t1 = r1 ** (1 - n) / (r1 - 1)
        t2 = r2 ** (1 - n) / (r2 - 1)
        eps = 8 * Mmax * (t1 * S2 + S1 * t2)
        B = arb((S_nodes + eps).upper())
        return dict(ie=ie, isub=isub, iq=iq, ok=bool(ok), B=B.str(radius=True), S_nodes=S_nodes.str(radius=True),
                    M=Mmax.str(radius=True), eps=eps.str(radius=True), glo=glo.str(radius=True),
                    maxnode=maxnode)


def stage_cells(case, ck, t_start):
    R = Rig(case, ck)
    done = ck.setdefault('cells', {})
    P = R.P
    todo = [(ie, i, j) for ie in range(len(P.elems)) for i in range(VER['ns']) for j in range(VER['nq'])]
    for (ie, i, j) in todo:
        key = f'{ie}_{i}_{j}'
        if key in done:
            continue
        tc = time.time()
        res = R.cell(ie, i, j)
        res['time'] = time.time() - tc
        done[key] = res
        print(f'cell {key}: B={float(arb(res["B"]).mid()):.3e} nodes={float(arb(res["S_nodes"]).mid()):.3e} '
              f'M={float(arb(res["M"]).mid()):.3e} glo={float(arb(res["glo"]).mid()):.3f} ok={res["ok"]} '
              f'({res["time"]:.1f}s)', flush=True)
        if time.time() - t_start > TIME_BUDGET:
            return ck, False
    return ck, True


# ----------------------------------------------------------------------------------------------
# stage 3: interfaces, throat tail, assembly
# ----------------------------------------------------------------------------------------------

def eta_poly_bounds(R, coeffs, sub=16):
    """for an eta-polynomial (list of arb Chebyshev coeffs): lower bound of min and upper bound of sup |.|,
    by local re-expansion on `sub` sub-intervals (c0 -/+ sum |c_l|)."""
    nb = len(coeffs)
    cm = arb_mat(nb, 1)
    for k in range(nb):
        cm[k, 0] = coeffs[k]
    lo, hi = None, arb(0)
    for i in range(sub):
        h = arb(1) / sub
        m = arb(-1) + arb(2 * i + 1) / sub
        A = local_reexp_matrix(nb, m, h)
        cl = A * cm
        rest = arb(0)
        for k in range(1, nb):
            rest += abs(cl[k, 0])
        l_ = cl[0, 0] - rest
        u_ = abs(cl[0, 0]) + rest
        lo = l_.lower() if lo is None else min(lo, l_.lower())
        hi = max(hi, u_.upper())
    return arb(lo), arb(hi)


def throat_bounds(R):
    """t >= T+: F = F_T(eta) (constant extension).  |Phi| <= A + B e^{-t}."""
    P = R.P
    C = R.coef[-1]
    nb, nbq = C.nrows(), C.ncols()
    FT = [sum((C[k, l] for k in range(nb)), arb(0)) for l in range(nbq)]
    Dq = R.dmat(nbq)
    FTm = arb_mat(nbq, 1)
    for l in range(nbq):
        FTm[l, 0] = FT[l]
    case = R.case
    p = CASES[case]['p']
    kamp = arb(CASES[case]['kamp_num']) / CASES[case]['kamp_den']
    n, rho, nar, sub = VER['n1d'], VER['rho1d'], VER['arcs1d'], VER['nthroat']
    Cn = dct_first_kind(n)
    Bq = first_kind_basis(n, nbq)
    Aall, Fmin, Bmax = arb(0), None, arb(0)
    okall = True
    for isub in range(sub):
        h = arb(1) / sub
        m = arb(-1) + arb(2 * isub + 1) / sub
        A = local_reexp_matrix(nbq, m, h)
        c0 = A * FTm
        c1 = (Dq * c0) * (1 / h)
        c2 = (Dq * (Dq * c0)) * (1 / (h * h))

        def ev(eta, F, Fq, Fqq):
            LF = 4 * (1 - eta * eta) * Fqq - 2 * (1 + 3 * eta) * Fq
            if p == 7:
                K = kamp * (1 - eta) / 2
                src = 2 * K / (F * F)
            else:
                src = 2 * kamp
            EV = F * LF - 2 * (1 - eta * eta) * Fq * Fq - F * F / 2 + src
            return EV / (2 * F * F.sqrt()), LF
        # real nodes
        vals = [Bq * cm for cm in (c0, c1, c2)]
        fv = arb_mat(n, 1)
        for i in range(n):
            tau = (arb(2 * i + 1) / (2 * n)).cos_pi()
            phv, _ = ev(m + h * tau, vals[0][i, 0], vals[1][i, 0], vals[2][i, 0])
            fv[i, 0] = phv
        cc = Cn * fv
        Snodes = mat_abs_sum(cc)
        # complex ellipse (half circle by conjugate symmetry)
        W = arc_balls(rho, nar, True)
        Mmax = arb(0)
        for Wb in W:
            T, sg = cheb_T_arc(Wb, nbq)
            Fz = [sum((acb(cm[k, 0]) * T[k] for k in range(nbq)), acb(0)) for cm in (c0, c1, c2)]
            ph, _ = ev(m + h * sg, *Fz)
            if not (Fz[0].real.lower() > 0):
                okall = False
            if ph.is_finite():
                Mmax = max(Mmax, abs(ph).upper())
            else:
                okall = False
        r = arb(rho)
        eps = 4 * Mmax * r ** (1 - n) / (r - 1)
        Aall = max(Aall, arb((Snodes + eps).upper()))
        # crude bounds for Fmin, sup|LF|
        rest = arb(0)
        for k in range(1, nbq):
            rest += abs(c0[k, 0])
        fl = (c0[0, 0] - rest).lower()
        Fmin = fl if Fmin is None else min(Fmin, fl)
        # |LF| <= 4 sup|F_etaeta| + 8 sup|F_eta|  (|1-eta^2| <= 1, |2(1+3 eta)| <= 8 on [-1,1])
        s1 = sum((abs(c1[k, 0]) for k in range(nbq)), arb(0))
        s2 = sum((abs(c2[k, 0]) for k in range(nbq)), arb(0))
        Bmax = max(Bmax, arb((4 * s2 + 8 * s1).upper()))
    Fmin = arb(Fmin)
    assert Fmin > 0
    if p == 7:
        Kmax = kamp
        Bcoef = (Bmax + 4 * Kmax / (Fmin ** 3)) / (2 * Fmin * Fmin.sqrt())
    else:
        Bcoef = Bmax / (2 * Fmin * Fmin.sqrt())
    return dict(A=Aall, B=arb(Bcoef.upper()), Fmin=Fmin, ok=okall)


def stage_assemble(case, ck):
    R = Rig(case, ck)
    P = R.P
    cells = ck['cells']
    ns, nq = VER['ns'], VER['nq']
    delta = arb(0)
    parts = {}
    allok = all(c['ok'] for c in cells.values())
    glo_all = None
    # cells
    dcell = arb(0)
    for ie, el in enumerate(P.elems):
        for i in range(ns):
            Bmax = max(arb(cells[f'{ie}_{i}_{j}']['B']) for j in range(nq))
            for j in range(nq):
                gl = arb(cells[f'{ie}_{i}_{j}']['glo'])
                if el['kind'] == 't':
                    glo_all = gl if glo_all is None else min(glo_all, gl)
            sa = arb(-1) + arb(2 * i) / ns
            sb = sa + arb(2) / ns
            if el['kind'] == 't':
                ta = R.ec[ie] + R.eh[ie] * sa
                tb_ = R.ec[ie] + R.eh[ie] * sb
                dcell += 2 * Bmax * ((-ta / 2).exp() - (-tb_ / 2).exp())
            else:
                x1 = arb(P.x1)
                xa = x1 * (1 + sa) / 2
                xb = x1 * (1 + sb) / 2
                dcell += Bmax * (xb - xa)
    parts['cells'] = dcell
    # interface kinks: contribution r^{1/2} sup|[F_t]| / (2 sqrt(g_lo))
    dshell = arb(0)
    nbq = P.nE
    for ie in range(len(P.elems)):
        el = P.elems[ie]
        C = R.coef[ie]
        nb = C.nrows()
        DC = R.dmat(nb) * C
        right = [sum((DC[k, l] for k in range(nb)), arb(0)) for l in range(nbq)]
        if el['kind'] == 'x':
            right = [2 * v for v in right]  # F_t = x F_x = x1 (2/x1) F_s
        else:
            right = [v / R.eh[ie] for v in right]
        if ie + 1 < len(P.elems):
            el2 = P.elems[ie + 1]
            C2 = R.coef[ie + 1]
            nb2 = C2.nrows()
            DC2 = R.dmat(nb2) * C2
            left = [sum((DC2[k, l] * (1 if k % 2 == 0 else -1) for k in range(nb2)), arb(0)) / R.eh[ie + 1]
                    for l in range(nbq)]
            tk = R.tbx[ie]
            gl = min(arb(cells[f'{ie+1}_0_{j}']['glo']) for j in range(nq))
        else:
            left = [arb(0)] * nbq
            tk = R.tbx[-1]
            gl = min(arb(cells[f'{ie}_{ns-1}_{j}']['glo']) for j in range(nq))
        jump = [right[l] - left[l] for l in range(nbq)]
        sj = sum((abs(v) for v in jump), arb(0))
        assert gl > 0
        dshell += (-tk / 2).exp() * sj / (2 * gl.sqrt())
    parts['shells'] = dshell
    # throat tail
    th = throat_bounds(R)
    Tp = R.tbx[-1]
    rp = (-Tp).exp()
    dthroat = 2 * th['A'] * rp.sqrt() + arb(2) / 3 * th['B'] * rp * rp.sqrt()
    parts['throat'] = dthroat
    delta = dcell + dshell + dthroat
    # sub-solution positivity (psi_a - h > 0) checks
    # r <= r+ : 4A + (2B r+ + delta/r+) sqrt(r+) < sqrt(FTmin)
    c1 = 4 * th['A'] + (2 * th['B'] * rp + delta / rp) * rp.sqrt() < th['Fmin'].sqrt()
    # r >= r+ (t-region): delta/r < sqrt(g_lo/r)  <=  r > delta^2/g_lo ; need delta^2/g_lo < r+
    c2 = delta * delta / glo_all < rp
    # far region (r >= 1/x1): h <= delta/r <= delta x1 < psi_a (>= sqrt of Re-lower of 1+xF)
    farlo = min(arb(cells[f'0_{i}_{j}']['glo']) for i in range(ns) for j in range(nq))
    c3 = delta * arb(P.x1) < farlo.sqrt()
    Ma = R.Ma
    Mlo = Ma - 2 * delta
    Mhi = Ma + 2 * delta
    out = dict(case=case, Ma=Ma.str(radius=True), delta=delta.str(radius=True),
               parts={k: v.str(radius=True) for k, v in parts.items()},
               M_lo=arb(Mlo.lower()).str(30), M_hi=arb(Mhi.upper()).str(30),
               M_lo_exact=str(Mlo.lower().mid()), M_hi_exact=str(Mhi.upper().mid()),
               cells_ok=allok, throat_ok=th['ok'], throat_A=th['A'].str(radius=True), throat_B=th['B'].str(radius=True),
               throat_Fmin=th['Fmin'].str(radius=True), glo_t=glo_all.str(radius=True), far_lo=farlo.str(radius=True),
               subsol_checks=[bool(c1), bool(c2), bool(c3)], f0_correction=R.f0_correction)
    return out


def main():
    cmd = sys.argv[1]
    global DISC, VER
    if cmd == 'run':
        case = sys.argv[2]
        if 'disc' in CASES[case]:
            DISC = CASES[case]['disc']
        if 'ver' in CASES[case]:
            VER = CASES[case]['ver']
        t_start = time.time()
        ck = load_ckpt(case)
        ck.setdefault('meta', dict(DISC=DISC, VER=VER, PREC=PREC))
        if 'solve' not in ck:
            ck = stage_solve(case, ck)
            save_ckpt(case, ck)
            print('solve done', ck['solve']['M_float'], flush=True)
        if not ck.get('cells_done'):
            ck, fin = stage_cells(case, ck, t_start)
            ck['cells_done'] = fin
            save_ckpt(case, ck)
            if not fin:
                print('budget reached; rerun to resume', flush=True)
                sys.exit(3)
        if 'final' not in ck:
            ck['final'] = stage_assemble(case, ck)
            save_ckpt(case, ck)
        print(json.dumps(ck['final'], indent=1))
        print('DONE', case)
        sys.exit(0)
    elif cmd == 'assemble':
        assemble()


def arb_to_frac(x):
    from fractions import Fraction
    m, e = x.mid().man_exp()
    m, e = int(m), int(e)
    return Fraction(m) * (Fraction(2) ** e)


def dec_floor(fr, nd):
    import math
    q = 10 ** nd
    v = math.floor(fr * q)
    sgn = '-' if v < 0 else ''
    v = abs(v)
    return sgn + str(v // q) + '.' + str(v % q).zfill(nd)


def dec_ceil(fr, nd):
    import math
    q = 10 ** nd
    v = math.ceil(fr * q)
    sgn = '-' if v < 0 else ''
    v = abs(v)
    return sgn + str(v // q) + '.' + str(v % q).zfill(nd)


def assemble():
    from fractions import Fraction
    ctx.prec = PREC
    ND = 15
    out = dict(seat='REP7', object='extreme (mu=0) Bowen-York, J=1: Delta u = -S(1+u)^-7, S = 9 sin^2/(4 r^6)',
               rounding='outward: lower ends floored, upper ends ceiled, %d decimals, from exact binary endpoints' % ND)
    win = {}
    for case in CASES:
        ck = load_ckpt(case)
        if 'final' not in ck:
            continue
        fin = ck['final']
        lo = arb_to_frac(arb(fin['M_lo_exact']).lower())
        hi = arb_to_frac(arb(fin['M_hi_exact']).upper())
        win[case] = (lo, hi)
        d = dict(M_lo=dec_floor(lo, ND), M_hi=dec_ceil(hi, ND), width=dec_ceil(hi - lo, 15),
                 M_a=fin['Ma'], delta=fin['delta'], parts=fin['parts'], cells_ok=fin['cells_ok'],
                 throat_ok=fin['throat_ok'], subsol_checks=fin['subsol_checks'],
                 newton_time_s=ck['solve']['time'], disc=ck['meta']['DISC'], ver=ck['meta']['VER'])
        if CASES[case]['p'] == 7:
            if case in ('BY', 'BY2', 'BYP'):
                d['chi_lo'] = dec_floor(1 / hi ** 2, ND)
                d['chi_hi'] = dec_ceil(1 / lo ** 2, ND)
        out[case] = d
    if 'BY' in win:
        lo, hi = win['BY']
        out['M0_lo'] = dec_floor(lo, ND)
        out['M0_hi'] = dec_ceil(hi, ND)
        out['M0_width'] = dec_ceil(hi - lo, 15)
        out['chi0_lo'] = dec_floor(1 / hi ** 2, ND)
        out['chi0_hi'] = dec_ceil(1 / lo ** 2, ND)
        out['chi0_width'] = dec_ceil(1 / lo ** 2 - 1 / hi ** 2, 15)
        out['width_target_1e-8_met'] = bool(hi - lo <= Fraction(1, 10 ** 8))
    if 'RN' in win:
        lo, hi = win['RN']
        out['RN_M_lo'] = dec_floor(lo, ND)
        out['RN_M_hi'] = dec_ceil(hi, ND)
        out['RN_contains_1'] = bool(lo <= 1 <= hi)
    if 'KN' in win:
        lo, hi = win['KN']
        out['KN_M_lo'] = dec_floor(lo, ND)
        out['KN_M_hi'] = dec_ceil(hi, ND)
        if 'BY' in win:
            blo, bhi = win['BY']
            out['KN_disjoint_from_BY'] = bool(lo > bhi or hi < blo)
            # scaling law (r -> b r, b^4 = J^2): M scales as kamp^(1/4); compare in ball arithmetic
            ctx.prec = 256
            bl = arb(blo.numerator) / arb(blo.denominator)
            bh = arb(bhi.numerator) / arb(bhi.denominator)
            f = (arb(1001) / 1000) ** (arb(1) / 4)
            plo, phi = (bl * f).lower(), (bh * f).upper()
            out['KN_predicted_from_BY_by_scaling'] = [dec_floor(arb_to_frac(plo), ND), dec_ceil(arb_to_frac(phi), ND)]
            plo_f, phi_f = arb_to_frac(plo), arb_to_frac(phi)
            out['KN_scaling_overlap'] = bool(not (plo_f > hi or phi_f < lo))
            # KN has J^2 = 1.001: chi = sqrt(1.001)/M^2 (must equal chi(0) by scaling)
            sj = (arb(1001) / 1000).sqrt()
            kl = arb(lo.numerator) / arb(lo.denominator)
            kh = arb(hi.numerator) / arb(hi.denominator)
            out['KN_chi'] = [dec_floor(arb_to_frac((sj / (kh * kh)).lower()), ND),
                             dec_ceil(arb_to_frac((sj / (kl * kl)).upper()), ND)]
            ctx.prec = PREC
    if 'BYP' in win and 'BY' in win:
        a_, b_ = win['BY'], win['BYP']
        out['BYP_contains_BY'] = bool(b_[0] <= a_[0] and a_[1] <= b_[1])
    if 'BY2' in win and 'BY' in win:
        a_, b_ = win['BY'], win['BY2']
        out['BY_BY2_overlap'] = bool(not (a_[0] > b_[1] or b_[0] > a_[1]))
        out['BY_BY2_intersection'] = [dec_floor(max(a_[0], b_[0]), ND), dec_ceil(min(a_[1], b_[1]), ND)]
    out = {k: v for k, v in out.items() if v is not None}
    with open(os.path.join(HERE, 'result_rep7.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
