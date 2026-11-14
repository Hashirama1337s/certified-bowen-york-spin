"""REP6 numerics (float64): Newton/continuation solver for the Bowen-York spinning puncture
in the log-tanh formulation.  See RESULT.md for the derivation.

Coordinates: sigma = log(2 r) (mu = 1), xi = tanh(sigma/8) in (-1, 1), c = cos(theta).
  A = 2 (1 + 6 xi^2 + xi^4) = (1+xi)^4 + (1-xi)^4,  w = 1 - xi^2.
Unknown: u = 2 G / (2 r + 1),  G = g0 + w^4 H(xi, c)   (G even in xi and in c for the main problem).
Normalised residual Rhat = r^2 (Delta u + S Phi^-7) / (e^{-sigma/2} sech^3(sigma/2)):
  Rhat = (A^2/256)[w^2 H'' - 18 xi w H' + (64 xi^2 - 8 w) H] - (xi (1+xi^2) A/4)(w H' - 8 xi H)
         + (A^2/4) Lap_S H - G/2 + CS (1 - c^2) A^10 w^8 / E^7,
  E = A^2 + 2 w^4 G,  CS = (9/32) j^2  (j = 4 J / mu^2, mu = 1).
"""
import numpy as np
from numpy.polynomial import chebyshev as C
from numpy.polynomial import legendre as Lg


def cheb_even_vals(xi, N):
    """T_{2n}(xi), first and second xi-derivatives, n = 0..N-1; arrays (len(xi), N)."""
    xi = np.asarray(xi, float)
    V = np.zeros((xi.size, N)); V1 = np.zeros_like(V); V2 = np.zeros_like(V)
    for n in range(N):
        cf = np.zeros(2 * n + 1); cf[2 * n] = 1.0
        V[:, n] = C.chebval(xi, cf)
        V1[:, n] = C.chebval(xi, C.chebder(cf, 1)) if n > 0 else 0.0
        V2[:, n] = C.chebval(xi, C.chebder(cf, 2)) if n > 0 else 0.0
    return V, V1, V2


def cheb_all_vals(xi, N):
    """T_n(xi), n = 0..N-1 with derivatives (for the non-even tangent problem)."""
    xi = np.asarray(xi, float)
    V = np.zeros((xi.size, N)); V1 = np.zeros_like(V); V2 = np.zeros_like(V)
    for n in range(N):
        cf = np.zeros(n + 1); cf[n] = 1.0
        V[:, n] = C.chebval(xi, cf)
        V1[:, n] = C.chebval(xi, C.chebder(cf, 1)) if n > 0 else 0.0
        V2[:, n] = C.chebval(xi, C.chebder(cf, 2)) if n > 1 else 0.0
    return V, V1, V2


def leg_even_vals(c, K):
    c = np.asarray(c, float)
    P = np.zeros((c.size, K + 1)); L = np.zeros_like(P)
    for k in range(K + 1):
        cf = np.zeros(2 * k + 1); cf[2 * k] = 1.0
        P[:, k] = Lg.legval(c, cf)
        L[:, k] = -(2 * k) * (2 * k + 1) * P[:, k]
    return P, L


class Grid:
    def __init__(self, xi, c, N, K, even=True):
        XI, CC = np.meshgrid(xi, c, indexing='ij')
        self.xi = XI.ravel(); self.c = CC.ravel()
        if even:
            V, V1, V2 = cheb_even_vals(xi, N)
        else:
            V, V1, V2 = cheb_all_vals(xi, N)
        P, LP = leg_even_vals(c, K)
        # basis tensors: column index n*(K+1)+k
        nb = V.shape[1] * (K + 1)
        self.B = np.einsum('in,jk->ijnk', V, P).reshape(-1, nb)
        self.B1 = np.einsum('in,jk->ijnk', V1, P).reshape(-1, nb)
        self.B2 = np.einsum('in,jk->ijnk', V2, P).reshape(-1, nb)
        self.BL = np.einsum('in,jk->ijnk', V, LP).reshape(-1, nb)
        x = self.xi
        self.A = 2 * (1 + 6 * x ** 2 + x ** 4)
        self.w = 1 - x ** 2
        A, w = self.A, self.w
        # linear operator on H coefficients (excluding the -G/2 term's w^4 H part, added below)
        self.Lin = ((A ** 2 / 256)[:, None] * (w[:, None] ** 2 * self.B2 - (18 * x * w)[:, None] * self.B1
                                                 + (64 * x ** 2 - 8 * w)[:, None] * self.B)
                    - (x * (1 + x ** 2) * A / 4)[:, None] * (w[:, None] * self.B1 - (8 * x)[:, None] * self.B)
                    + (A ** 2 / 4)[:, None] * self.BL
                    - 0.5 * (w ** 4)[:, None] * self.B)


def residual(grid, g0, h, CS):
    H = grid.B @ h
    G = g0 + grid.w ** 4 * H
    E = grid.A ** 2 + 2 * grid.w ** 4 * G
    src = CS * (1 - grid.c ** 2) * grid.A ** 10 * grid.w ** 8 / E ** 7
    R = grid.Lin @ h - 0.5 * g0 + src
    return R, G, E, src


def newton(grid, g0, h, CS, iters=30, tol=1e-13, verbose=False):
    for it in range(iters):
        R, G, E, src = residual(grid, g0, h, CS)
        if np.any(E <= 0):
            raise RuntimeError('E <= 0 in Newton')
        dsrc = -7 * src * 2 * grid.w ** 4 / E       # d src / d G
        J = np.empty((R.size, 1 + h.size))
        J[:, 0] = -0.5 + dsrc
        J[:, 1:] = grid.Lin + (dsrc * grid.w ** 4)[:, None] * grid.B
        dx = np.linalg.lstsq(J, -R, rcond=None)[0]
        # damping
        lam = 1.0
        for _ in range(30):
            g1 = g0 + lam * dx[0]; h1 = h + lam * dx[1:]
            R1, G1, E1, _ = residual(grid, g1, h1, CS)
            if np.all(E1 > 0) and np.linalg.norm(R1) < np.linalg.norm(R) * (1 - 0.3 * lam) + 1e-300:
                break
            lam *= 0.5
        g0, h = g1, h1
        nr = np.max(np.abs(R1))
        if verbose:
            print(f'   it {it} lam {lam:.3g} maxres {nr:.3e} g0 {g0:.12g}')
        if lam * np.max(np.abs(dx)) < tol * (1 + abs(g0)):
            break
    return g0, h


def colloc_grid(N, K, extra=8):
    nx = N + extra
    i = np.arange(nx)
    xi = np.cos(np.pi * (i + 0.5) / (2 * nx))  # positive nodes of T_{2nx}
    xi = np.concatenate([[1.0], xi])
    nc = K + 1 + 3
    cn, _ = Lg.leggauss(2 * nc)
    cc = np.sort(cn[cn > 0])
    cc = np.concatenate([cc, [1.0]])
    return xi, cc


def solve(jtarget, N, K, CSfac=1.0, j0=2.0 ** -4, step=2.0, verbose=False, start=None):
    xi, cc = colloc_grid(N, K)
    grid = Grid(xi, cc, N, K)
    if start is None:
        g0 = 0.0; h = np.zeros(N * (K + 1)); j = j0
    else:
        g0, h, j = start
    while True:
        j = min(j, jtarget)
        CS = CSfac * (9.0 / 32.0) * j ** 2
        g0, h = newton(grid, g0, h, CS, verbose=verbose)
        if verbose:
            print(f'j = {j:.6g}  delta = {2*g0:.12g}')
        if j >= jtarget:
            break
        j *= step
    return g0, h, grid


def check_grid(g0, h, N, K, CS, nx=4001, nc=201):
    xi = np.linspace(0, 1, nx); c = np.linspace(0, 1, nc)
    g = Grid(xi, c, N, K)
    R, G, E, src = residual(g, g0, h, CS)
    return R, G, E, g
