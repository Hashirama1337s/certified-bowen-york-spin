"""BHSLICE FLUXT instrument: the throat flux E(sigma) on the independent 2-D extreme-data float solution (BAR-FLUXT.md).

Identity (scout/B5-DIGITS.md sec. 1; REFEREE-FLUX items 1-2): Omega := r^(1/2) Phi solves the autonomous equation
  EQ := Omega_ss + L_S Omega - Omega/4 + (9/4) sin^2(theta) Omega^-7 = 0   on R x S^2   (J = 1),
and E(sigma) = int_{S^2} [Omega_s^2/2 - |grad_S Omega|^2/2 - Omega^2/8 - (3/8) sin^2 Omega^-6] dOmega has
dE/dsigma = int_{S^2} Omega_s EQ dOmega (= 0 on an exact solution).

The 2-D solution: b3t.py float stage (seal 10a65e940a; imported, NOT edited), MAIN_CFG N = 112, K = 12, re-run here.
Ansatz (BUILD-B3T.md sec. 1): sigma = log r, x = tanh(sigma/8), a := sqrt r = e^(sigma/2) = ((1+x)/(1-x))^2,
  w = r^-1/2 G/(1 + a), G = g0 + (1-x)^2 H(x, c), H = sum h_nk T_n(x) P_2k(c).
So Omega = a + b with b := G/(1 + a), and by exact differentiation of the ansatz (x' = dx/dsigma = (1 - x^2)/8):
  Omega_s = a/2 + b_s,   b_s = G_x x'/(1 + a) - G a/(2 (1 + a)^2),   G_x = -2(1-x) H + (1-x)^2 H_x,
  Omega_c = G_c/(1 + a), G_c = (1-x)^2 H_c,   |grad_S Omega|^2 = (1 - c^2) Omega_c^2  (axisymmetric).
The e^sigma pieces cancel analytically (no round-off loss at sigma = 6):
  Omega_s^2/2 - Omega^2/8 = a b_s/2 + b_s^2/2 - a b/4 - b^2/8.
Density used (generic model: source coef ang r^-(p+5)/2 Phi^-p, so the potential is coef ang Omega^(1-p)/(p-1)):
  e = a b_s/2 + b_s^2/2 - a b/4 - b^2/8 - gcoef (1-c^2) Omega_c^2 - coef ang Omega^(1-p)/(p-1),
  gcoef = 1/2 (the identity), gcoef = 1/4 (F3 null).   BY: coef 9/4, p 7, ang 1-c^2 -> 3/8 sin^2 Omega^-6.
E(sigma) = 2 pi int_{-1}^{1} e dc by nc-point Gauss-Legendre in c = cos theta (nc = 64; control: 128).

Tolerance (stated BEFORE any slice is evaluated; derivation by this build seat, for the lead to pin):
  The Omega-equation residual of the ansatz is EQ = r^(5/2) (Delta w + S (1+w)^-7) = E_G/(1 + a)
  = Q (1 + 3a)/(4 (1 + a)^3)  (BUILD-B3T.md sec. 2 normalisation, Q = the b3t Q-normalised residual). Hence for any
  sigma_1 < sigma_2 (including +-oo)
    |E(sigma_2) - E(sigma_1)| <= sup|Q| * W,   W := int_R int_{S^2} |Omega_s| (1 + 3a)/(4 (1 + a)^3) dOmega dsigma.
  tol := sup|Q| * W, with sup|Q| = max(sampled supQ, -sampled infQ) from b3t.float_check (2001 x 41 grid; a float,
  NOT rigorous, sup) and W by Gauss-Legendre in x (dsigma = 8 dx/(1 - x^2)) on nW x nc nodes (W also at 2 nW, reported).
  On the float solution E(+oo) = -pi (2 g0) exactly (the ansatz has lim r w = g0), so the same tol bounds
  |E(sigma) + 2 pi g0| for every sigma.

Scoring (BAR-FLUXT.md, as sealed):
  F1: spread_i E(sigma_i) <= 10 tol.   F2: max_i |E(sigma_i)/(-pi) - M_0| <= 10 tol, M_0 = 1.03795639891365137036534431609.
  F3 (null, gcoef 1/4): spread_i E_1/4(sigma_i) > 100 tol.
  Controls (any failure = VOID): RN: analytic Omega^2 = e^s + q (q = 1 and q = 0.37) through the same density code,
  max_i |E/(-pi) - q| <= RN_TOL = 1e-13 (double-precision floor, fixed here); QUAD: max_i |E_128 - E_64| < tol.
  Verdict: HOLD iff F1, F2, F3 and both controls pass; F1 failing = STOP.
  Extra diagnostic (NOT a bar control, not scored): RN through the ansatz path (b3t float solve, model RN, the b3t
  RN control configuration N = 64, K = 0) must give E = -pi within its own tol.

Usage (every run: py -3.11 <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 fluxt.py ..., from <LAB>/BHSLICE):
  fluxt.py --rn        RN control only (writes seats/build-fluxt/rn_control.json)
  fluxt.py --dry-run   float solve + tolerance only, no slice evaluated (writes seats/build-fluxt/dryrun.json)
  fluxt.py             the full run: refuses unless GO-FLUXT exists and the lab ledger has a 'go' row naming fluxt.py
                       with this file's sha256; writes result_fluxt.json and FLUXT-REPORT.md.
Exit codes: 0 finished (any verdict), 1 control failed (VOID), 2 no GO-FLUXT, 6 self-sha not pinned.
"""
import sys
sys.dont_write_bytecode = True
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
import json, time, hashlib, argparse
from fractions import Fraction
import numpy as np
from numpy.polynomial import legendre as Lg

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import b3t                                   # sealed (10a65e940a); used read-only: float_solve, float_check, tbasis

SCRATCH = os.path.join(HERE, 'seats', 'build-fluxt')
GO_FILE = os.path.join(HERE, 'GO-FLUXT')
INSTRUMENT_NAME = 'fluxt.py'

M0_FLUX = '1.03795639891365137036534431609'   # scout/b5 b5_flux_coef.py value (BAR-FLUXT F2)
SIGMAS = np.linspace(-6.0, 6.0, 10)          # 10 slices in [-6, 6]
NC = 64                                      # Gauss-Legendre nodes in c on [-1, 1]
NC_DOUBLE = 128                              # quadrature-doubling control
NW = 2000                                    # Gauss-Legendre nodes in x for W
RN_TOL = 1e-13
RN_QS = (1.0, 0.37)
CFG = dict(b3t.MAIN_CFG)                     # model BY, fac 1, N 112, K 12
RN_ANSATZ_CFG = dict(model='RN', fac='1', N=64, K=0)


def sha_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        h.update(f.read())
    return h.hexdigest()


def sha_self():
    return sha_file(os.path.abspath(__file__))


# ================================================================== density (shared by every path)
def density(a, b, bs, Oc, c, mp, gcoef):
    """e(sigma, c) with Omega = a + b, Omega_s = a/2 + bs, Omega_c = Oc; arrays broadcast (rows sigma, cols c)."""
    p = mp['p']
    coef = float(mp['coef'])
    ang = (1 - c * c) if mp['ang'] == '1-c2' else np.ones_like(c)
    Om = a + b
    return (a * bs / 2 + bs * bs / 2 - a * b / 4 - b * b / 8
            - gcoef * (1 - c * c) * Oc * Oc - coef * ang * Om ** (1 - p) / (p - 1))


def E_quad(dens, wc):
    """2 pi int_{-1}^{1} e dc (axisymmetric S^2 integral), rows = slices."""
    return 2 * np.pi * (dens @ wc)


# ================================================================== the 2-D ansatz, exact derivatives
def pbasis_d(c, K):
    """P_2k(c) and d/dc P_2k(c), k = 0..K."""
    V = np.zeros((len(c), K + 1)); D = np.zeros_like(V)
    for k in range(K + 1):
        cf = np.zeros(2 * k + 1); cf[2 * k] = 1
        V[:, k] = Lg.legval(c, cf)
        D[:, k] = Lg.legval(c, Lg.legder(cf)) if k > 0 else 0.0
    return V, D


class Ansatz:
    def __init__(self, z, N, K):
        self.g0 = float(z[0]); self.N, self.K = N, K
        self.h = np.asarray(z[1:], dtype=float).reshape(N + 1, K + 1)   # column order of np.kron(Tv, Pv)

    def fields_x(self, x, c):
        """(a, b, b_s, Omega_c) on the grid x (rows) times c (cols)."""
        Tv, T1, _ = b3t.tbasis(x, self.N)
        Pv, Pc = pbasis_d(c, self.K)
        H = Tv @ self.h @ Pv.T; Hx = T1 @ self.h @ Pv.T; Hc = Tv @ self.h @ Pc.T
        X = x[:, None]
        w1 = 1 - X
        a = ((1 + X) / w1) ** 2                      # = e^(sigma/2)
        G = self.g0 + w1 ** 2 * H
        Gx = -2 * w1 * H + w1 ** 2 * Hx
        Gc = w1 ** 2 * Hc
        xs = (1 - X * X) / 8
        b = G / (1 + a)
        bs = Gx * xs / (1 + a) - G * a / (2 * (1 + a) ** 2)
        Oc = Gc / (1 + a)
        return a, b, bs, Oc, G

    def fields_sigma(self, sig, c):
        return self.fields_x(np.tanh(np.asarray(sig, dtype=float) / 8), c)


def W_integral(ans, nw, nc):
    """W = int_R int_{S^2} |Omega_s| (1+3a)/(4(1+a)^3) dOmega dsigma, in x with dsigma = 8 dx/(1-x^2)."""
    xq, wx = Lg.leggauss(nw)
    cq, wc = Lg.leggauss(nc)
    tot = 0.0
    for j0 in range(0, nw, 250):
        x = xq[j0:j0 + 250]
        a, b, bs, Oc, G = ans.fields_x(x, cq)
        Os = a / 2 + bs
        rho = (1 + 3 * a) / (4 * (1 + a) ** 3)
        integ = np.abs(Os) * rho * (8 / (1 - x * x))[:, None]
        tot += float(wx[j0:j0 + 250] @ (integ @ wc)) * 2 * np.pi
    return tot


# ================================================================== RN analytic control
def rn_analytic(sigmas, q, nc):
    """Omega^2 = e^s + q (extreme RN, p = 3, C = q^2/4): a = e^(s/2), b = q/(Omega + a), b_s = -a b/(2 Omega)."""
    cq, wc = Lg.leggauss(nc)
    a = np.exp(np.asarray(sigmas)[:, None] / 2) * np.ones((1, nc))
    Om = np.sqrt(np.exp(np.asarray(sigmas)[:, None]) + q) * np.ones((1, nc))
    b = q / (Om + a)
    bs = -a * b / (2 * Om)
    mp = dict(p=3, coef=Fraction(q * q / 4), ang='1')
    E = E_quad(density(a, b, bs, np.zeros_like(a), cq[None, :], mp, 0.5), wc)
    return E


def run_rn(log):
    out = dict(criterion=f'max_i |E(sigma_i)/(-pi) - q| <= {RN_TOL:g}', sigmas=SIGMAS.tolist(), cases=[])
    ok = True
    for q in RN_QS:
        for nc in (NC, NC_DOUBLE):
            E = rn_analytic(SIGMAS, q, nc)
            dev = float(np.max(np.abs(E / (-np.pi) - q)))
            out['cases'].append(dict(q=q, nc=nc, E=E.tolist(), max_dev=dev, ok=dev <= RN_TOL))
            ok &= dev <= RN_TOL
            log(f"RN analytic q={q} nc={nc}: E/(-pi) - q max |dev| = {dev:.3e}  {'ok' if dev <= RN_TOL else 'FAIL'}")
    out['ok'] = bool(ok)
    return out


def run_rn_ansatz(log):
    """diagnostic, not a bar control: b3t RN float solve, E via the ansatz path; expect -pi within sup|Q| W."""
    mp = b3t.model_params('RN', Fraction(RN_ANSATZ_CFG['fac']))
    N, K = RN_ANSATZ_CFG['N'], RN_ANSATZ_CFG['K']
    z, info = b3t.float_solve(mp, N, K, log=log)
    ck = b3t.float_check(mp, N, K, z)
    supabs = max(ck['supQ'], -ck['infQ'])
    ans = Ansatz(z, N, K)
    W = W_integral(ans, NW, NC)
    tol = supabs * W
    cq, wc = Lg.leggauss(NC)
    a, b, bs, Oc, _ = ans.fields_sigma(SIGMAS, cq)
    E = E_quad(density(a, b, bs, Oc, cq[None, :], mp, 0.5), wc)
    dev = float(np.max(np.abs(E / (-np.pi) - 1.0)))
    devE = float(np.max(np.abs(E + np.pi)))
    log(f"RN ansatz (diagnostic) N={N} K={K}: sup|Q| {supabs:.3e} W {W:.6f} tol {tol:.3e}; "
        f"max|E + pi| = {devE:.3e} ({devE / tol if tol > 0 else float('inf'):.3f} tol); 2 g0 = {2 * z[0]:.16f}")
    return dict(cfg=RN_ANSATZ_CFG, supabsQ=supabs, W=W, tol=tol, E=E.tolist(), max_abs_E_plus_pi=devE,
                max_dev_over_pi=dev, within_tol=bool(devE <= tol), two_g0=float(2 * z[0]), float_info=info)


# ================================================================== tolerance (before any slice)
def solve_and_tol(log):
    mp = b3t.model_params(CFG['model'], Fraction(CFG['fac']))
    N, K = CFG['N'], CFG['K']
    t0 = time.time()
    z, info = b3t.float_solve(mp, N, K, log=log)
    ck = b3t.float_check(mp, N, K, z)
    supabs = max(ck['supQ'], -ck['infQ'])
    ans = Ansatz(z, N, K)
    W = W_integral(ans, NW, NC)
    W2 = W_integral(ans, 2 * NW, NC_DOUBLE)
    tol = supabs * W
    t = dict(cfg=CFG, float_info=dict(info, **ck), g0=float(z[0]), two_g0=float(2 * z[0]), supabsQ=supabs,
             W=W, W_doubled=W2, tol=tol, solve_seconds=round(time.time() - t0, 1),
             rule='tol = max(sampled supQ, -sampled infQ) * W, W = int int |Omega_s| (1+3a)/(4(1+a)^3) dOmega dsigma')
    log(f"float solve BY N={N} K={K}: {info['iterations']} it, sampled supQ {ck['supQ']:.4e} infQ {ck['infQ']:.4e}, "
        f"Gmin {ck['Gmin']:.6f}, 2 g0 = {2 * z[0]:.16f}")
    log(f"W = {W:.12f} (nW {NW}, nc {NC});  W at doubled nodes = {W2:.12f}")
    log(f"TOLERANCE (stated before any slice): tol = {supabs:.6e} * {W:.6f} = {tol:.6e}")
    log(f"  F1 bar: spread <= 10 tol = {10 * tol:.6e};  F2 bar: |E/(-pi) - M_0| <= {10 * tol:.6e};  "
        f"F3 bar: spread > 100 tol = {100 * tol:.6e}")
    return ans, mp, t


# ================================================================== slices, scoring
def evaluate(ans, mp, nc):
    cq, wc = Lg.leggauss(nc)
    a, b, bs, Oc, _ = ans.fields_sigma(SIGMAS, cq)
    E = E_quad(density(a, b, bs, Oc, cq[None, :], mp, 0.5), wc)
    E4 = E_quad(density(a, b, bs, Oc, cq[None, :], mp, 0.25), wc)
    return E, E4


def score(E, E4, Ed, E4d, tol, rn):
    M0 = float(M0_FLUX)
    spread = float(E.max() - E.min())
    f2dev = float(np.max(np.abs(E / (-np.pi) - M0)))
    spread4 = float(E4.max() - E4.min())
    qd = float(np.max(np.abs(Ed - E)))
    qd4 = float(np.max(np.abs(E4d - E4)))
    F1 = spread <= 10 * tol
    F2 = f2dev <= 10 * tol
    F3 = spread4 > 100 * tol
    QUAD = qd < tol
    controls_ok = bool(rn['ok'] and QUAD)
    if not controls_ok:
        verdict = 'VOID'
    elif not F1:
        verdict = 'FAILED (F1: stop)'
    elif F2 and F3:
        verdict = 'FLUXT HOLD'
    else:
        verdict = 'FAILED'
    return dict(F1=dict(spread=spread, spread_over_tol=spread / tol, bar='<= 10 tol', ok=bool(F1)),
                F2=dict(max_abs_E_over_minus_pi_minus_M0=f2dev, over_tol=f2dev / tol, M0=M0_FLUX, bar='<= 10 tol',
                        ok=bool(F2)),
                F3=dict(spread=spread4, spread_over_tol=spread4 / tol, bar='> 100 tol', drifts=bool(F3)),
                QUAD=dict(max_abs_change=qd, over_tol=qd / tol, F3_max_abs_change=qd4, bar='< tol', ok=bool(QUAD)),
                RN=dict(ok=bool(rn['ok'])), controls_ok=controls_ok, verdict=verdict)


def go_pinned(sha):
    from pathlib import Path
    p = str(Path(os.environ.get('LAB_DATA') or (Path.home() / '.pocketlab')) / 'ledger.jsonl')
    if not os.path.exists(p):
        return False, f'ledger {p} not found'
    for line in open(p, encoding='utf-8'):
        try:
            r = json.loads(line)
        except ValueError:
            continue
        if r.get('kind') == 'go':
            for i in r.get('instruments', []):
                if str(i.get('path', '')).lower().replace('\\', '/').endswith('/' + INSTRUMENT_NAME) and i.get('sha256') == sha:
                    return True, f"GO row {r.get('at')} seal {r.get('seal_id')}"
    return False, f'no GO row names {INSTRUMENT_NAME} with sha256 {sha}'


def write_report(path, R):
    s, t = R['score'], R['tolerance']
    L = ['# FLUXT-REPORT: throat flux E(sigma) on the 2-D float solution (written by fluxt.py)', '',
         f"- Instrument: fluxt.py sha256 `{R['instrument_sha256']}`; b3t.py sha256 `{R['b3t_sha256']}` (imported, float stage).",
         f"- Verdict: **{s['verdict']}**", '',
         '## Tolerance (stated before any slice)',
         f"- sup|Q| (sampled) = {t['supabsQ']:.6e}; W = {t['W']:.10f} (doubled nodes {t['W_doubled']:.10f}); tol = {t['tol']:.6e}.",
         f"- Float solution: 2 g0 = {t['two_g0']:.16f} (E(+oo) = -pi 2 g0 on this solution).", '',
         '## Slices (nc = %d)' % NC, '', '| sigma | E | E/(-pi) - M_0 | E (grad coef 1/4) |', '|---|---|---|---|']
    M0 = float(M0_FLUX)
    for sg, e, e4 in zip(R['sigmas'], R['E'], R['E_quarter']):
        L.append(f"| {sg:+.4f} | {e:.16f} | {e / (-np.pi) - M0:+.3e} | {e4:.16f} |")
    L += ['', '## Scorecard', '', '| test | value | bar | result |', '|---|---|---|---|',
          f"| F1 spread | {s['F1']['spread']:.3e} ({s['F1']['spread_over_tol']:.3f} tol) | <= 10 tol | {'HOLD' if s['F1']['ok'] else 'FAILED'} |",
          f"| F2 max abs(E/(-pi) - M_0) | {s['F2']['max_abs_E_over_minus_pi_minus_M0']:.3e} ({s['F2']['over_tol']:.3f} tol) | <= 10 tol | {'HOLD' if s['F2']['ok'] else 'FAILED'} |",
          f"| F3 null spread (1/4) | {s['F3']['spread']:.3e} ({s['F3']['spread_over_tol']:.3e} tol) | > 100 tol | {'DRIFTS' if s['F3']['drifts'] else 'FAILED (no drift)'} |",
          f"| QUAD doubling | {s['QUAD']['max_abs_change']:.3e} ({s['QUAD']['over_tol']:.3f} tol) | < tol | {'ok' if s['QUAD']['ok'] else 'VOID'} |",
          f"| RN analytic | max dev {max(c['max_dev'] for c in R['rn']['cases']):.3e} | <= {RN_TOL:g} | {'ok' if R['rn']['ok'] else 'VOID'} |",
          '', '## Diagnostic (not a bar control)',
          f"- RN through the ansatz path (b3t RN float solve N=64 K=0): max abs(E + pi) = {R['rn_ansatz']['max_abs_E_plus_pi']:.3e}, "
          f"its own tol {R['rn_ansatz']['tol']:.3e}, within: {R['rn_ansatz']['within_tol']}.",
          '', '## Notes',
          '- tol is a float tolerance from a sampled residual sup (not rigorous).',
          '- On the float solution E(+oo) = -pi 2 g0 exactly, so F2 tests 2 g0 (and the slices) against the 1-D flux value M_0.']
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(L) + '\n')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--rn', action='store_true')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()
    os.makedirs(SCRATCH, exist_ok=True)
    mode = 'rn' if args.rn else ('dryrun' if args.dry_run else 'main')
    logp = os.path.join(SCRATCH, f'{mode}.log') if mode != 'main' else os.path.join(HERE, 'result_fluxt.json.log')

    def log(s):
        line = f"{time.strftime('%H:%M:%S')} {s}"
        print(line, flush=True)
        with open(logp, 'a', encoding='utf-8') as f:
            f.write(line + '\n')

    sha = sha_self()
    log(f'--- fluxt.py {mode} invocation sha256 {sha}')
    if mode == 'rn':
        rn = run_rn(log)
        rn['ansatz_diagnostic'] = run_rn_ansatz(log)
        with open(os.path.join(SCRATCH, 'rn_control.json'), 'w', encoding='utf-8') as f:
            json.dump(dict(instrument_sha256=sha, **rn), f, indent=1)
        sys.exit(0 if rn['ok'] else 1)
    if mode == 'dryrun':
        _, _, t = solve_and_tol(log)
        with open(os.path.join(SCRATCH, 'dryrun.json'), 'w', encoding='utf-8') as f:
            json.dump(dict(instrument_sha256=sha, tolerance=t, note='dry run: no slice evaluated'), f, indent=1)
        sys.exit(0)
    # ---- full run
    if not os.path.exists(GO_FILE):
        log(f'REFUSED: {GO_FILE} does not exist'); sys.exit(2)
    ok, why = go_pinned(sha)
    if not ok:
        log(f'REFUSED: {why}'); sys.exit(6)
    log('GO: ' + why)
    ans, mp, t = solve_and_tol(log)                         # tolerance printed and logged first
    rn = run_rn(log)
    rna = run_rn_ansatz(log)
    E, E4 = evaluate(ans, mp, NC)
    Ed, E4d = evaluate(ans, mp, NC_DOUBLE)
    for sg, e, e4 in zip(SIGMAS, E, E4):
        log(f"sigma {sg:+.4f}: E = {e:.16f}  E/(-pi) - M0 = {e / (-np.pi) - float(M0_FLUX):+.3e}  E_1/4 = {e4:.16f}")
    s = score(E, E4, Ed, E4d, t['tol'], rn)
    log(f"F1 spread {s['F1']['spread']:.3e} ({s['F1']['spread_over_tol']:.3f} tol)  F2 {s['F2']['over_tol']:.3f} tol  "
        f"F3 spread {s['F3']['spread_over_tol']:.3e} tol  QUAD {s['QUAD']['over_tol']:.3f} tol  RN {rn['ok']}")
    log('VERDICT: ' + s['verdict'])
    R = dict(bar='BAR-FLUXT.md', instrument_sha256=sha, b3t_sha256=sha_file(os.path.join(HERE, 'b3t.py')),
             go=why, tolerance=t, sigmas=SIGMAS.tolist(), nc=NC, nc_double=NC_DOUBLE, E=E.tolist(), E_quarter=E4.tolist(),
             E_nc_double=Ed.tolist(), E_quarter_nc_double=E4d.tolist(), rn=rn, rn_ansatz=rna, score=s)
    with open(os.path.join(HERE, 'result_fluxt.json'), 'w', encoding='utf-8') as f:
        json.dump(R, f, indent=1)
    write_report(os.path.join(HERE, 'FLUXT-REPORT.md'), R)
    sys.exit(0 if s['controls_ok'] else 1)


if __name__ == '__main__':
    main()
