"""BHSLICE P2b instrument: certified Bowen-York spinning puncture, mu = 1, J = j/4, j in [10, 1000].

A thin runner over the shared engine <LAB>/POCKETLAB/gravity_engine/lichnerowicz_axisym.py (imported, never copied,
never edited). The engine holds the whole certificate (R1-R5: exact residual Q, abs-Chebyshev bound S, margin eta,
exact Bernstein positivity of U_-/s, slope route for delta = M_ADM - 1), the Bowen-York problem object, the K1 cell
record (arb, outward) and the K2 gap check (build_result). This file adds only:
  - the P2b degree schedule (N, K) rising with j (justified in BUILD-P2b.md),
  - the P2b grid: j_{k+1} = (21/20) j_k from j_0 = 10 to the first node >= 1000, plus j = 30, 100, 300, 1000,
  - a resumable, chunked runner (exact "p/q" strings in JSON checkpoints; no pickle),
  - the controls (DESIGN-P2 sec 5 adapted to P2b; brief).
Design: DESIGN-P2.md sections 4.3 and 5. Blind: BAR-P2b.md and TARGETS-P2.md were not opened by the build.

Usage (every run under memcap):
  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 p2b.py --probe --j 1000 [--N 112 --K 14]
  ... p2b.py --controls controls_p2b.json     (writes controls_p2b.log; < 10 min)
  ... p2b.py --main result_p2b.json           (refuses without GO-P2b; resumable; <= 520 s per invocation)
"""
import sys
sys.dont_write_bytecode = True
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # single-thread LAPACK => deterministic fits
    os.environ[_v] = '1'
ENGINE_DIR = '<LAB>/POCKETLAB'
sys.path.insert(0, ENGINE_DIR)
from gravity_engine import lichnerowicz_axisym as la          # import BEFORE numpy (deterministic fit)
import json, time, math, hashlib, argparse, subprocess, shutil, random
from fractions import Fraction
from flint import arb, fmpq, fmpq_poly, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-P2b')
ENGINE_FILE = os.path.join(ENGINE_DIR, 'gravity_engine', 'lichnerowicz_axisym.py')
ENGINE_TESTS = os.path.join(ENGINE_DIR, 'gravity_engine', 'tests')
SCRATCH = os.path.join(HERE, 'seats', 'build-p2b')
TIME_BUDGET = 520.0
J0, RATIO, J_END = Fraction(10), Fraction(21, 20), Fraction(1000)
EXTRA_NODES = (30, 100, 300, 1000)
GRID_LINE = 'j_{k+1} = (21/20) j_k from 10 to the first node >= 1000, plus j = 30, 100, 300, 1000'

# Degree schedule: (upper end of the band in j, N, K). N = Chebyshev degree of p, q_l in X = 2s - 1; Legendre modes
# l = 2..2K. Measured at the band tops (scratch probes, BUILD-P2b sec 2): every band top certified with relative width
# of delta <= 8e-8. K is NOT the limiting degree (top Legendre mode <= 1e-15 relative at K = 12); N is.
SCHEDULE = ((Fraction(20), 36, 10), (Fraction(40), 40, 12), (Fraction(80), 48, 12), (Fraction(150), 60, 12),
            (Fraction(300), 80, 12), (Fraction(600), 96, 14), (None, 112, 14))


def degrees(j):
    j = Fraction(j)
    for top, N, K in SCHEDULE:
        if top is None or j <= top:
            return N, K


class P2bBowenYork(la.BowenYorkSpinningPuncture):
    """The engine's Bowen-York problem (same source, tau, observables, brackets, toy, cells) with the P2b degree
    schedule. Fit-operator caches are cleared whenever (N, K) changes, so memory stays at one degree's worth."""
    EXTRA_NODES = EXTRA_NODES

    def degrees(self, j):
        return degrees(j)

    def float_ops(self, N, K):
        if (N, K) not in self._fops:
            self._fops.clear()
            self._aops.clear()
        return la.AxisymLichnerowiczProblem.float_ops(self, N, K)

    def arb_ops(self, N, K):
        if (N, K, la.PREC) not in self._aops:
            self._aops.clear()
        return la.AxisymLichnerowiczProblem.arb_ops(self, N, K)


PROBLEM = P2bBowenYork()


def grid_nodes():
    return PROBLEM.grid_nodes(j0=J0, ratio=RATIO, j_end=J_END, extra=EXTRA_NODES)


# ================================================================== small helpers
def sha_file(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def sha_self():
    return sha_file(os.path.abspath(__file__))


def atomic_json(path, obj, sort=True):
    tmp = path + '.tmp'
    with open(tmp, 'w') as f:
        json.dump(obj, f, indent=1, sort_keys=sort)
    os.replace(tmp, path)


def peak_memory():
    """peak working set and peak commit (private bytes) of THIS process, GB (Windows psapi); None elsewhere."""
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


def certify(j, log, keep=False, N=None, K=None, coef_override=None):
    return la.certify_node(PROBLEM, Fraction(j), N=N, K=K, log=log, keep=keep, coef_override=coef_override)


def finite_record(rec):
    """NaN/inf guard on a certified record: every float field finite, every exact endpoint parses."""
    if rec.get('status') != 'CERTIFIED':
        return True
    try:
        for k in ('delta_lo', 'delta_hi', 'chi_lo', 'chi_hi', 'S', 'eta'):
            la.qparse(rec[k])
        return all(math.isfinite(rec[k]) for k in ('rel_width', 'chi_rel_width', 'S_float', 'eta_float'))
    except Exception:
        return False


def build_result(node_recs, sha):
    """engine K1/K2 result + the P2b summary keys."""
    res = PROBLEM.build_result(node_recs, sha, j_first=J0, j_last=J_END, grid_line=GRID_LINE)
    nodes = res['nodes']
    rels = [r['rel_width'] for r in nodes if r['status'] == 'CERTIFIED']
    res['max_node_rel_width'] = max(rels) if rels else None
    res['max_node_rel_width_j'] = (max(nodes, key=lambda r: r.get('rel_width', -1))['j_dec'] if rels else None)
    res['nodes_finite'] = all(finite_record(r) for r in nodes)
    cert = [r for r in nodes if r['status'] == 'CERTIFIED']
    if cert:
        best = max(cert, key=lambda r: la.qparse(r['chi_lo']))
        res['chi_sup_lower_bound'] = dict(j=best['j_dec'], chi_lo=best['chi_lo'], chi_lo_dec=best['chi_lo_dec'],
                                          note='X3: the supremum of chi over the family is >= this certified value')
    if res['cells']:
        top = max(res['cells'], key=lambda c: la.qparse(c['chi_max']))
        res['max_cell_chi_max'] = dict(j_lo=top['j_lo_dec'], j_hi=top['j_hi_dec'], chi_max_dec=top['chi_max_dec'])
    res['degree_schedule'] = [dict(j_upto=(None if t is None else str(t)), N=N, K=K) for t, N, K in SCHEDULE]
    res['engine'] = dict(file=ENGINE_FILE, sha256=sha_file(ENGINE_FILE), seal=la.SEAL_ID)
    res['peak_memory_note'] = ('per-invocation peak memory is logged in the run log and <out>.runinfo.json (kept out of '
                               'this file so that it is byte-reproducible); the MEM control in controls_p2b.json gives '
                               'the peak of the largest-degree node (N, K) = (112, 14)')
    ef = []
    if not res['no_gaps'] and res['all_nodes_certified'] and res['nodes'] and \
            Fraction(*map(int, res['nodes'][-1]['j'].split('/'))) < J_END:
        ef.append('no_gaps: sub-grid run (does not reach j = 1000); false by design for the RESUME control only')
    res['expected_false'] = ef
    return res


# ================================================================== resumable grid runner (main and RESUME control)
class RunState:
    """checkpoint = JSON with exact "p/q" strings; refuses a checkpoint from another instrument or engine sha."""

    def __init__(self, out, nodes):
        self.out = out
        self.ck = out + '.ckpt'
        self.sha = sha_self()
        self.esha = sha_file(ENGINE_FILE)
        self.spec = [la.qstr(la.fq(j)) for j in nodes]
        self.st = None

    def load(self, log):
        if os.path.exists(self.ck):
            st = json.load(open(self.ck))
            if st.get('instrument_sha256') != self.sha or st.get('engine_sha256') != self.esha:
                log(f"REFUSED: checkpoint {self.ck} was written by instrument {st.get('instrument_sha256')} / engine "
                    f"{st.get('engine_sha256')}; this is {self.sha} / {self.esha}.")
                return False
            if st.get('nodes_spec') != self.spec:
                log("REFUSED: checkpoint node list differs from this run's node list.")
                return False
        else:
            st = dict(instrument_sha256=self.sha, engine_sha256=self.esha, nodes_spec=self.spec, done={}, invocations=0)
        st['invocations'] += 1
        self.st = st
        self.save()
        return True

    def save(self):
        atomic_json(self.ck, self.st)


def run_grid(out, nodes, budget, log):
    """Exit codes: 0 = finished or budget reached (rerun to resume), 3 = checkpoint refused."""
    T0 = time.time()
    rs = RunState(out, nodes)
    if not rs.load(log):
        return 3
    info_p = out + '.runinfo.json'
    info = json.load(open(info_p)) if os.path.exists(info_p) else dict(timing={}, invocations=[])
    worst = 5.0
    did = 0
    for j in nodes:
        key = la.qstr(la.fq(j))
        if key in rs.st['done']:
            continue
        N, K = degrees(j)
        # time guard: predicted node time from this invocation's worst at the same or lower degree, scaled by N^2 K
        if did >= 1 and time.time() - T0 + 1.5 * worst * pred_scale(N, K) > budget - 20:
            log(f"time budget reached: {len(rs.st['done'])}/{len(nodes)} nodes done; rerun the same command to resume")
            info['invocations'].append(dict(ended='budget', seconds=time.time() - T0, peak=peak_memory()))
            atomic_json(info_p, info)
            return 0
        t = time.time()
        rec, _ = certify(j, log)
        dt = time.time() - t
        worst = max(worst, dt / pred_scale(N, K))
        info['timing'][key] = dict(N=N, K=K, total=dt, parts=rec.pop('_timing', None))
        rs.st['done'][key] = rec
        rs.save()
        atomic_json(info_p, info)
        did += 1
    res = build_result(list(rs.st['done'].values()), rs.sha)
    atomic_json(out, res)
    info['invocations'].append(dict(ended='finished', seconds=time.time() - T0, peak=peak_memory()))
    atomic_json(info_p, info)
    log(f"grid finished: {len(nodes)} nodes, all certified = {res['all_nodes_certified']}, no gaps = {res['no_gaps']}, "
        f"max node rel width = {res['max_node_rel_width']}, max cell rel width = {res['max_cell_rel_width']}, "
        f"invocations = {rs.st['invocations']}, peak memory this invocation = {peak_memory()}")
    return 0


def pred_scale(N, K):
    """relative cost model for the time guard (measured: the exact S step dominates; ~ N^2.5 K fits 1.5 s at (36, 10) and 34 s at (112, 14))."""
    return (N / 36.0) ** 2.5 * (K / 10.0)


# ================================================================== independent reference checks (REF)
def _cheb_int(n, x):
    """T_k(x) q^k as python ints for x = p/q (own recurrence, k = 0..n)."""
    p, q = int(x.p), int(x.q)
    a = [1, p]
    for k in range(1, n):
        a.append(2 * p * a[k] - q * q * a[k - 1])
    return a[:n + 1], q


def ref_chebyshev(fn, rec, npts, seed, log):
    """exact Chebyshev reconstruction of Q at random rational points vs Q straight from the ansatz terms."""
    t = time.time()
    S, Q = fn.residual_S(want_Q=True)
    out = dict(S_matches_record=bool(la.qstr(S) == rec['S']))
    C, sc = la.chebmat(Q, fn.W)
    rows, cols = C.nrows(), C.ncols()
    Ci = [[int(C[i, k]) for k in range(cols)] for i in range(rows)]
    out['abs_sum_equals_S'] = bool(fmpq(sum(abs(v) for r in Ci for v in r)) * sc == S)
    rng = random.Random(seed)
    pts = []
    for _ in range(npts):
        s0 = fmpq(rng.randint(1, 996), 997)
        c0 = fmpq(rng.randint(-990, 990), 991)
        As, qs = _cheb_int(rows - 1, 2 * s0 - 1)
        Ac, qc = _cheb_int(cols - 1, c0)
        ws = [As[i] * qs ** (rows - 1 - i) for i in range(rows)]
        tot = 0
        for k in range(cols):
            inner = 0
            for i in range(rows):
                v = Ci[i][k]
                if v:
                    inner += v * ws[i]
            tot += inner * Ac[k] * qc ** (cols - 1 - k)
        series = fmpq(tot) / (fmpq(qs) ** (rows - 1) * fmpq(qc) ** (cols - 1)) * sc
        u = sum((f(s0) * g(c0) for f, g in fn.Ut), fmpq(0))
        Lu = sum((f(s0) * g(c0) for f, g in fn.Lt), fmpq(0))
        D = 1 + (1 - s0) * u
        src = fmpq(9, 4) * (1 - c0 ** 2) * s0 ** 2 * (1 - s0) ** 2
        direct = D ** 7 * Lu / fn.tau + src
        pts.append(dict(s=la.qstr(s0), c=la.qstr(c0), equal=bool(series == direct), within_S=bool(abs(series) <= S),
                        Q_over_S=la.qfloat(series / S)))
    out['points'] = pts
    out['ok'] = bool(out['S_matches_record'] and out['abs_sum_equals_S'] and all(p['equal'] and p['within_S'] for p in pts))
    out['seconds'] = time.time() - t
    log(f"REF Chebyshev j={rec['j_dec'][:10]}: S matches {out['S_matches_record']}, sum|C| sc == S {out['abs_sum_equals_S']}, "
        f"{sum(p['equal'] for p in pts)}/{npts} points series == direct, all |Q| <= S: {all(p['within_S'] for p in pts)} "
        f"({out['seconds']:.1f}s)")
    return out


def ref_spherical(fn, rec, j, seed, log):
    """residual signs from the flat Laplacian in spherical coordinates (mpmath, numerical differentiation; no use of
    the engine's operator Lh or its exact polynomials): E[u] = u_rr + 2u_r/r + (u_thth + cot u_th)/r^2
    + (9/4) J^2 sin^2 r^-6 Phi^-7 must be < 0 for u_+ and > 0 for u_-; modes l > 2 must be present."""
    import mpmath as mp
    t = time.time()
    mp.mp.dps = 60
    K = len(fn.qcs)
    lmodes = sum(1 for k in range(1, K) if any(c != 0 for c in fn.qcs[k]))

    def mq(q):
        return mp.mpf(int(q.p)) / int(q.q)
    pcm = [mq(c) for c in fn.pc]
    qcm = [[mq(c) for c in q] for q in fn.qcs]
    taum = mp.mpf(Fraction(j).numerator) ** 2 / mp.mpf(Fraction(j).denominator) ** 2
    etam = mq(la.qparse(rec['eta']))
    Sm = mq(la.qparse(rec['S']))
    J = mp.mpf(Fraction(j).numerator) / Fraction(j).denominator / 4

    def cheb(cs, X):
        th = mp.acos(X)
        return mp.fsum(c * mp.cos(n * th) for n, c in enumerate(cs))

    def uN(r, th):
        s = 1 / (2 * r + 1)
        X = 2 * s - 1
        c = mp.cos(th)
        val = s * cheb(pcm, X)
        for k, q in enumerate(qcm, 1):
            val += s ** 2 * (1 - s) * cheb(q, X) * mp.legendre(2 * k, c)
        return val

    rng = random.Random(seed)
    pts, teeth_hits = [], []
    for _ in range(6):
        r = mp.mpf(rng.choice([1, 2, 3, 5, 7])) / rng.choice([1, 2, 4, 8, 16]) * rng.choice([1, 1, 3, 10])
        th = mp.pi * rng.randint(1, 15) / 32
        ur = mp.diff(uN, (r, th), (1, 0))
        urr = mp.diff(uN, (r, th), (2, 0))
        ut = mp.diff(uN, (r, th), (0, 1))
        utt = mp.diff(uN, (r, th), (0, 2))
        lapN = urr + 2 * ur / r + (utt + mp.cot(th) * ut) / r ** 2
        u0 = uN(r, th)
        s = 1 / (2 * r + 1)
        lap_s = 8 / (2 * r + 1) ** 3 - 4 / ((2 * r + 1) ** 2 * r)     # flat Laplacian of the radial function s(r)
        src = mp.mpf(9) / 4 * J ** 2 * mp.sin(th) ** 2 * r ** -6

        def E(sgn, fac=1):
            u = u0 + sgn * taum * etam * s
            return lapN + sgn * taum * etam * lap_s + fac * src * (1 + 1 / (2 * r) + u) ** -7
        Ep, Em, E0 = E(+1), E(-1), E(0)
        fac = 4 * s ** 4 / (1 - s) * taum
        Q0 = E0 / fac * (1 + (1 - s) * u0) ** 7
        # teeth: scan a relative source error 10^-eps, eps = 12..1; report the smallest error that breaks the super
        # sign. The super margin relative to the source term src Phi^-7 is large at big j (Phi^-7 suppresses the
        # source near the throat, DESIGN-P2 B2), so the bite appears at a j-dependent decade; the ratio is reported.
        srcterm = src * (1 + 1 / (2 * r) + u0 + taum * etam * s) ** -7
        bite = None
        for eps in range(12, 0, -1):
            if E(+1, 1 + mp.mpf(10) ** -eps) >= 0:
                bite = eps
                break
        teeth_hits.append(bite)
        pts.append(dict(r=float(r), theta=float(th), E_plus=float(Ep), E_minus=float(Em),
                        margin_over_source_term=float(-Ep / srcterm),
                        signs_ok=bool(Ep < 0 < Em), Q_within_S=bool(abs(Q0) <= Sm)))
    out = dict(points=pts, l_gt_2_modes_present=lmodes, teeth_largest_eps_that_bites_per_point=teeth_hits,
               teeth_ok=any(b is not None for b in teeth_hits), seconds=time.time() - t)
    out['ok'] = bool(lmodes >= 1 and all(p['signs_ok'] and p['Q_within_S'] for p in pts) and out['teeth_ok'])
    log(f"REF spherical j={rec['j_dec'][:10]}: signs E(u+) < 0 < E(u-) at {sum(p['signs_ok'] for p in pts)}/6 points, "
        f"|Q| <= S {all(p['Q_within_S'] for p in pts)}, modes l >= 4 nonzero: {lmodes}, teeth (10^-eps source error "
        f"breaks a sign) {teeth_hits} ({out['seconds']:.1f}s)")
    return out


def ref_bernstein(fn, rec, seed, log):
    """Bernstein reconstruction at random rational points on [0,1] and on a sub-piece (exact)."""
    eta = la.qparse(rec['eta'])
    T = fn.V_terms(eta)
    rng = random.Random(seed)
    res = []
    for a, b in ((fmpq(0), fmpq(1)), (fmpq(1, 4), fmpq(3, 8))):
        B, n, m = la.bernstein_matrix(T, a, b)
        for _ in range(4):
            tt = fmpq(rng.randint(0, 97), 97)
            w = fmpq(rng.randint(0, 89), 89)
            bt = [math.comb(n, i) * tt ** i * (1 - tt) ** (n - i) for i in range(n + 1)]
            bw = [math.comb(m, k) * w ** k * (1 - w) ** (m - k) for k in range(m + 1)]
            lhs = sum((B[i, k] * bt[i] * bw[k] for i in range(n + 1) for k in range(m + 1)), fmpq(0))
            s0, c0 = a + (b - a) * tt, 2 * w - 1
            rhs = sum((f(s0) * g(c0) for f, g in T), fmpq(0))
            res.append(bool(lhs == rhs))
    out = dict(points=len(res), all_equal=all(res), degrees=[n, m], ok=all(res) and m == 2 * len(fn.qcs))
    log(f"REF Bernstein j={rec['j_dec'][:10]}: {sum(res)}/{len(res)} reconstructions exact (n = {n}, m = {m})")
    return out


# ================================================================== controls
def run_controls(out, log):
    T0 = time.time()
    sha = sha_self()
    R = dict(instrument_sha256=sha, engine_sha256=sha_file(ENGINE_FILE), started=time.strftime('%Y-%m-%d %H:%M:%S'),
             degree_schedule=[dict(j_upto=(None if t is None else str(t)), N=N, K=K) for t, N, K in SCHEDULE])
    oks, expected_false = {}, []
    nodes = grid_nodes()
    R['grid'] = dict(line=GRID_LINE, n_nodes=len(nodes), first=la.qstr(la.fq(nodes[0])), last_dec=float(nodes[-1]))

    # --- MEM + control nodes. The largest degree (112, 14) runs FIRST so the process peak is that node's peak.
    R['MEM_baseline'] = peak_memory()
    recs, fns, secs = {}, {}, {}
    order = [Fraction(1000), Fraction(10), Fraction(30), Fraction(100), Fraction(300)]
    for jc in order:
        t = time.time()
        rec, fnc = certify(jc, log, keep=True)
        secs[jc] = time.time() - t
        recs[jc], fns[jc] = rec, fnc
        if jc == 1000:
            R['MEM'] = dict(node='1000/1', degrees=list(degrees(jc)), peak_after_node=peak_memory(),
                            baseline_before=R['MEM_baseline'], seconds=secs[jc], cap_gb=8)
            log(f"MEM: largest-degree node j=1000 (N, K) = {degrees(jc)}: process peak {R['MEM']['peak_after_node']} "
                f"(baseline {R['MEM_baseline']}), {secs[jc]:.1f}s")
    cnodes = sorted(order)
    R['control_nodes'] = [{k: v for k, v in recs[jc].items() if not k.startswith('_')} for jc in cnodes]
    R['node_seconds'] = {str(jc): round(secs[jc], 2) for jc in cnodes}
    R['node_timing_parts'] = {str(jc): recs[jc].get('_timing') for jc in cnodes}
    oks['control_nodes_certified'] = all(recs[jc]['status'] == 'CERTIFIED' for jc in cnodes)
    memp = R["MEM"]["peak_after_node"] or {}
    oks["MEM_under_cap"] = bool(memp and memp["peak_commit_gb"] < 8 and memp["peak_working_set_gb"] < 8)

    # --- T0c: stage-0 bracket, RN bound, chi < 1 at every control node
    R['T0c'] = {str(jc): recs[jc].get('t0c') for jc in cnodes}
    oks['T0c'] = all(recs[jc].get('t0c', {}).get('ok', False) for jc in cnodes)
    log(f"T0c (stage-0 bracket, RN, chi < 1) at control nodes {[str(x) for x in cnodes]}: {oks['T0c']}")

    # --- NaN / inf guards: (i) engine guard on absurd coefficients, (ii) finite-record check on every control node
    fn = fns[Fraction(100)]
    N, K = degrees(100)
    base = [la.qfloat(c) for c in fn.pc] + [la.qfloat(c) for q in fn.qcs for c in q]
    g = {}
    for bad in (float('nan'), float('inf'), -float('inf'), 1e300):
        cf = list(base)
        cf[N + 1 + 3] = bad
        r, _ = certify(Fraction(100), log, N=N, K=K, coef_override=cf)
        g[str(bad)] = r['status']
    R['nan_inf_guard'] = dict(engine=g, records_finite={str(jc): finite_record(recs[jc]) for jc in cnodes})
    oks['nan_inf_guard'] = all(v == 'INVALID' for v in g.values()) and all(R['nan_inf_guard']['records_finite'].values())
    log(f"guard: {g}; control records finite: {R['nan_inf_guard']['records_finite']}")

    # --- T0: toy J = 0 refused/trivial with S == 0 (and p(0) = 1/40 exactly)
    r0, fn0 = certify(Fraction(0), log, keep=True, N=36, K=10)
    S0 = la.qparse(r0['S'])
    R['T0'] = dict(status=r0['status'], S=la.qstr(S0), S_is_zero=bool(S0 == 0), p0=la.qstr(fn0.p0),
                   certified=bool(r0['status'] == 'CERTIFIED'))
    oks['T0'] = bool(r0['status'] == 'REFUSED/TRIVIAL' and S0 == 0 and fn0.p0 == fmpq(1, 40))
    expected_false.append('T0.certified (toy J = 0 must be refused as trivial, never certified)')
    log(f"T0 j=0 (U = v1): status {r0['status']}, S == 0: {S0 == 0}, p(0) == 1/40: {fn0.p0 == fmpq(1, 40)}")

    # --- N1: u_+ - 1e-6 j^2 s (U-scale: U_N + (eta - 1e-6) s) at j = 100 and j = 1000 -> certified VIOLATION
    for jc in (Fraction(100), Fraction(1000)):
        eta = la.qparse(recs[jc]['eta'])
        r = la.certify_pair(fns[jc], eta, e_sup=eta - fmpq(1, 10 ** 6), log=log, tag=f"N1 null j={jc}")
        R[f'null_j{jc}'] = r
        oks[f'null_j{jc}'] = bool(r['status'] == 'VIOLATION')

    # --- N4: eta' = 0.25 x needed margin (S/8) on the super side at j = 100 -> FAIL through the bound path
    jc = Fraction(100)
    S = la.qparse(recs[jc]['S'])
    eta = la.qparse(recs[jc]['eta'])
    r4 = la.certify_pair(fns[jc], eta, e_sup=S / 8, log=log, tag='N4 j=100')
    R['N4_j100'] = r4
    oks['N4_j100'] = bool(r4['status'] == 'FAIL' and r4['path'][0].startswith('point checks: sup OK')
                          and any('bound path' in p and 'sup' in p and 'FAILS' in p for p in r4['path']))
    expected_false.append('N4_j100 bound path "sup ... ok" (a quarter margin must FAIL the bound path)')
    rp = la.certify_pair(fns[jc], eta, log=log, tag='pair+ j=100 (true margin)')
    R['pair_positive_j100'] = rp['status']
    oks['pair_positive_j100'] = bool(rp['status'] == 'PASS')

    # --- REF: independent reconstructions
    R['REF_chebyshev_j1000'] = ref_chebyshev(fns[Fraction(1000)], recs[Fraction(1000)], 2, 2026, log)
    R['REF_spherical_j300'] = ref_spherical(fns[Fraction(300)], recs[Fraction(300)], Fraction(300), 377, log)
    R['REF_spherical_j1000'] = ref_spherical(fns[Fraction(1000)], recs[Fraction(1000)], Fraction(1000), 1077, log)
    R['REF_bernstein_j1000'] = ref_bernstein(fns[Fraction(1000)], recs[Fraction(1000)], 4243, log)
    for k in ('REF_chebyshev_j1000', 'REF_spherical_j300', 'REF_spherical_j1000', 'REF_bernstein_j1000'):
        oks[k] = bool(R[k]['ok'])
    del fns, fn, fn0

    # --- ENGINE: the engine's own tests (subprocess; inherits the memcap job)
    t = time.time()
    rE = subprocess.run([sys.executable, '-B', '-m', 'pytest', ENGINE_TESTS, '-q', '-p', 'no:cacheprovider'],
                        capture_output=True, text=True)
    tail = [ln for ln in rE.stdout.strip().splitlines() if ln.strip()][-1:] or ['']
    R['ENGINE'] = dict(rc=rE.returncode, summary=tail[0], seconds=time.time() - t)
    oks['ENGINE'] = bool(rE.returncode == 0 and 'passed' in tail[0] and 'failed' not in tail[0])
    log(f"ENGINE tests: rc={rE.returncode} '{tail[0]}' ({R['ENGINE']['seconds']:.1f}s)")

    # --- RESUME: uninterrupted vs killed-by-PID + resumed vs budget-chunked, byte-identical; changed sha refused
    scratch = os.path.join(SCRATCH, 'resume')
    if os.path.isdir(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)
    me = [sys.executable, os.path.abspath(__file__)]
    oa, ob, oc, od = (os.path.join(scratch, x) for x in ('A.json', 'B.json', 'C.json', 'D.json'))
    t = time.time()
    rA = subprocess.run(me + ['--_grid', oa, '--_n', '6'], capture_output=True, text=True)
    tA = time.time() - t
    pB = subprocess.Popen(me + ['--_grid', ob, '--_n', '6'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    killed_at, killed_pid = None, pB.pid
    while pB.poll() is None:
        try:
            nd = len(json.load(open(ob + '.ckpt'))['done'])
        except Exception:
            nd = 0
        if nd >= 3:
            pB.kill()              # our own child, by its PID (TerminateProcess on this handle)
            pB.wait()
            killed_at = nd
            break
        time.sleep(0.05)
    rB2 = subprocess.run(me + ['--_grid', ob, '--_n', '6'], capture_output=True, text=True)
    invC = 0
    while not os.path.exists(oc) and invC < 12:
        subprocess.run(me + ['--_grid', oc, '--_n', '6', '--_budget', '1'], capture_output=True, text=True)
        invC += 1
    bA = open(oa, 'rb').read() if os.path.exists(oa) else b'A'
    bB = open(ob, 'rb').read() if os.path.exists(ob) else b'B'
    bC = open(oc, 'rb').read() if os.path.exists(oc) else b'C'
    stD = json.load(open(oa + '.ckpt'))
    stD['instrument_sha256'] = '0' * 64
    stD['done'].pop(list(stD['done'].keys())[-1])
    atomic_json(od + '.ckpt', stD)
    rD = subprocess.run(me + ['--_grid', od, '--_n', '6'], capture_output=True, text=True)
    resA = json.loads(bA) if bA != b'A' else {}
    R['RESUME'] = dict(sub_grid=[la.qstr(la.fq(x)) for x in nodes[:6]], uninterrupted_rc=rA.returncode,
                       uninterrupted_seconds=tA, killed_pid=killed_pid, killed_after_nodes=killed_at,
                       resumed_rc=rB2.returncode, B_invocations=json.load(open(ob + '.ckpt'))['invocations'],
                       C_invocations=json.load(open(oc + '.ckpt'))['invocations'] if os.path.exists(oc + '.ckpt') else None,
                       A_eq_B=bool(bA == bB), A_eq_C=bool(bA == bC), sha_changed_rc=rD.returncode,
                       sha_changed_output_exists=os.path.exists(od), A_sha256=hashlib.sha256(bA).hexdigest(),
                       A_all_certified=resA.get('all_nodes_certified'), A_no_gaps=resA.get('no_gaps'),
                       A_expected_false=resA.get('expected_false'), no_pickle='checkpoints are JSON with p/q strings')
    RS = R['RESUME']
    oks['RESUME'] = bool(rA.returncode == 0 and killed_at is not None and killed_at < 6 and RS['B_invocations'] >= 2
                         and bA == bB and bA == bC and (RS['C_invocations'] or 0) >= 2 and rD.returncode == 3
                         and not os.path.exists(od) and RS['A_all_certified'] is True)
    expected_false += ['RESUME.A_no_gaps (6-node sub-grid does not reach j = 1000)',
                       'RESUME.sha_changed_output_exists (changed-sha checkpoint must be refused)']
    log(f"RESUME: A==B (killed own PID {killed_pid} after {killed_at} nodes, resumed) {bA == bB}; A==C (budget-chunked, "
        f"{RS['C_invocations']} invocations) {bA == bC}; changed-sha checkpoint refused rc={rD.returncode}")

    # --- --main refuses without GO-P2b
    if os.path.exists(GO_FILE):
        R['main_refusal'] = 'not exercised: GO-P2b exists (the controls never start --main)'
        oks['main_refusal'] = True
    else:
        mo = os.path.join(scratch, 'M.json')
        rM = subprocess.run(me + ['--main', mo], capture_output=True, text=True)
        R['main_refusal'] = dict(rc=rM.returncode, out=rM.stdout.strip()[-200:], output_exists=os.path.exists(mo))
        oks['main_refusal'] = bool(rM.returncode == 2 and not os.path.exists(mo))
        expected_false.append('main_refusal.output_exists (no output without GO-P2b)')
    log(f"--main without GO-P2b: {R['main_refusal']}")

    R['oks'] = oks
    R['expected_false'] = expected_false
    R['resume_tested'] = bool(oks['RESUME'])
    R['all_ok'] = all(oks.values())
    R['seconds'] = time.time() - T0
    atomic_json(out, R, sort=False)
    log(f"controls done in {R['seconds']:.1f}s; all_ok = {R['all_ok']}; {oks}; instrument sha256 = {sha}")
    return 0 if R['all_ok'] else 1


def run_probe(args, log):
    js = [Fraction(x) for x in args.j.split(',')] if args.j else [Fraction(10), Fraction(100)]
    for j in js:
        N, K = (args.N, args.K) if args.N else degrees(j)
        t = time.time()
        rec, _ = certify(j, log, N=N, K=K)
        log(f"probe j={j} N={N} K={K}: {rec['status']} rel {rec.get('rel_width')} S {rec.get('S_float')} "
            f"total {time.time()-t:.1f}s timing {rec.get('_timing')} peak {peak_memory()}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--probe', action='store_true')
    ap.add_argument('--j')
    ap.add_argument('--N', type=int)
    ap.add_argument('--K', type=int)
    ap.add_argument('--controls')
    ap.add_argument('--main')
    ap.add_argument('--_grid')
    ap.add_argument('--_n', type=int)
    ap.add_argument('--_budget', type=float, default=TIME_BUDGET)
    args = ap.parse_args()
    logf = None

    def log(msg):
        print(msg, flush=True)
        if logf:
            logf.write(msg + '\n')
            logf.flush()

    if la.NUMPY_PREIMPORTED:
        log("WARNING: numpy was imported before the engine; byte-identical fits are not guaranteed")
    if args.controls:
        logf = open(os.path.splitext(args.controls)[0] + '.log', 'w')
        sys.exit(run_controls(args.controls, log))
    elif args.main:
        if not os.path.exists(GO_FILE):
            print("REFUSED: <LAB>/BHSLICE/GO-P2b does not exist (the lead arms it with lab.py go after sealing).")
            sys.exit(2)
        logf = open(os.path.splitext(args.main)[0] + '.log', 'a')
        log(f"--- main invocation {time.strftime('%Y-%m-%d %H:%M:%S')} sha256 {sha_self()} engine {sha_file(ENGINE_FILE)}")
        sys.exit(run_grid(args.main, grid_nodes(), TIME_BUDGET, log))
    elif args._grid:
        nodes = grid_nodes()[:args._n] if args._n else grid_nodes()
        sys.exit(run_grid(args._grid, nodes, args._budget, log))
    elif args.probe:
        run_probe(args, log)
