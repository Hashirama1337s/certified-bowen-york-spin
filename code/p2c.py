"""BHSLICE P2c instrument: certified uniform spin ceiling of the Bowen-York spinning-puncture family (mu > 0, every J).

A thin runner over the shared engine <LAB>/POCKETLAB/gravity_engine/lichnerowicz_axisym.py, VERSION 2 (imported, never
copied, never edited; fit_ops 'auto' = exact-midpoint fit operators above N = 100, cache_ops=False). The engine holds the
whole point certificate (R1-R5), the K1 cell record (arb, outward) and the K2 gap check. This file adds only:
  - the P2c degree schedule (N rising with j, K = 12; BUILD-P2c.md sec 2),
  - the P2c grid: j_{k+1} = (5/4) j_k from j_0 = P2b's last node (exact string from result_p2b.json) to the first node
    >= 10000, plus j = 3000 and j = 10000,
  - the tail X1 at the top node j1: chi(j) <= j1 / (4 delta_-(j1)^2) for every j >= j1 (exact rational),
  - a resumable runner (exact "p/q" strings in JSON checkpoints; no pickle; lab exit codes 0 finished / 3 rerun),
  - a ceiling assembler (--ceiling) that reads result_p1/p2a/p2b and the P2c result BY SHA256 and computes
        c = max(P1 range bound 1/4 on j in [0, 1], every cell chi_max on [1, j1] from P2a/P2b/P2c, X1)
    exactly, with the gap-free chain [0, 1] u [1, j1] u [j1, oo) checked across the four files,
  - the controls (BUILD-P2c.md sec 5).
Design: DESIGN-P2.md sec 2, 4.2-4.4 (X1-X3); REFEREE-P2DESIGN item 4; REFEREE-P2b; P2C-SCOUT. Blind: BAR-P2c.md and
TARGETS-P2.md were not opened by the build.

Usage (every run under memcap):
  python <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 p2c.py --probe --j 10000 [--N 160 --K 12]
  ... p2c.py --controls controls_p2c.json --part node     (then --part rest; each part < 10 min)
  ... p2c.py --main result_p2c.json      (refuses without GO-P2c; resumable; <= 520 s per invocation; exit 0 finished,
                                          3 = rerun the same command to resume, anything else = error)
  ... p2c.py --ceiling ceiling_p2c.json --p2c-sha256 <sha256 of result_p2c.json>     (after --main has finished)
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
from flint import arb, fmpq, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-P2c')
ENGINE_FILE = os.path.join(ENGINE_DIR, 'gravity_engine', 'lichnerowicz_axisym.py')
ENGINE_TESTS = os.path.join(ENGINE_DIR, 'gravity_engine', 'tests')
SCRATCH = os.path.join(HERE, 'seats', 'build-p2c')
TIME_BUDGET = 520.0
EXIT_FINISHED, EXIT_RESUME, EXIT_NO_GO, EXIT_CKPT_REFUSED, EXIT_REFUSED_INPUT, EXIT_CONTROLS_FAILED = 0, 3, 2, 4, 5, 1
assert la.ENGINE_VERSION == 2, 'p2c.py is built on gravity_engine lichnerowicz_axisym version 2'

# ---- inputs, pinned by sha256 (the sealed results this ceiling chains to)
RESULT_P1 = os.path.join(HERE, 'result_p1.json')
RESULT_P2A = os.path.join(HERE, 'result_p2a.json')
RESULT_P2B = os.path.join(HERE, 'result_p2b.json')
PIN_SHA256 = {'p1': '927a72d7f44aa3cbef1454565cb2a37e3bd22970ca8a9538c7d943c4c2ee77f0',
              'p2a': '6fcc3a94833ef9205fe3bca87d64f1ee92a7225b22f643c59370548cbdcacd80',
              'p2b': '9f4fa988dfb66bb6978f7e41fe6b26d32c0c152b908061ffe91076bc154776a3'}
# P2b's last node (VERDICT-P2b: j = 1030.3467644609451531), exact; re-read from result_p2b.json and compared at start.
P2B_LAST_J = ('408162404503791007909934673028597120680242890131500424287099038285814437270363261494520596976526022833642'
              '193972454699977107901/39614081257132168796771975168000000000000000000000000000000000000000000000000000000'
              '0000000000000000000000000000000000000000')

RATIO, J_END = Fraction(5, 4), Fraction(10000)
EXTRA_NODES = (3000, 10000)
GRID_LINE = ("j_{k+1} = (5/4) j_k from j_0 = P2b's last node (1030.3467644609451531, exact) to the first node >= 10000, "
             "plus j = 3000 and j = 10000")

# Degree schedule (upper end of the band in j, N, K). K = 12 everywhere (P2C-SCOUT: K = 12 gives the same S and width as
# K = 16-30 at j = 300, 1000, 8000). N from the probes in BUILD-P2c.md sec 2. N > 100 => engine 'mid' fit operators.
SCHEDULE = ((Fraction(1300), 112, 12), (Fraction(2100), 120, 12), (Fraction(3200), 128, 12), (Fraction(5000), 140, 12),
            (Fraction(8000), 152, 12), (None, 192, 12))


def degrees(j):
    j = Fraction(j)
    for top, N, K in SCHEDULE:
        if top is None or j <= top:
            return N, K


class P2cBowenYork(la.BowenYorkSpinningPuncture):
    """The engine's Bowen-York problem (same source, tau, observables, brackets, toy, cells) with the P2c degree
    schedule, engine-v2 defaults fit_ops='auto' and cache_ops=False (only the current degree's fit operators are held,
    and they are freed before the exact proof part)."""
    EXTRA_NODES = EXTRA_NODES

    def degrees(self, j):
        return degrees(j)


PROBLEM = P2cBowenYork(cache_ops=False)


def p2b_last_j_from_file(path=RESULT_P2B):
    return json.load(open(path))['nodes'][-1]['j']


def j0():
    return Fraction(*map(int, P2B_LAST_J.split('/')))


def grid_nodes():
    return PROBLEM.grid_nodes(j0=j0(), ratio=RATIO, j_end=J_END, extra=EXTRA_NODES)


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


def frac(pq):
    a, b = pq.split('/')
    return Fraction(int(a), int(b))


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
                    peak_commit_gb=round(c.PeakPagefileUsage / 1024 ** 3, 3),
                    commit_now_gb=round(c.PagefileUsage / 1024 ** 3, 3))
    except Exception:
        return None


def certify(j, log, keep=False, N=None, K=None, coef_override=None, problem=None):
    return la.certify_node(problem or PROBLEM, Fraction(j), N=N, K=K, log=log, keep=keep, coef_override=coef_override)


def strip_private(rec):
    """record without '_' keys (fit diagnostics and timing are kept out of the reproducible result)."""
    return {k: v for k, v in rec.items() if not k.startswith('_')}


def fit_summary(rec):
    f = rec.get('_fit') or {}
    h = f.get('arb_hist') or []
    return dict(fit_ops=f.get('fit_ops'), reached_stop=f.get('reached_stop'), iters=len(h),
                last_residual=(h[-1] if h else None), float_residual=f.get('float_residual'))


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


def tail_record(top):
    """X1 at the top node j1 (DESIGN-P2 sec 4.2; REFEREE-P2DESIGN item 4, CORRECT): for every j >= j1,
    chi(j) = j/(4(1+delta)^2) < j/(4 delta^2) <= j1/(4 delta_-(j1)^2), since delta(j)^2 >= (j/j1) delta_-(j1)^2
    (delta(j)/sqrt(j) nondecreasing: DGC Lemma 3.1 + scaling). Exact rational: no rounding needed."""
    j1 = la.qparse(top['j'])
    dlo = la.qparse(top['delta_lo'])
    if not (top.get('status') == 'CERTIFIED' and dlo > 0):
        return dict(ok=False, why='top node not certified or delta_lo <= 0')
    X1 = j1 / (4 * dlo ** 2)
    return dict(ok=True, j1=top['j'], j1_dec=top['j_dec'], delta_lo=top['delta_lo'], X1=la.qstr(X1),
                X1_dec_up=la.dec_ceil(X1, 20),
                statement='for every j >= j1: chi(j) < j1 / (4 delta_lo(j1)^2) = X1 (exact rational)')


def build_result(node_recs, sha, mock=False):
    """engine K1/K2 result + tail X1 + the P2c summary keys."""
    res = PROBLEM.build_result(node_recs, sha, j_first=j0(), j_last=J_END, grid_line=GRID_LINE)
    nodes = res['nodes']
    js = [la.qstr(la.fq(frac(r['j']))) for r in nodes]
    res['grid_complete'] = bool(js == [la.qstr(la.fq(x)) for x in grid_nodes()])
    res['mock_for_controls'] = bool(mock)
    rels = [r['rel_width'] for r in nodes if r['status'] == 'CERTIFIED']
    res['max_node_rel_width'] = max(rels) if rels else None
    res['nodes_finite'] = all(finite_record(r) for r in nodes)
    cert = [r for r in nodes if r['status'] == 'CERTIFIED']
    if cert:
        best = max(cert, key=lambda r: la.qparse(r['chi_lo']))
        res['chi_sup_lower_bound'] = dict(j=best['j_dec'], chi_lo=best['chi_lo'], chi_lo_dec=best['chi_lo_dec'],
                                          note='X3: the supremum of chi over the family is >= this certified value')
    if res['cells']:
        top = max(res['cells'], key=lambda c: la.qparse(c['chi_max']))
        res['max_cell_chi_max'] = dict(j_lo=top['j_lo_dec'], j_hi=top['j_hi_dec'], chi_max=top['chi_max'],
                                       chi_max_dec=top['chi_max_dec'])
    res['tail_X1'] = tail_record(nodes[-1]) if (nodes and res['all_nodes_certified']) else dict(ok=False)
    res['p2b_chain'] = dict(p2b_last_node=P2B_LAST_J, first_node_equals_p2b_last=bool(nodes and nodes[0]['j'] == P2B_LAST_J))
    res['degree_schedule'] = [dict(j_upto=(None if t is None else str(t)), N=N, K=K) for t, N, K in SCHEDULE]
    res['engine'] = dict(file='gravity_engine/lichnerowicz_axisym.py', sha256=sha_file(ENGINE_FILE),
                         version=la.ENGINE_VERSION, fit_ops='auto', cache_ops=False, seal=la.SEAL_ID)
    res['peak_memory_note'] = ('per-node timing and process peak memory are in <out>.runinfo.json (kept out of this file '
                               'so that it is byte-reproducible); controls_p2c.json MEM gives the top-degree node')
    ef = []
    if not res['grid_complete']:
        ef.append('grid_complete (sub-grid or mock run for the controls; the ceiling refuses it)')
    if not res['no_gaps'] and res['all_nodes_certified'] and nodes and frac(nodes[-1]['j']) < J_END:
        ef.append('no_gaps (sub-grid run that does not reach j = 10000; false by design for the RESUME control only)')
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


def pred_scale(N, K):
    """relative node cost (the exact Chebyshev conversion dominates and grows ~N^3; linear in K; P2C-SCOUT (B))."""
    return (N / 160.0) ** 3 * (K / 12.0)


def run_grid(out, nodes, budget, log, degfn=degrees, prior=140.0):
    """Lab exit codes: 0 = finished (result written), 3 = budget reached (rerun the same command to resume),
    4 = checkpoint refused (other instrument/engine sha or node list)."""
    T0 = time.time()
    rs = RunState(out, nodes)
    if not rs.load(log):
        return EXIT_CKPT_REFUSED
    info_p = out + '.runinfo.json'
    info = json.load(open(info_p)) if os.path.exists(info_p) else dict(timing={}, invocations=[])
    worst = prior                       # seconds per unit pred_scale (prior from the probes; replaced by measurements)
    did = 0
    for j in nodes:
        key = la.qstr(la.fq(j))
        if key in rs.st['done']:
            continue
        N, K = degfn(j)
        if did >= 1 and time.time() - T0 + 1.3 * worst * pred_scale(N, K) > budget - 15:
            log(f"time budget reached: {len(rs.st['done'])}/{len(nodes)} nodes done; exit {EXIT_RESUME} = rerun the same "
                f"command to resume")
            info['invocations'].append(dict(ended='budget', seconds=time.time() - T0, peak=peak_memory()))
            atomic_json(info_p, info)
            return EXIT_RESUME
        t = time.time()
        rec, _ = certify(j, log, N=N, K=K)
        dt = time.time() - t
        worst = dt / pred_scale(N, K) if did == 0 else max(worst, dt / pred_scale(N, K))
        info['timing'][key] = dict(j_dec=rec.get('j_dec'), N=N, K=K, total=dt, parts=rec.get('_timing'),
                                   fit=fit_summary(rec), process_peak_after_node=peak_memory())
        rs.st['done'][key] = strip_private(rec)
        rs.save()
        atomic_json(info_p, info)
        did += 1
    res = build_result(list(rs.st['done'].values()), rs.sha)
    atomic_json(out, res)
    info['invocations'].append(dict(ended='finished', seconds=time.time() - T0, peak=peak_memory()))
    atomic_json(info_p, info)
    tx = res.get('tail_X1', {})
    log(f"grid finished: {len(nodes)} nodes, all certified = {res['all_nodes_certified']}, no gaps = {res['no_gaps']}, "
        f"grid complete = {res['grid_complete']}, max node rel width = {res['max_node_rel_width']}, "
        f"max cell chi_max = {res.get('max_cell_chi_max', {}).get('chi_max_dec')}, X1 = {tx.get('X1_dec_up')}, "
        f"invocations = {rs.st['invocations']}, peak memory this invocation = {peak_memory()}")
    return EXIT_FINISHED


# ================================================================== ceiling assembler (X2, X3)
class Refusal(Exception):
    pass


def _check_cells(name, res):
    """K1/K2 consistency of one node+cell result file; returns its cell candidates [(label, chi_max fmpq)]."""
    nodes, cells = res['nodes'], res['cells']
    if not (res.get('all_nodes_certified') is True and res.get('no_gaps') is True):
        raise Refusal(f'{name}: all_nodes_certified / no_gaps not both True')
    if len(cells) != len(nodes) - 1 or len(nodes) < 2:
        raise Refusal(f'{name}: {len(cells)} cells for {len(nodes)} nodes')
    js = [frac(n['j']) for n in nodes]
    if not all(a < b for a, b in zip(js[:-1], js[1:])):
        raise Refusal(f'{name}: nodes not strictly increasing')
    out = []
    ctx.prec = 256
    for i, c in enumerate(cells):
        na, nb = nodes[i], nodes[i + 1]
        if na['status'] != 'CERTIFIED' or nb['status'] != 'CERTIFIED':
            raise Refusal(f'{name}: cell {i} touches a non-certified node')
        if not (c['j_lo'] == na['j'] and c['j_hi'] == nb['j'] and c['delta_lo'] == na['delta_lo']
                and c['delta_hi'] == nb['delta_hi'] and c['nonempty'] is True):
            raise Refusal(f'{name}: cell {i} does not match its node records')
        ja, jb, dla = la.qparse(na['j']), la.qparse(nb['j']), la.qparse(na['delta_lo'])
        if not dla > 0:
            raise Refusal(f'{name}: cell {i} delta_lo <= 0')
        # independent arb evaluation of chi_max = j_b / (4 (1 + sqrt(j_b/j_a) delta_lo(j_a))^2) (REFEREE-P2DESIGN 3a):
        # the stored exact chi_max must be >= the upper end of this ball
        val = arb(jb) / (4 * (1 + (arb(jb) / arb(ja)).sqrt() * arb(dla)) ** 2)
        stored = la.qparse(c['chi_max'])
        mine = la.arb_to_fmpq(val.upper())
        if stored < la.arb_to_fmpq(val.lower()):
            raise Refusal(f'{name}: cell {i} stored chi_max is below the recomputed value (not an upper bound)')
        # candidate = the larger of the stored bound and this recomputation's upper end (both are upper bounds when the
        # stored one is valid; if the stored one were a few ulps low, the recomputed upper end still bounds chi)
        out.append((f"{name} cell [{c['j_lo_dec']}, {c['j_hi_dec']}]", stored if stored >= mine else mine))
    return out


def assemble(paths, p2c_sha256, allow_mock=False):
    """c = max(P1 range bound, every cell chi_max on [1, j1], X1) in exact rationals; X3; gap-free chain. Raises
    Refusal on any sha mismatch or inconsistency."""
    want = dict(PIN_SHA256, p2c=p2c_sha256)
    shas = {}
    for k in ('p1', 'p2a', 'p2b', 'p2c'):
        if not want[k] or not os.path.exists(paths[k]):
            raise Refusal(f'{k}: missing file or sha')
        shas[k] = sha_file(paths[k])
        if shas[k] != want[k]:
            raise Refusal(f'{k}: sha256 {shas[k]} of {paths[k]} != expected {want[k]}')
    P1 = json.load(open(paths['p1']))
    A, B, C = (json.load(open(paths[k])) for k in ('p2a', 'p2b', 'p2c'))
    chain = []
    # ---- P1: J in [0, 1/4], i.e. j in [0, 1]. What it certifies (BUILD-P1 sec 1, VERDICT-P1): for EVERY j in each of
    # 19 slabs, u_- <= u <= u_+ with V = U_-/s > 0 on the square, so u >= u_- >= 0 and delta = 2 lim r u >= 0. Hence
    # chi = j/(4 (1 + delta)^2) <= j/4 <= 1/4 on [0, 1]. Only all_slabs_ok, the slab statuses and the edges are used.
    if not (P1.get('all_slabs_ok') is True and P1.get('finished') is True):
        raise Refusal('p1: all_slabs_ok / finished not True')
    slabs = sorted(P1['slabs'].values(), key=lambda s: s['k'])
    for s in slabs:
        if not (s['ok'] is True and all(s['res'][t]['status'] == 'PASS' for t in ('sup', 'sub', 'pos'))):
            raise Refusal(f"p1: slab {s['k']} not certified (sup/sub/pos)")
    edges_ok = (slabs[0]['j'][0] == 0.0 and slabs[-1]['j'][1] == 1.0
                and all(a['j'][1] == b['j'][0] for a, b in zip(slabs[:-1], slabs[1:])))
    if not edges_ok:
        raise Refusal('p1: slab j-edges do not chain from 0 to 1')
    chain.append(dict(piece='P1', j_from='0', j_to='1', how=f'{len(slabs)} slabs, float edges shared, 0.0 -> 1.0'))
    cands = [('P1 range bound on j in [0, 1]: chi <= j/4 <= 1/4 (delta >= 0 from P1 positivity of U_-/s)', fmpq(1, 4))]
    # ---- P2a, P2b, P2c cells
    for name, R in (('P2a', A), ('P2b', B), ('P2c', C)):
        cands += _check_cells(name, R)
    if not (C.get('instrument_sha256') == sha_self()):
        raise Refusal('p2c: result written by another instrument version')
    if not C.get('grid_complete') is True and not allow_mock:
        raise Refusal('p2c: grid_complete is not True (sub-grid or mock result)')
    # ---- chain: [0,1] (P1) u [1, a_last] (P2a) u [b_first, b_last] (P2b) u [c_first, j1] (P2c), each next first <=
    # previous last (closed intervals), so the union is [0, j1] without gaps.
    fj = lambda R, i: frac(R['nodes'][i]['j'])
    links = [('P1 -> P2a', fj(A, 0) <= 1 and fj(A, 0) >= 0), ('P2a -> P2b', fj(B, 0) <= fj(A, -1)),
             ('P2b -> P2c', fj(C, 0) <= fj(B, -1))]
    for nm, okl in links:
        if not okl:
            raise Refusal(f'chain gap at {nm}')
    if C['nodes'][0]['j'] != B['nodes'][-1]['j']:
        raise Refusal("P2c's first node is not P2b's last node (exact string)")
    for nm, R in (('P2a', A), ('P2b', B), ('P2c', C)):
        chain.append(dict(piece=nm, j_from=R['nodes'][0]['j_dec'], j_to=R['nodes'][-1]['j_dec'],
                          how=f"{len(R['cells'])} K1 cells, consecutive, no_gaps"))
    chain_links = {nm: bool(okl) for nm, okl in links}
    # ---- tail X1 at j1 (recomputed here from the top node record, not read from the file's tail_X1)
    top = C['nodes'][-1]
    tx = tail_record(top)
    if not tx['ok']:
        raise Refusal('tail: ' + tx.get('why', ''))
    if C.get('tail_X1', {}).get('X1') != tx['X1']:
        raise Refusal("tail: the P2c file's stored X1 differs from the recomputation")
    cands.append((f"X1 tail: every j >= j1 = {top['j_dec']}", la.qparse(tx['X1'])))
    chain.append(dict(piece='X1 tail', j_from=top['j_dec'], j_to='infinity', how='delta(j)/sqrt(j) nondecreasing'))
    # ---- c and the binding piece
    c_lab, c = max(cands, key=lambda t: t[1])
    nbind = sum(1 for _, v in cands if v == c)
    X3 = la.qparse(top['chi_lo'])
    allnodes = A['nodes'] + B['nodes'] + C['nodes']
    best = max(allnodes, key=lambda r: la.qparse(r['chi_lo']))
    per_piece = {}
    for lab, v in cands:
        p = lab.split(' ')[0]
        if p not in per_piece or v > per_piece[p][1]:
            per_piece[p] = (lab, v)
    out = dict(
        object='Bowen-York spinning puncture, every bare mass mu > 0 and every J (chi = J/M_ADM^2 depends on j = 4J/mu^2 '
               'only, by scale invariance; J enters squared)',
        inputs={k: dict(path=os.path.basename(paths[k]), sha256=shas[k]) for k in shas},
        instrument_sha256=sha_self(), engine_sha256=sha_file(ENGINE_FILE),
        mock_p2c=bool(C.get('mock_for_controls')),
        c=la.qstr(c), c_dec_up=la.dec_ceil(c, 20), c_float_up=la.qfloat_up(c), c_lt_1=bool(c < 1),
        binding_piece=c_lab, n_pieces_attaining_c=nbind, n_candidates=len(cands),
        per_piece_max={p: dict(label=l, value=la.qstr(v), value_dec_up=la.dec_ceil(v, 20)) for p, (l, v) in
                       per_piece.items()},
        tail_X1=tx,
        X3_top_node=dict(j=top['j_dec'], chi_lo=top['chi_lo'], chi_lo_dec=la.dec_floor(X3, 20),
                         note='sup chi >= this certified member value'),
        X3_max_over_nodes=dict(j=best['j_dec'], chi_lo_dec=best['chi_lo_dec']),
        c_minus_X3=la.qstr(c - X3), c_minus_X3_dec_up=la.dec_ceil(c - X3, 12),
        chain=chain, chain_links=chain_links, chain_ok=True,
        statement=('For every Bowen-York spinning puncture with mu > 0: J/M_ADM^2 <= c (exact rational above). '
                   'Covers [0, 1] by P1 (chi <= j/4), [1, j1] by the K1 cells of P2a, P2b, P2c, and [j1, oo) by X1. '
                   'Rests on: the comparison lemma C1 (refereed, P1), DGC Lemma 3.1 + scaling (delta(j)/sqrt(j) '
                   'nondecreasing; REFEREE-P2DESIGN), the exact engine certificates. mu = 0 is not covered.'))
    return out


def run_ceiling(args, log):
    paths = dict(p1=args.p1 or RESULT_P1, p2a=args.p2a or RESULT_P2A, p2b=args.p2b or RESULT_P2B,
                 p2c=args.p2c or os.path.join(HERE, 'result_p2c.json'))
    try:
        out = assemble(paths, args.p2c_sha256, allow_mock=args._allow_mock)
    except Refusal as e:
        log(f'REFUSED: {e}')
        return EXIT_REFUSED_INPUT
    atomic_json(args.ceiling, out)
    log(f"ceiling: c = {out['c_dec_up']} (c < 1: {out['c_lt_1']}), binding = {out['binding_piece']}; "
        f"X1 = {out['tail_X1']['X1_dec_up']}; X3(top) = {out['X3_top_node']['chi_lo_dec']}; "
        f"c - X3 = {out['c_minus_X3_dec_up']}; mock = {out['mock_p2c']}")
    return EXIT_FINISHED


# ================================================================== independent reference checks (REF)
# ref_chebyshev / ref_spherical / ref_bernstein follow BHSLICE p2b.py (project code, not engine; read-only original),
# with the spherical check's points extended toward the throat for j = 1e4.
def _cheb_int(n, x):
    """T_k(x) q^k as python ints for x = p/q (own recurrence, k = 0..n)."""
    p, q = int(x.p), int(x.q)
    a = [1, p]
    for k in range(1, n):
        a.append(2 * p * a[k] - q * q * a[k - 1])
    return a[:n + 1], q


def ref_chebyshev(fn, rec, npts, seed, log):
    """exact Chebyshev reconstruction of Q at random rational points vs Q straight from the ansatz terms. S is
    recomputed from scratch by a different route than the engine's residual_S (which sums |cheb(A) + tau cheb(src)|
    with A = D^7 Lh u): here Q = D^7 Lh u / tau + src is formed first and converted to Chebyshev ONCE; sum |C| * scale
    must equal the record's S exactly."""
    t = time.time()
    W = fn.W
    srcP = la.pack(fn.problem.src_terms, W)
    DP = la.pack([(la.ONE, la.ONE)] + [(fn.basis.D_factor * f, g) for f, g in fn.Ut], W)
    Q = DP ** fn.P * la.pack(fn.Lt, W) / fn.tau + srcP
    C, sc = la.chebmat(Q, W)
    rows, cols = C.nrows(), C.ncols()
    Ci = [[int(C[i, k]) for k in range(cols)] for i in range(rows)]
    S = fmpq(sum(abs(v) for r in Ci for v in r)) * sc
    out = dict(S_matches_record=bool(la.qstr(S) == rec['S']))
    out['abs_sum_equals_S'] = out['S_matches_record']
    out['Q_degrees'] = [rows - 1, cols - 1]
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
    log(f"REF Chebyshev j={rec['j_dec'][:10]}: S recomputed == record {out['S_matches_record']}, sum|C| sc == S "
        f"{out['abs_sum_equals_S']}, {sum(p['equal'] for p in pts)}/{npts} points series == direct, all |Q| <= S: "
        f"{all(p['within_S'] for p in pts)} ({out['seconds']:.1f}s)")
    return out


def ref_spherical(fn, rec, j, seed, log):
    """residual signs from the flat Laplacian in spherical coordinates (mpmath, numerical differentiation; no use of
    the engine's operator Lh or its exact polynomials): E[u] = u_rr + 2u_r/r + (u_thth + cot u_th)/r^2
    + (9/4) J^2 sin^2 r^-6 Phi^-7 must be < 0 for u_+ and > 0 for u_-; modes l >= 4 must be present."""
    import mpmath as mp
    t = time.time()
    mp.mp.dps = 60
    K = len(fn.qcs)
    lmodes = sum(1 for k in range(1, K) if any(c != 0 for c in fn.qcs[k]))

    def mq(q):
        return mp.mpf(int(q.p)) / int(q.q)
    pcm = [mq(c) for c in fn.pc]
    qcm = [[mq(c) for c in q] for q in fn.qcs]
    jF = Fraction(j)
    taum = mp.mpf(jF.numerator) ** 2 / mp.mpf(jF.denominator) ** 2
    etam = mq(la.qparse(rec['eta']))
    Sm = mq(la.qparse(rec['S']))
    J = mp.mpf(jF.numerator) / jF.denominator / 4

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
    pts, teeth_hits, qteeth_hits = [], [], []
    for _ in range(6):
        r = mp.mpf(rng.choice([1, 2, 3, 5, 7])) / rng.choice([1, 2, 4, 8, 16, 64]) * rng.choice([1, 1, 3, 10])
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
        srcterm = src * (1 + 1 / (2 * r) + u0 + taum * etam * s) ** -7
        bite = None
        for eps in range(12, 0, -1):
            if E(+1, 1 + mp.mpf(10) ** -eps) >= 0:
                bite = eps
                break
        teeth_hits.append(bite)
        # sharper teeth (REFEREE-P2b item 5): the |Q| <= S branch. Smallest relative source error 10^-eps that makes
        # the reconstructed |Q| exceed S at this point.
        qbite = None
        for eps in range(16, 0, -1):
            if abs(E(0, 1 + mp.mpf(10) ** -eps) / fac * (1 + (1 - s) * u0) ** 7) > Sm:
                qbite = eps
                break
        qteeth_hits.append(qbite)
        pts.append(dict(r=float(r), theta=float(th), E_plus=float(Ep), E_minus=float(Em),
                        margin_over_source_term=float(-Ep / srcterm), Q_over_S=float(Q0 / Sm),
                        signs_ok=bool(Ep < 0 < Em), Q_within_S=bool(abs(Q0) <= Sm)))
    out = dict(points=pts, l_ge_4_modes_present=lmodes, teeth_largest_eps_that_bites_per_point=teeth_hits,
               Q_teeth_largest_eps_that_bites_per_point=qteeth_hits,
               teeth_ok=any(b is not None for b in teeth_hits) and any(b is not None for b in qteeth_hits),
               seconds=time.time() - t)
    out['ok'] = bool(lmodes >= 1 and all(p['signs_ok'] and p['Q_within_S'] for p in pts) and out['teeth_ok'])
    log(f"REF spherical j={rec['j_dec'][:10]}: signs E(u+) < 0 < E(u-) at {sum(p['signs_ok'] for p in pts)}/6 points, "
        f"|Q| <= S {all(p['Q_within_S'] for p in pts)} (max |Q|/S {max(abs(p['Q_over_S']) for p in pts):.3e}), "
        f"modes l >= 4 nonzero: {lmodes}, sign teeth (10^-eps source error flips E(u+)) {teeth_hits}, |Q| <= S teeth "
        f"{qteeth_hits} ({out['seconds']:.1f}s)")
    return out


def ref_bernstein(fn, rec, seed, log):
    """Bernstein reconstruction at random rational points on [0,1] and on a sub-piece (exact)."""
    t = time.time()
    eta = la.qparse(rec['eta'])
    T = fn.V_terms(eta)
    rng = random.Random(seed)
    res = []
    for a, b in ((fmpq(0), fmpq(1)), (fmpq(1, 4), fmpq(3, 8))):
        B, n, m = la.bernstein_matrix(T, a, b)
        for _ in range(3):
            tt = fmpq(rng.randint(0, 97), 97)
            w = fmpq(rng.randint(0, 89), 89)
            bt = [math.comb(n, i) * tt ** i * (1 - tt) ** (n - i) for i in range(n + 1)]
            bw = [math.comb(m, k) * w ** k * (1 - w) ** (m - k) for k in range(m + 1)]
            lhs = sum((B[i, k] * bt[i] * bw[k] for i in range(n + 1) for k in range(m + 1)), fmpq(0))
            s0, c0 = a + (b - a) * tt, 2 * w - 1
            rhs = sum((f(s0) * g(c0) for f, g in T), fmpq(0))
            res.append(bool(lhs == rhs))
    out = dict(points=len(res), all_equal=all(res), degrees=[n, m], ok=all(res) and m == 2 * len(fn.qcs),
               seconds=time.time() - t)
    log(f"REF Bernstein j={rec['j_dec'][:10]}: {sum(res)}/{len(res)} reconstructions exact (n = {n}, m = {m}) "
        f"({out['seconds']:.1f}s)")
    return out


# ================================================================== independent mpmath recomputations (TAIL, CEIL)
def mp_tail(top, dps=120):
    import mpmath as mp
    mp.mp.dps = dps
    q = lambda s: mp.mpf(int(s.split('/')[0])) / mp.mpf(int(s.split('/')[1]))
    return q(top['j']) / (4 * q(top['delta_lo']) ** 2)


def mp_ceiling(files_nodes, top, dps=120):
    """c from the stored NODE records only (cells re-derived here): max(1/4, chi_max of every consecutive node pair,
    j1/(4 delta_lo(j1)^2)); mpmath at dps digits."""
    import mpmath as mp
    mp.mp.dps = dps
    q = lambda s: mp.mpf(int(s.split('/')[0])) / mp.mpf(int(s.split('/')[1]))
    best = (mp.mpf(1) / 4, 'P1')
    ncell = 0
    for name, nodes in files_nodes:
        for a, b in zip(nodes[:-1], nodes[1:]):
            ja, jb, dla = q(a['j']), q(b['j']), q(a['delta_lo'])
            v = jb / (4 * (1 + mp.sqrt(jb / ja) * dla) ** 2)
            ncell += 1
            if v > best[0]:
                best = (v, f'{name} cell [{a["j_dec"]}, {b["j_dec"]}]')
    x1 = mp_tail(top, dps)
    if x1 > best[0]:
        best = (x1, 'X1')
    return best, ncell


# ================================================================== controls
def _load_controls(out, sha):
    if os.path.exists(out):
        try:
            R = json.load(open(out))
            if R.get('instrument_sha256') == sha and R.get('engine_sha256') == sha_file(ENGINE_FILE):
                return R
        except Exception:
            pass
    return dict(instrument_sha256=sha, engine_sha256=sha_file(ENGINE_FILE), parts={}, oks={}, expected_false=[])


def _finish_controls(out, R, log):
    oks = R['oks']
    need = {'node', 'rest'}
    R['parts_done'] = sorted(R['parts'].keys())
    R['resume_tested'] = bool(oks.get('RESUME', False))
    R['all_ok'] = bool(need <= set(R['parts'].keys()) and all(oks.values()))
    R['degree_schedule'] = [dict(j_upto=(None if t is None else str(t)), N=N, K=K) for t, N, K in SCHEDULE]
    nodes = grid_nodes()
    R['grid'] = dict(line=GRID_LINE, n_nodes=len(nodes), first=la.qstr(la.fq(nodes[0])),
                     nodes_dec=[la.dec_floor(la.fq(x), 12) for x in nodes],
                     degrees=[list(degrees(x)) for x in nodes])
    atomic_json(out, R, sort=False)
    log(f"controls parts {R['parts_done']}; all_ok = {R['all_ok']}; oks = {oks}; instrument sha256 = {R['instrument_sha256']}")
    return 0 if all(oks.values()) else EXIT_CONTROLS_FAILED


def controls_node(R, log):
    """top-degree node j = 10000: MEM, REF (Chebyshev with S recomputed, spherical, Bernstein), N1, N4, T0c, finite."""
    T0 = time.time()
    oks, ef = R['oks'], R['expected_false']
    P = R['parts'].setdefault('node', {})
    jc = Fraction(10000)
    N, K = degrees(jc)
    base = peak_memory()
    t = time.time()
    rec, fn = certify(jc, log, keep=True)
    sec = time.time() - t
    P['MEM'] = dict(node='10000/1', degrees=[N, K], baseline_before=base, peak_after_node=peak_memory(), seconds=sec,
                    cap_gb=8, timing_parts=rec.get('_timing'), fit=fit_summary(rec))
    memp = P['MEM']['peak_after_node'] or {}
    oks['MEM_under_cap'] = bool(memp and memp['peak_commit_gb'] < 8 and memp['peak_working_set_gb'] < 8)
    log(f"MEM: top-degree node j=10000 (N, K) = ({N}, {K}): process peak {memp} (baseline {base}), {sec:.1f}s, "
        f"fit {fit_summary(rec)}")
    P['control_node_j10000'] = strip_private(rec)
    oks['control_node_certified'] = bool(rec['status'] == 'CERTIFIED')
    oks['fit_reached_stop'] = bool(fit_summary(rec)['reached_stop'])
    oks['T0c_j10000'] = bool(rec.get('t0c', {}).get('ok', False))
    oks['finite_j10000'] = finite_record(rec)
    if rec['status'] != 'CERTIFIED':
        log('control node did not certify; remaining node controls skipped')
        return
    S, eta = la.qparse(rec['S']), la.qparse(rec['eta'])
    # REF Chebyshev FIRST: it recomputes S from scratch and checks it equals the record's S; the pair checks below then
    # reuse that exact S (memoised on this ExactFn) instead of a third exact pass at the top degree.
    P['REF_chebyshev_j10000'] = ref_chebyshev(fn, rec, 2, 2027, log)
    oks['REF_chebyshev_j10000'] = bool(P['REF_chebyshev_j10000']['ok'])
    if P['REF_chebyshev_j10000']['S_matches_record']:
        fn.residual_S = lambda want_Q=False, _S=S: (_S, None)
        P['S_memo'] = 'pair checks reuse the exact S (recomputed and matched by REF_chebyshev_j10000)'
    # N1: u_+ - 1e-6 j^2 s, i.e. U-scale e_sup = eta - 1e-6 -> certified VIOLATION
    r1 = la.certify_pair(fn, eta, e_sup=eta - fmpq(1, 10 ** 6), log=log, tag='N1 null j=10000')
    P['N1_j10000'] = r1
    oks['N1_j10000_violation'] = bool(r1['status'] == 'VIOLATION')
    # N4: near-threshold super margin S/2 (1 - 2^-10) at the top degree -> FAIL through the bound path (point checks OK)
    e4 = S / 2 * (1 - fmpq(1, 2 ** 10))
    r4 = la.certify_pair(fn, eta, e_sup=e4, log=log, tag='N4 j=10000 (margin S/2 (1 - 2^-10))')
    P['N4_j10000'] = r4
    oks['N4_j10000_fail_bound_path'] = bool(r4['status'] == 'FAIL' and r4['path'][0].startswith('point checks: sup OK')
                                            and any('bound path' in p and 'sup' in p and 'FAILS' in p for p in r4['path']))
    rp = la.certify_pair(fn, eta, log=log, tag='pair+ j=10000 (true margin)')
    P['pair_positive_j10000'] = rp
    oks['pair_positive_j10000'] = bool(rp['status'] == 'PASS')
    P['REF_spherical_j10000'] = ref_spherical(fn, rec, jc, 10077, log)
    oks['REF_spherical_j10000'] = bool(P['REF_spherical_j10000']['ok'])
    P['REF_bernstein_j10000'] = ref_bernstein(fn, rec, 4244, log)
    oks['REF_bernstein_j10000'] = bool(P['REF_bernstein_j10000']['ok'])
    ef.append('parts.node.N4_j10000.status == FAIL (a margin S/2 (1 - 2^-10) must NOT pass the bound path)')
    ef.append('parts.node.N1_j10000.status == VIOLATION (u_+ - 1e-6 j^2 s must be caught)')
    P['seconds'] = time.time() - T0


def controls_rest(R, log):
    """ENGINE tests, T0, NaN/inf, CHAIN (+ sha refusal), TAIL, CEIL (mock P2c), RESUME (+ exit codes), --main refusal."""
    T0 = time.time()
    oks, ef = R['oks'], R['expected_false']
    P = R['parts'].setdefault('rest', {})
    nodeP = R['parts'].get('node', {})
    me = [sys.executable, os.path.abspath(__file__)]
    scratch = os.path.join(SCRATCH, 'controls')
    if os.path.isdir(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)

    # --- T0: toy J = 0 refused/trivial with S == 0 (and p(0) = 1/40 exactly)
    r0, fn0 = certify(Fraction(0), log, keep=True, N=36, K=10)
    S0 = la.qparse(r0['S'])
    P['T0'] = dict(status=r0['status'], S=la.qstr(S0), S_is_zero=bool(S0 == 0), p0=la.qstr(fn0.p0),
                   certified=bool(r0['status'] == 'CERTIFIED'))
    oks['T0'] = bool(r0['status'] == 'REFUSED/TRIVIAL' and S0 == 0 and fn0.p0 == fmpq(1, 40))
    ef.append('parts.rest.T0.certified (toy J = 0 must be refused as trivial, never certified)')
    log(f"T0 j=0 (U = v1): status {r0['status']}, S == 0: {S0 == 0}, p(0) == 1/40: {fn0.p0 == fmpq(1, 40)}")

    # --- NaN / inf guard at the top degree (the guard runs before any exact work)
    N, K = degrees(10000)
    g = {}
    for bad in (float('nan'), float('inf'), -float('inf'), 1e300):
        cf = [0.0] * ((N + 1) * (K + 1))
        cf[0] = 0.5
        cf[N + 1 + 3] = bad
        r, _ = certify(Fraction(10000), log, N=N, K=K, coef_override=cf)
        g[str(bad)] = r['status']
    P['nan_inf_guard'] = g
    oks['nan_inf_guard'] = all(v == 'INVALID' for v in g.values())
    log(f"NaN/inf guard at (N, K) = ({N}, {K}): {g}")

    # --- CHAIN: first grid node == P2b's last node (exact string, file and constant); sha refusal on a modified copy
    nodes = grid_nodes()
    f_last = p2b_last_j_from_file()
    P['CHAIN'] = dict(first_node=la.qstr(la.fq(nodes[0])), p2b_last_from_file=f_last, pinned=P2B_LAST_J,
                      equal=bool(la.qstr(la.fq(nodes[0])) == f_last == P2B_LAST_J),
                      p2b_sha_matches_pin=bool(sha_file(RESULT_P2B) == PIN_SHA256['p2b']),
                      top_node_dec=float(nodes[-1]), n_nodes=len(nodes))
    # modified copy of result_p2b.json (one byte of whitespace added): the ceiling must refuse it (exit 5)
    mod = os.path.join(scratch, 'result_p2b_modified.json')
    with open(RESULT_P2B, 'rb') as fsrc:
        data = fsrc.read()
    with open(mod, 'wb') as fdst:
        fdst.write(data + b' ')
    rS = subprocess.run(me + ['--ceiling', os.path.join(scratch, 'ceil_mod.json'), '--p2b', mod,
                              '--p2c-sha256', '0' * 64], capture_output=True, text=True)
    P['CHAIN']['modified_p2b_rc'] = rS.returncode
    P['CHAIN']['modified_p2b_msg'] = rS.stdout.strip()[-240:]
    P['CHAIN']['modified_p2b_output_exists'] = os.path.exists(os.path.join(scratch, 'ceil_mod.json'))
    oks['CHAIN'] = bool(P['CHAIN']['equal'] and P['CHAIN']['p2b_sha_matches_pin'] and rS.returncode == EXIT_REFUSED_INPUT
                        and 'p2b: sha256' in rS.stdout and not P['CHAIN']['modified_p2b_output_exists'])
    ef.append('parts.rest.CHAIN.modified_p2b_output_exists (a modified input copy must be refused)')
    log(f"CHAIN: first node == P2b last (file, pinned) {P['CHAIN']['equal']}; modified p2b copy refused rc={rS.returncode}")

    # --- mock P2c result (controls only): P2b's last node record + the j = 10000 control node (part 'node')
    rec10k = nodeP.get('control_node_j10000')
    if rec10k is None or rec10k.get('status') != 'CERTIFIED':
        oks['TAIL'] = oks['CEIL'] = False
        log('TAIL/CEIL: no certified j = 10000 record from part node; run --part node first')
    else:
        lastB = json.load(open(RESULT_P2B))['nodes'][-1]
        mock = build_result([lastB, rec10k], sha_self(), mock=True)
        mock_p = os.path.join(scratch, 'result_p2c_mock.json')
        atomic_json(mock_p, mock)
        msha = sha_file(mock_p)
        # TAIL: X1 recomputed independently in mpmath at 120 digits
        import mpmath as mp
        x1mp = mp_tail(rec10k, 120)
        tx = mock['tail_X1']
        x1ex = mp.mpf(int(tx['X1'].split('/')[0])) / int(tx['X1'].split('/')[1])
        P['TAIL'] = dict(X1_exact=tx['X1'], X1_dec_up=tx['X1_dec_up'], X1_mpmath_120=mp.nstr(x1mp, 40),
                         rel_diff=float(abs(x1ex - x1mp) / x1mp), dec_up_ge_mp=bool(mp.mpf(tx['X1_dec_up']) >= x1mp))
        oks['TAIL'] = bool(P['TAIL']['rel_diff'] < 1e-100 and P['TAIL']['dec_up_ge_mp'])
        log(f"TAIL: X1(j1 = 10000, mock top) exact {tx['X1_dec_up']} vs mpmath-120 {mp.nstr(x1mp, 25)}: rel diff "
            f"{P['TAIL']['rel_diff']:.1e}")
        # CEIL: assembler (CLI, by sha, mock allowed) vs independent mpmath recomputation from all stored node records
        ce_p = os.path.join(scratch, 'ceiling_mock.json')
        rC = subprocess.run(me + ['--ceiling', ce_p, '--p2c', mock_p, '--p2c-sha256', msha, '--_allow_mock'],
                            capture_output=True, text=True)
        rW = subprocess.run(me + ['--ceiling', os.path.join(scratch, 'ceil_wrongsha.json'), '--p2c', mock_p,
                                  '--p2c-sha256', 'f' * 64, '--_allow_mock'], capture_output=True, text=True)
        rG = subprocess.run(me + ['--ceiling', os.path.join(scratch, 'ceil_nomock.json'), '--p2c', mock_p,
                                  '--p2c-sha256', msha], capture_output=True, text=True)
        ok_cli = rC.returncode == 0 and os.path.exists(ce_p)
        P['CEIL'] = dict(cli_rc=rC.returncode, cli_out=rC.stdout.strip()[-300:], wrong_p2c_sha_rc=rW.returncode,
                         incomplete_grid_without_flag_rc=rG.returncode, mock_sha256=msha)
        if ok_cli:
            ce = json.load(open(ce_p))
            A, B = json.load(open(RESULT_P2A)), json.load(open(RESULT_P2B))
            (cmp_, lab), ncell = mp_ceiling([('P2a', A['nodes']), ('P2b', B['nodes']), ('P2c', mock['nodes'])],
                                            rec10k, 120)
            cex = mp.mpf(int(ce['c'].split('/')[0])) / int(ce['c'].split('/')[1])
            P['CEIL'].update(c_assembler=ce['c_dec_up'], binding_assembler=ce['binding_piece'],
                             c_mpmath_120=mp.nstr(cmp_, 30), binding_mpmath=lab, cells_recomputed=ncell,
                             n_candidates_assembler=ce['n_candidates'], rel_diff=float(abs(cex - cmp_) / cmp_),
                             assembler_ge_mp=bool(cex >= cmp_ * (1 - mp.mpf(10) ** -100)),
                             X3_top=ce['X3_top_node']['chi_lo_dec'], c_minus_X3=ce['c_minus_X3_dec_up'],
                             chain=ce['chain'], chain_links=ce['chain_links'], mock_flag_in_output=ce['mock_p2c'])
            same_bind = (lab == 'X1' and ce['binding_piece'].startswith('X1')) or (lab != 'X1' and lab in ce['binding_piece'])
            oks['CEIL'] = bool(P['CEIL']['rel_diff'] < 1e-60 and P['CEIL']['assembler_ge_mp'] and same_bind
                               and ncell + 2 == ce['n_candidates'] and ce['mock_p2c'] is True
                               and rW.returncode == EXIT_REFUSED_INPUT and rG.returncode == EXIT_REFUSED_INPUT)
            log(f"CEIL (mock P2c [1030.35, 10000]): assembler c = {ce['c_dec_up']} ({ce['binding_piece']}) vs mpmath-120 "
                f"from {ncell} node pairs {mp.nstr(cmp_, 22)} ({lab}); rel diff {P['CEIL']['rel_diff']:.1e}; wrong P2c "
                f"sha rc={rW.returncode}; mock without flag rc={rG.returncode}")
        else:
            oks['CEIL'] = False
            log(f"CEIL: assembler CLI failed rc={rC.returncode}: {rC.stdout[-400:]} {rC.stderr[-400:]}")
        ef.append('parts.rest.CEIL.mock_flag_in_output is True and the mock c is not the P2c value (mock grid = two '
                  'nodes 1030.35 and 10000; one coarse cell)')

    # --- ENGINE: the engine's own tests (subprocess; inherits the memcap job)
    t = time.time()
    rE = subprocess.run([sys.executable, '-B', '-m', 'pytest', ENGINE_TESTS, '-q', '-p', 'no:cacheprovider'],
                        capture_output=True, text=True)
    tail = [ln for ln in rE.stdout.strip().splitlines() if ln.strip()][-1:] or ['']
    P['ENGINE'] = dict(rc=rE.returncode, summary=tail[0], seconds=time.time() - t)
    oks['ENGINE'] = bool(rE.returncode == 0 and 'passed' in tail[0] and 'failed' not in tail[0])
    log(f"ENGINE tests: rc={rE.returncode} '{tail[0]}' ({P['ENGINE']['seconds']:.1f}s)")

    # --- RESUME on a cheap test grid (same runner path; degrees fixed at (36, 10)):
    # A uninterrupted; B killed by its own PID after 3 nodes, then resumed; C budget-chunked (exit 3 ... 3, 0);
    # D checkpoint with a changed instrument sha -> refused (exit 4), no output. A == B == C byte for byte.
    rd = os.path.join(scratch, 'resume')
    os.makedirs(rd)
    oa, ob, oc, od = (os.path.join(rd, x) for x in ('A.json', 'B.json', 'C.json', 'D.json'))
    t = time.time()
    rA = subprocess.run(me + ['--_grid', oa], capture_output=True, text=True)
    tA = time.time() - t
    pB = subprocess.Popen(me + ['--_grid', ob], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
    rB2 = subprocess.run(me + ['--_grid', ob], capture_output=True, text=True)
    codesC = []
    while not os.path.exists(oc) and len(codesC) < 12:
        rc = subprocess.run(me + ['--_grid', oc, '--_budget', '1'], capture_output=True, text=True).returncode
        codesC.append(rc)
    bA = open(oa, 'rb').read() if os.path.exists(oa) else b'A'
    bB = open(ob, 'rb').read() if os.path.exists(ob) else b'B'
    bC = open(oc, 'rb').read() if os.path.exists(oc) else b'C'
    stD = json.load(open(oa + '.ckpt'))
    stD['instrument_sha256'] = '0' * 64
    stD['done'].pop(list(stD['done'].keys())[-1])
    atomic_json(od + '.ckpt', stD)
    rD = subprocess.run(me + ['--_grid', od], capture_output=True, text=True)
    rA2 = subprocess.run(me + ['--_grid', oa], capture_output=True, text=True)     # rerun after finish: idempotent
    bA2 = open(oa, 'rb').read()
    resA = json.loads(bA) if bA != b'A' else {}
    P['RESUME'] = dict(test_grid=[la.qstr(la.fq(x)) for x in test_grid()], uninterrupted_rc=rA.returncode,
                       uninterrupted_seconds=tA, killed_pid=killed_pid, killed_after_nodes=killed_at,
                       resumed_rc=rB2.returncode, B_invocations=json.load(open(ob + '.ckpt'))['invocations'],
                       C_exit_codes=codesC, A_eq_B=bool(bA == bB), A_eq_C=bool(bA == bC),
                       rerun_after_finish_rc=rA2.returncode, rerun_identical=bool(bA == bA2),
                       sha_changed_rc=rD.returncode, sha_changed_output_exists=os.path.exists(od),
                       A_sha256=hashlib.sha256(bA).hexdigest(), A_all_certified=resA.get('all_nodes_certified'),
                       A_no_gaps=resA.get('no_gaps'), no_pickle='checkpoints are JSON with p/q strings')
    RS = P['RESUME']
    oks['RESUME'] = bool(rA.returncode == 0 and killed_at is not None and killed_at < 6 and RS['B_invocations'] >= 2
                         and rB2.returncode == 0 and bA == bB and bA == bC and len(codesC) >= 2
                         and all(c == EXIT_RESUME for c in codesC[:-1]) and codesC[-1] == EXIT_FINISHED
                         and rD.returncode == EXIT_CKPT_REFUSED and not os.path.exists(od)
                         and rA2.returncode == 0 and bA == bA2 and RS['A_all_certified'] is True)
    ef.append('parts.rest.RESUME.sha_changed_output_exists (changed-sha checkpoint must be refused)')
    log(f"RESUME: A==B (killed own PID {killed_pid} after {killed_at} nodes, resumed rc {rB2.returncode}) {bA == bB}; "
        f"A==C (budget-chunked, exit codes {codesC}) {bA == bC}; changed-sha checkpoint rc={rD.returncode}; "
        f"rerun after finish rc={rA2.returncode} identical {bA == bA2}")

    # --- --main refuses without GO-P2c
    if os.path.exists(GO_FILE):
        P['main_refusal'] = 'not exercised: GO-P2c exists (the controls never start --main)'
        oks['main_refusal'] = True
    else:
        mo = os.path.join(scratch, 'M.json')
        rM = subprocess.run(me + ['--main', mo], capture_output=True, text=True)
        P['main_refusal'] = dict(rc=rM.returncode, out=rM.stdout.strip()[-200:], output_exists=os.path.exists(mo))
        oks['main_refusal'] = bool(rM.returncode == EXIT_NO_GO and not os.path.exists(mo))
        ef.append('parts.rest.main_refusal.output_exists (no output without GO-P2c)')
    log(f"--main without GO-P2c: {P['main_refusal']}")
    P['seconds'] = time.time() - T0


def test_grid():
    """cheap RESUME grid: 6 nodes j = 10 (5/4)^k at fixed degrees (36, 10); exercises the same runner path."""
    return [Fraction(10) * Fraction(5, 4) ** k for k in range(6)]


def run_controls(out, part, log):
    sha = sha_self()
    R = _load_controls(out, sha)
    R['oks'] = {k: v for k, v in R['oks'].items()}
    # a part rerun replaces that part's entries
    prefix = {'node': ('MEM', 'control_node', 'fit_reached', 'T0c_j', 'finite_j', 'REF_', 'N1_', 'N4_', 'pair_'),
              'rest': ('T0', 'nan_inf', 'CHAIN', 'TAIL', 'CEIL', 'ENGINE', 'RESUME', 'main_refusal')}
    parts = ['node', 'rest'] if part == 'all' else [part]
    for p in parts:
        R['oks'] = {k: v for k, v in R['oks'].items() if not (k.startswith(prefix[p]) and not
                                                               (p == 'rest' and k.startswith('T0c_j')))}
        R['expected_false'] = [e for e in R['expected_false'] if not e.startswith(f'parts.{p}.')]
        R['parts'].pop(p, None)
        R.setdefault('started', {})[p] = time.strftime('%Y-%m-%d %H:%M:%S')
        (controls_node if p == 'node' else controls_rest)(R, log)
        atomic_json(out, R, sort=False)
    return _finish_controls(out, R, log)


def run_probe(args, log):
    gn = grid_nodes()
    js = ([gn[int(x[4:])] if x.startswith('node') else Fraction(x) for x in args.j.split(',')] if args.j
          else [Fraction(3000)])
    for j in js:
        N, K = (args.N, args.K or 12) if args.N else degrees(j)
        t = time.time()
        rec, _ = certify(j, log, N=N, K=K)
        log(f"probe j={float(j):.4f} N={N} K={K}: {rec['status']} rel {rec.get('rel_width')} S {rec.get('S_float')} "
            f"Bernstein {rec.get('positivity', {}).get('pieces')} total {time.time()-t:.1f}s timing {rec.get('_timing')} "
            f"fit {fit_summary(rec)} peak {peak_memory()}")


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--probe', action='store_true')
    ap.add_argument('--j')
    ap.add_argument('--N', type=int)
    ap.add_argument('--K', type=int)
    ap.add_argument('--controls')
    ap.add_argument('--part', default='all', choices=['node', 'rest', 'all'])
    ap.add_argument('--main')
    ap.add_argument('--ceiling')
    ap.add_argument('--p2c-sha256', dest='p2c_sha256')
    ap.add_argument('--p1')
    ap.add_argument('--p2a')
    ap.add_argument('--p2b')
    ap.add_argument('--p2c')
    ap.add_argument('--_allow_mock', action='store_true')
    ap.add_argument('--_grid')
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
        logf = open(os.path.splitext(args.controls)[0] + '.log', 'a')
        log(f"--- controls part {args.part} {time.strftime('%Y-%m-%d %H:%M:%S')} sha256 {sha_self()}")
        sys.exit(run_controls(args.controls, args.part, log))
    elif args.main:
        if not os.path.exists(GO_FILE):
            print("REFUSED: <LAB>/BHSLICE/GO-P2c does not exist (the lead arms it with lab.py go after sealing).")
            sys.exit(EXIT_NO_GO)
        if p2b_last_j_from_file() != P2B_LAST_J or sha_file(RESULT_P2B) != PIN_SHA256['p2b']:
            print("REFUSED: result_p2b.json does not match the pinned sha256 / last node.")
            sys.exit(EXIT_REFUSED_INPUT)
        logf = open(os.path.splitext(args.main)[0] + '.log', 'a')
        log(f"--- main invocation {time.strftime('%Y-%m-%d %H:%M:%S')} sha256 {sha_self()} engine {sha_file(ENGINE_FILE)}")
        sys.exit(run_grid(args.main, grid_nodes(), TIME_BUDGET, log))
    elif args.ceiling:
        if not args.p2c_sha256:
            print("REFUSED: --ceiling needs --p2c-sha256 (the sha256 of the finished P2c result).")
            sys.exit(EXIT_REFUSED_INPUT)
        logf = open(os.path.splitext(args.ceiling)[0] + '.log', 'a')
        sys.exit(run_ceiling(args, log))
    elif args._grid:
        sys.exit(run_grid(args._grid, test_grid(), args._budget, log, degfn=lambda j: (36, 10), prior=0.05))
    elif args.probe:
        run_probe(args, log)
