"""BHSLICE P2e instrument: decide the sixth digit of the Bowen-York maximal spin (BAR-P2e, seal 560847f6a4).

Cloned from p2d.py (P2d, HOLD 9677804a51; sealed, read-only). Same object, same certified route:
single, non-boosted, conformally flat, maximal Bowen-York spinning puncture (puncture BC, J = J_ADM), mu > 0. At nodes
j_k = 2^k, k = 0..K_TOP (K_TOP = 21; mu = 1 units, j = 4J):
  (a) point certificate: delta(j_k) = M_ADM - 1 in [delta_lo, delta_hi] (sqrt-map sandwich, REP4's 14 fixed pieces);
  (b) tangent certificate: psi_+ = Y^2 (Z_N + eta_Z) with TL-H, slope bound 2L_hi exact (Lemma TL, LEMMA-REP4.md sec. 6).
Both from the shared engine <LAB>/POCKETLAB/gravity_engine/bowen_york_sqrt.py, module version 2 (imported, never copied;
sha256 pinned below; NO engine change for P2e). Pieces: REP4's 14 fixed Y-pieces at every node; AMEND-P2d-2's targeted
refinement is NOT used at any node (recorded per node as pieces_rule).
Changes relative to p2d.py (BUILD-P2e.md):
  - nodes k = 0..21 and the degree schedule SCHEDULE below (fixed here before any main run; chosen from the build seat's
    double and rigorous probes, BUILD-P2e.md sec. 2): the top bands use N = 180 (k = 18, 19) and N = 220 (k = 20, 21),
    because at large j the tangent and point residuals are set by the radial degree N (K = 12..14 is already enough);
  - TL cells [2^k, 2^(k+1)] from node k (k = 0..20), tail [2^21, oo) from node 21; P1 piece chi <= 1/4 on [0, 1];
  - predictions E1-E4 of BAR-P2e evaluated in the result (E4 against result_p2d.json, read by sha256);
  - CEIL and TAIL independent exact-Fractions recomputations run inside build_result on EVERY result, so the REAL result
    of --main stores them (result keys CEIL_independent, TAIL_independent; REFEREE-P2d 5, BAR-P2e CEIL);
  - controls split into parts that each stay under ~520 s (engine_a, engine_b, node_a, node_b, nulls_a, nulls_b, rest);
    the CEIL control also exercises the cell-binding and P1-binding branches on synthetic variants of the mock result;
  - GO-P2e and the self-sha check against the lab ledger's go row for seal 560847f6a4 naming p2e.py.
Kept from p2d.py: exact assembly, TC and N-TAN at every node pair, K1 / K1T at j = 2^-64 (consistency info), the resumable
runner (JSON checkpoint with exact "p/q" strings, no pickle; lab exit codes 0 finished / 3 rerun / other = error; at most
TIME_BUDGET = 520 s per invocation), structural asserts, NaN/inf guards.
Blind: TARGETS-*.md were not opened by the build seat. Expected values come from the build seat's own probes (BUILD-P2e.md).

Usage (every run under memcap, one process at a time):
  py -3.11 <LAB>/POCKETLAB/lab/tools/memcap.py 8 py -3.11 p2e.py --controls controls_p2e.json --part engine_a
  ... --part engine_b | node_a | node_b | nulls_a | nulls_b | rest      (each part < 520 s; node_a, then node_b, then
                                                                          nulls_a / nulls_b; rest last)
  ... p2e.py --main result_p2e.json   (refuses without GO-P2e and a ledger GO row naming this file's sha256; resumable;
                                       exit 0 finished, 3 = rerun the same command, anything else = error)
"""
import sys
sys.dont_write_bytecode = True
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):   # single-thread LAPACK => deterministic fits
    os.environ[_v] = '1'
ENGINE_DIR = '<LAB>/POCKETLAB'
sys.path.insert(0, ENGINE_DIR)
from gravity_engine import bowen_york_sqrt as by          # import BEFORE numpy users (deterministic fit)
import json, time, math, hashlib, argparse, subprocess, shutil
from fractions import Fraction
from flint import fmpq

HERE = os.path.dirname(os.path.abspath(__file__))
GO_FILE = os.path.join(HERE, 'GO-P2e')
SEAL_ID = '560847f6a4'                                    # BAR-P2e claim
INSTRUMENT_NAME = 'p2e.py'
ENGINE_FILE = os.path.join(ENGINE_DIR, 'gravity_engine', 'bowen_york_sqrt.py')
ENGINE_TESTS = os.path.join(ENGINE_DIR, 'gravity_engine', 'tests')
ENGINE_SHA256 = '6f3a40993bf6ca4b082ee7699913fd49dec3beecf15365dfa12e730ccc7bc9b3'   # module v2, POCKETLAB 6a4a597
SCRATCH = os.path.join(HERE, 'seats', 'build-p2e')
TIME_BUDGET = 520.0
EXIT_FINISHED, EXIT_RESUME, EXIT_NO_GO, EXIT_CKPT_REFUSED, EXIT_REFUSED_INPUT, EXIT_NOT_PINNED, EXIT_CONTROLS_FAILED = \
    0, 3, 2, 4, 5, 6, 1
assert by.MODULE_VERSION == 2, 'p2e.py is built on gravity_engine bowen_york_sqrt module version 2'

# ---- earlier certified results, pinned by sha256: P2c (psi = 0 control), P2d (E4 consistency)
INPUTS = {'p2c': os.path.join(HERE, 'result_p2c.json'), 'p2d': os.path.join(HERE, 'result_p2d.json')}
PIN_SHA256 = {'p2c': '4b2dbf7b95742ce309a361fa693074779fb4f895f1f8afaecff5d874e60c33e7',     # as in p2d.py
              'p2d': '5620dd7a941afceb7d843a485eb11f5c94127eab06e522441eeebf1b6ab1b639'}     # VERDICT-P2d: 5620dd7a...

# ---- grid, schedule, constants (fixed before any main run)
K_TOP = 21
KS = tuple(range(K_TOP + 1))                              # nodes j = 2^k, k = 0..21 (BAR-P2e: K_top in [20, 24])
assert 20 <= K_TOP <= 24
N_PIECES = 8                                              # REP4 make_pieces(8): 14 fixed Y-pieces (AMEND-P2d-2 item 1)
PIECES_RULE = 'REP4 fixed 14 Y-pieces x [0, 1]; AMEND-P2d-2 targeted refinement NOT used'
# (k upper end, N, K, prec bits). k <= 17: P2d's schedule (same degrees and code => the same certificates as
# result_p2d.json); k = 18, 19 at N = 180 and k = 20, 21 at N = 220 (BUILD-P2e.md sec. 2: double + rigorous probes).
SCHEDULE = ((6, 64, 10, 192), (11, 90, 12, 192), (15, 110, 12, 192), (17, 140, 14, 224), (19, 180, 14, 224),
            (21, 220, 14, 224))
assert SCHEDULE[-1][0] == K_TOP
K1_RUNS = (('K1T_192', 64, 8, 192), ('K1T_256', 64, 8, 256))   # j = 2^-64, 8 start pieces (AMEND-P2d-1 item 5)
J_K1 = fmpq(1, 2 ** 64)
K1T_TOL = fmpq(1, 10 ** 9)                                # DESIGN-P2d sec. 4.4 K1T tolerance
S_HALF_LO_STR = '0.707106781186547524400844362104'        # sqrt(1/2) truncated at 30 digits: a LOWER bound
S_HALF_LO = fmpq(707106781186547524400844362104, 10 ** 30)
assert S_HALF_LO * S_HALF_LO * 2 <= 1 < (S_HALF_LO + fmpq(1, 10 ** 30)) ** 2 * 2
NTAN_SHIFT = fmpq(1, 10 ** 4)
# BUILD-P2e.md sec. 4 (written before any main run; REPRODUCIBILITY, not blind - REFEREE-P2d 5): pairs (k, k+1) where
# N-TAN is expected to bite (expected certified gap below (1 - 1/sqrt2) 1e-4 = 2.93e-5); every other pair is listed there
# as "cannot bite".
NTAN_EXPECTED_BITES = ((14, 15), (15, 16), (16, 17), (17, 18), (18, 19), (19, 20), (20, 21))
# BAR-P2e thresholds (exact)
E2_BAR = fmpq(1, 10 ** 7)                                 # E2: c - X3 < 1e-7
E3_BAR = fmpq(9282005, 10 ** 7)                           # E3: c < 0.9282005 -> HOLD; X3 >= 0.9282005 -> FAILED
GRID_LINE = ('j = 2^k, k = 0..21; cells [2^k, 2^(k+1)] from node k (TL, ratio 2); tail [2^21, oo) from node 21; '
             '[0, 1] by P1')


def degrees(k):
    for top, N, K, prec in SCHEDULE:
        if k <= top:
            return N, K, prec
    raise ValueError(k)


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


Q = by.parse_q
qs = by.qstr


def F(pq):
    a, b = pq.split('/')
    return Fraction(int(a), int(b))


def jq(k):
    return fmpq(2) ** k


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


def load_inputs(paths=None):
    """earlier results by sha256 (P2c: psi = 0 control; P2d: E4); Refusal on any mismatch."""
    paths = dict(INPUTS, **(paths or {}))
    out = {}
    for k, p in paths.items():
        if not os.path.exists(p):
            raise Refusal(f'{k}: {p} missing')
        s = sha_file(p)
        if s != PIN_SHA256[k]:
            raise Refusal(f'{k}: sha256 {s} of {p} != pinned {PIN_SHA256[k]}')
        out[k] = json.load(open(p))
    return out


# ================================================================== per-step work
def point_step(j, N, K, prec, log):
    """point certificate; returns (record without timing keys, private dict with the exact ansatz for the tangent step)."""
    rec, (lo, hi), c, wq, CS = by.certify_point(j, N, K, prec, n_pieces=N_PIECES, adaptive=False, keep=True, log=log)
    secs = {k: rec.pop(k) for k in list(rec) if k.startswith('secs')}
    priv = dict(wq=[qs(x) for x in wq], eta=qs(c['eta']), CS=qs(CS)) if rec['status'] == 'CERTIFIED' else None
    return rec, priv, secs


def tangent_step(priv, N, K, prec, log, **kw):
    wq = [Q(x) for x in priv['wq']]
    t = by.certify_tangent(wq, N, K, Q(priv['CS']), Q(priv['eta']), prec=prec, n_pieces=N_PIECES, **kw)
    r = by.tangent_record(t)
    if log:
        log(f"  tangent (N,K,prec)=({N},{K},{prec}): {r['status']} sup|H|<={r['sup_H_le']} eta_Z<={r['eta_Z_up']} "
            f"E_->={r['E_minus_min_ge']} 2L_hi<={r['twoL_hi_up']} ({t['secs']}s)")
    return r, t['secs']


def k1_step(label, N, K, prec, log):
    rec, priv, secs = point_step(J_K1, N, K, prec, log)
    out = dict(label=label, j=qs(J_K1), N=N, K=K, prec_bits=prec, point=rec)
    t2 = J_K1 * J_K1
    if priv is None:
        out.update(K1_contains_1_40=False, tangent=None, K1T_holds=False)
        return out, secs
    lo, hi = Q(rec['delta_lo_exact']), Q(rec['delta_hi_exact'])
    out['K1_contains_1_40'] = bool(lo / t2 <= fmpq(1, 40) <= hi / t2)
    tr, ts = tangent_step(priv, N, K, prec, log)
    out['tangent'] = tr
    if tr['status'] == 'CERTIFIED':
        x = (1 - Q(tr['twoL_hi_exact'])) / t2
        out['K1T_value'] = by.dec_floor(x, 25)
        out['K1T_minus_3_40'] = by.sci_up(x - fmpq(3, 40), 4) if x != fmpq(3, 40) else '0'
        out['K1T_holds'] = bool(abs(x - fmpq(3, 40)) <= K1T_TOL and x <= fmpq(3, 40) + K1T_TOL)
    else:
        out['K1T_holds'] = False
    return out, dict(secs, tangent=ts)


# ================================================================== assembly (exact)
def chi_of(j, d):
    return j / (4 * (1 + d) ** 2)


def le_x_times(P, Qv):
    """exact test P <= Qv / sqrt(2)."""
    if Qv >= 0:
        return P <= 0 or 2 * P * P <= Qv * Qv
    return P < 0 and 2 * P * P >= Qv * Qv


def tc_pair(ka, ra, kb, rb, shift=fmpq(0)):
    """tangent consistency (DESIGN-P2d sec. 4.4 TC) for nodes a = 2^ka, b = 2^kb = 2 j_a: TL's lower bound on M at mu' =
    1/sqrt2, 1 + delta_lo(a) - (1 - 1/sqrt2) (2L_hi(a) - shift), must not exceed the certified upper end
    (1/sqrt2)(1 + delta_hi(b)). Exact. Gap = (1/sqrt2)(B - t) - P (decimal, 30 digits, display)."""
    import mpmath as mp
    t = Q(ra['tangent']['twoL_hi_exact']) - shift
    P = 1 + Q(ra['point']['delta_lo_exact']) - t
    B = 1 + Q(rb['point']['delta_hi_exact'])
    holds = le_x_times(P, B - t)
    with mp.workdps(60):
        g = (mp.mpf(int((B - t).p)) / int((B - t).q)) / mp.sqrt(2) - mp.mpf(int(P.p)) / int(P.q)
        gs = mp.nstr(g, 6)
    return holds, gs


def build_result(node_recs, k1_recs, sha, ks=KS, mock=False, inputs=None):
    """node_recs: {k: {'point': rec, 'tangent': trec}}; k1_recs: {label: rec}. Exact assembly of c, X3, cells, tail;
    BAR-P2e predictions E1-E4; independent CEIL / TAIL recomputations (exact Fractions) stored in the result."""
    res = dict(object='BHSLICE P2e: Bowen-York spinning puncture (single, non-boosted, conformally flat, maximal, '
                      'puncture BC, J = J_ADM, mu > 0); sixth digit of the maximal spin via the tangent lemma',
               grid_line=GRID_LINE, instrument_sha256=sha, mock_for_controls=bool(mock),
               engine=dict(file='gravity_engine/bowen_york_sqrt.py', sha256=sha_file(ENGINE_FILE),
                           module_version=by.MODULE_VERSION, seal=by.SEAL_ID),
               degree_schedule=[dict(k_upto=t, N=N, K=K, prec_bits=p) for t, N, K, p in SCHEDULE],
               pieces_rule=PIECES_RULE, s_half_lo=S_HALF_LO_STR, ntan_shift='1/10000')
    nodes = []
    for k in ks:
        r = node_recs.get(k) or node_recs.get(str(k))
        if r is None:
            continue
        p, t = r['point'], r.get('tangent')
        j = jq(k)
        n = dict(k=k, j=qs(j), point=p, tangent=t, point_certified=bool(p['status'] == 'CERTIFIED'),
                 tangent_certified=bool(t is not None and t['status'] == 'CERTIFIED'),
                 pieces_rule=PIECES_RULE, targeted_refinement_used=False)
        if n['point_certified']:
            lo, hi = Q(p['delta_lo_exact']), Q(p['delta_hi_exact'])
            n['chi_lo'] = qs(chi_of(j, hi))
            n['chi_hi'] = qs(chi_of(j, lo))
            n['chi_lo_dec'] = by.dec_floor(chi_of(j, hi), 20)
            n['chi_hi_dec'] = by.dec_ceil(chi_of(j, lo), 20)
        nodes.append(n)
    res['nodes'] = nodes
    res['n_nodes'] = len(nodes)
    res['grid_complete'] = bool([n['k'] for n in nodes] == list(KS))
    res['all_point_certified'] = bool(nodes and all(n['point_certified'] for n in nodes))
    res['all_tangent_certified'] = bool(nodes and all(n['tangent_certified'] for n in nodes))
    res['k1'] = [k1_recs[l] for l, *_ in K1_RUNS if l in k1_recs]
    res['k1_all_hold'] = bool(len(res['k1']) == len(K1_RUNS) and all(x.get('K1_contains_1_40') and x.get('K1T_holds')
                                                                      for x in res['k1']))
    ef = []
    if not res['grid_complete']:
        ef.append('grid_complete (sub-grid or mock run for the controls)')
    if len(res['k1']) != len(K1_RUNS):
        ef.append('k1_all_hold (K1 runs not part of this sub-grid run)')
    pieces, ok_all = [], res['all_point_certified'] and res['all_tangent_certified']
    if ok_all:
        pieces.append(dict(piece='P1', j_lo=qs(fmpq(0)), j_hi=qs(fmpq(1)), value=qs(fmpq(1, 4)),
                           label='P1 range: chi = j/(4(1+delta)^2) <= j/4 <= 1/4 on j in [0, 1] (delta >= 0: E2)'))
        for a, b in zip(nodes[:-1], nodes[1:]):
            if b['k'] != a['k'] + 1:
                raise AssertionError('node list has a gap')
            ja, jb = jq(a['k']), jq(b['k'])
            B = by.tl_cell_bound(ja, jb, S_HALF_LO, Q(a['point']['delta_lo_exact']), Q(a['tangent']['twoL_hi_exact']))
            pieces.append(dict(piece='cell', k=a['k'], j_lo=qs(ja), j_hi=qs(jb), value=qs(B),
                               label=f'TL cell [2^{a["k"]}, 2^{b["k"]}] from node 2^{a["k"]}'))
        top = nodes[-1]
        T = by.tl_tail_bound(jq(top['k']), Q(top['point']['delta_lo_exact']), Q(top['tangent']['twoL_hi_exact']))
        pieces.append(dict(piece='tail', k=top['k'], j_lo=qs(jq(top['k'])), j_hi='oo', value=qs(T),
                           label=f'TL tail [2^{top["k"]}, oo) from node 2^{top["k"]}'))
        # exact cover [0, 1] u [1, 2] u ... u [2^K, oo): consecutive, starts at 0, ends at oo
        chain_ok = bool(Q(pieces[0]['j_lo']) == 0 and pieces[-1]['j_hi'] == 'oo' and nodes[0]['k'] == 0
                        and all(Q(pieces[i]['j_hi']) == Q(pieces[i + 1]['j_lo']) for i in range(len(pieces) - 1)))
        if not chain_ok:
            raise AssertionError('pieces do not cover [0, oo) exactly')
        for x in pieces:
            x['value_dec_up'] = by.dec_ceil(Q(x['value']), 20)
        vals = [Q(x['value']) for x in pieces]
        c = max(vals)
        bind = [x['label'] for x, v in zip(pieces, vals) if v == c]
        best = max(nodes, key=lambda n: Q(n['chi_lo']))
        X3 = Q(best['chi_lo'])
        X3top = Q(nodes[-1]['chi_lo'])
        res.update(pieces=pieces, cover_exact=chain_ok, c=qs(c), c_dec_up=by.dec_ceil(c, 20), binding_piece=bind[0],
                   n_pieces_attaining_c=len(bind), X3=qs(X3), X3_dec_down=by.dec_floor(X3, 20), X3_node_k=best['k'],
                   X3_top=qs(X3top), X3_top_dec_down=by.dec_floor(X3top, 20), c_minus_X3=qs(c - X3),
                   c_minus_X3_dec_up=by.sci_up(c - X3, 6), c_lt_1=bool(c < 1),
                   statement='for every j > 0 (every J > 0 at mu > 0): chi = J/M_ADM^2 <= c; the supremum of chi over '
                             'the family is >= X3 (certified lower bound at a node)')
        # TC + N-TAN at every pair
        pairs = []
        for a, b in zip(nodes[:-1], nodes[1:]):
            h0, g0 = tc_pair(a['k'], a, b['k'], b)
            h1, g1 = tc_pair(a['k'], a, b['k'], b, NTAN_SHIFT)
            pairs.append(dict(pair=[a['k'], b['k']], tc_holds=h0, tc_gap=g0, ntan_still_consistent=h1, ntan_gap=g1,
                              ntan_predicted_to_bite=bool((a['k'], b['k']) in NTAN_EXPECTED_BITES)))
        res['tc_pairs'] = pairs
        res['tc_all_hold'] = bool(all(p['tc_holds'] for p in pairs))
        pred = [p for p in pairs if p['ntan_predicted_to_bite']]
        res['ntan_bites_where_predicted'] = bool(all(not p['ntan_still_consistent'] for p in pred)) if pred else None
        res['ntan_bites_at'] = [p['pair'] for p in pairs if not p['ntan_still_consistent']]
        res['ntan_bites_only_where_predicted'] = bool(all(p['ntan_predicted_to_bite'] for p in pairs
                                                          if not p['ntan_still_consistent']))
        if not pred:
            ef.append('ntan_bites_where_predicted (no predicted pair in this sub-grid)')
        for i, p in enumerate(pairs):
            if p['ntan_predicted_to_bite']:
                ef.append(f'tc_pairs[{i}].ntan_still_consistent (pair {p["pair"][0]}-{p["pair"][1]}: 2L_hi - 1e-4 must '
                          f'break TC here, BUILD-P2e sec. 4)')
        # BAR-P2e predictions (E4 below, with the P2d input)
        res['E1_all_certify'] = bool(res['all_point_certified'] and res['all_tangent_certified'] and res['grid_complete'])
        res['E2_c_minus_X3_lt_1e-7'] = bool(0 <= c - X3 < E2_BAR)
        res['E3_c_lt_0.9282005'] = bool(c < E3_BAR)
        res['E3_X3_ge_0.9282005'] = bool(X3 >= E3_BAR)
        res['E3_scoring'] = ('HOLD (c < 0.9282005: the supremum rounds to 0.928200)' if c < E3_BAR else
                             'FAILED (X3 >= 0.9282005: the supremum rounds to 0.928201)' if X3 >= E3_BAR else
                             'UNDECIDED-FAILED (the window straddles 0.9282005)')
        if not res['E1_all_certify']:
            ef.append('E1_all_certify (sub-grid or mock run)')
        if mock:
            for key in ('E2_c_minus_X3_lt_1e-7', 'E3_c_lt_0.9282005', 'E3_X3_ge_0.9282005'):
                if not res[key]:
                    ef.append(f'{key} (mock sub-grid for the controls: coarse low-degree tail, not P2e)')
        elif not res['E3_X3_ge_0.9282005']:
            ef.append('E3_X3_ge_0.9282005 (the FAILED side of E3; false unless the opposite rounding is certified)')
        # independent recomputations on THIS result (REFEREE-P2d 5 / BAR-P2e CEIL): stored with the result
        res['CEIL_independent'] = ceil_compare(res)
        res['TAIL_independent'] = tail_compare(res)
    else:
        res.update(pieces=[], c=None, binding_piece=None)
    res['nodes_finite'] = bool(all(finite_node(n) for n in nodes))
    if inputs is not None:
        res['E4'] = e4_checks(nodes, res, inputs)
        if mock and res['E4'].get('ineq_all_hold') is False:
            ef.append('E4.ineq_all_hold (mock sub-grid for the controls: its c is a coarse low-degree tail, not P2e)')
    res['expected_false'] = ef
    return res


def finite_node(n):
    """NaN/inf guard on a stored node: every exact string parses, every decimal string is finite."""
    try:
        for part in (n['point'], n.get('tangent') or {}):
            for k, v in part.items():
                if isinstance(v, str) and k.endswith('_exact'):
                    F(v)
                if isinstance(v, str) and k in ('delta_lo', 'delta_hi', 'twoL_hi_up', 'eta_Z_up', 'sup_H_le',
                                                'sup_Gh_N_le') and not math.isfinite(float(v)):
                    return False
                if isinstance(v, float) and not math.isfinite(v):
                    return False
        return True
    except Exception:
        return False


# ================================================================== E4 consistency with P2d (BAR-P2e)
def e4_checks(nodes, res, inp):
    """E4 against result_p2d.json (sha256-pinned): (i) at every shared node k the delta intervals overlap (exact strings);
    (ii) c <= P2d's c and X3 >= P2d's X3 (exact rationals stored by P2d)."""
    out = dict(shared=[])
    p2d = inp['p2d']
    earlier = {}
    for r in p2d['nodes']:
        if r.get('point_certified'):
            earlier[int(r['k'])] = (F(r['point']['delta_lo_exact']), F(r['point']['delta_hi_exact']))
    for n in nodes:
        if not n['point_certified'] or n['k'] not in earlier:
            continue
        lo, hi = F(n['point']['delta_lo_exact']), F(n['point']['delta_hi_exact'])
        a, b = earlier[n['k']]
        out['shared'].append(dict(k=n['k'], overlap=bool(max(lo, a) <= min(hi, b)),
                                  identical_to_P2d=bool(lo == a and hi == b), nested_in_P2d=bool(a <= lo and hi <= b)))
    out['n_shared'] = len(out['shared'])
    out['shared_all_overlap'] = bool(out['shared'] and all(x['overlap'] for x in out['shared']))
    cP2d, X3P2d = F(p2d['c']), F(p2d['X3'])
    out['P2d_c_dec_up'] = p2d['c_dec_up']
    out['P2d_X3_dec_down'] = p2d['X3_dec_down']
    if res.get('c') is not None:
        c, X3 = F(res['c']), F(res['X3'])
        out['c_le_P2d_c'] = bool(c <= cP2d)
        out['X3_ge_P2d_X3'] = bool(X3 >= X3P2d)
        out['ineq_all_hold'] = bool(out['c_le_P2d_c'] and out['X3_ge_P2d_X3'])
        out['E4_holds'] = bool(out['ineq_all_hold'] and out['shared_all_overlap'])
    return out


# ================================================================== resumable runner
def main_steps(ks=KS, with_k1=True):
    st = []
    if with_k1:
        st += [('K1', l) for l, *_ in K1_RUNS]
    for k in ks:
        st += [('point', k), ('tangent', k)]
    return st


def step_key(s):
    return f'{s[0]}:{s[1]}'


def est_seconds(step, degfn):
    if step[0] == 'K1':
        return 25.0
    N, K, prec = degfn(step[1])
    base = 45.0 if step[0] == 'point' else 65.0          # P2d main at (140, 14, 224); P2e probes: (180): 84 s / 94 s,
    return base * (N / 140.0) ** 2 * (K / 14.0) * (prec / 224.0)   # (220): 90 s / 136 s (BUILD-P2e sec. 5); budget only


class RunState:
    """checkpoint = JSON with exact "p/q" strings; refuses a checkpoint from another instrument or engine sha."""

    def __init__(self, out, steps):
        self.out = out
        self.ck = out + '.ckpt'
        self.sha = sha_self()
        self.esha = sha_file(ENGINE_FILE)
        self.spec = [step_key(s) for s in steps]
        self.st = None

    def load(self, log):
        if os.path.exists(self.ck):
            st = json.load(open(self.ck))
            if st.get('instrument_sha256') != self.sha or st.get('engine_sha256') != self.esha:
                log(f"REFUSED: checkpoint {self.ck} was written by instrument {st.get('instrument_sha256')} / engine "
                    f"{st.get('engine_sha256')}; this is {self.sha} / {self.esha}.")
                return False
            if st.get('steps_spec') != self.spec:
                log('REFUSED: checkpoint step list differs from this run.')
                return False
        else:
            st = dict(instrument_sha256=self.sha, engine_sha256=self.esha, steps_spec=self.spec, done={}, private={},
                      invocations=0)
        st['invocations'] += 1
        self.st = st
        self.save()
        return True

    def save(self):
        atomic_json(self.ck, self.st)


def run_steps(out, steps, budget, log, degfn=degrees, ks=KS, est_scale=1.0, inputs=None):
    """Lab exit codes: 0 finished (result written), 3 budget reached (rerun the same command), 4 checkpoint refused."""
    T0 = time.time()
    rs = RunState(out, steps)
    if not rs.load(log):
        return EXIT_CKPT_REFUSED
    info_p = out + '.runinfo.json'
    info = json.load(open(info_p)) if os.path.exists(info_p) else dict(timing={}, invocations=[])
    ratio = 1.0
    did = 0
    for s in steps:
        key = step_key(s)
        if key in rs.st['done']:
            continue
        est = est_seconds(s, degfn) * est_scale * ratio
        if did >= 1 and time.time() - T0 + 1.3 * est > budget - 15:
            log(f"time budget reached: {len(rs.st['done'])}/{len(steps)} steps done; exit {EXIT_RESUME} = rerun the same "
                f"command to resume")
            info['invocations'].append(dict(ended='budget', seconds=round(time.time() - T0, 1), peak=peak_memory()))
            atomic_json(info_p, info)
            return EXIT_RESUME
        t = time.time()
        if s[0] == 'K1':
            N, K, prec = next((n, kk, p) for l, n, kk, p in K1_RUNS if l == s[1])
            rec, secs = k1_step(s[1], N, K, prec, log)
            rs.st['done'][key] = rec
        elif s[0] == 'point':
            N, K, prec = degfn(s[1])
            rec, priv, secs = point_step(jq(s[1]), N, K, prec, log)
            rs.st['done'][key] = rec
            if priv is not None:
                rs.st['private'][str(s[1])] = priv
        else:
            N, K, prec = degfn(s[1])
            priv = rs.st['private'].get(str(s[1]))
            if priv is None:
                rec, secs = dict(status='NOT RUN', certified=False, fail_reason='point certificate not certified'), 0
            else:
                rec, secs = tangent_step(priv, N, K, prec, log)
            rs.st['done'][key] = rec
            rs.st['private'].pop(str(s[1]), None)           # the exact ansatz is not needed any more
        dt = time.time() - t
        ratio = max(ratio if did else 0.0, dt / max(1e-3, est_seconds(s, degfn) * est_scale))
        info['timing'][key] = dict(seconds=round(dt, 1), engine_secs=secs, process_peak_after_step=peak_memory())
        rs.save()
        atomic_json(info_p, info)
        did += 1
    node_recs = {k: dict(point=rs.st['done'][f'point:{k}'], tangent=rs.st['done'][f'tangent:{k}']) for k in ks}
    k1 = {s[1]: rs.st['done'][step_key(s)] for s in steps if s[0] == 'K1'}
    res = build_result(node_recs, k1, rs.sha, ks=ks, mock=(tuple(ks) != KS), inputs=inputs)
    atomic_json(out, res)
    info['invocations'].append(dict(ended='finished', seconds=round(time.time() - T0, 1), peak=peak_memory()))
    atomic_json(info_p, info)
    log(f"finished: {len(steps)} steps, points certified {res['all_point_certified']}, tangents certified "
        f"{res['all_tangent_certified']}, c = {res.get('c_dec_up')} ({res.get('binding_piece')}), X3 = "
        f"{res.get('X3_dec_down')}, c - X3 = {res.get('c_minus_X3_dec_up')}, TC all hold {res.get('tc_all_hold')}, "
        f"N-TAN bites at {res.get('ntan_bites_at')}, K1/K1T {res.get('k1_all_hold')}, invocations {rs.st['invocations']}, "
        f"peak this invocation {peak_memory()}")
    if res.get('c') is not None:
        log(f"BAR-P2e: E1 {res.get('E1_all_certify')}; E2 (c - X3 < 1e-7) {res.get('E2_c_minus_X3_lt_1e-7')}; E3 "
            f"{res.get('E3_scoring')}; E4 {(res.get('E4') or {}).get('E4_holds')}; CEIL independent identical "
            f"{res['CEIL_independent'].get('ok')}; TAIL independent identical {res['TAIL_independent'].get('ok')}")
    return EXIT_FINISHED


# ================================================================== GO gate and self-sha check (REFEREE-P2c 5a)
def ledger_path():
    from pathlib import Path
    return str(Path(os.environ.get('LAB_DATA') or (Path.home() / '.pocketlab')) / 'ledger.jsonl')


def go_pinned(sha, ledger=None):
    """True iff the lab ledger has a 'go' row for SEAL_ID naming an instrument p2e.py with exactly this sha256."""
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
    """returns (exit code or None, message). Order: GO file, self-sha in the ledger, engine pin, inputs pins."""
    if not os.path.exists(go_file):
        return EXIT_NO_GO, f'REFUSED: {go_file} does not exist (the lead arms it with lab.py go after sealing).'
    ok, why = go_pinned(sha_self(), ledger)
    if not ok:
        return EXIT_NOT_PINNED, f'REFUSED: self-sha check failed: {why}'
    if sha_file(ENGINE_FILE) != ENGINE_SHA256:
        return EXIT_REFUSED_INPUT, f'REFUSED: engine sha256 {sha_file(ENGINE_FILE)} != pinned {ENGINE_SHA256}'
    try:
        load_inputs()
    except Refusal as e:
        return EXIT_REFUSED_INPUT, f'REFUSED: {e}'
    return None, 'GO: ' + why


# ================================================================== CEIL (independent exact Fractions) and TAIL
def isqrt_half_lo(digits=30):
    """floor(sqrt(1/2) 10^d) / 10^d, by integer square root (independent of the literal S_HALF_LO_STR)."""
    return Fraction(math.isqrt(10 ** (2 * digits) // 2), 10 ** digits)


def ceil_independent(res):
    """recompute every piece, c, the binding piece, X3 and c - X3 from the stored node records with Python Fractions only
    (no engine code; sqrt(1/2) lower bound by isqrt)."""
    s = isqrt_half_lo()
    nodes = res['nodes']
    vals = [(Fraction(1, 4), 'P1')]
    for a, b in zip(nodes[:-1], nodes[1:]):
        ja = Fraction(2) ** a['k']
        dl = F(a['point']['delta_lo_exact'])
        t = max(F(a['tangent']['twoL_hi_exact']), Fraction(0))
        br = 1 + dl - t * (1 - s)
        if not br > 0:
            raise AssertionError('independent CEIL: cell bracket <= 0')
        vals.append((ja / 4 / br ** 2, f'cell {a["k"]}'))
    top = nodes[-1]
    jt = Fraction(2) ** top['k']
    br = 1 + F(top['point']['delta_lo_exact']) - max(F(top['tangent']['twoL_hi_exact']), Fraction(0))
    if not br > 0:
        raise AssertionError('independent CEIL: tail bracket <= 0')
    vals.append((jt / 4 / br ** 2, 'tail'))
    c = max(v for v, _ in vals)
    bind = [l for v, l in vals if v == c][0]
    X3 = max(Fraction(2) ** n['k'] / (4 * (1 + F(n['point']['delta_hi_exact'])) ** 2) for n in nodes)
    return dict(c=c, binding=bind, X3=X3, c_minus_X3=c - X3, s_equal=bool(s == F(qs(S_HALF_LO))), n=len(vals))


def tail_independent(res):
    top = res['nodes'][-1]
    jt = Fraction(2) ** top['k']
    dl = F(top['point']['delta_lo_exact'])
    t = F(top['tangent']['twoL_hi_exact'])
    t = t if t > 0 else Fraction(0)
    return jt / (4 * (1 + dl - t) * (1 + dl - t)), bool(1 + dl - t > 0)


def ceil_compare(res):
    """run ceil_independent on an assembled result and compare with the instrument's own values (exact). Stored in every
    result (key CEIL_independent), so the REAL main result carries its own independent recomputation."""
    try:
        ci = ceil_independent(res)
    except Exception as e:                                # an independent failure is recorded, never hidden
        return dict(ok=False, error=repr(e)[:300])
    bi, bp = ci['binding'], res['binding_piece']
    same_bind = bool((bi == 'tail' and bp.startswith('TL tail')) or (bi == 'P1' and bp.startswith('P1'))
                     or (bi.startswith('cell') and bp.startswith('TL cell [2^%s,' % bi[5:])))
    out = dict(method='Python Fractions only (no engine code); sqrt(1/2) lower bound by integer isqrt',
               c_instrument_dec_up=res['c_dec_up'],
               c_independent_dec_up=by.dec_ceil(fmpq(ci['c'].numerator, ci['c'].denominator), 20),
               identical=bool(F(res['c']) == ci['c']), binding_instrument=bp, binding_independent=bi,
               same_binding=same_bind, X3_identical=bool(F(res['X3']) == ci['X3']),
               c_minus_X3_identical=bool(F(res['c_minus_X3']) == ci['c_minus_X3']),
               s_half_isqrt_equals_literal=ci['s_equal'], n_pieces=ci['n'], n_pieces_instrument=len(res['pieces']))
    out['ok'] = bool(out['identical'] and out['same_binding'] and out['X3_identical'] and out['c_minus_X3_identical']
                     and out['s_half_isqrt_equals_literal'] and out['n_pieces'] == out['n_pieces_instrument'])
    return out


def tail_compare(res):
    """independent tail value (Fractions) vs the instrument's tail piece (exact). Stored in every result."""
    try:
        tv, tb_ok = tail_independent(res)
    except Exception as e:
        return dict(ok=False, error=repr(e)[:300])
    tl = next(x for x in res['pieces'] if x['piece'] == 'tail')
    out = dict(identical=bool(F(tl['value']) == tv), bracket_positive_independent=tb_ok,
               tail_dec_up=by.dec_ceil(fmpq(tv.numerator, tv.denominator), 20))
    out['ok'] = bool(out['identical'] and tb_ok)
    return out


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


PARTS = ('engine_a', 'engine_b', 'node_a', 'node_b', 'nulls_a', 'nulls_b', 'rest')
PREFIX = {'engine_a': ('ENGINE_A',), 'engine_b': ('ENGINE_B',),
          'node_a': ('NODEA_',),
          'node_b': ('NODEB_', 'MEM', 'TC_TOP', 'NTAN_TOP', 'nan_inf', 'finite_nodes'),
          'nulls_a': ('N1_', 'N4_'),
          'nulls_b': ('N4T_', 'NBUMP_T_'),
          'rest': ('psi0', 'STRUCT', 'CEIL', 'TAIL', 'RESUME', 'main_refusal', 'selfsha', 'E4_mock', 'finite_guard')}


def _finish_controls(out, R, log):
    oks = R['oks']
    R['parts_done'] = sorted(R['parts'].keys())
    R['resume_tested'] = bool(oks.get('RESUME', False))
    R['all_ok'] = bool(set(PARTS) <= set(R['parts'].keys()) and all(oks.values()))
    R['degree_schedule'] = [dict(k_upto=t, N=N, K=K, prec_bits=p) for t, N, K, p in SCHEDULE]
    R['ntan_expected_bites'] = [list(p) for p in NTAN_EXPECTED_BITES]
    atomic_json(out, R, sort=False)
    log(f"controls parts {R['parts_done']}; all_ok = {R['all_ok']}; oks = {oks}; instrument sha256 = {R['instrument_sha256']}")
    return 0 if all(oks.values()) else EXIT_CONTROLS_FAILED


def _pytest(files):
    t = time.time()
    r = subprocess.run([sys.executable, '-B', '-m', 'pytest'] + files + ['-q', '-p', 'no:cacheprovider'],
                       capture_output=True, text=True)
    tail = [ln for ln in r.stdout.strip().splitlines() if ln.strip()][-1:] or ['']
    return dict(rc=r.returncode, summary=tail[0], seconds=round(time.time() - t, 1), files=[os.path.basename(f) for f in files])


ENGINE_A_FILES = ['test_bowen_york_sqrt.py', 'test_bowen_york_tangent.py']


def controls_engine(R, log, which):
    P = R['parts'].setdefault(which, {})
    allf = sorted(f for f in os.listdir(ENGINE_TESTS) if f.startswith('test_') and f.endswith('.py'))
    files = ENGINE_A_FILES if which == 'engine_a' else [f for f in allf if f not in ENGINE_A_FILES]
    r = _pytest([os.path.join(ENGINE_TESTS, f) for f in files])
    P['ENGINE'] = r
    key = 'ENGINE_A' if which == 'engine_a' else 'ENGINE_B'
    R['oks'][key] = bool(r['rc'] == 0 and 'passed' in r['summary'] and 'failed' not in r['summary']
                         and 'error' not in r['summary'])
    log(f"{key} tests {r['files']}: rc={r['rc']} '{r['summary']}' ({r['seconds']}s)")


def _node_control(R, log, part, k, prefix):
    """one real node at the scheduled degrees: point + tangent; returns (record dict, private ansatz or None)."""
    T0 = time.time()
    oks = R['oks']
    P = R['parts'].setdefault(part, {})
    base = peak_memory()
    N, K, prec = degrees(k)
    rec, priv, secs = point_step(jq(k), N, K, prec, log)
    tr, ts = tangent_step(priv, N, K, prec, log) if priv else (dict(status='NOT RUN'), 0)
    P['node'] = dict(k=k, degrees=[N, K, prec], seconds=round(time.time() - T0, 1), point_secs=secs, tangent_secs=ts,
                     point=rec, tangent=tr, pieces_rule=PIECES_RULE, targeted_refinement_used=False)
    oks[f'{prefix}point_certified_k{k}'] = bool(rec['status'] == 'CERTIFIED')
    oks[f'{prefix}tangent_certified_k{k}'] = bool(tr.get('status') == 'CERTIFIED')
    P['peak'] = dict(baseline_before=base, after_node=peak_memory())
    return dict(point=rec, tangent=tr, k=k), priv


def controls_node_a(R, log):
    """node k = K_TOP - 1 at its scheduled degrees (the lower node of the top pair)."""
    T0 = time.time()
    _node_control(R, log, 'node_a', K_TOP - 1, 'NODEA_')
    R['parts']['node_a']['seconds'] = round(time.time() - T0, 1)


def controls_node_b(R, log):
    """top node k = K_TOP at the scheduled degrees: point + tangent; MEM; TC and N-TAN on the real top pair
    (K_TOP - 1, K_TOP) (node_a's record); NaN/inf guard of the tangent fit; finite-record guard."""
    T0 = time.time()
    oks = R['oks']
    rb, priv = _node_control(R, log, 'node_b', K_TOP, 'NODEB_')
    P = R['parts']['node_b']
    P['top_private'] = priv                               # exact ansatz of the top node for nulls_a/b (p/q strings)
    m = P['peak']['after_node'] or {}
    P['MEM'] = dict(peak_after_top_node=m, cap_gb=8, target_gb=6)
    oks['MEM_under_cap'] = bool(m and m['peak_commit_gb'] < 8 and m['peak_working_set_gb'] < 8)
    log(f"MEM: top node k = {K_TOP}: process peak {m}")
    na = (R['parts'].get('node_a') or {}).get('node')
    pair = (K_TOP - 1, K_TOP)
    if na and na['tangent'].get('status') == 'CERTIFIED' and rb['tangent'].get('status') == 'CERTIFIED':
        a = dict(point=na['point'], tangent=na['tangent'], k=K_TOP - 1)
        h0, g0 = tc_pair(K_TOP - 1, a, K_TOP, rb)
        h1, g1 = tc_pair(K_TOP - 1, a, K_TOP, rb, NTAN_SHIFT)
        pred = bool(pair in NTAN_EXPECTED_BITES)
        P['TC_TOP'] = dict(pair=list(pair), holds=h0, gap=g0)
        P['NTAN_TOP'] = dict(pair=list(pair), still_consistent=h1, gap=g1, shift='1e-4', predicted_to_bite=pred)
        oks['TC_TOP_holds'] = bool(h0)
        oks['NTAN_TOP_as_predicted'] = bool((not h1) if pred else h1)
        if pred:
            R['expected_false'].append('parts.node_b.NTAN_TOP.still_consistent (2L_hi lowered by 1e-4 must break TC here)')
        log(f"TC pair {pair[0]}-{pair[1]}: holds {h0}, gap {g0}; N-TAN (2L_hi - 1e-4): still consistent {h1} (gap {g1}); "
            f"predicted to bite {pred}")
    else:
        oks['TC_TOP_holds'] = oks['NTAN_TOP_as_predicted'] = False
        log('TC/N-TAN top pair: run --part node_a first (or a top tangent did not certify)')
    # NaN / inf guard (tangent fit replaced before any arb work)
    g = {}
    N, K, prec = degrees(K_TOP)
    if priv:
        for bad in (float('nan'), float('inf'), -float('inf'), 1e300):
            z = [0.0] * ((N + 1) * (K + 1))
            z[N + 4] = bad
            r, _ = tangent_step(priv, N, K, prec, None, z_override=z)
            g[str(bad)] = r['status']
    P['nan_inf_guard'] = g
    oks['nan_inf_guard'] = bool(len(g) == 4 and all(v == 'INVALID' for v in g.values()))
    log(f"NaN/inf guard (tangent fit at the top degree): {g}")
    oks['finite_nodes'] = bool(finite_node(rb) and (na is None or finite_node(na)))
    P['seconds'] = round(time.time() - T0, 1)


def _top(R):
    nd = R['parts'].get('node_b', {})
    node = nd.get('node') or {}
    return nd.get('top_private'), node.get('point'), node.get('tangent')


def controls_nulls_a(R, log):
    """top node (k = K_TOP): N1 wrong ansatz -> FAIL; N4 near-threshold point margin -> FAIL through the bound path."""
    T0 = time.time()
    oks, ef = R['oks'], R['expected_false']
    P = R['parts'].setdefault('nulls_a', {})
    priv, rec, tr = _top(R)
    if not priv or not rec or rec.get('status') != 'CERTIFIED' or (tr or {}).get('status') != 'CERTIFIED':
        oks['N1_top_fail'] = oks['N4_top_fail_bound_path'] = False
        log('nulls_a: run --part node_b first (top node not certified there)')
        return
    N, K, prec = degrees(K_TOP)
    wq = [Q(x) for x in priv['wq']]
    CS, eta = Q(priv['CS']), Q(priv['eta'])
    sup = Q(rec['sup_Gh_N_le_exact'])
    # N1: wrong ansatz (P_0 raised by 1e-6 of itself) with the honest eta forced
    wb = list(wq)
    wb[0] = wb[0] + wq[0] / 10 ** 6
    c1 = by.certify(wb, N, K, CS, n_pieces=N_PIECES, prec=prec, adaptive=False, eta=eta)
    s1 = by.sample_signs(wb, N, K, CS, eta, n_y=24, n_v=6)
    P['N1_top'] = dict(k=K_TOP, status=c1['status'], strict=c1['strict'], sup_Gh_le=by.sci_up(c1['sup_Gh_le']),
                       honest_sup=by.sci_up(sup), sample_signs_ok=s1['signs_ok'])
    oks['N1_top_fail'] = bool(c1['status'] == 'FAIL' and not s1['signs_ok'])
    ef.append('parts.nulls_a.N1_top.sample_signs_ok (the wrong ansatz violates the exact sample signs)')
    log(f"N1 (P_0 + 1e-6 P_0, honest eta): {c1['status']}, exact sample signs ok {s1['signs_ok']}")
    # N4: near-threshold point margin (2/3) sup (1 - 2^-10): FAIL through the bound path, sample signs still fine
    e4 = fmpq(2, 3) * sup * (1 - fmpq(1, 2 ** 10))
    c4 = by.certify(wq, N, K, CS, n_pieces=N_PIECES, prec=prec, adaptive=False, eta=e4)
    s4 = by.sample_signs(wq, N, K, CS, e4, n_y=24, n_v=6)
    P['N4_top'] = dict(k=K_TOP, status=c4['status'], strict=c4['strict'], sup_equal_honest=bool(c4['sup_Gh_le'] == sup),
                       sample_signs_ok=s4['signs_ok'], eta_forced=by.sci_up(e4))
    oks['N4_top_fail_bound_path'] = bool(c4['status'] == 'FAIL' and not c4['strict'] and c4['sup_Gh_le'] == sup
                                         and s4['signs_ok'])
    ef.append('parts.nulls_a.N4_top.strict (a margin (2/3) sup (1 - 2^-10) must not pass the bound path)')
    log(f"N4 (eta = (2/3) sup (1 - 2^-10)): {c4['status']}, strict {c4['strict']}, sample signs ok {s4['signs_ok']}")
    P['seconds'] = round(time.time() - T0, 1)


def controls_nulls_b(R, log):
    """top node (k = K_TOP): N4-T near-threshold eta_Z -> FAIL; tangent bump with the honest eta_Z -> FAIL."""
    T0 = time.time()
    oks, ef = R['oks'], R['expected_false']
    P = R['parts'].setdefault('nulls_b', {})
    priv, rec, tr = _top(R)
    if not priv or not rec or rec.get('status') != 'CERTIFIED' or (tr or {}).get('status') != 'CERTIFIED':
        oks['N4T_top_fail'] = oks['NBUMP_T_top_fail'] = False
        log('nulls_b: run --part node_b first (top node not certified there)')
        return
    N, K, prec = degrees(K_TOP)
    wq = [Q(x) for x in priv['wq']]
    CS, eta = Q(priv['CS']), Q(priv['eta'])
    supH = Q(tr['sup_H_le_exact'])
    eZ = fmpq(2, 3) * supH * (1 - fmpq(1, 2 ** 10))
    t4 = by.certify_tangent(wq, N, K, CS, eta, prec=prec, n_pieces=N_PIECES, eta_Z=eZ)
    P['N4T_top'] = dict(k=K_TOP, status=t4['status'], strict=t4['strict'], sup_equal_honest=bool(t4['sup_H_le'] == supH),
                        Z0_equal_honest=bool(qs(t4['Z0']) == tr['Z0_exact']))
    oks['N4T_top_fail'] = bool(t4['status'] == 'FAIL' and not t4['strict'] and t4['sup_H_le'] == supH
                               and qs(t4['Z0']) == tr['Z0_exact'])
    ef.append('parts.nulls_b.N4T_top.strict (eta_Z = (2/3) sup (1 - 2^-10) must FAIL)')
    log(f"N4-T (eta_Z = (2/3) sup (1 - 2^-10)): {t4['status']}")
    # tangent bump: Z_N with P_Z0 raised by 1e-6 max|Z|, honest eta_Z forced
    zq = list(t4['zq'])
    mx = max(abs(x) for x in zq)
    zq[0] = zq[0] + mx / 10 ** 6
    tb = by.certify_tangent(wq, N, K, CS, eta, prec=prec, n_pieces=N_PIECES, zq_override=zq,
                            eta_Z=Q(tr['eta_Z_exact']))
    P['NBUMP_T_top'] = dict(k=K_TOP, status=tb['status'], strict=tb['strict'], sup_H_le=by.sci_up(tb['sup_H_le']),
                            honest_sup=by.sci_up(supH))
    oks['NBUMP_T_top_fail'] = bool(tb['status'] == 'FAIL')
    ef.append('parts.nulls_b.NBUMP_T_top.strict (a bumped tangent with the honest eta_Z must FAIL)')
    log(f"tangent bump (Z_0 + 1e-6 max|Z|, honest eta_Z): {tb['status']}")
    P['seconds'] = round(time.time() - T0, 1)


def test_ks():
    return (0, 1, 2, 3, 4)


def test_degrees(k):
    return 24, 4, 128


def _variant(resA, how):
    """synthetic variants of the mock result A (exercise the CEIL branches the real data may not reach):
    'cell': node k = 2 gets a large 2L_hi so that its TL cell binds; 'P1': every delta set to 10^6 so that P1 binds."""
    recs = {}
    for n in resA['nodes']:
        p = json.loads(json.dumps(n['point']))
        t = json.loads(json.dumps(n['tangent']))
        if how == 'P1':
            p['delta_lo_exact'] = p['delta_hi_exact'] = '1000000/1'
        if how == 'cell' and n['k'] == 2:
            dl = Q(p['delta_lo_exact'])
            t['twoL_hi_exact'] = qs((1 + dl) / 2 / (1 - S_HALF_LO))       # cell bracket = (1 + delta_lo)/2 > 0
        recs[n['k']] = dict(point=p, tangent=t)
    return build_result(recs, {}, resA['instrument_sha256'], ks=test_ks(), mock=True, inputs=None)


def controls_rest(R, log):
    """psi = 0, structural asserts, RESUME (+ exit codes), CEIL (mock + cell/P1-binding variants), TAIL, E4 on the mock,
    finite guard, --main refusal, self-sha check."""
    T0 = time.time()
    oks, ef = R['oks'], R['expected_false']
    P = R['parts'].setdefault('rest', {})
    priv_top, rec_top, _ = _top(R)
    me = [sys.executable, os.path.abspath(__file__)]
    scratch = os.path.join(SCRATCH, 'controls')
    if os.path.isdir(scratch):
        shutil.rmtree(scratch)
    os.makedirs(scratch)
    inp = load_inputs()

    # --- psi = 0: exact reduction to P2c's plain cells and X1 (REFEREE-TL 1)
    rP = {}
    top = inp['p2c']['nodes'][-1]
    j1, dl = Q(top['j']), Q(top['delta_lo'])
    X1_stored = Q(inp['p2c']['tail_X1']['X1'])
    rP['tail_with_2L_1_equals_P2c_X1'] = bool(by.tl_tail_bound(j1, dl, fmpq(1)) == X1_stored)
    eq = []
    for r in inp['p2c']['nodes'][:5]:
        ja, d = Q(r['j']), Q(r['delta_lo'])
        eq.append(by.tl_cell_bound(ja, 4 * ja, fmpq(1, 2), d, fmpq(1)) == (4 * ja) / (4 * (1 + 2 * d) ** 2))
    rP['cells_with_2L_1_equal_P2c_form'] = bool(all(eq))
    rec8, priv8, _ = point_step(fmpq(8), 40, 6, 192, log)
    t0, _ = tangent_step(priv8, 40, 6, 192, log, psi_zero=True)
    rP['psi0_certificate'] = dict(status=t0['status'], twoL_hi=t0['twoL_hi_exact'], eta_Z=t0['eta_Z_exact'])
    rP['psi0_2L_is_1'] = bool(t0['status'] == 'CERTIFIED' and t0['twoL_hi_exact'] == '1/1' and t0['eta_Z_exact'] == '0/1')
    if rec_top and rec_top.get('status') == 'CERTIFIED':
        dtop = Q(rec_top['delta_lo_exact'])
        rP['top_node_tail_psi0_equals_X1_form'] = bool(
            by.tl_tail_bound(jq(K_TOP), dtop, fmpq(1)) == jq(K_TOP) / (4 * dtop ** 2))
    else:
        rP['top_node_tail_psi0_equals_X1_form'] = False
    P['psi0'] = rP
    oks['psi0'] = bool(all(v for v in rP.values() if isinstance(v, bool)))
    log(f"psi = 0: {rP}")

    # --- structural asserts (each must fire)
    S = {}

    def fires(fn):
        try:
            fn()
            return False
        except AssertionError:
            return True
    w0 = [fmpq(0)] * 11
    S['K_eq_0_refused'] = fires(lambda: by.certify_tangent(w0, 10, 0, fmpq(1), fmpq(0), prec=128))
    S['tangent_degrees_K0_refused'] = fires(lambda: by.tangent_degrees(10, 0))
    wq8 = [Q(x) for x in priv8['wq']]
    S['Z_shape_refused'] = fires(lambda: by.certify_tangent(wq8, 40, 6, Q(priv8['CS']), Q(priv8['eta']), prec=128,
                                                            zq_override=[fmpq(0)] * 7))
    S['Ansatz_shape_refused'] = fires(lambda: by.Ansatz([fmpq(0)] * 5, 2, 1))
    S['cover_gap_refused'] = fires(lambda: by.check_cover([(fmpq(0), fmpq(1, 2), fmpq(0), fmpq(1))]))
    S['cover_overlap_refused'] = fires(lambda: by.check_cover([(fmpq(0), fmpq(3, 4), fmpq(0), fmpq(1)),
                                                               (fmpq(1, 2), fmpq(1), fmpq(0), fmpq(1))]))
    d8 = Q(rec8['delta_lo_exact'])
    S['tail_bracket_refused'] = fires(lambda: by.tl_tail_bound(fmpq(8), d8, 1 + d8))
    S['cell_bracket_refused'] = fires(lambda: by.tl_cell_bound(fmpq(8), fmpq(16), S_HALF_LO, d8, 10 * (1 + d8)))
    S['cell_s_not_lower_refused'] = fires(lambda: by.tl_cell_bound(fmpq(8), fmpq(16), fmpq(71, 100), d8, fmpq(1, 2)))
    fake = {0: dict(point=rec8, tangent=t0), 2: dict(point=rec8, tangent=t0)}
    S['assembler_gap_refused'] = fires(lambda: build_result(fake, {}, 'x', ks=(0, 2), mock=True))
    P['STRUCT'] = S
    oks['STRUCT'] = bool(all(S.values()))
    log(f"structural asserts fire: {S}")

    # --- RESUME on a cheap test grid (k = 0..4 at (24, 4, 128), no K1 steps): same runner path; exit codes 0 / 3 / 4
    rd = os.path.join(scratch, 'resume')
    os.makedirs(rd)
    oa, ob, oc, od = (os.path.join(rd, x) for x in ('A.json', 'B.json', 'C.json', 'D.json'))
    nsteps = len(main_steps(test_ks(), with_k1=False))
    t = time.time()
    rA = subprocess.run(me + ['--_grid', oa], capture_output=True, text=True)
    tA = time.time() - t
    pB = subprocess.Popen(me + ['--_grid', ob], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    killed_at, killed_pid = None, pB.pid
    while pB.poll() is None:
        try:
            nd_ = len(json.load(open(ob + '.ckpt'))['done'])
        except Exception:
            nd_ = 0
        if nd_ >= 3:
            pB.kill()                                    # our own child, by its handle
            pB.wait()
            killed_at = nd_
            break
        time.sleep(0.05)
    rB2 = subprocess.run(me + ['--_grid', ob], capture_output=True, text=True)
    codesC = []
    while not os.path.exists(oc) and len(codesC) < nsteps + 4:
        codesC.append(subprocess.run(me + ['--_grid', oc, '--_budget', '1'], capture_output=True, text=True).returncode)
    rd_ = lambda p, d: open(p, 'rb').read() if os.path.exists(p) else d
    bA, bB, bC = rd_(oa, b'A'), rd_(ob, b'B'), rd_(oc, b'C')
    stD = json.load(open(oa + '.ckpt'))
    stD['instrument_sha256'] = '0' * 64
    stD['done'].pop(list(stD['done'].keys())[-1])
    atomic_json(od + '.ckpt', stD)
    rD = subprocess.run(me + ['--_grid', od], capture_output=True, text=True)
    rA2 = subprocess.run(me + ['--_grid', oa], capture_output=True, text=True)
    bA2 = rd_(oa, b'A2')
    resA = json.loads(bA) if bA != b'A' else {}
    P['RESUME'] = dict(test_grid=list(test_ks()), degrees=list(test_degrees(0)), steps=nsteps,
                       uninterrupted_rc=rA.returncode, uninterrupted_seconds=round(tA, 1), killed_pid=killed_pid,
                       killed_after_steps=killed_at, resumed_rc=rB2.returncode,
                       B_invocations=json.load(open(ob + '.ckpt'))['invocations'], C_exit_codes=codesC,
                       A_eq_B=bool(bA == bB), A_eq_C=bool(bA == bC), rerun_after_finish_rc=rA2.returncode,
                       rerun_identical=bool(bA == bA2), sha_changed_rc=rD.returncode,
                       sha_changed_output_exists=os.path.exists(od), A_sha256=hashlib.sha256(bA).hexdigest(),
                       A_points=resA.get('all_point_certified'), A_tangents=resA.get('all_tangent_certified'),
                       no_pickle='checkpoints are JSON with p/q strings', stderr_A=rA.stderr[-300:])
    RS = P['RESUME']
    oks['RESUME'] = bool(rA.returncode == 0 and killed_at is not None and killed_at < nsteps and RS['B_invocations'] >= 2
                         and rB2.returncode == 0 and bA == bB and bA == bC and len(codesC) >= 2
                         and all(c == EXIT_RESUME for c in codesC[:-1]) and codesC[-1] == EXIT_FINISHED
                         and rD.returncode == EXIT_CKPT_REFUSED and not os.path.exists(od)
                         and rA2.returncode == 0 and bA == bA2 and RS['A_points'] is True and RS['A_tangents'] is True)
    ef.append('parts.rest.RESUME.sha_changed_output_exists (a changed-sha checkpoint must be refused)')
    log(f"RESUME: A==B (killed own PID {killed_pid} after {killed_at} steps, resumed rc {rB2.returncode}) {bA == bB}; "
        f"A==C (budget-chunked, exit codes {codesC}) {bA == bC}; changed-sha checkpoint rc={rD.returncode}; rerun after "
        f"finish rc={rA2.returncode} identical {bA == bA2}")

    # --- CEIL: the block stored in the mock result by build_result (the same code path --main uses on the REAL result),
    #     a fresh recomputation, and synthetic cell-binding / P1-binding variants (the branches the mock may not reach)
    if resA.get('c'):
        stored = resA.get('CEIL_independent') or {}
        fresh = ceil_compare(resA)
        vC = _variant(resA, 'cell')
        vP = _variant(resA, 'P1')
        cC, cP = vC.get('CEIL_independent') or {}, vP.get('CEIL_independent') or {}
        P['CEIL'] = dict(mock_stored=stored, mock_fresh_equals_stored=bool(fresh == stored),
                         variant_cell=dict(binding_instrument=vC.get('binding_piece'), block=cC),
                         variant_P1=dict(binding_instrument=vP.get('binding_piece'), block=cP),
                         note='CEIL also runs on the REAL result at the end of --main (result key CEIL_independent)')
        oks['CEIL'] = bool(stored.get('ok') and fresh == stored and stored.get('binding_independent') == 'tail'
                           and cC.get('ok') and cC.get('binding_independent') == 'cell 2'
                           and str(vC.get('binding_piece')).startswith('TL cell [2^2,')
                           and cP.get('ok') and cP.get('binding_independent') == 'P1'
                           and str(vP.get('binding_piece')).startswith('P1'))
        st_t = resA.get('TAIL_independent') or {}
        P['TAIL'] = dict(mock_stored=st_t, mock_fresh_equals_stored=bool(tail_compare(resA) == st_t),
                         variant_cell_ok=(vC.get('TAIL_independent') or {}).get('ok'),
                         variant_P1_ok=(vP.get('TAIL_independent') or {}).get('ok'))
        oks['TAIL'] = bool(st_t.get('ok') and P['TAIL']['mock_fresh_equals_stored'] and P['TAIL']['variant_cell_ok']
                           and P['TAIL']['variant_P1_ok'])
        E = resA.get('E4') or {}
        P['E4_mock'] = dict(shared=E.get('shared'), n_shared=E.get('n_shared'),
                            shared_all_overlap=E.get('shared_all_overlap'), ineq_all_hold=E.get('ineq_all_hold'))
        oks['E4_mock'] = bool(E.get('shared_all_overlap') and E.get('n_shared') == len(test_ks()))
        log(f"CEIL: mock {stored.get('c_instrument_dec_up')} ({stored.get('binding_instrument')}) identical "
            f"{stored.get('identical')}; cell variant binds {cC.get('binding_independent')} ok {cC.get('ok')}; P1 variant "
            f"binds {cP.get('binding_independent')} ok {cP.get('ok')}; TAIL {P['TAIL']}; E4 mock shared overlap "
            f"{E.get('shared_all_overlap')} ({E.get('n_shared')})")
    else:
        oks['CEIL'] = oks['TAIL'] = oks['E4_mock'] = False
        log(f"CEIL/TAIL: mock result A has no c ({rA.stdout[-300:]} {rA.stderr[-300:]})")

    # --- finite guard: an injected NaN in a stored record is caught
    if resA.get('nodes'):
        bad = json.loads(json.dumps(resA['nodes'][0]))
        bad['point']['delta_hi'] = 'nan'
        bad2 = json.loads(json.dumps(resA['nodes'][0]))
        bad2['tangent']['twoL_hi_exact'] = 'inf/1'
        P['finite_guard'] = dict(good=finite_node(resA['nodes'][0]), nan_caught=not finite_node(bad),
                                 inf_caught=not finite_node(bad2))
        oks['finite_guard'] = bool(P['finite_guard']['good'] and P['finite_guard']['nan_caught']
                                   and P['finite_guard']['inf_caught'])
    else:
        oks['finite_guard'] = False

    # --- --main refusal without GO-P2e; self-sha check against a scratch ledger
    if os.path.exists(GO_FILE):
        P['main_refusal'] = 'not exercised: GO-P2e exists (the controls never start --main)'
        oks['main_refusal'] = True
    else:
        mo = os.path.join(scratch, 'M.json')
        rM = subprocess.run(me + ['--main', mo], capture_output=True, text=True)
        P['main_refusal'] = dict(rc=rM.returncode, out=rM.stdout.strip()[-200:], output_exists=os.path.exists(mo))
        oks['main_refusal'] = bool(rM.returncode == EXIT_NO_GO and not os.path.exists(mo))
        ef.append('parts.rest.main_refusal.output_exists (no output without GO-P2e)')
    fake_go = os.path.join(scratch, 'GO-fake')
    open(fake_go, 'w').write('controls only\n')
    led_bad = os.path.join(scratch, 'ledger_bad.jsonl')
    led_good = os.path.join(scratch, 'ledger_good.jsonl')
    led_other = os.path.join(scratch, 'ledger_other_seal.jsonl')
    with open(led_bad, 'w') as f:
        f.write(json.dumps(dict(kind='go', seal_id=SEAL_ID, instruments=[dict(path='<lab>/bhslice/p2e.py',
                                                                               sha256='0' * 64)])) + '\n')
    with open(led_good, 'w') as f:
        f.write(json.dumps(dict(kind='go', seal_id=SEAL_ID, at='controls', instruments=[dict(
            path='<lab>/bhslice/p2e.py', sha256=sha_self())])) + '\n')
    with open(led_other, 'w') as f:                     # right sha but the P2d seal: must not count
        f.write(json.dumps(dict(kind='go', seal_id='9677804a51', at='controls', instruments=[dict(
            path='<lab>/bhslice/p2e.py', sha256=sha_self())])) + '\n')
    gb, mb = main_gate(fake_go, led_bad)
    gg, mg = main_gate(fake_go, led_good)
    go_, mo_ = main_gate(fake_go, led_other)
    eng_ok = sha_file(ENGINE_FILE) == ENGINE_SHA256
    P['selfsha'] = dict(wrong_sha_rc=gb, wrong_sha_msg=mb, right_sha_rc=gg, right_sha_msg=mg, other_seal_rc=go_,
                        other_seal_msg=mo_, engine_pin_matches=eng_ok)
    oks['selfsha'] = bool(gb == EXIT_NOT_PINNED and gg is None and go_ == EXIT_NOT_PINNED and eng_ok)
    log(f"self-sha gate: wrong sha -> {gb}; other seal -> {go_}; right sha -> {gg} ({mg[:60]}); engine pin matches {eng_ok}")
    P['seconds'] = round(time.time() - T0, 1)


def run_controls(out, part, log):
    sha = sha_self()
    R = _load_controls(out, sha)
    parts = list(PARTS) if part == 'all' else [part]
    for p in parts:
        R['oks'] = {k: v for k, v in R['oks'].items() if not k.startswith(PREFIX[p])}
        R['expected_false'] = [e for e in R['expected_false'] if not e.startswith(f'parts.{p}.')]
        R['parts'].pop(p, None)
        R.setdefault('started', {})[p] = time.strftime('%Y-%m-%d %H:%M:%S')
        dict(engine_a=lambda R_, l_: controls_engine(R_, l_, 'engine_a'),
             engine_b=lambda R_, l_: controls_engine(R_, l_, 'engine_b'),
             node_a=controls_node_a, node_b=controls_node_b, nulls_a=controls_nulls_a, nulls_b=controls_nulls_b,
             rest=controls_rest)[p](R, log)
        atomic_json(out, R, sort=False)
    return _finish_controls(out, R, log)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--controls')
    ap.add_argument('--part', default='all', choices=list(PARTS) + ['all'])
    ap.add_argument('--main')
    ap.add_argument('--_grid')
    ap.add_argument('--_budget', type=float, default=TIME_BUDGET)
    args = ap.parse_args()
    logf = None

    def log(msg):
        print(msg, flush=True)
        if logf:
            logf.write(msg + '\n')
            logf.flush()

    if args.controls:
        logf = open(os.path.splitext(args.controls)[0] + '.log', 'a')
        log(f"--- controls part {args.part} {time.strftime('%Y-%m-%d %H:%M:%S')} sha256 {sha_self()}")
        sys.exit(run_controls(args.controls, args.part, log))
    elif args.main:
        rc, msg = main_gate()
        if rc is not None:
            print(msg)
            sys.exit(rc)
        logf = open(os.path.splitext(args.main)[0] + '.log', 'a')
        log(f"--- main invocation {time.strftime('%Y-%m-%d %H:%M:%S')} sha256 {sha_self()} engine {sha_file(ENGINE_FILE)}; {msg}")
        sys.exit(run_steps(args.main, main_steps(), TIME_BUDGET, log, inputs=load_inputs()))
    elif args._grid:
        sys.exit(run_steps(args._grid, main_steps(test_ks(), with_k1=False), args._budget, log, degfn=test_degrees,
                           ks=test_ks(), est_scale=0.05, inputs=load_inputs()))
    else:
        ap.print_help()
