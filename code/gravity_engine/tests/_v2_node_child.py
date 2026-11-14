"""Child process for tests/test_engine_v2.py: certify ONE Bowen-York node in a fresh process (so the peak memory is
this node's alone) and print one JSON line. Usage: py -3.11 _v2_node_child.py <j> <N> <K> <fit_ops> <cache_ops 0|1>"""
import sys
import os
for _v in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[_v] = '1'
sys.dont_write_bytecode = True
sys.path.insert(0, '<LAB>/POCKETLAB')
import json
import math
import time
from fractions import Fraction

from gravity_engine import lichnerowicz_axisym as la


def peak_commit_gb():
    """peak private commit of THIS process, GB (Windows psapi); None elsewhere."""
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
        return c.PeakPagefileUsage / 1024 ** 3
    except Exception:
        return None


def main():
    j, N, K, mode, cache = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5] == '1'
    P = la.BowenYorkSpinningPuncture(fit_ops=mode, cache_ops=cache)
    t0 = time.time()
    rec, _ = la.certify_node(P, Fraction(j), N=N, K=K, log=lambda m: None)
    dt = time.time() - t0
    fit = rec.get('_fit', {})
    hist = fit.get('arb_hist', [])
    out = dict(status=rec['status'], prec=la.PREC, engine_version=la.ENGINE_VERSION, fit_ops=fit.get('fit_ops'),
               arb_hist=[h if math.isfinite(h) else repr(h) for h in hist],
               any_nan=any(not math.isfinite(h) for h in hist), reached_stop=fit.get('reached_stop'),
               delta_lo_dec=rec.get('delta_lo_dec'), delta_hi_dec=rec.get('delta_hi_dec'),
               rel_width=rec.get('rel_width'), bern_pieces=rec.get('positivity', {}).get('pieces'),
               caches_empty=(not P._fops and not P._aops), time_s=round(dt, 1), peak_commit_gb=peak_commit_gb())
    print('RESULT ' + json.dumps(out))


if __name__ == '__main__':
    main()
