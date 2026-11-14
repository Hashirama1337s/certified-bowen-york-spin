"""REP6b (lead seal 3e929bb723): tighter tangent certificate at j = 2^21.
Reuses the REP6 instrument functions from rep6.py UNCHANGED (formulas, interpolation certificate, sandwich logic);
changes only: (i) a higher-resolution sub-solution u_- (N = 96) so V1 is closer to the true potential and eta_m is
~30x smaller; (ii) a larger tangent basis (Ne = 200); (iii) an alternating split policy; (iv) KT and KZ run through the
machinery.  Usage: rep6b.py solve | cert <case> <target> | tangent | kz | final   (exit 0 done, 3 resume)."""
import sys, os, json, time
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import rep6 as R
import rep6num as RN
from flint import arb

PREC, NSAMP = R.PREC, R.NSAMP
J21 = 2.0 ** 21
PATCHES_FROM = None


def sf(case):
    return os.path.join(HERE, f'certb_{case}_p{PREC}_n{NSAMP}.json')


def run_cert(case, fam, xlo, nx_init, target, onesided, budget=480):
    """Same acceptance logic as rep6.run_cert; split policy alternates (split the longer side, c counted 1:1)."""
    f = sf(case)
    if os.path.exists(f):
        st = json.load(open(f))
    else:
        pts = [Fraction(xlo) + Fraction(1 - xlo, nx_init) * i for i in range(nx_init + 1)]
        st = dict(todo=[[str(pts[i]), str(pts[i + 1]), '0', '1'] for i in range(nx_init)], done=[], failed=[])
        if PATCHES_FROM:
            st['todo'] = [p['patch'] for p in json.load(open(sf(PATCHES_FROM)))['done']]
    t0 = time.time()
    while st['todo'] and time.time() - t0 < budget:
        p = st['todo'].pop(0)
        xa, xb, ca, cb = [Fraction(s) for s in p]
        t1 = time.time()
        rec = R.certify_patch(fam, xa, xb, ca, cb)
        good = rec is not None and float(rec['sup']) <= target and (onesided or float(rec['inf']) >= -target)
        if good:
            st['done'].append(dict(patch=p, rho=rec['rho'], M=rec['M'], err=rec['err'],
                                   sup=str(R._ball_upper(rec['sup'])), inf=str(R._ball_lower(rec['inf'])),
                                   Einf=str(R._ball_lower(rec['Einf'])), sec=round(time.time() - t1, 2)))
        else:
            if max(xb - xa, cb - ca) < Fraction(1, 2048):
                st['failed'].append(dict(patch=p, sup=None if rec is None else float(rec['sup'])))
                print('FAIL', p, flush=True)
            elif (xb - xa) >= (cb - ca):
                m = (xa + xb) / 2
                st['todo'][0:0] = [[str(xa), str(m), p[2], p[3]], [str(m), str(xb), p[2], p[3]]]
            else:
                m = (ca + cb) / 2
                st['todo'][0:0] = [[p[0], p[1], str(ca), str(m)], [p[0], p[1], str(m), str(cb)]]
        json.dump(st, open(f, 'w'))
    print(case, 'done', len(st['done']), 'todo', len(st['todo']), 'failed', len(st['failed']), flush=True)
    return not st['todo']


def state(case):
    st = json.load(open(sf(case)))
    assert not st['todo'], case
    return st


def main_fam(name):
    d = R.load_sol(name)
    return R.MainF(R.main_sol(d), R.CS_of(d))


def sandwich(name):
    st = state(name)
    assert not st['failed']
    R.cover_check(st, True)
    d = R.load_sol(name)
    g0 = Fraction(R.unhex(d['g0']))
    sup = max(Fraction(p['sup']) for p in st['done']); inf = min(Fraction(p['inf']) for p in st['done'])
    Einf = min(Fraction(p['Einf']) for p in st['done'])
    eta_p = R.dyadic_up(2 * sup * (1 + Fraction(1, 2 ** 20)) + Fraction(1, 2 ** 200), 96)
    eta_m = R.dyadic_up(-2 * inf * (1 + Fraction(1, 2 ** 20)) + Fraction(1, 2 ** 200), 96)
    assert eta_p > 2 * sup and eta_m > -2 * inf and Einf - 2 * eta_m > 0
    return dict(g0=g0, eta_p=eta_p, eta_m=eta_m, sup=sup, inf=inf, Einf=Einf, patches=len(st['done']),
                delta_lo=2 * (g0 - eta_m), delta_hi=2 * (g0 + eta_p), N=d['N'], K=d['K'])


def tan_fam(info, zeta_arb=None):
    d = R.load_sol(info['main'])
    sol = R.main_sol(d)
    es = R.Sol(0.0, np.array([R.unhex(s) for s in info['he']]).reshape(info['Ne'], info['Ke'] + 1), even=info['even'])
    j = arb(R.unhex(d['j']))
    z = arb(R.unhex(info['zeta0'])) if zeta_arb is None else zeta_arb
    return R.TanF(sol, Fraction(info['eta_m']), es, z, j * j, info['dscale'])


def tan_result(case, info, extra=Fraction(0)):
    st = state(case)
    R.cover_check(st, False)
    supT = max(Fraction(p['sup']) for p in st['done'])
    sw = sandwich(info['main'])
    kappa = Fraction(0) if supT <= 0 else R.dyadic_up(2 * supT * (1 + Fraction(1, 2 ** 20)) + Fraction(1, 2 ** 200), 96)
    zeta0 = Fraction(R.unhex(info['zeta0'])) + extra
    Zinf = (sw['g0'] if info['dscale'] else 0) + zeta0 + kappa
    return dict(supT=supT, kappa=kappa, zeta0=zeta0, twoL_hi=1 + 2 * Zinf, patches=len(st['done']),
                failed=len(st['failed']), Eminus=min(Fraction(p['Einf']) for p in st['done']))


class KTF:
    """KT through the machinery: psi = 0 gives That(0) = -(63/32) j^2 (1-c^2) A^11 w^8 (1-xi)^4 / E_-^8
    = q * F with q = (1-c^2)(1-xi)^4 w^8 >= 0 on the real domain; the machinery certifies sup F < 0 (=> That(0) <= 0,
    kappa = 0, 2L = 1).  The factorisation is checked against rep6.tan_point at sample points (ktf_check)."""

    def __init__(self, sol, eta_m, j2):
        self.sol = sol; self.eta_m = eta_m; self.j2 = j2

    def grid(self, xa, xb, ca, cb, xs, cs, cplx):
        loc = self.sol.local(xa, xb, ca, cb)
        Hg, _, _, _, xi, cc = R.fields(self.sol, loc, xs, cs, cplx)
        g0m = arb(self.sol.g0) - R.arb_exact(self.eta_m)
        k = arb(63) / 32 * self.j2
        Rr = []; Eg = []; Aa = []
        for i, x in enumerate(xi):
            A = 2 * (1 + 6 * x * x + x ** 4); w = 1 - x * x; w4 = w ** 4
            rowR = []; rowE = []
            for l, c in enumerate(cc):
                Em = A * A + 2 * w4 * (g0m + w4 * Hg[i, l])
                rowR.append(-k * A ** 11 / Em ** 8); rowE.append(Em)
            Rr.append(rowR); Eg.append(rowE); Aa.append(A)
        return {'R': Rr, 'E': Eg, 'A': Aa}


def ktf_check():
    """F * q == That(psi = 0) from rep6.tan_point at random real points (arb)."""
    import random
    random.seed(1)
    d = R.load_sol('main96'); sol = R.main_sol(d); sw = sandwich('main96')
    j2 = arb(R.unhex(d['j'])) ** 2
    fam = KTF(sol, sw['eta_m'], j2)
    worst = 0.0
    for _ in range(6):
        xa = Fraction(random.randint(-15, 14), 16); ca = Fraction(random.randint(0, 14), 16)
        x0 = [arb(random.uniform(-1, 1))]; y0 = [arb(random.uniform(-1, 1))]
        g = fam.grid(xa, xa + Fraction(1, 16), ca, ca + Fraction(1, 16), x0, y0, False)
        loc = sol.local(xa, xa + Fraction(1, 16), ca, ca + Fraction(1, 16))
        Hg, Hg1, Hg2, LHg, xi, cc, Hg3, LHg1 = R.fields(sol, loc, x0, y0, False, extra=True)
        x, c = xi[0], cc[0]
        z = arb(0)
        v, _, _ = R.tan_point(x, c, arb(sol.g0), arb(sol.g0) - R.arb_exact(sw['eta_m']), z, arb(63) / 16 * j2,
                              arb(63) / 32 * j2, Hg[0, 0], Hg1[0, 0], Hg2[0, 0], Hg3[0, 0], LHg[0, 0], LHg1[0, 0],
                              z, z, z, z, 0)
        q = (1 - c * c) * (1 - x) ** 4 * (1 - x * x) ** 8
        worst = max(worst, float(abs(v - q * g['R'][0][0]).upper() / (abs(v).upper() + 1e-300)))
    return worst


def info_file(case):
    return os.path.join(HERE, f'tanb_{case}.json')


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'solve':
        g0, h, _ = RN.solve(J21, 96, 11)
        R.save_sol('main96', dict(j=R.fhex(J21), fac_num='1', N=96, K=11, g0=R.fhex(g0), h=[R.fhex(v) for v in h],
                                  float_delta=2 * g0))
        print('main96 delta(float)', 2 * g0)
    elif cmd == 'cert':
        case, target = sys.argv[2], float(sys.argv[3])
        if case.startswith('main'):
            done = run_cert(case, main_fam(case), 0, 16, target, False)
        elif case == 'KTF':
            d = R.load_sol('main96'); sw = sandwich('main96')
            done = run_cert(case, KTF(R.main_sol(d), sw['eta_m'], arb(R.unhex(d['j'])) ** 2), -1, 32, target, True)
        else:
            info = json.load(open(info_file(case)))
            zeta = None
            if 'patches_from' in info:
                PATCHES_FROM = info['patches_from']
            if 'zeta_extra' in info:
                zeta = arb(R.unhex(info['zeta0'])) + R.arb_exact(Fraction(info['zeta_extra']))
            done = run_cert(case, tan_fam(info, zeta), -1, 32, target, True)
        sys.exit(0 if done else 3)
    elif cmd == 'tangent':     # args: case main Ne Ke [dscale]
        case, mainname, Ne, Ke = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
        dscale = int(sys.argv[6]) if len(sys.argv) > 6 else 1
        sw = sandwich(mainname)
        if dscale:
            z0, he, res = R.solve_tangent(mainname, float(sw['eta_m']), Ne, Ke, even=False)
        else:
            z0, he, res = 0.0, np.zeros((Ne, Ke + 1)), 0.0
        info = dict(main=mainname, eta_m=str(sw['eta_m']), Ne=Ne, Ke=Ke, even=False, dscale=dscale,
                    zeta0=R.fhex(z0), he=[R.fhex(v) for v in np.asarray(he).ravel()], float_colloc_res=float(res))
        json.dump(info, open(info_file(case), 'w'))
        print(case, 'zeta0', z0, 'float res', res, '2L float', 1 + 2 * (float(sw['g0']) + z0) if dscale else 1)
    elif cmd == 'kz':          # args: source_case  -> builds KZ info with kappa_Z = 2 sup (1 - 2^-10)
        src = sys.argv[2]
        info = json.load(open(info_file(src)))
        tr = tan_result(src, info)
        kz = 2 * tr['supT'] * (1 - Fraction(1, 2 ** 10))
        info['zeta_extra'] = str(kz)
        info['patches_from'] = src
        json.dump(info, open(info_file('KZ'), 'w'))
        print('kappa_Z', float(kz), 'supT', float(tr['supT']))
    elif cmd == 'show':
        case = sys.argv[2]
        if case.startswith('main'):
            print({k: float(v) if isinstance(v, Fraction) else v for k, v in sandwich(case).items()})
        else:
            info = json.load(open(info_file(case)))
            extra = Fraction(info.get('zeta_extra', '0'))
            print({k: float(v) if isinstance(v, Fraction) else v for k, v in tan_result(case, info, extra).items()})


def final():
    D = R.dec_round
    out = dict(seat='REP6b', seal='3e929bb723', prec_bits=PREC, nsamp=NSAMP,
               disclosure='brief named the decision threshold chi_tail < 0.9282005 and the prior score; computation '
                          'choices (N, Ne, targets) were driven by the residual levels, the threshold was only compared at the end')
    sw = sandwich('main96')
    out['delta_recertified'] = dict(N=96, K=11, delta_lo=D(sw['delta_lo'], 13, False), delta_hi=D(sw['delta_hi'], 13, True),
                                    delta_lo_exact=R.fstr(sw['delta_lo']), delta_hi_exact=R.fstr(sw['delta_hi']),
                                    eta_plus=R.fstr(sw['eta_p']), eta_minus=R.fstr(sw['eta_m']),
                                    Rhat_sup_upper=D(sw['sup'], 6, True), Rhat_inf_lower=D(sw['inf'], 6, False),
                                    E_minus_inf_lower=D(sw['Einf'] - 2 * sw['eta_m'], 8, False), patches=sw['patches'])
    info = json.load(open(info_file('tan200')))
    tr = tan_result('tan200', info)
    assert tr['failed'] == 0 and tr['Eminus'] > 0
    out['tangent'] = dict(sub_solution='u_- of main96 (eta_minus above)', Ne=200, Ke=11,
                          twoL_hi=D(tr['twoL_hi'], 12, True), twoL_hi_exact=R.fstr(tr['twoL_hi']),
                          That_sup_upper=D(tr['supT'], 6, True), kappa=R.fstr(tr['kappa']), zeta0=R.fstr(tr['zeta0']),
                          float_twoL_estimate=1 + 2 * (float(sw['g0']) + R.unhex(info['zeta0'])),
                          float_colloc_residual=info['float_colloc_res'], patches=tr['patches'], TL_H_certified=True,
                          class_K='structural, as in REP6 (same ansatz form, Ne changed only)')
    j = Fraction(2) ** 21
    den = 1 + sw['delta_lo'] - max(tr['twoL_hi'], Fraction(0))
    chi = (j / 4) / den ** 2
    out['chi_tail'] = D(chi, 15, True)
    out['chi_tail_exact'] = R.fstr(chi)
    out['chi_tail_below_0.9282005'] = bool(chi < Fraction(9282005, 10 ** 7))
    # controls
    stk = state('KTF'); R.cover_check(stk, False)
    supF = max(Fraction(p['sup']) for p in stk['done'])
    assert supF <= 0 and not stk['failed']
    out['KT'] = dict(method='machinery on the cofactor F = That(0)/q, q = (1-c^2)(1-xi)^4 w^8 >= 0 (factorisation checked '
                            'against rep6.tan_point, worst rel. diff %.1e)' % ktf_check(),
                     F_sup_upper=D(supF, 6, True), patches=len(stk['done']), kappa='0', twoL_hi='1',
                     note='direct interpolation of That(0) itself was abandoned: That(0) attains its max 0 on c = 1 and '
                          'xi = +-1, so the coefficient-sum bound can never reach <= 0 and patches kept splitting (840 accepted, 35 pending when the seat stopped it)')
    ikz = json.load(open(info_file('KZ')))
    stz = json.load(open(sf('KZ')))
    sups = [Fraction(p['sup']) for p in stz['done']]
    out['KZ'] = dict(kappa_Z=R.fstr(Fraction(ikz['zeta_extra'])), kappa_Z_rule='2 sup That (1 - 2^-10)  [= (2/3) sup (1-2^-10) in the Lin[1] = -3/2 normalisation]',
                     patches=len(sups), patches_with_sup_bound_positive=sum(1 for s in sups if s > 0),
                     max_sup_bound=D(max(sups), 6, True), certificate_issued=bool(max(sups) <= 0), expected='FAIL')
    out['convergence'] = {}
    t160 = tan_result('tan160', json.load(open(info_file('tan160'))))
    out['convergence']['tan160_twoL_hi'] = D(t160['twoL_hi'], 12, True)
    out['convergence']['rep6_main_N84_delta'] = ['750.5607142022', '750.5607145913']
    json.dump(out, open(os.path.join(HERE, 'result_rep6b.json'), 'w'), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__' and sys.argv[1] == 'ktfcheck':
    print('KTF factorisation worst relative difference', ktf_check())
if __name__ == '__main__' and sys.argv[1] == 'final':
    final()
