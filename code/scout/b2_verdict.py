"""B2 verdict instrument (for the lead to seal; NOT run on real data by the proof seat).
Input: a P2e-format JSON: {"nodes": [{"j": "2^k as rational string", "point": {"delta_lo_exact", "delta_hi_exact"}}, ...]}.
Units: m(x) = M/sqrt(J), x = mu/sqrt(J) = 2/sqrt(j), m = x (1 + delta(j)).  All bounds outward (mpmath iv, 80 digits).
Tests (theorems of B2-PROOF.md):
  T1  every chord slope between consecutive nodes lies in (0, 1)        [Thm A: 0 < (M1-M2)/(mu1-mu2) < 1]
  T2  chord slopes strictly decrease as j increases                     [Thm B: M convex in mu]
  T3  chi(j) = j/(4(1+delta)^2) strictly increases over the nodes       [Cor C]
  T4  convexity-only tail from the top two nodes (no tangent certificate):
        node b = top (largest j), node a = next.  In units mu_b = 1 at J = j_b/4: mu_a = sqrt(j_b/j_a).
        slope bound s_b <= chord(a,b) (convexity), so for every mu < 1: M(mu) >= M(1) - chord  (mu -> 0 worst),
        chi_sup <= (j_b/4) / (1 + delta_b_lo - chord_hi)^2.  Prints the bound; bar compares it with 0.9282005.
  T5  local j-exponent of chord slopes e_k = log(chord_k/chord_{k+1}) / log(j_{k+1}/j_k) (POINT estimate from interval mids)
      predicted limit nu - 1/2 = 0.88201913609 (b2_throat.log); wrong-exponent null (p = 6) 0.79642947556.
Usage: py -3.11 b2_verdict.py file.json
"""
import json, sys
from fractions import Fraction as F
from mpmath import mp, iv, log

mp.dps = 80; iv.dps = 80
NU_HALF = mp.mpf('0.88201913609')
NU_HALF_NULL = mp.mpf('0.79642947556')

def load(path):
    d = json.load(open(path))
    rows = []
    for n in d['nodes']:
        j = F(n['j'])
        p = n['point']
        rows.append((j, F(p['delta_lo_exact']), F(p['delta_hi_exact'])))
    rows.sort()
    return rows

def fr(x):
    x = F(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)

def ivl(lo, hi):
    return iv.mpf([fr(lo).a, fr(hi).b])

def analyse(rows, bar=None):
    out = {}
    X = [2 / iv.sqrt(fr(j)) for j, _, _ in rows]
    Mx = [X[i] * (1 + ivl(rows[i][1], rows[i][2])) for i in range(len(rows))]
    ch = [(Mx[i] - Mx[i + 1]) / (X[i] - X[i + 1]) for i in range(len(rows) - 1)]
    out['T1_pos'] = all(c.a > 0 for c in ch)
    out['T1_lt1'] = all(c.b < 1 for c in ch)
    t2 = [ch[i].a > ch[i + 1].b for i in range(len(ch) - 1)]
    out['T2_decided'] = f"{sum(t2)}/{len(t2)}"
    out['T2_fail_idx'] = [i for i, ok in enumerate(t2) if not ok]
    chi = [fr(j) / (4 * (1 + ivl(lo, hi)) ** 2) for j, lo, hi in rows]
    t3 = [chi[i].b < chi[i + 1].a for i in range(len(chi) - 1)]
    out['T3_decided'] = f"{sum(t3)}/{len(t3)}"
    # T4
    ja, alo, ahi = rows[-2]; jb, blo, bhi = rows[-1]
    mua = iv.sqrt(fr(jb) / fr(ja))
    Ma = mua * (1 + ivl(alo, ahi)); Mb = 1 + ivl(blo, bhi)
    chord = (Ma - Mb) / (mua - 1)
    base = 1 + fr(blo) - iv.mpf(chord.b)
    out['T4_bracket_positive'] = bool(base.a > 0)
    tail = fr(jb) / 4 / base ** 2
    out['T4_tail_chi_up'] = mp.nstr(tail.b, 22)
    out['T4_chord_slope'] = [mp.nstr(chord.a, 15), mp.nstr(chord.b, 15)]
    if bar is not None:
        out['T4_below_bar'] = bool(tail.b < mp.mpf(bar))
    # T5
    # interval exponents e_k = log(ch_k/ch_{k+1}) / log(sqrt(j_{k+2}/j_k)); "decided" if width < 2e-4
    es = []
    for i in range(len(ch) - 1):
        if ch[i + 1].a <= 0:
            es.append(None); continue
        e = iv.log(ch[i] / ch[i + 1]) / iv.log(iv.sqrt(fr(rows[i + 2][0]) / fr(rows[i][0])))
        es.append(e)
    out["T5_exponents"] = [None if e is None else (mp.nstr(e.mid, 9), mp.nstr(e.delta, 2)) for e in es]
    dec = [i for i, e in enumerate(es) if e is not None and e.delta < mp.mpf("2e-4")]
    if dec:
        top = dec[-1]; e = es[top]
        out["T5_top_decided_index"] = top
        out["T5_top_minus_pred"] = mp.nstr(e.mid - NU_HALF, 6)
        out["T5_top_minus_null"] = mp.nstr(e.mid - NU_HALF_NULL, 6)
        first = [i for i in dec if i >= 6]
        if first:
            out["T5_first_ge6_minus_pred"] = mp.nstr(es[first[0]].mid - NU_HALF, 6)
    out['top_chi_lo'] = mp.nstr(chi[-1].a, 22)
    return out

if __name__ == '__main__':
    if sys.argv[1] == '--selftest':
        # Positive control: m(x) = sqrt(x^2 + 1) -> delta(j) = sqrt(1 + j/4) - 1 (convex, increasing, slope < 1,
        # chi_sup = 1, slope ~ x ~ 2 j^-1/2 so exponent -> 0.5).  Null: add 0.02 sin(6x) x (breaks convexity).
        import math
        def mk(fun, ks):
            rows = []
            for k in ks:
                j = F(2) ** k
                v = mp.mpf(fun(mp.mpf(j.numerator) / j.denominator))
                lo = F(mp.nstr(v - mp.mpf('1e-40'), 60)); hi = F(mp.nstr(v + mp.mpf('1e-40'), 60))
                rows.append((j, lo, hi))
            return rows
        pos = mk(lambda j: mp.sqrt(1 + j / 4) - 1, range(0, 22))
        print('POSITIVE CONTROL', analyse(pos, bar='1.0000001'))
        def nul(j):
            x = 2 / mp.sqrt(j); m = mp.sqrt(x * x + 1) + mp.mpf('0.02') * x * mp.sin(6 * x)
            return m / x - 1
        print('NULL (nonconvex)', analyse(mk(nul, range(-4, 8))))
    else:
        bar = sys.argv[2] if len(sys.argv) > 2 else None
        print(json.dumps(analyse(load(sys.argv[1]), bar), indent=1))
