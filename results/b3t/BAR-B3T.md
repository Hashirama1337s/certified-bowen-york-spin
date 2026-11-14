# BAR-B3T: certified spin of the EXTREME (mu = 0) Bowen-York data, directly on the throat

Written 2026-11-14 by the lead before b3t.py exists. Raid Wing 1, boss B3. Design: scout/B3-EXTREME.md (route B; unrefereed - referee R runs in parallel on B3-1..B3-5 and B2-PROOF). Gate: GATE-WING1.md (certified mu = 0 value OPEN; numerical 0.928 at 3 digits by Hannam, Husa, O Murchadha 0908.1063).

instruments: b3t.py

## Object
The extreme solution u_0 of Delta u = -S (1 + u)^-7 (mu = 0), J fixed (scale: chi(0) = J / M_0^2 is J-independent). A certified super-solution w in B3-EXTREME's class (C^2 off 0, 1 + w > 0, w -> 0 at infinity, w -> +infinity at 0, Delta w <= -S (1 + w)^-7), built in the throat coordinates (sigma = log r, rational ansatz w = r^-1/2 G / (1 + sqrt r), xi = tanh(sigma / 8)), giving
  M_0 <= 2 lim r w   and   M_0 >= (1 / (2 pi)) Int S (1 + w)^-7 dx (rigorous quadrature),
hence a certified window for chi(0) = J / M_0^2. Comparison basis: B3-5 (super-solution above every u_mu, hence above u_0), B2/B3 monotone-limit results (M_0 = lim M(mu)), to be refereed by R before the verdict is final.

## Predictions
- T1: the certificate holds (super-solution inequality certified on an exact cover; positivity of 1 + w).
- T2: the certified chi(0) window has width <= 1e-9.
- T3: the window lies inside P2e's certified window [0.9282004538, 0.9282004687] (consistency with the family supremum = chi(0), Corollary D3).
- T4: the window contains or lies within 1e-9 of the design probe value 0.9282004587652 (double precision, B3-EXTREME Q4) - stated as the design expectation, before any certified number.

## Controls (any failure = VOID)
- RN: the same machinery on extreme Reissner-Nordstrom (exponent 3, exact M = 1) gives a window containing 1.
- NULL: source 9/4 -> 9/4 (1 + 1e-3) moves the window off P2e's window.
- Margin null: forcing the shift below the certified residual bound FAILS.
- Cover asserted exact; outward rounding; RESUME; exit codes 0 / 3 / other; memcap; no stdlib-named scripts.

## Verdict rule
B3T CERTIFIED if T1-T3 hold and controls pass; T4 scored as stated. A T3 failure = contradiction with P2e or with Corollary D3: stop and investigate.

## Seal fields
- formula_ref: scout/B3-EXTREME.md:1-220
- calibrated_on: 560847f6a4, a580a8ad0f
- depends_on: 560847f6a4, a580a8ad0f, 3e929bb723
- null_not_symmetry: NULL changes the source; RN is an exactly known case; margin null forces failure.
- target_form: T1 boolean; T2 width in [0, 1e-9]; T3 containment boolean; T4 distance in [0, 1e-9].
- ceiling: P2e window [0.9282004538 (down), 0.9282004687 (up)]; literature: Hannam, Husa, O Murchadha (2009) give a 3-digit numerical value at mu = 0, Lovelace et al. (2008) a 6-digit extrapolated one (exact quotes and locators in scout/GATE-WING1.md).
- toy_model: extreme RN (exact).
- tolerances: from need; 8 GB; one process; chunks <= 520 s.
- parameters: mu = 0, J = 1 (scale-free)
- files: BAR-B3T.md, b3t.py
- command: memcap 8 python b3t.py --main result_b3t.json
