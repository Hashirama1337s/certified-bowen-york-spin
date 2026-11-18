# How fast can a Bowen-York black hole spin?

**A widely used recipe for spinning black holes on a computer has a spin ceiling. We proved, with a computer that
carries a guaranteed error bar, that its J/M² always stays below 0.9282004588 and, for a small enough bare mass,
rises above 0.9282004587.**

[![DOI 10.5281/zenodo.24233955](brand/doi.svg)](https://doi.org/10.5281/zenodo.24233955) ![Licence: text and data CC BY 4.0, code MIT](brand/licence.svg) [![Moki & Julio lab map](brand/lab.svg)](https://github.com/Hashirama1337s)

Moki & Julio · 2026 · version 1.0.1

**Read** [paper.md](paper.md) · **Re-run** [tools/reproduce.md](tools/reproduce.md) · **Check** [CLAIMS.md](CLAIMS.md) and [SCORECARD.md](SCORECARD.md)

**The magic trick.** The fastest-spinning member of this family does not exist: every member can be beaten. Yet the
ceiling is a real object, an endless throat, and one certified function drawn on that throat pins its spin from both
sides to ten digits.

![Top: certified spin values climbing toward the ceiling of the endless throat as j grows; bottom: their distance below the ceiling falling to about 5 billionths, never zero](figures/spin-ceiling.svg)

*A spin ceiling, proved. Top: the spin J/M² at j = 1, 2, 4, …, 2²¹ (green, each certified to within 8.2 × 10⁻⁹) and a certified upper bound for every j in each stretch (blue steps; the last covers all j beyond 2²¹). Orange: the ceiling, the spin of the endless throat at μ = 0, certified in [0.92820045875619, 0.92820045878266]. Dashed: 1, the most a settled black hole can spin, for comparison. Bottom: distance below the ceiling (log scale), from 0.69 at j = 1 to at least 4.7 × 10⁻⁹ at j = 2²¹; no setting reaches the ceiling.*

## In plain words

A spinning black hole has a mass M and an angular momentum J. Einstein's theory allows J/M² up to 1 for a black hole
that has settled down. Computer simulations need a starting snapshot of the black hole, and a common recipe for a
spinning one, due to Bowen and York, cannot reach 1. Numerical work since 1990 pointed to a ceiling near 0.928 (first
printed in 2002), and a 2008 study gave 0.928200.

We certified the ceiling. For every setting up to j ≈ 12,000 we trapped the mass between guaranteed bounds, and beyond
that we bounded the spin from above, so every bound is a proof, not an estimate. The spin ratio grows as one parameter, the "bare mass", shrinks,
but no setting reaches the top. In the limit the black hole becomes an infinitely long throat. We measured that limit
directly: J/M² lies between 0.92820045875619 and 0.92820045878266.

![Three nested certified windows for the spin ceiling: P2d, P2e inside it, and the endless throat's window inside that](figures/pinning-the-ceiling.svg)

*Pinning the ceiling. Three certified windows for the same number, each inside the one before (J/M² − 0.9282, in billionths): the earlier certificate P2d, [0.9282003665, 0.9282009626]; the final certificate P2e, [0.9282004538, 0.9282004687], enough to round the ceiling to 0.928200; the endless throat certified directly (B3T), [0.92820045875619, 0.92820045878266], 2.6 × 10⁻¹¹ wide. Ends rounded outward.*

Along the way we proved four facts by hand. The mass is a convex curve with no kinks in the bare mass. The mass of the
limiting throat is fixed by one curve on a sphere. For any finite bare mass the black hole's surface sits exactly at
the throat, with no other candidate. The limiting throat has no horizon at all, a question Dain and Gabach Clement left
open in 2008 (numerics already pointed that way).

**The idea in one picture.** A ladder whose rungs climb toward a ceiling. We weigh every rung with a certified scale,
prove that the rungs only go up, and then measure the ceiling itself. **Where the picture breaks.** The ceiling is not
a top rung. It is a different object, the infinite throat, and it needs its own certificate.

## At a glance

- **Claim (computer-assisted).** For the single, non-boosted, conformally flat, maximal Bowen-York spinning puncture
  (puncture boundary conditions, J = J_ADM, bare mass μ > 0): sup J/M_ADM² is in [0.9282004538, 0.9282004687], so it
  rounds to 0.928200; it is not attained, and it equals the spin of the extreme (μ = 0) data, certified in
  [0.92820045875619, 0.92820045878266].
- **Theorems.** At fixed J ≠ 0 the mass M(μ) is strictly increasing with slopes in (0, 1), strictly convex and C^1. The
  extreme data's mass is M_0 = (1/(6π)) ∫ (4|∇V|² + V²) dΩ = (2/3) Σ (2l+1) a_l² from the throat profile V. For μ > 0
  the throat is the only closed connected MOTS and minimal surface. The extreme data contain no closed minimal surface, MOTS or
  apparent horizon (two independent proofs).
- **Not new.** Monotonicity, slopes in [0, 1] and continuity at μ = 0 are corollaries of Dain and Gabach Clement's
  lemmas (we record them). The flux behind the mass formula is the radial Pohozaev invariant of Korevaar, Mazzeo, Pacard
  and Schoen. The sliding method is Berestycki and Nirenberg's. A ≥ 8π|J| is Dain and Reiris's. The values 0.928 and
  0.928200 are in print as numerical results.
- **Checked.** Every certificate run was a sealed test (16 HOLD, 1 FAILED); the probe values are unsealed and labelled.
  Masses, the supremum decision and the extreme spin were each reproduced by blind seats with independent code.
  Fifteen internal hostile referee reports (not external peer review), plus one on this release text. A clean-room
  reproduction from the repository alone (re-running the P2c ceiling assembly, B2V, B3T, FLUXT and the full P2e
  certificate; the other stages ship their recorded outputs) ended with "ALL CHECKS PASS". Lean checks one cell inequality.
- **Conditions.** Ball arithmetic (arb) is trusted. Existence of the solution for μ > 0 is cited. A few standard
  analysis facts are used without opening a textbook; the maximum principles are proved in the paper.

## What we do not claim

- No certified area-spin constant. The probe value A_0/(8πJ) ≈ 1.1991 is not certified; A ≥ 8π|J| is Dain and Reiris's
  theorem for the μ > 0 throats (it does not directly cover the μ = 0 data), and the probe constant is not proved sharp.
- No certified digits from the mass formula beyond the windows above; its long digit strings are probe values.
- No proof of the throat exponent (Conjecture E); a sealed numerical test agrees with it.
- One blind tangent certificate (REP6) FAILED as sealed; its tightened rerun (REP6b) held, with a thin margin.
- Nothing about boosted data, binaries, Kerr, or the evolution (event horizons).
- The final text has had one internal referee round only, and the novelty searches (to 2026-11-14) must be re-run before
  release.

## How it works

### Step 1 - the idea (for anyone)

```mermaid
flowchart LR
  A[Bowen-York recipe: bare mass mu, spin J] --> B[Certified sandwich: lower and upper solution]
  B --> C[Certified mass M and spin J/M^2 for every setting]
  C --> D[Proof: J/M^2 only grows as mu shrinks]
  D --> E[Ceiling = the mu = 0 limit, an endless throat]
  E --> F[Certified directly on the throat]
```

The true black-hole snapshot solves one equation, the Hamiltonian constraint. We cannot write its solution down, so
we trap it. We build two explicit functions, one provably below the solution and one provably above. Their tails at
large distance give a lower and an upper bound on the mass.

### Step 2 - how it works (for a curious student)

Everything depends only on one number, j = 4J/μ². The computer builds a polynomial approximation, works out its error
exactly in rational arithmetic, and shifts it up and down by more than that error. A comparison lemma then guarantees
that the true solution sits in between.

A worked example from the repository: at j = 10 the certified mass excess δ = M/μ − 1 lies in
[0.738555349918561234, 0.738555349918583462], so J/M² lies in [0.82710941638178, 0.82710941638181]. A blind seat that
never saw this code put δ(10) inside that interval, 21 digits deep.

To reach the ceiling, a "tangent lemma" turns one certified point into a straight line that the mass can never fall
below. Twenty-two points at j = 1, 2, 4, ..., 2^21 then bound the spin for every j at once.

### Step 3 - the real thing (for specialists)

The paper states nine results with their exact hypotheses (results 1-9 in [paper.md](paper.md), Section 1). The proof
appendix writes out the derivative decay at μ = 0, the C^1 convergence at the throat, a self-contained strong maximum
principle, and a repair of one step in a published proof. [CLAIMS.md](CLAIMS.md) lists every claim with its status and
every number with its source file.

## Use it as a tool

- **Questions it answers.** What is the mass of a Bowen-York spinning puncture at a given spin, with a guarantee? How
  close to the ceiling is a given setting? Does a candidate surface in these data form a horizon?
- **Reproduce.** Instructions are in [tools/reproduce.md](tools/reproduce.md). The certificate runs, their result files and the
  sealed verdicts are in `code/` and `results/`.
- **Reuse.** The engine in `code/gravity_engine/` certifies sub/super-solution pairs for axisymmetric Lichnerowicz
  equations; the tangent certificate gives one-sided global mass bounds from one node. The throat formula turns the
  extreme mass into a one-dimensional problem.

## Words you will meet

- **Bowen-York data**: a standard formula for a spinning black hole's snapshot, used to start simulations.
- **Puncture**: the point r = 0 of the snapshot, which stands for a second, mirror-image region of space (for μ > 0).
- **Bare mass μ**: the recipe's mass parameter; the physical mass M is larger.
- **ADM mass**: the total mass, read off far away from the black hole.
- **Sub/super-solution**: two explicit functions proved to lie below and above the true solution.
- **Ball arithmetic**: every number is stored with a guaranteed error radius, so results are bounds.
- **MOTS / apparent horizon**: a closed surface from which outgoing light just fails to escape; the usual marker of a
  black hole's surface in a snapshot.
- **Sealed test**: a prediction and its pass bar, hashed before the run and scored as written.
- **Seat**: a separate internal working session of the lab with its own files. Referee seats and blind seats are
  internal checks, not external peer review; a blind seat did not see our code or the values it was asked to reproduce; partial context is disclosed per test in the paper (Section 9).
- **Lean**: an interactive proof assistant whose kernel checks every step of a formal proof.

## Open doors

What this work leaves open, and what we noticed along the way but did not follow up, for anyone who wants to take the next step.

- **Nudge the extreme throat: does a horizon always appear.** We proved that the extreme (μ = 0) data have no horizon.
  Dain and Gabach Clement also expected the second half of the picture: that a small change of extreme data, in
  particular of the extreme Bowen-York data, always produces a horizon, and they asked for a proof or a disproof. Along
  the Bowen-York family itself it holds: for every bare mass μ > 0 the throat r = μ/2 is a minimal surface. For general
  small changes it is not settled. In 2010 they also asked whether small deformations of extreme Kerr data contain
  trapped surfaces. *Status:* open according to Dain and Gabach Clement 2009 (arXiv:0806.2180, Section 6; no later answer
  found in our search of 17 November 2026). *Start here:* the two no-horizon proofs in paper.md Sections 8.2-8.3; an idea
  we did not test is that the sliding argument needs no explicit spacetime, which is the situation of deformed data.
  *Ceiling:* a theorem that extreme spinning data sit exactly on the edge between "no horizon" and "horizon".
- **How fast the spin approaches its ceiling.** A heuristic matched-asymptotics argument at the throat (Conjecture E)
  suggests M(μ) − M_0 ~ μ^{2ν} with ν ≈ 1.38201913608948921114, set by the lowest eigenvalue of one operator on the
  sphere. It is a heuristic, not a proof; one sealed numerical check agrees (local exponent 0.882059 against the predicted
  0.882019). Whether Lovelace et al.'s fitted power law is a two-digit view of the same exponent is not established.
  *Status:* UNTESTED observation from this work (a heuristic with one sealed numerical check). *Start here:* paper.md
  Section 5 (Conjecture E and the B2V test) and the throat profile equation of Section 7. *Ceiling:* a proved law for how
  the spin of these black holes nears 0.9282.
- **A certified area-spin constant.** Every Bowen-York throat with bare mass μ > 0 obeys A ≥ 8π|J| (Dain and Reiris).
  Our probe value for the best constant in this family is A_0/(8πJ) ≈ 1.1991, the value at the extreme throat; numerical
  work already showed Bowen-York horizons staying about 20% away from the extreme-Kerr equality. We proved only one
  direction (every μ > 0 throat has area above A_0); that the constant is sharp, and its digits, are not certified.
  *Status:* open as far as our dated search found (17 November 2026; arXiv, INSPIRE and Zenodo: no proved or certified
  family-wide constant). *Start here:* paper.md Section 8.4 (area corollaries and probe values) and the one-dimensional
  throat problem of Section 7, where A_0 = ∫ V⁴ dΩ. *Ceiling:* a sharp, certified bound on how close Bowen-York black
  holes can come to extreme Kerr.
- **A one-line mass for related throats.** The flux identity of Section 7 turns the mass of the extreme throat into an
  integral over one sphere. We noticed that the same argument seems to apply to other single throats whose equation does
  not change along the throat, such as the maximally sliced spinning "trumpet" data of Hannam, Husa and Ó Murchadha
  (2009); boosted data do not qualify. This was a quick derivation, not refereed and not tested. *Status:* UNTESTED
  observation from this work. *Start here:* paper.md Section 7 (the general sphere equation L_S V − V/4 + C V^-p = 0 and
  its extreme Reissner-Nordström control) and Appendices A-B for the two limit steps. *Ceiling:* the mass of a family that
  simulations already use, as a one-dimensional integral with a proof.

## How to cite

```
Moki & Julio (2026). Certified spin, mass and horizons of the Bowen-York spinning puncture family. Version 1.0.1. Zenodo. https://doi.org/10.5281/zenodo.24233955
```

The DOI always resolves to the newest archived version. Version 1.0 is 10.5281/zenodo.24233956. Machine-readable: [CITATION.cff](CITATION.cff). Text and data: CC BY 4.0 (`LICENSE`). Code: MIT (`LICENSE-CODE`).

v1.0.1: adds Open doors and figures; corrections dated 2026-11-17; no change to the certified results. Correction 2026-11-17: the "Re-run" and "Reproduce" links pointed to `code/RUN.md`, which never existed in this repository; they now point to [tools/reproduce.md](tools/reproduce.md).

## Related projects

- [certified-page-time](https://github.com/Hashirama1337s/certified-page-time): certified Page evaporation constants (entropy ratio and Page time) of a Schwarzschild black hole emitting photons and gravitons
- [certified-evans-coleman](https://github.com/Hashirama1337s/certified-evans-coleman): a certified Evans-Coleman critical solution of radiation collapse, with exactly one growing spherical mode
- [certified-kummer-k3](https://github.com/Hashirama1337s/certified-kummer-k3): a certified closeness bound for a Kummer K3 surface
- [certified-chen-lebrun-weber](https://github.com/Hashirama1337s/certified-chen-lebrun-weber): a certified approximation of the Chen-LeBrun-Weber Einstein metric
- All projects: [the Moki & Julio lab map](https://github.com/Hashirama1337s)

---

*Prepared with AI assistance; what was checked, and by what, is stated above.*
