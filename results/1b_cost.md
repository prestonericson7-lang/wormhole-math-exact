# 1b, a reading of the paper's wormhole at every size: what it costs, 2026-11-12

**First line:** this computes a geometry on paper; it does not create one. The table below is Maldacena and Milekhin's own construction
(arXiv:2008.06618) evaluated at mouth sizes other than the one they chose, with our exact constant in place of their estimate. Every row is
COMPUTED on a PC from their formula chain; nothing on this page is a new solution of the equations, and nothing on it is within reach of
any technology, as the columns show.

## Why the size is a free choice

The paper's section 3.3 fixes the mouth radius by one condition only: a human must survive the tidal acceleration, about 20 g over a
0.5 m body, which gives re = 1.5e7 m. Nothing else in the construction sets the size. The construction itself (section 3.2) needs the mouth
far larger than the fifth dimension's size R5 = 50 micrometres, and the transition scale l_B small against the mouth (lhat = l_B/re small),
so a mouth that nothing living passes through can be smaller, down to about a hundred times R5. This page asks what each size costs.

## The chain (the paper's, with our constant)

With R5 = 5e-5 m, G4 = l_p^2 (l_p = 1.616255e-35 m), G5 = G4 R5/2 (the paper's 3.23), and the coupling g4 = e (alpha = 1/137.036):

- ln_eff = ln(L_IR/R5) solves ln_eff = (1/2) ln(re/R5) + (1/4) ln(6 ln_eff) + kappa_flat, with kappa_flat = 0.128210970524 (`1a_flow.md`; the paper
  took the first term alone, its 3.25);
- Q2/R5^2 = (re/R5)^2/ln_eff; the fixed point w = (3 + sqrt(9 + 96 Q2))/(4 Q2), S* = R5/sqrt(w), L3 = R5/sqrt(3 + w/4); c2 = 6 pi L3 S*^2/G5;
- l_B = sqrt(re R5) (6 ln_eff)^(1/4) and lhat = l_B/re;
- the throat length ell = 16 re^3/(G4 c2) (the paper's 2.12), the binding energy E_bin = -c2/(16 ell) (2.13), gamma = ell/re (2.14);
- the integer magnetic charge q = re g4/(sqrt(pi) l_p) (2.5), in Dirac units;
- the field at the mouth, the extremal horizon's, B = c^2 sqrt(mu0/(4 pi G))/re, which is the paper's (2.6) in SI and does not depend on g4;
- the tidal acceleration on a 0.5 m body, c^2 (0.5 m)/re^2 (3.26); the proper time to cross, pi re/c; the time to cross as seen from outside, pi ell/c.

At the paper's size this chain returns the numbers of `1a_flow.md`: ln_eff 14.4574779, l_B 83.578 m, c2 6.709e72, ell 3,257 light-years,
gamma 2.054e12, E_bin -4.788e9 kg, and the tidal acceleration 20.4 g, the paper's 20 g.

## The table

| mouth radius re | ln_eff | lhat | charge q (Dirac units) | field at the mouth | c2 | throat length | time to cross, from outside | proper time to cross | binding energy | tides on a 0.5 m body | what bites |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 mm (100 R5, the construction's floor) | 3.167 | 0.21 | 5.3e31 | 7.0e20 T | 4.7e63 | 1.61 m | 17 ns | 0.05 ns | -6.5e19 kg | 1.8e20 g | lhat 0.21, the small-lhat premise is marginal; no human; B is 1.6e11 times the Schwinger field |
| 5 cm | 4.400 | 0.072 | 5.3e32 | 7.0e19 T | 4.1e64 | 189 m | 2.0 us | 0.5 ns | -4.7e18 kg | 1.8e18 g | no human; B is 1.6e10 times the Schwinger field |
| 1 m | 5.975 | 0.017 | 1.1e34 | 3.5e18 T | 7.0e65 | 88 km | 0.92 ms | 10 ns | -1.7e17 kg | 4.6e15 g | no human; B is 7.9e8 times the Schwinger field |
| 1 km | 9.546 | 6.2e-4 | 1.1e37 | 3.5e15 T | 5.5e68 | 1.1e8 km | 19 min | 10 us | -1.1e14 kg | 4.6e9 g | no human; B is 7.9e5 times the Schwinger field |
| 1,000 km | 13.08 | 2.1e-5 | 1.1e40 | 3.5e12 T | 4.7e71 | 13.8 light-years | 43 years | 10 ms | -7.9e10 kg | 4.6e3 g | no human; B is 790 times the Schwinger field |
| 15,000 km (the paper's) | 14.46 | 5.6e-6 | 1.6e41 | 2.3e11 T | 6.7e72 | 3,257 light-years | 10,230 years | 0.16 s | -4.8e9 kg | 20 g | the human limit; B is 53 times the Schwinger field |
| 1 astronomical unit | 19.13 | 6.0e-8 | 1.6e45 | 2.3e7 T | 5.8e76 | 3.7e11 light-years | 1.2e12 years | 26 min | -3.6e5 kg | 2e-7 g | the field alone |

The scaling laws behind the columns, for small lhat: q grows like re; B falls like 1/re; c2 grows like re/sqrt(ln_eff); the throat length
grows like re^2 sqrt(ln_eff)/R5, and with it the outside crossing time; the binding energy's size falls like 1/(re ln_eff); the tides fall
like 1/re^2. The Schwinger field, 4.4e9 T, is quoted only as a scale of the field's own energy density; a static magnetic field does not
make pairs by itself.

## What the table says

- There is no cheap size. Shrinking the mouth shortens the throat and the crossing time but raises the field at the mouth and the negative
  energy the throat needs; a 5 mm mouth wants minus 6.5e19 kg of binding energy, the paper's wants minus 4.8e9 kg, and only a mouth the size
  of the Earth's orbit gets the negative energy down to tonnes, at the price of a crossing time longer than the age of the universe.
- Every size needs a magnetic monopole charge of 1e31 to 1e45 Dirac units. No magnetic monopole has ever been observed.
- Every size needs a field at the mouth of at least 2e7 T. The strongest steady fields made are below 50 T, the strongest pulsed ones about
  1,000 T, and the strongest known in nature, on magnetars, about 1e11 T.
- Even the paper's wormhole is not a shortcut for the outside world: the traveller ages 0.16 s while the outside ages 10,230 years, because
  the time inside runs slow by gamma = 2e12. By construction nothing beats light through it, so it is not a time machine; it is a one-way
  trip into the future.
- The negative energy has no known source at any of these scales. The construction gets it classically from the fifth dimension's geometry;
  that fifth dimension, 50 micrometres in size, is itself at the edge of what gravity experiments allow.

What would change the table is not engineering but physics: different matter, a different source of negative energy, or a different
construction. That is the open question named in `PLAN.md`, and this page is the honest yardstick for it: any proposal has to beat these
columns.
