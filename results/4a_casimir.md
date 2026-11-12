# Result 4a: the negative energy itself, the Casimir sum on a lattice, exact, 2026-11-11 (v4 of the program, the evening)

**First line:** the Casimir effect is real, measured since 1997, and its negative energy density is of the order of nanojoules per
cubic metre at micron spacings; nothing here makes it larger. This result computes the textbook case for one scalar field as a finite sum
with no divergent subtraction at all, and finds, in exact arithmetic, the law by which the lattice sum approaches the continuum:
1 + (5 pi^2/24)/N^2 + O(1/N^4).

## What was computed

A massless field between two perfectly reflecting plates a distance d apart has less zero-point energy than the open vacuum. For the
electromagnetic field (two polarizations) E_C / A = -pi^2 hbar c / (720 d^3) per unit area (Casimir, 1948); for one massless scalar with
Dirichlet plates it is half that, -pi^2 hbar c / (1440 d^3). That deficit is the "negative energy" that every traversable-wormhole and
warp-drive paper in the record leans on; the Maldacena-Milekhin glue of computation 1 is its fermionic version along the field lines.

On a cubic lattice with spacing a = 1 (hbar = c = a = 1), the zero-point energy of one scalar field in a slab of N spacings between the
plates (Dirichlet in x) and a periodic transverse box of M x M sites is the finite sum

    E(N) = (1/2) SUM_{n=1}^{N-1} SUM_{m=0}^{M-1} SUM_{l=0}^{M-1} omega,   omega^2 = 4 sin^2(pi n / 2N) + 4 sin^2(pi m / M) + 4 sin^2(pi l / M).

E(N) = e_bulk N + 2 e_surf + E_C(N) + ...: the bulk and surface pieces are the same for every N, so the second difference

    F2(N) = E(N+1) - 2 E(N) + E(N-1)

is the Casimir energy's own second difference, with nothing divergent subtracted. In the continuum, with A = M^2,

    F2(N) / M^2  ->  -(pi^2 / 1440) [ (N+1)^-3 - 2 N^-3 + (N-1)^-3 ]  ~  -pi^2 / (120 N^5),

and the lattice's deviation from that is itself exact and printed as a ratio. The transverse box is kept at M >= 6 N so that its own
finite-size term, which goes as e^{-pi M/N}, is below the arithmetic (M = 128 and 256 agree to 8e-5 in (ratio - 1) N^2 at N = 20, where
M/N is 6.4 and 12.8).

**What "exact" had to mean here, learned the hard way.** F2 at N = 32 is one part in 10^10 of E(N), so every mode frequency, each used
M^2 times coherently, must carry far more than the printing precision. The first three versions of the program computed the
frequencies in Q.40 and lost digits past N = 16: version 2, which ran on both machines at 7:39 PM and was published for two hours, printed
1.00262 at N = 32 where the true value is 1.002017; an independent review (another AI session with its own float and long-double runs and
an analytic derivation of the coefficient) caught it the same evening. Version 4 computes 4 sin^2(k/2) = 2(1 - cos k) with a cosine in Q.60
(128-bit products, pi as floor(pi 2^60) = 3622009729038561421), omega^2 as an exact Q.60 sum, omega as the integer square root of a
124-bit radicand, the sums in 128-bit accumulators, the second difference on the raw sums, the continuum's second difference as the exact
rational Nn / D with D = N^3 (N+1)^3 (N-1)^3, and the ratio as one 128-bit division: its only roundings are pi^2's own (1e-13) and F2's
single rounding to Q.40 for printing. A self-test prints 1 + cos 0, 1 + cos pi, 4 sin^2(pi/4) and sqrt 4 in Q.60.

- Program: `code/wh_casimir.c` (v4) on `code/whq.h` (v2); binary md5 4340931e (ARM build).
- Grid: M = 128 for N = 4..20 and M = 256 for N = 4..32 (M >= 6 N throughout).
- Run: the pair, both machines at once, 9:33 to 9:36 PM on 11-11 (2 minutes 15 seconds), the two outputs identical line for line and
  equal, digit for digit, to the PC build of the same source; chain 98ce8074ad037b40; the raw output `results/4a_casimir_board1.txt`.

## The numbers

M = 256 (A = 65536), MEASURED on the pair (the M = 128 rows for N = 4..20 are in the raw output and agree to 2e-9 in the ratio):

| N | E(N)/M^2 | F2(N), the whole box | F2(N)/M^2 | continuum, per unit area | ratio lattice/continuum | (ratio - 1) N^2 |
|---|---|---|---|---|---|---|
| 4 | 3.594731192354 | -7.363033365061 | -0.000112350972 | -0.000094494862 | 1.188963828182 | 3.023421 |
| 6 | 5.982414094820 | -0.792898293997 | -0.000012098667 | -0.000011351365 | 1.065833706021 | 2.370013 |
| 8 | 8.370035259392 | -0.177039419271 | -0.000002701407 | -0.000002610949 | 1.034645673856 | 2.217323 |
| 12 | 13.145249318011 | -0.022367806908 | -0.000000341306 | -0.000000336344 | 1.014749695393 | 2.123956 |
| 16 | 17.920456114246 | -0.005233439645 | -0.000000079856 | -0.000000079207 | 1.008178201236 | 2.093620 |
| 20 | 22.695661420178 | -0.001703803092 | -0.000000025997 | -0.000000025863 | 1.005199857676 | 2.079943 |
| 24 | 27.470866267894 | -0.000682321580 | -0.000000010410 | -0.000000010373 | 1.003598277148 | 2.072608 |
| 32 | 37.021275526506 | -0.000161355786 | -0.000000002462 | -0.000000002457 | 1.002016970339 | 2.065378 |

**The law.** The review derives the coefficient of the 1/N^2 term analytically from the lattice dispersion's |k|^5 term:
c = 5 pi^2 / 24 = 2.05617 (COMPUTED by the reviewer; its verification and the 1/N^4 coefficient, fitted so far as about 9.6, are to
follow and will be added here with their labels). The exact values above read 2.1240, 2.0936, 2.0799, 2.0726, 2.0654 at N = 12, 16, 20, 24,
32, falling toward c as c + 9.6/N^2 would (2.123, 2.094, 2.080, 2.073, 2.066), and they agree with the reviewer's own exact evaluations at
N = 12, 16, 24 and 32 (2.1240, 2.0936, 2.0724, 2.0654) to 3e-7 in the ratio. The lattice's dispersion omega = 2 sin(k/2) softens every
mode relative to omega = k, so the lattice sum carries slightly more negative energy than the continuum at any finite N, by c/N^2. The bulk
energy per site, the slope of E(N)/M^2 in N, is 1.19384 from the table, and the surface term is what remains after it: both cancel
exactly in F2, which is why no subtraction is needed.

## Error bound, on paper

Each omega is exact to one unit of 2^-60 (the integer square root's floor); the sums have up to 33 x 65536 terms, so the worst-case
accumulated floor bias on E(N) is 2e6 x 2^-60 = 2e-12 and it cancels in the second difference to the extent the counts match (it is
linear in N: exactly); the random part is below 1e-15. The frequencies carry the cosine's twelve floor roundings in Q.60, 1e-17 each. F2
is therefore exact to its printed twelve places at every N here, and the ratio to pi^2's rounding. The two machines' bit-equality checks
the arithmetic, not this bound.

## In SI, for the record (DERIVED, not computed here)

Electromagnetic: E_C / A = -pi^2 hbar c / (720 d^3) with hbar c = 3.16153e-26 J m: at d = 1 um, -4.33e-10 J/m^2 (a pressure of 1.3 mPa);
at d = 10 nm, -4.33e-4 J/m^2 (1.3 kPa). The one scalar computed here gives half of each. The energy density between the plates,
-pi^2 hbar c / (720 d^4) electromagnetic: -4.3e-4 J/m^3 at 1 um. A traversable wormhole of a metre needs negative energy of planetary-mass
scale (computation 1's own papers say so); this is the gap, stated in numbers.

## Lineage

v1 (7:35 PM, never run): the electromagnetic coefficient pi^2/720 for a scalar sum, and each E(N) truncated to Q.40 before the second
difference. v2 (ran 7:39 PM, both machines identical, published 7:58 PM): pi^2/1440, the second difference on the raw sums, M >= 6N; its
Q.40 frequencies limited it to about five digits of F2 at N = 16 and three at N = 32, which the review quantified. v3 (never run): the
roots with their remainders in Q.60 sums and the exact continuum rational; the frequencies still Q.40. v4: the frequencies in Q.60. The
published v2 table is superseded by the one above; the chain hash and the raw output file are v4's, from the pair's run at 9:36 PM.
