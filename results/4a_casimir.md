# Result 4a: the negative energy itself, the Casimir sum on a lattice, exact, 2026-11-11 (v5 of the program, the evening)

**First line:** the Casimir effect is real, measured since 1997, and its negative energy density is of the order of nanojoules per
cubic metre at micron spacings; nothing here makes it larger. This result computes the textbook case for one scalar field as a finite sum
with no divergent subtraction at all, and finds, in exact arithmetic, the law by which the lattice sum approaches the continuum:
1 + (5 pi^2/24)/N^2 + [pi^2 (2600 + 291 pi^2)/5760]/N^4 + 52.5/N^6 + ... = 1 + 2.0561676/N^2 + 9.3762182/N^4 + ..., the coefficients
derived in review and matched by the exact sums to 8e-18 at N = 32.

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
124-bit radicand, the sums in 128-bit accumulators and the second difference on the raw sums. It ran on both machines at 9:36 PM and was
right to the law's seven digits, but the same review then certified the exact values (three independent routes agreeing to 25 digits:
N = 32, M = 256, ratio 1.0020169686478780) and found v4 still +1.69e-9 high at N = 32 and +2.50e-9 at N = 24, one unit in the ninth
decimal, from three roundings it located: the mode arguments formed as a floored quotient times the index (which rescales k differently
for N-1, N and N+1), the second difference rounded to Q.40 before the division, and pi^2/1440 held in Q.40. Version 5 forms every
argument as round(pi 2^60 n / N) in 128 bits, keeps the second difference unrounded in the ratio's one 128-bit division, and takes the
continuum's second difference for the whole box as the exact rational (pi^2/1440) (12 N^4 - 6 N^2 + 2) M^2 / (N^3 (N^2 - 1)^3) in Q.60
with pi^2/1440 as round(pi^2/1440 2^60) = 7901999413873821 (the self-test re-derives it from pi 2^60): its ratio at N = 32 is
1.002016968639, within 8.4e-12 of the certified value, and within 4.1e-11 at every N, the residual being the frequencies' and roots' own
floor roundings in Q.60 (ten decimals at N = 20 to 32, twelve at N = 12 and 16). A self-test prints
1 + cos 0, 1 + cos pi, 4 sin^2(pi/4), sqrt 4 and pi^2/1440 in Q.60.

- Program: `code/wh_casimir.c` (v5) on `code/whq.h` (v2); binary md5 b38a166c (ARM build; v4 was 4340931e).
- Grid: M = 128 for N = 4..20 and M = 256 for N = 4..32 (M >= 6 N throughout).
- Run: v5 on the pair, both machines at once, 10:18:15 to 10:20:25 PM on 11-11 (2 minutes 10 seconds), the two outputs identical line for
  line and equal, digit for digit, to the PC build of the same source; chain 708f274a7b9e6043; the raw output `results/4a_casimir_board1.txt`.
  (v4 had run 9:33 to 9:36 PM, chain 98ce8074ad037b40.)

## The numbers

M = 256 (A = 65536), MEASURED on the pair (the M = 128 rows for N = 4..20 are in the raw output and agree to 2e-9 in the ratio):

| N | E(N)/M^2 | F2(N), the whole box | F2(N)/M^2 | continuum, per unit area | ratio lattice/continuum | (ratio - 1) N^2 |
|---|---|---|---|---|---|---|
| 4 | 3.594731192354 | -7.363033365061 | -0.000112350972 | -0.000094494862 | 1.188963828140 | 3.023421250240 |
| 6 | 5.982414094820 | -0.792898293997 | -0.000012098667 | -0.000011351365 | 1.065833705984 | 2.370013415424 |
| 8 | 8.370035259392 | -0.177039419271 | -0.000002701407 | -0.000002610949 | 1.034645673821 | 2.217323124583 |
| 12 | 13.145249318011 | -0.022367806908 | -0.000000341306 | -0.000000336344 | 1.014749695375 | 2.123956134004 |
| 16 | 17.920456114246 | -0.005233439645 | -0.000000079856 | -0.000000079207 | 1.008178201417 | 2.093619562918 |
| 20 | 22.695661420178 | -0.001703803091 | -0.000000025997 | -0.000000025863 | 1.005199857059 | 2.079942823911 |
| 24 | 27.470866267894 | -0.000682321579 | -0.000000010410 | -0.000000010373 | 1.003598274684 | 2.072606218513 |
| 32 | 37.021275526506 | -0.000161355786 | -0.000000002462 | -0.000000002457 | 1.002016968639 | 2.065375886857 |

(The ratio column is good to ten decimals at N = 20 to 32 and twelve at N = 12 and 16, the last column to eight at N = 32; the error bound
below says why. The certified exact ratio at N = 32, M = 256 is 1.0020169686478780.)

**The law.** The review derives it (COMPUTED; a second, independent derivation agrees): per x-mode, the transverse-integrated zero-point
energy has the non-analytic part -(1/12 pi)|k|^3 + (7/480 pi)|k|^5, where 7/480 = 1/96 + 1/240, the 1/96 from the x dispersion
2 sin(k/2) and the 1/240 from the slope of the two-dimensional density of states; the Dirichlet sum turns |k|^s into zeta(-s) (pi/N)^s, so
the energy itself is E_C/A = -(pi^2/1440) N^-3 [1 + (pi^2/12) N^-2 + ...], and the second difference multiplies the N^-2 term by 5/2:

    (ratio - 1) N^2 = 5 pi^2/24 + [pi^2 (2600 + 291 pi^2)/5760]/N^2 + 52.54/N^4 + ... = 2.0561676 + 9.3762182/N^2 + 52.54/N^4 + ...

(An earlier draft of this page fitted the second coefficient as "about 9.6": a two-term fit that was absorbing the 52.5/N^4.) The exact
values at M -> infinity are 2.12396, 2.09362, 2.07994, 2.07261, 2.06537 at N = 12, 16, 20, 24, 32 (the review, certified); the table's
M = 256 values agree with them to the box's own term and the arithmetic's residual. The lattice's dispersion omega = 2 sin(k/2) softens
every mode relative to omega = k, so the lattice sum carries slightly more negative energy than the continuum at any finite N: the energy's
own excess is (pi^2/12)/N^2 = 0.822/N^2, and the second-difference ratio's is the 2.056/N^2 above. The bulk energy per site, the slope of
E(N)/M^2 in N, is 1.19384 from the table, and the surface term is what remains after it: both cancel exactly in F2, which is why no
subtraction is needed.

## Error bound, on paper

Each mode argument is rounded once to half a unit of 2^-60; each 4 sin^2(k/2) carries the cosine's twelve floor roundings in Q.60 and
each omega the integer square root's floor, so a frequency is good to about 15 units of 2^-60 (1.3e-17). The sums have up to 33 x 65536
terms and are exact; the floor biases are one-sided and cancel in the second difference to the extent the counts match (they are linear in
N), which leaves the residual the review measured against the certified values for this arithmetic: under 1e-12 in the ratio at N = 12
and 16, +3.3e-11 at N = 20, +4.1e-11 at N = 24 and -8.4e-12 at N = 32 (random, the Q.60 root and frequency floors), where version 4's
three extra roundings had left 1.7e-9 to 2.5e-9. So the ratio is good to twelve decimals at N = 12 and 16 and to ten at N = 20 to 32, and
(ratio - 1) N^2 to eight at N = 32; eleven decimals at N = 32 would need a 2^-72 frequency table and Q.62 roots (the reviewer's route, which
gave +1.1e-12 there). F2's printed twelve places are its own single rounding to Q.40 for printing, with the ratio taken from the unrounded
value. The two machines' bit-equality checks the arithmetic, not this bound.

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
published v2 table is superseded by the one above. v4 (ran 9:36 PM, chain 98ce8074ad037b40): the Q.60 frequencies, but three roundings
left one unit in the ninth decimal of the ratio at N = 24 and 32 (above). v5 (ran 10:20 PM, chain 708f274a7b9e6043): those three removed;
its table and chain are the ones on this page.
