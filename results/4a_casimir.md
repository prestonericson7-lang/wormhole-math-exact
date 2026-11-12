# Result 4a: the negative energy itself, the Casimir sum on a lattice, exact, 2026-11-11

**First line:** the Casimir effect is real, measured since 1997, and its negative energy density is of the order of nanojoules per
cubic metre at micron spacings; nothing here makes it larger. This result computes the textbook case for one scalar field as a finite sum
with no divergent subtraction at all, and finds the exact law by which the lattice sum approaches the continuum: 1 + 2.1/N^2.

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
finite-size term, which goes as e^{-pi M/N}, is below the arithmetic.

- Program: `code/wh_casimir.c` on `code/whq.h` (v2); binary md5 5e191971 (ARM build).
- Arithmetic: Q.40, 128-bit intermediates, round to nearest; each omega is the exact-to-the-floor square root of a sum of three
  2(1 - cos k) terms with cos from the shared series; the sums run in 128-bit accumulators, the second difference is taken on the raw
  128-bit sums before any shift, and the whole-box F2 is rounded once; the ratio divides whole-box numbers (full precision).
- Grid: M = 128 for N = 4..20 and M = 256 for N = 4..32 (M >= 6 N throughout).
- Run: 7:37-7:39 PM, both machines at once; 20 output lines each; **identical line for line**, chain 74d672bf99243ecf. The raw output:
  `4a_casimir_board1.txt` (board 2's lines are byte-identical). A native PC build of the same source gives the same 20 lines byte for
  byte: the arithmetic does not depend on the machine.

## The numbers

M = 256 (A = 65536):

| N | E(N)/M^2 | F2(N), the whole box | F2(N)/M^2 | continuum, per unit area | ratio lattice/continuum | (ratio - 1) N^2 |
|---|---|---|---|---|---|---|
| 4 | 3.594731192354 | -7.363033362708 | -0.000112350972 | -0.000094494862 | 1.188963823640 | 3.02 |
| 6 | 5.982414094820 | -0.792898382356 | -0.000012098668 | -0.000011351365 | 1.065833842791 | 2.37 |
| 8 | 8.370035259392 | -0.177039442934 | -0.000002701407 | -0.000002610949 | 1.034645760296 | 2.22 |
| 12 | 13.145249318011 | -0.022367810508 | -0.000000341306 | -0.000000336344 | 1.014749505424 | 2.12 |
| 16 | 17.920456114245 | -0.005233460006 | -0.000000079856 | -0.000000079207 | 1.008185658075 | 2.10 |
| 20 | 22.695661420177 | -0.001703843066 | -0.000000025998 | -0.000000025863 | 1.005230620276 | 2.09 |
| 24 | 27.470866267894 | -0.000682324644 | -0.000000010411 | -0.000000010373 | 1.003639131559 | 2.10 |
| 32 | 37.021275526505 | -0.000161473702 | -0.000000002463 | -0.000000002457 | 1.002619982647 | 2.68 |

M = 128 gives the same ratios to 6e-7 for N <= 20 (its lines are in the raw output), which measures the box's finite-size term: below
1e-6 at M = 6.4 N, as e^{-pi M/N} says.

**The law.** The ratio approaches 1 from above as 1 + c/N^2 with c settling at 2.09-2.12 from N = 12 to 24 (the N = 32 entry, 2.68, is
where the rounding floor of the five-digit F2 begins to show: see the bound). The lattice's dispersion omega = 2 sin(k/2) softens every
mode relative to omega = k, so the lattice sum carries slightly more negative energy than the continuum at any finite N, by about 2.1/N^2.
The bulk energy per site, the slope of E(N)/M^2 in N, is 1.19384 from the table (37.021276 - 27.470866 over 8), and the surface term
2 e_surf is what remains after it: both cancel exactly in F2, which is why no subtraction is needed.

## Error bound, on paper

sqrtq is a floor, so every term carries a bias of about -1/2 unit of 2^-40; the bias is the same per term for every N and cancels in the
second difference. What remains is the random part: 0.29 units times the square root of the number of terms. At N = 32, M = 256 the three
sums have 31, 32 and 33 times 65536 terms, so F2 carries about 0.29 x sqrt(6 x 65536 x 32) = 1000 units = 9e-10 against F2 =
1.6e-4: six digits at N = 32, more at smaller N; the (ratio - 1) N^2 column at N = 32 moves by 0.3 for a 3e-7 change in the ratio, which
is why its last entry is the first to wander. The two machines' bit-equality checks the arithmetic, not this bound.

## In SI, for the record (DERIVED, not computed here)

Electromagnetic: E_C / A = -pi^2 hbar c / (720 d^3) with hbar c = 3.16153e-26 J m: at d = 1 um, -4.33e-10 J/m^2 (a pressure of 1.3 mPa);
at d = 10 nm, -4.33e-4 J/m^2 (1.3 kPa). The one scalar computed here gives half of each. The energy density between the plates,
-pi^2 hbar c / (720 d^4) electromagnetic: -4.3e-4 J/m^3 at 1 um. A traversable wormhole of a metre needs negative energy of planetary-mass
scale (computation 1's own papers say so); this is the gap, stated in numbers.
