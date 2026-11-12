# Result 2, second run: the SYK teleportation protocol exact for one realization, the full density matrix, N = 8 to 16, 2026-11-11

**First line:** this is the same finite-N arithmetic as the first run, with one correction from review: the paper's density-matrix
formulas are averaged over the random couplings, and for a single realization they are off by up to 14 percent at N = 8 and 4 percent at
N = 12. The exact per-realization matrix is computed here. The traversable-wormhole sign still wins over the wrong sign at every N, and
the signal still shrinks as N grows at these parameters. Nothing here builds anything; it is numbers others can check.

## What was computed

Gao and Jafferis (arXiv:1911.07416) write the reduced density matrix of the reference qubit R (half of an EPR pair, swapped into the
left Dirac mode chi_l = (psi_l^1 + i psi_l^2)/sqrt2 at time -t') and the right readout mode T (chi_r at time t, after the coupling
U = e^{i mu V}) in two forms: the six independent elements (2.15)-(2.20), and the two-parameter form (2.63)-(2.66) that follows from the
SO(N) symmetry of the ensemble average (rho22 = rho33 = 1 - rho11, rho23 = 0). The first run on this page's sister (`2_syk.md`) used the
two-parameter form per seed; an independent review (another AI session, with its own float implementation) showed that for one seed the
two forms differ at the percent level. This run computes the six elements.

- The elements, in the basis |a c> (a = R, c = T; index k = 2a + c), with M_ab = |a><b| on the left mode and N_cc' = |c><c'| on the right:
  rho_(ac),(a'c') = (1/2) SUM_b < M_ba'(-t') U^dag N_c'c(t) U M_ab(-t') > over the thermofield double; nonzero only for a + a' + c + c'
  even: rho00, rho11, rho22, rho33, rho03, rho12. Hermitian with trace 1 by construction. Two identities hold exactly and are printed
  as checks: rho00 + rho11 = rho22 + rho33 = 1/2 (the reference is never touched after the swap). The readout's marginal, rho00 + rho22,
  is not 1/2 for one seed when mu is not 0.
- I_RT = S(R) + S(T) - S(RT), with S(RT) from the two 2 x 2 blocks {0, 3} and {1, 2}; perfect teleportation is 2 ln 2.
- The layout: the left Majoranas on one half of the qubits and the right on the other, so H = H_L + H_R, the state is a D x D matrix
  with D = 2^(N/2), and e^{-iHt} Psi = U_L Psi U_R^T with U(1) a Taylor matrix exponential (128-bit accumulation, one rounding per
  entry) and U(t) = U(1)^t. The entangled state |I> in this layout has amplitude i^e(L) on |L>|~L> with e(L) = SUM over set bits k of L of
  (3 if k even, 1 if k odd) + 2 C(|L|, 2) mod 4, found numerically and proved by the review in a few lines; the program checks
  (psi_l^j + i psi_r^j)|I> = 0 exactly at every N it runs.
- The cost: the four states Z_ab = e^{-iHt} U e^{-iHt'} (2 M_ab) e^{iHt'} |TFD> serve all twelve correlators of the six elements
  (5t' + 4t evolution steps); 2M and 2N are sums of two Pauli strings (2 chi chi^dag = 1 - i M1 M2, 2 chi = M1 + i M2, and so on),
  applied exactly; each element is (1/16) SUM_b <Z_a'b | 2N_c'c | Z_ab>, accumulated in 128-bit and rounded once. The averaged forms are
  computed beside, on the same seed, by the first run's code (16t steps), so each line carries both and their difference.
- The couplings J ~ N(0, 6/N^3) from the same generator and seeds as the first run (the printed first coefficient is the same number),
  beta = 4, mu = -12, +12 and 0, t = t' = 1..6, three seeds at N = 8, 10, 12; N = 14 and 16 as the machines finish them tonight.
- Checks: the identities above (to 4e-11); mu = 0 gives I_RT = 0 to 5e-11; an independent double-precision implementation of the six
  elements on the same couplings agrees with every element to 1e-12 and with I_RT to 4e-12 at N = 8; the averaged-form columns equal the
  first run's lines digit for digit.
- Program: `code/wh_syk3.c` on `code/whq.h` (v2); binary md5 1c51c4a5 (ARM build).
- Runs on the pair, both machines at once, each identical line for line and equal to the PC build's chain: N = 8 at 9:45 PM (28 s,
  chain 9cb64d3bbcbbf2b3), N = 10 at 9:47 PM (81 s, chain 8cceeccb931b6c19), N = 12 at 09:55:32 PM (9 min, chain 3db1596b490e2767); N = 14, seed 1, at 10:16:52 PM (21 min, chain d4f859fe7106985d); N = 14 seeds 2-3 and N = 16
  follow and are added here when they land. The raw output of one machine: `results/2_syk_exact_board1.txt`.

## Error bound, on paper

Each D x D product rounds every entry once to Q.40, an error matrix with entries at most 2^-41, so of norm at most D 2^-41 (5.8e-11 at
D = 128); products of unitaries do not amplify it, so after about 100 products an expectation value moves at most 1.2e-8 in the worst
case and about 1e-10 if the roundings are independent (the review's bound, COMPUTED). The identities printed on every line (the reference's
marginal, the imaginary part of the diagonal) sit at 4e-11 or below at N = 8 to 12, and the double-precision check agrees to 1e-12 at
N = 8: of the twelve printed places the first ten are the protocol's. The two machines' bit-equality checks the arithmetic, not this bound.

## What the exact numbers say

- The averaged forms are a good approximation to one realization only on average: over 36 points at N = 8 their error in I_RT/(2 ln 2)
  runs from -0.019 to +0.011 (up to 14 percent of the value where the value is above 0.01), at N = 10 from -0.009 to +0.006, at N = 12
  from -0.004 to +0.004 (4 percent). The departures from the averaged structure are visible in the elements themselves: rho00 - rho33 up
  to 0.008 and |rho12| up to 0.07 at N = 8, both falling with N.
- The sign asymmetry is exact and grows relatively with N: the mean over three seeds of the peak exact I_RT/(2 ln 2) is 0.180 for the
  wormhole sign (mu = -12) against 0.125 for the wrong sign at N = 8; 0.159 against 0.101 at N = 10; 0.144 against 0.091 at N = 12.
  The peaks move later in t as N grows (t = 3-4 at N = 8, 4 at N = 10 for the wormhole sign).
- mu = 0 is the exact zero at every N and seed.

## The numbers (each table from the machines' own lines; the tables are regenerated by `syk3_table.py` from any log of the run)

### N = 8: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 3 | 0.146969053887 | 0.146132801740 | +0.000836252147 | -0.002134192 | 0.006268142 |
| 1 | +12 | 6 | 0.115218913771 | 0.116014285693 | -0.000795371922 | +0.008315894 | 0.001402536 |
| 2 | -12 | 4 | 0.217788340764 | 0.212457505041 | +0.005330835723 | +0.005129828 | 0.070712431 |
| 2 | +12 | 4 | 0.156039246864 | 0.156817382270 | -0.000778135406 | +0.000428150 | 0.044984349 |
| 3 | -12 | 3 | 0.174185823175 | 0.170045474534 | +0.004140348641 | -0.000669403 | 0.007264842 |
| 3 | +12 | 6 | 0.103983083720 | 0.096001428623 | +0.007981655097 | +0.001046484 | 0.004581906 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.179648; the averaged forms at the same points: 0.176212
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.125080; the averaged forms at the same points: 0.122944
exact - averaged over all 36 points: from -0.018696 (seed 2, mu -12, t 6) to +0.010516 (seed 2, mu -12, t 3); the largest relative departure where I > 0.01: 13.9%
mu = 0: the exact I_RT/(2 ln 2) is at most 3.4e-11 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 2.7e-11 over every point

### N = 10: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.158524590199 | 0.158534205653 | -0.000009615454 | +0.002549185 | 0.005754578 |
| 1 | +12 | 6 | 0.099490708534 | 0.100576400542 | -0.001085692008 | -0.000814569 | 0.011295968 |
| 2 | -12 | 4 | 0.171120669214 | 0.165888550538 | +0.005232118676 | -0.000640565 | 0.029056145 |
| 2 | +12 | 4 | 0.103499517612 | 0.103510044533 | -0.000010526921 | -0.001274533 | 0.019436651 |
| 3 | -12 | 4 | 0.147321127403 | 0.151541877187 | -0.004220749784 | +0.003942014 | 0.015266956 |
| 3 | +12 | 6 | 0.101040449401 | 0.109903307432 | -0.008862858031 | -0.002027956 | 0.009449319 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.158989; the averaged forms at the same points: 0.158655
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.101344; the averaged forms at the same points: 0.104663
exact - averaged over all 36 points: from -0.008863 (seed 3, mu +12, t 6) to +0.005924 (seed 2, mu -12, t 5); the largest relative departure where I > 0.01: 8.8%
mu = 0: the exact I_RT/(2 ln 2) is at most 4.0e-11 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 4.1e-11 over every point

### N = 12: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.146953107241 | 0.148357534011 | -0.001404426770 | -0.000728107 | 0.006622384 |
| 1 | +12 | 5 | 0.092173545097 | 0.095484204402 | -0.003310659305 | -0.002310609 | 0.005744434 |
| 2 | -12 | 4 | 0.140981138551 | 0.139841139117 | +0.001139999434 | +0.000600403 | 0.003179067 |
| 2 | +12 | 6 | 0.092471120442 | 0.092211175660 | +0.000259944782 | -0.000978394 | 0.001723603 |
| 3 | -12 | 4 | 0.143299680336 | 0.142190589434 | +0.001109090902 | -0.002077237 | 0.007944345 |
| 3 | +12 | 6 | 0.088439394066 | 0.084911137320 | +0.003528256746 | -0.000094734 | 0.004124797 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.143745; the averaged forms at the same points: 0.143463
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.091028; the averaged forms at the same points: 0.090869
exact - averaged over all 36 points: from -0.003774 (seed 1, mu +12, t 6) to +0.003528 (seed 3, mu +12, t 6); the largest relative departure where I > 0.01: 4.1%
mu = 0: the exact I_RT/(2 ln 2) is at most 5.0e-11 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 3.6e-11 over every point

### N = 14: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.141559308799 | 0.142215050231 | -0.000655741432 | -0.000839398 | 0.010612170 |
| 1 | +12 | 6 | 0.086986690628 | 0.087356583623 | -0.000369892995 | +0.000561942 | 0.005864685 |
mean over 1 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.141559; the averaged forms at the same points: 0.142215
mean over 1 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.086987; the averaged forms at the same points: 0.087357
exact - averaged over all 12 points: from -0.001732 (seed 1, mu -12, t 6) to +0.000062 (seed 1, mu +12, t 2); the largest relative departure where I > 0.01: 1.6%
mu = 0: the exact I_RT/(2 ln 2) is at most 6.2e-11 over 6 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 5.8e-11 over every point

## Version 4, 2026-11-11 night: the same six elements in Q.60 end to end, and six times faster

The reviewer's bit-exact emulation of version 3's arithmetic (10:20 PM) showed that its twelve printed decimals were good to ten: the Q.40
chain of products U(1)^t left up to 21 units of 2^-40 in the evolved states at t = 6, and the Q.40 cosine and sine of the coupling (whose
squares summed to 1 minus 2.4 units) moved K by up to 27 units. Version 4 (`code/wh_syk4.c`) follows the reviewer's route and rounds nothing
to Q.40 at all:

- every matrix and state carries 60 fraction bits (int64 entries, 128-bit products, one rounding per entry per product);
- one set of Hamiltonian powers per seed gives both the time step W(1) = e^{-iH} (summed with the exact phases (-i)^k) and the thermal
  factor e^{-beta H/4} (exact for beta = 4), about 18-19 terms to below one unit of 2^-60; W(t) = W(1) W(t-1) once per seed, so an
  evolution by any t is one pair of products, and e^{+iHt} is the conjugate transpose of W(t) with no new rounding;
- the coupling's cosine and sine, the eigenvalues' square root and the entropies' logarithm are all Q.60 (ln 2 is computed the same way and
  printed as a self-test); the entangled state |I> is built normalized exactly, with a global phase (1 + i)/sqrt2 when N/2 is odd, so that
  no sum of squares can approach the 128-bit limit (the first build had kept |I> unnormalized, whose squared norm is exactly 2^127 at
  D = 128: caught at N = 14 before any run);
- products skip exact zeros, which is the fermion-parity block structure of H, W and e^{-beta H/4} without bookkeeping.

**Checks (MEASURED, PC build):** at N = 8 and N = 10 every one of the twelve printed decimals of every element and of I_RT equals the
independent double-precision implementation of the same elements on the same couplings (version 3 had differed in the last one or two);
the identities rho00 + rho11 = rho22 + rho33 = 1/2 hold to 0-2 units of 2^-60 at every point; c^2 + s^2 - 1 is 0 to -2 units; the
unitarity defect of W(1), max |W^dag W - 1|, is 5-7 units; e^{-beta H_L/4} Psi_I = Psi_I (e^{-beta H_R/4})^T to 1-3 units; the thermal
norm Z_beta/Z_inf = 1.419913512654 at N = 8, seed 1 agrees with two double-precision routes (1.419913512653848 by both); mu = 0 prints exact
zeros. The 18-point N = 14 set takes 8 s on one PC core (version 3: 49 s, with the averaged forms beside), N = 16 about a minute.

**What the extra digits are for:** nothing in the physics changes at the tenth decimal; they let the two machines' bit-equality and the
comparison with any independent integration be stated to the full printed width, and they make N = 16 a half-hour run on the pair. The
pair's version-4 runs at N = 8, 10, 12, 14 (three seeds) and 16 (seed 1) are below, each both machines line for line.

**Toward N = 18 (a development note, not a result).** Version 4.1 derives the right-hand operators from the left-hand ones exactly (the
right Hamiltonian is a phased re-indexing of the left, and with the left-hand series multiplied as P H and the chain as W(t-1) W(1) the
identity holds bit for bit, so no right-hand series or storage is needed), and version 5 (`code/wh_syk5.c`) lays the same arithmetic out
in a pool of seventeen matrices that fits N = 18 (D = 512, 68 MB) inside the machines' memory; it reproduces version 4.1 bit for bit at
N = 8, 10 and 12. A development run of it on one PC core at N = 18, seed 1 (601 s; the identities at 0 units, the step operator's unitarity
defect 7 units of 2^-60, mu = 0 exact zeros) gives a peak exact I_RT/(2 ln 2) of 0.1233 at t = 5 for the wormhole sign against 0.0824 at
t = 6 for the wrong one; with the pair's N = 16 value 0.1315 (t = 4, below), the peak keeps falling slowly with N (0.180, 0.159, 0.144, 0.142 on the
pair at N = 8 to 14 for seed 1) and the sign asymmetry stays at 1.4 to 1.6. These two numbers are the PC build's and are labelled so until
the pair runs them (N = 16 ran on the pair at 12:08 AM, below, and printed the PC's digits; N = 18 another night, about four hours a seed on one of the machines' cores at this precision).

**Version 7, the same bits in a quarter of the instructions (11-12, 1:17 AM).** The reviewer's session wrote and proved a three-limb
kernel for the Q.60 products (`code/limb_q60.c`): each 64-bit operand split into 24-bit limbs, the complex product by three real ones, the
128-bit sums rebuilt exactly from int64 accumulators, one rounding per entry by the same rule, so the result is bit-equal to the 128-bit
kernel for every input (22,059 random cases, and the whole program's output at N = 8 to 16). Version 7 (`code/wh_syk7.c`) is version 6 with
every matrix product through it: on the PC the printed output is byte-identical to version 6 at N = 8, 10, 12 and 14 (three seeds, every
check on) and the N = 16 chains are equal; the A9 instruction count per product falls four- to fivefold (the reviewer's qemu counts), so
N = 18 is about 50 minutes a seed on one of the machines' cores instead of four hours. The machines run it tonight, first at N = 8 to 14
against version 6's chains, then at N = 18 for three seeds.

## Version 4 on the pair (2026-11-11 night; each step both machines line for line, each chain equal to the PC build's)

Runs (the time each landed, both machines line for line, every chain equal to the PC build's): N = 8 at 11:59:50 PM (24 s, chain
f5b7631e8db77784); N = 10 at 12:00:19 AM (29 s, c6fd6b0164a9e5c8); N = 12 at 12:01:31 AM (72 s, e947d9359ae0dd57; version 3 had taken
nine minutes for the same points); N = 14, three seeds, at 12:08:12 AM (6.7 min, 67c068cb529dac63; version 3 had needed about half an hour
for one seed); N = 16, seed 1, at 12:24:34 AM (16.4 min, c6032abd10ef74d6); the 2022 experiment's model at t = 2.8 at 12:25:02 AM (28 s,
a9379694959a987b). The raw output of one machine: `results/2_syk_exact_v4_board1.txt`. The numbers below are the exact matrix's; the
"averaged forms" columns are empty for version 4, which no longer computes them (version 3's tables above have them).

### N = 8: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 3 | 0.146969053882 | n/a | n/a | -0.002134192 | 0.006268142 |
| 1 | +12 | 6 | 0.115218913751 | n/a | n/a | +0.008315894 | 0.001402536 |
| 2 | -12 | 4 | 0.217788340772 | n/a | n/a | +0.005129828 | 0.070712431 |
| 2 | +12 | 4 | 0.156039246871 | n/a | n/a | +0.000428150 | 0.044984349 |
| 3 | -12 | 3 | 0.174185823166 | n/a | n/a | -0.000669403 | 0.007264842 |
| 3 | +12 | 6 | 0.103983083688 | n/a | n/a | +0.001046484 | 0.004581906 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.179648
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.125080
mu = 0: the exact I_RT/(2 ln 2) is at most 0.0e+00 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 0.0e+00 over every point

### N = 10: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.158524590182 | n/a | n/a | +0.002549185 | 0.005754578 |
| 1 | +12 | 6 | 0.099490708516 | n/a | n/a | -0.000814569 | 0.011295968 |
| 2 | -12 | 4 | 0.171120669180 | n/a | n/a | -0.000640565 | 0.029056145 |
| 2 | +12 | 4 | 0.103499517574 | n/a | n/a | -0.001274533 | 0.019436651 |
| 3 | -12 | 4 | 0.147321127373 | n/a | n/a | +0.003942014 | 0.015266956 |
| 3 | +12 | 6 | 0.101040449351 | n/a | n/a | -0.002027956 | 0.009449319 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.158989
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.101344
mu = 0: the exact I_RT/(2 ln 2) is at most 0.0e+00 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 0.0e+00 over every point

### N = 12: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.146953107224 | n/a | n/a | -0.000728107 | 0.006622384 |
| 1 | +12 | 5 | 0.092173545074 | n/a | n/a | -0.002310609 | 0.005744434 |
| 2 | -12 | 4 | 0.140981138537 | n/a | n/a | +0.000600403 | 0.003179067 |
| 2 | +12 | 6 | 0.092471120418 | n/a | n/a | -0.000978394 | 0.001723603 |
| 3 | -12 | 4 | 0.143299680320 | n/a | n/a | -0.002077237 | 0.007944345 |
| 3 | +12 | 6 | 0.088439394037 | n/a | n/a | -0.000094734 | 0.004124797 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.143745
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.091028
mu = 0: the exact I_RT/(2 ln 2) is at most 0.0e+00 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 0.0e+00 over every point

### N = 14: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.141559308749 | n/a | n/a | -0.000839398 | 0.010612170 |
| 1 | +12 | 6 | 0.086986690547 | n/a | n/a | +0.000561942 | 0.005864685 |
| 2 | -12 | 4 | 0.136228568414 | n/a | n/a | +0.000295675 | 0.008760569 |
| 2 | +12 | 6 | 0.088662992114 | n/a | n/a | +0.002470893 | 0.003045576 |
| 3 | -12 | 4 | 0.135315758588 | n/a | n/a | +0.000341361 | 0.005906064 |
| 3 | +12 | 6 | 0.087609856281 | n/a | n/a | +0.000773986 | 0.000358937 |
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.137701
mean over 3 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.087753
mu = 0: the exact I_RT/(2 ln 2) is at most 0.0e+00 over 18 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 0.0e+00 over every point

### N = 16: the peak over t of the exact I_RT/(2 ln 2), per seed and sign of mu

| seed | mu | t at the peak | exact I_RT/(2 ln 2) | averaged forms (2.63)-(2.66) | exact - averaged | rho00 - rho33 | abs(rho12) |
|---|---|---|---|---|---|---|---|
| 1 | -12 | 4 | 0.131525302861 | n/a | n/a | +0.000030040 | 0.007056335 |
| 1 | +12 | 6 | 0.085268129072 | n/a | n/a | +0.000597384 | 0.002465491 |
mean over 1 seeds of the peak exact I_RT/(2 ln 2) at mu = -12: 0.131525
mean over 1 seeds of the peak exact I_RT/(2 ln 2) at mu = +12: 0.085268
mu = 0: the exact I_RT/(2 ln 2) is at most 0.0e+00 over 6 points (the exact zero to the rounding)
the reference's marginal: |rho00 + rho11 - 1/2| and |rho22 + rho33 - 1/2| at most 0.0e+00 over every point

### The 2022 experiment's model at t = 2.8 (version 4 with `model=m1 dt=2.8`, both machines line for line, chain a9379694959a987b)

The learned Hamiltonian of the 2022 experiment, as defined on `2_syk.md`, at beta = 4 and the experiment's t = t' = 2.8, now in Q.60 end to
end with the exact per-realization matrix (the model has no random couplings, so this is the one exact answer for it). The exact
I_RT/(2 ln 2): 0.120532639717 at mu = -12 and 0.050104226326 at mu = +12; version 1's averaged forms had given 0.123572 and 0.050165 for the
same points, so the averaged formulas are off by 2.5 percent for the right sign here. The six elements at mu = -12: rho00 = rho33 =
0.330084437435, rho11 = rho22 = 0.169915562565, rho03 = -0.185893876352, rho12 = -0.025344460252, every imaginary part exactly zero;
K = -0.675272946464. At mu = +12: rho00 = rho33 = 0.243452031122, rho11 = rho22 = 0.256547968878, rho03 = 0.124934526964,
rho12 = 0.021723232433; K = 0.529263917433. The reference's marginal holds to 23 units of 2^-60 and the step operator's unitarity to the
usual few units; the time unit dt = 2.8 is parsed through the Q.40 input format, so the point is t = 2.8 + 1.8e-13 (a review found this;
it moves K by about 3e-14, below the printed width, and version 6 parses dt straight into Q.60).
