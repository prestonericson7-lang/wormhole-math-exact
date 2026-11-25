# Result 2, first run: the SYK teleportation protocol's averaged formulas in exact arithmetic at N = 8, and the 2022 experiment's model, 2026-11-11

**In plain words.** A paper by Gao and Jafferis (arXiv:1911.07416) described a recipe called wormhole teleportation: two entangled
(quantum-linked) copies of the SYK model (a model of many interacting particles that physicists use as a stand-in for a wormhole) are
prepared, a quantum bit (qubit) is put into one copy, a brief coupling, which can be given a plus or a minus sign, joins the copies,
and part of the qubit's information can then be read out of the other copy, which in the matching gravity picture (the same system
described as a wormhole) means it went through the wormhole. The paper's formulas for this are averages over the model's random
interaction strengths, and its best case, 0.540 of a perfect teleportation (scored by how much information the read-out qubit shares
with a reference qubit kept back, entangled with the one sent in), is for very many particles at low temperature; in 2022 a team ran a
small version on nine qubits of a quantum computer, with a machine-learned model in place of SYK, and whether that experiment shows
wormhole physics is disputed. This first run evaluates the paper's formulas exactly, in fixed-point integer arithmetic (every
calculation done with whole numbers in a set format, so any computer prints the same digits), for three random draws of the full SYK
model (seeds 1, 2 and 3) with 8 particles a side (N = 8 Majoranas in the text below) and for the 2022 model; two small ARM boards
(called the machines, or the pair, below) printed identical lines, with the same chain hash (a fingerprint of the run's results, taken
at full precision, so two runs can be compared exactly). The arithmetic is exact, but for any one draw the averaged formulas themselves
are only an approximation; the exact per-draw values are on `2_syk_exact.md`. Here the sign that the paper says teleports does better
than the other sign, 0.172 against 0.116 of a perfect teleportation at the peak over time of the three-draw average, and the 2022 model
favours that sign more strongly (2.5 times the other sign's score, against 1.5 times for SYK; all COMPUTED from the machines' lines),
so that lopsidedness alone cannot tell the learned model from the full SYK model at this size. Nothing here builds anything.

**First line:** this is arithmetic about a quantum protocol, not a wormhole; the protocol is a map between a toy quantum system's
teleportation and a gravitational picture, and the 2022 experiment that ran it on nine qubits is disputed on exactly that point. This
result evaluates the paper's density-matrix formulas exactly (in integer arithmetic) for three realizations of the full SYK model at
N = 8 Majoranas a side, and for the experiment's learned five-term model at its own parameters; the exact per-realization matrix, which
those formulas approximate, is on `2_syk_exact.md`. It is the first run of a program built to go to N = 14 and beyond.

## What was computed

Gao and Jafferis, "A traversable wormhole teleportation protocol in the SYK model" (arXiv:1911.07416), sections 2.1-2.3. Two SYK systems
l and r, N Majoranas each ({psi_i, psi_j} = delta_ij), each with H = -SUM_{i<j<k<l} J_{ijkl} psi_i psi_j psi_k psi_l and Gaussian couplings
of variance 3!/N^3 (J = 1), the same couplings on both sides (their eq. 2.37 for q = 4). The state |I> is annihilated by every
psi_l^j + i psi_r^j (2.35) and the thermofield double is e^{-beta H/4}|I> normalized (2.36). The protocol inserts a reference-entangled qubit
into the Dirac fermion chi_l^1 = (psi_l^1 + i psi_l^2)/sqrt2 at time -t', applies the two-sided coupling U = e^{i mu V} with
V = (i/qN) SUM_j psi_l^j psi_r^j (2.50) at time 0, and reads chi_r^1 at +t. The reduced state of the reference R and the read qubit T is
their (2.13) with (2.63)-(2.66):

    rho_TR = (1/2) [[rho11, 0, 0, rho14], [0, 1-rho11, 0, 0], [0, 0, 1-rho11, 0], [rho14*, 0, 0, rho11]]
    rho11 = (1 - <{A1, [A2, B12]}>)/2,   rho14 = <{A1, B1}>/2 - <{A2, A1 B2 A1}>,

with A_i = psi_l^i(-t'), B_i = U^dag psi_r^i(t) U, B12 = B1 B2, all in the thermofield double. The mutual information
I_RT = S(R) + S(T) - S(RT) = 2 ln 2 + SUM lambda ln lambda over the four eigenvalues (rho11 +- |rho14|)/2 and (1 - rho11)/2 twice; perfect
teleportation is I_RT = 2 ln 2. K = <{A1, B1}> is their large-N variable, in which rho11 = (1 + K^2)/2, rho14 = K and their eq. (4.9)
gives I_RT; their Table 1 puts the large-N, low-temperature (J -> infinity) maximum of I_RT/(2 ln 2) at 0.540 for q = 4 (their eqs. 4.19
and 4.9). The teleporting sign of the coupling in
their convention is mu < 0.

The 2022 experiment (Jafferis et al., Nature 612, 51) could not run SYK on hardware and ran a machine-learned five-term Hamiltonian on
7 Majoranas a side instead; Kobrin, Schuster and Yao (arXiv:2302.07897) give it as their eq. (1), "Model 1", and argue it does not carry
the gravitational physics. The same protocol is computed here on Model 1 exactly as they state it: H = -0.36 psi1psi2psi4psi5
+ 0.19 psi1psi3psi4psi7 - 0.71 psi1psi3psi5psi6 + 0.22 psi2psi3psi4psi6 + 0.49 psi2psi3psi5psi7 (its coefficients declared in Q.40),
beta = 4, t = t' = 2.8, mu = -12 and +12 with V normalized as if N were 10 (their note on the experiment's convention).

- Program: `code/wh_syk.c` on `code/whq.h` (v2); binary md5 c1c5ada6 (ARM build, the binary that ran; later input checks change no
  printed line).
- Arithmetic: complex Q.40 amplitudes (2^-40 = 9.1e-13), 2^N of them; the Majoranas are the paper's Pauli strings (their eq. 2.1-2.2),
  applied as exact permutations with phases; H v is accumulated in signed 128-bit over every term and rounded once; e^{sH} by Taylor in
  steps of 1 (the norm of H_l + H_r is 0.79 to 1.06 at N = 8 over the three seeds, COMPUTED), terms rounded to nearest, the sum in
  128-bit, stopped below one unit; U applied as the exact product
  over j of cos(mu/2qN) + i sin(mu/2qN) W_j with W_j = 2i psi_l^j psi_r^j (W_j^2 = 1); the couplings from xorshift64* with declared seeds
  1, 2, 3 (twelve 40-bit uniforms per Gaussian). Self-tests printed and exact: (psi_l^j + i psi_r^j)|I> = 0 for every j;
  {bare_a, bare_b} = 2 delta on |I>; (H_l - H_r)|I> = 0; and mu = 0 gives rho11 = 1/2, rho14 = 0, I_RT = 0 to 1e-11 at all 18 points (the
  two sides' operators anticommute, so the protocol without the coupling teleports nothing: an exact zero of the algebra, reached here
  through 18 time evolutions per point).
- Cross-check: a double-precision implementation with dense 2^8 x 2^8 matrices and exact matrix exponentials, on the same couplings bit
  for bit, agrees with the integer program at every printed digit but the last one or two (1e-11) at N = 8, seed 1, at the points it ran.
  A comparison of every point with the reviewer's __float128 reference on the same couplings (11-11, 10:20 PM; COMPUTED) puts the error
  at up to 4.8e-11 in K, from version 1's Q.40 arithmetic: the rounding of every Taylor step on the 2^N amplitudes, and the Q.40 cosine
  and sine of the coupling (their squares sum to 1 minus 2.4 units); the factorized program carries both in Q.60 from version 4 on
  (`2_syk_exact.md`).
- Runs: 7:40-7:46 PM (the full SYK, 54 points) and 7:46 PM (Model 1, 2 points), both machines at once; 64 and 8 output lines each;
  **identical line for line**, chains 2c70424b9fa1aca3 and 9899c884ff17ebb5. The raw output: `2_syk_board1.txt` (board 2 byte-identical).

**A correction from review, 11-11 evening, before the numbers.** The formulas (2.63)-(2.66) above are the paper's ensemble-averaged
forms: the SO(N) symmetry of the disorder average sets rho22 = rho33 = 1 - rho11 and rho23 = 0, and the paper derives them for the
average over couplings. For one realization of the couplings, which is what every line on this page is, the exact rho_TR has six
independent elements, the paper's (2.15)-(2.20): the four diagonal ones, rho14 and rho23. An independent review (another AI session,
with its own float implementation) finds that I_RT from the averaged forms differs from the exact per-seed value by -12.9% to +5.2% over
four seeds of its own at N = 8 (COMPUTED by the reviewer). Every number below is therefore "(2.63)-(2.66) evaluated per seed", exact as
arithmetic and labelled as the quantity it is; the sign comparison (mu < 0 against mu > 0) and the identity checks stand as stated.
Version 3 of the program computes the exact six-element matrix on the same seeds, side by side with these forms, and its results land
on their own page with the same labels.

## The numbers

I_RT / (2 ln 2) for the full q = 4 SYK, N = 8, beta = 4, t = t' (the mean of the three seeds, then each seed; MEASURED,
`2_syk_board1.txt`, the means COMPUTED from those lines):

| t | mu = -12 (mean) | mu = +12 (mean) | mu = -12: seeds 1, 2, 3 | mu = +12: seeds 1, 2, 3 |
|---|---|---|---|---|
| 1 | 0.1093 | 0.1048 | 0.1071, 0.1118, 0.1089 | 0.1050, 0.1087, 0.1008 |
| 2 | 0.1396 | 0.1119 | 0.1280, 0.1497, 0.1413 | 0.1133, 0.1287, 0.0936 |
| 3 | 0.1702 | 0.1142 | 0.1461, 0.1944, 0.1700 | 0.1120, 0.1469, 0.0838 |
| 4 | 0.1721 | 0.1161 | 0.1449, 0.2125, 0.1589 | 0.1054, 0.1568, 0.0861 |
| 5 | 0.1471 | 0.1171 | 0.1267, 0.1910, 0.1237 | 0.1064, 0.1522, 0.0926 |
| 6 | 0.1213 | 0.1150 | 0.1104, 0.1528, 0.1007 | 0.1160, 0.1330, 0.0960 |

mu = 0: 0.000000000000 at every t and seed (rho11 = 0.5 to 3e-12, |rho14| below 1e-11). The full twelve-digit lines, with rho11, rho14 and
K for every point, are in `2_syk_board1.txt`.

At the mean peak (t = 4) the teleporting sign gives 0.172 against 0.116 for the other sign: a ratio of 1.5 at N = 8 (COMPUTED from the
table). Seed by seed the peak ratio runs from 1.26 to 1.69, and seed 3's wrong-sign curve falls from t = 1 to t = 3 and recovers only to
0.096 by t = 6 (its peak is at t = 1), while its right-sign curve peaks at t = 3. The large-N, J -> infinity maximum is 0.540; N = 8
reaches 0.21 at best (seed 2, t = 4).

**What the finite-N numbers say that the large-N formula does not.** At seed 2's peak, K = -0.6345 but rho14 = -0.5256 and
rho11 = 0.6706. The large-N relations would put rho14 = K and rho11 = (1 + K^2)/2 = 0.701, and their eq. (4.9) applied to this K gives
I_RT/(2 ln 2) = 0.314 (COMPUTED); the averaged forms on this seed give 0.2125, and the exact per-realization matrix gives 0.2178
(`2_syk_exact.md`). At N = 8 the second term of rho14, <{A2, A1 B2 A1}>, is far from its large-N limit and costs about a third (31 percent
against the exact value) of the mutual information the large-N formula would promise for the same K. Any finite-N reading of the
protocol's "signal" through K alone overstates it.

Model 1 (the 2022 experiment's learned Hamiltonian), beta = 4, t = t' = 2.8, V normalized with N = 10 (MEASURED, `2_syk_board1.txt`:70-71):

| mu | rho11 | rho14 | K | I_RT / (2 ln 2) |
|---|---|---|---|---|
| -12 | 0.660168874864 | -0.384084462984 | -0.675272946432 | 0.123572174034 |
| +12 | 0.486904062238 | +0.253236405331 | +0.529263917417 | 0.050164564023 |

Model 1 shows a stronger asymmetry between the two signs (2.5) than the full SYK at N = 8 does on average (1.5), at a mutual information
in the same range (0.124 against 0.17). The asymmetry alone therefore does not separate the learned model from the real one at this size;
the next test, whether the full SYK's size distribution winds and how it compares with the learned model's, is on `2_syk_winding.md`.

## Error bound, on paper

Every evolution rounds each amplitude once per Taylor term (about 14 terms per unit of time) and once per step; with amplitudes of
order 2^-4 at N = 8 and 18 evolutions of up to 6 steps per point, the worst-case accumulated rounding on an inner product is below
2e-9 and the random-walk estimate below 1e-11; the printed norm drift of every real-time evolution stayed below 1e-10. The mu = 0 points,
which must be exactly zero, came out at 1e-11, and the double-precision cross-check agrees to 1e-11 at the points it ran, while the
reviewer's __float128 reference finds up to 4.8e-11 in K (above): the digits down to the 10th are the
protocol's, the last two belong to the declared arithmetic. The two machines' bit-equality checks the arithmetic, not this bound. The
disorder average over three seeds is a three-sample mean; its spread is visible in the table and is the real uncertainty of the "mean"
columns.

## What came next

Version 2 moved the program to a factorized layout; versions 3 to 7 (`2_syk_exact.md`) compute the exact per-realization density matrix
in it, on the pair to N = 16; version 6 computes the size-winding signature (`2_syk_winding.md`).

## The two layouts against each other at N = 10 (added 11-11, 11:56 PM)

The paper's interleaved layout (this program, version 1) and the factorized layout of versions 2 and later compute the same quantities
by different arithmetic paths. At N = 10, three seeds, mu = -12, +12 and 0, t = t' = 1..6 (54 points, 270 printed numbers), version 1's PC
build (chain fdb4636ff5ab97ee) and version 3's run on the pair (chain 8cceeccb931b6c19, whose averaged-form columns are version 2's
arithmetic) agree to 7.4e-11 in the worst number, K at seed 1,
mu = +12, t = 6, which is version 1's own Q.40 rounding on 2^N amplitudes as the review measured it. **On the pair (MEASURED, 11-12):**
version 1 at N = 10, the same 54 points, ran on the two machines from 12:33 to 1:52 AM (79 minutes; version 3 took 81 seconds for the
same points and version 4 29 seconds), both machines line for line, chain fdb4636ff5ab97ee, the PC build's exactly. So the comparison above now
stands between two runs on the machines: the paper's layout and ours agree to 7.4e-11 across the 270 numbers, and the chains are the
proof of which bits each produced. The raw lines are appended to `results/2_syk_board1.txt`.
