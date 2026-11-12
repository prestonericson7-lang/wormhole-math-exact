# 2, part 3: the size-winding signature of the SYK teleportation protocol, exact

**What this page is.** Kobrin, Schuster and Yao (2023), commenting on the 2022 "traversable wormhole on a quantum processor" experiment,
defined the signature that separates teleportation through a wormhole from ordinary scrambling: *size winding*. Expand the thermal
operator rho^{1/2} psi_i(t) in Majorana strings, sum the squared coefficients by the string's size l, and the resulting complex numbers
q(l) must have phases that climb linearly with l, at exactly the rate the coupling e^{i mu V} unwinds. The experiment's learned
Hamiltonian showed the asymmetry between the two signs of mu; the comment asked whether its size distribution actually winds. This page
computes q(l) and the size distribution P(l) for the full SYK model, per realization, in exact integer arithmetic on the two machines,
together with the quantity that ties the winding to the protocol, and reports KSY's own fit numbers from them.

Version 6 of the SYK program (`code/wh_syk6.c`) is version 5 (`2_syk_exact.md`) plus these lines. Its protocol lines and their chain are
version 5's bit for bit; the new lines carry their own chain. Labels as everywhere here: MEASURED (printed by a program, run named), COMPUTED
(derived, inputs named), ESTIMATE.

## Where the formulas come from

The reviewer (the independent AI session that checks every program here) derived the observable for this program's layout and verified the
derivation by having a second, adversarial build redo it from nothing: its own port of the coupling, the entangled state found numerically,
a direct Majorana-subset transform with no Walsh-Hadamard step and no size table, and KSY's own operator. Twelve small corrections came out
of that pass; none changed a conclusion. The route, as implemented:

- **The operator state.** The program holds A = e^{-iHt} g_i e^{iHt} T, the Majorana string g_i (Heisenberg-evolved to time -t) applied to
  the normalized thermofield double T; it is sqrt2 times the matrix of psi_l^i(-t)|tfd>, with Frobenius norm 1. Version 5 already computes
  it for the first operator (it is one of the states inside the correlator K); the second costs two matrix products.
- **Undo the entangled state's phases.** O' = A Pi^dag, i.e. O'[c][L] = A[c][L ^ (D-1)] i^{-e(L)}, where |I> = SUM_L i^{e(L)} |L>|~L> is the
  layout's maximally entangled state. This is a column permutation and quarter turns: exact. Without it q(l) comes out wrong by 0.3 to 0.8.
- **The Pauli transform.** For each x in 0..D-1, take the diagonal f_x(c) = O'[c][c ^ x] and its Walsh-Hadamard transform F_x(z) =
  SUM_c (-1)^{popcount(z & c)} f_x(c) (additions and subtractions only). Then T_{x,z} = i^{popcount(x & z)} F_x(z) = tr(sigma_{x,z} O'), with
  sigma_{x,z} the Hermitian Pauli string with X where x has a bit and z does not, Y where both do, Z where only z does. Under the
  Jordan-Wigner map every Hermitian Majorana string is plus or minus exactly one sigma_{x,z}, so no sign bookkeeping is needed for its square.
- **The size.** Run the qubits from the highest down with a parity bit p = 0: X or Y on a qubit counts 1; Z counts 2 if p = 0 and 0 if
  p = 1; the identity counts 2 if p = 1 and 0 if p = 0; p flips after each odd count. The result is the number of Majoranas in the string.
- **The distributions.** P(l) = SUM_{|P| = l} |T|^2 / D and q(l) = omega^-2 SUM_{|P| = l} T^2 / D, where omega is the global phase the
  layout gives its normalized |I> ((1 + i)/sqrt2 when the number of qubits a side is odd, so omega^-2 = -i; 1 when it is even). These are
  KSY's p(l) and q(l) for rho^{1/2} psi_i(t): the reviewer checked them against KSY's own operator to 1.9e-15 at N = 8 to 16.
  Only odd sizes occur (psi is fermion-odd), and only the diagonals with an odd number of bits in x carry anything; the others are exact
  zeros, and the program skips them (or, with `windcheck=1`, computes them, reports the largest, which must be 0, and sums them into the even sizes).
- **The link to the protocol, an exact identity.** The coupling e^{i mu V} multiplies a Majorana string of size l, applied to |I>, by exactly
  e^{i mu (l - N/2)/(4N)}, because each term W_j of V is +1 on strings containing Majorana j and -1 on the others. Therefore

      C(mu) = <tfd| psi_r^i(t) e^{i mu V} psi_l^i(-t) |tfd> = -(i/2) e^{-i mu/8} SUM_l q(l) e^{i mu l/(4N)}

  (for the 2022 model: N -> the number of Majoranas in the exponent's prefactor and the coupling's normalization in the denominator).
  So 2|C(mu)| is exactly KSY's "line-fit" function C(m) at m = -mu/(4N): the best any coupling of either sign can do is C* = max_m C(m), and
  the sign of mu that unwinds the phases of q(l) is the one that teleports. The protocol's own K (`2_syk_exact.md`) has the coupling on both
  sides of it and is **not** a function of q(l); it is checked separately, against the program's own exact value.

## The arithmetic, and what is checked for free

Each entry of O' is rounded once from Q.60 to Q.56 (the only rounding before the output). The Walsh-Hadamard sums are 64-bit integer
additions: by Cauchy-Schwarz |F_x(z)| <= sqrt(D) times the diagonal's 2-norm, and the diagonals' norms sum to 1, so no partial sum exceeds
2^61 at D = 512. T^2 and |T|^2 are exact 128-bit products summed in 128 bits per size; the per-size sums are rounded once to Q.60 at the
output (a shift of 52 + n bits). Error bound on paper, worst case: 2 D^{3/2} 2^-57 on any q(l) or P(l) (2e-14 at N = 14, 1.3e-13 at N = 18),
with the actual error far below it because the roundings do not align (see the identity defects below). C(mu) from q(l) uses the phases
e^{i mu l/(4N)} built as powers of the small rotation e^{i mu/(8N)}, one rounding per step, so at most a few dozen units of 2^-60.

Printed with every operator and time: SUM_l P(l) - 1 in units of 2^-60 (must be a few units), the largest even-size sum (0 by construction unless `windcheck=1`, which sums the even-popcount diagonals
into the even sizes and so makes it a real check; it was 0 in every check), W =
SUM_l q(l) (must be real and the same at every t: it is tr(psi rho^{1/2} psi rho^{1/2})), and, with `windcheck=1`, C(mu) computed a second
way, from the states (<B| e^{i mu V} |A>/2 with B = e^{iHt} g_i^R e^{-iHt} T), beside C(mu) from q(l), with their difference in units.
The two routes share nothing but the state A, so their agreement checks the permutation, the transform, the size rule, omega^-2 and the
coupling's phase formula together.

## Development runs on the PC (labelled; the machines' runs follow below as they land)

**The acceptance (MEASURED, PC build, 11-12 12:48 AM).** With the winding on, the protocol chain at N = 8 (three seeds, the standard
arguments) is version 5's, `6ff9e3ebd0fa0c7f`. At N = 14, seed 1, for the first two operators at KSY's times t = 2.8 (the dt = 2.8 series)
and t = 4.0: every q(l) and P(l) agrees with the reviewer's independent float64 computation to at most 4.9e-13 (the float64's own width for
sums of 8191 terms), the eight values of 2|C(mu)| equal the reviewer's printed digits, and W = 0.786902343281 (psi_1) and 0.797399709523
(psi_2), real to the unit and the same at both times. The identity between the two routes to C(mu) holds to 1-11 units of 2^-60 in all
108 checks at N = 8 and 20-44 units in the 8 checks at N = 14. Every even-popcount diagonal and every even size is exactly 0.

**KSY's fit numbers from the exact q(l)** (COMPUTED by `code/wind_summary.py` from the printed lines; beta = 4; the definitions are KSY's:
r_l = |q(l)|/P(l), r = SUM |q(l)|, rbar = (r - W)/(1 - W), C(m) = |SUM_l q(l) e^{-iml}|, C* = max C, m* = argmax, chi = (C* - L)/(r - L)
with L = max(W, M), M = 2|q(l1)| + 2|q(l2)| - r for the two largest |q|):

| N | seed | t | op | W | r | rbar | C* | m* | N m* | chi | 2abs(C(-12))/C* | 2abs(C(+12))/C* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 14 | 1 | 2.8 | 1 | 0.786902343 | 0.901104 | 0.5359 | 0.898487 | +0.2432 | +3.404 | 0.9771 | 0.9982 | 0.6157 |
| 14 | 1 | 2.8 | 2 | 0.797399710 | 0.909096 | 0.5513 | 0.906496 | +0.2419 | +3.386 | 0.9767 | 0.9984 | 0.6234 |
| 14 | 1 | 4.0 | 1 | 0.786902343 | 0.873485 | 0.4063 | 0.872773 | +0.2047 | +2.866 | 0.9918 | 0.9998 | 0.6385 |
| 14 | 1 | 4.0 | 2 | 0.797399710 | 0.881350 | 0.4144 | 0.880360 | +0.2049 | +2.868 | 0.9882 | 0.9998 | 0.6535 |

Per size, N = 14 seed 1, psi_1, t = 2.8 (the exact values rounded to five places):

| l | P(l) | abs q(l) | arg q(l) - arg q(1) | r_l |
|---|---|---|---|---|
| 1 | 0.23013 | 0.22369 | +0.000 | 0.9720 |
| 3 | 0.28298 | 0.26651 | +0.606 | 0.9418 |
| 5 | 0.33485 | 0.28327 | +0.933 | 0.8460 |
| 7 | 0.12528 | 0.10732 | +1.597 | 0.8566 |
| 9 | 0.02567 | 0.01946 | +1.846 | 0.7584 |
| 11 | 0.00109 | 0.00085 | +2.390 | 0.7765 |
| 13 | 0.00001 | 0.00000 | +3.100 | 0.2222 |

**What they say.** At N = 14 the phase of q(l) climbs about 0.24 radians per unit of size at t = 2.8, against the coupling's 12/(4N) = 0.214:
the winding sits on the coupling (N m* = 3.4 against mu/4 = 3). The teleporting sign mu = -12 reaches 99.8 percent of the best possible
line C*; the wrong sign reaches 62 percent, because it winds the phases the other way, past pi by l = 9 to 11. chi, KSY's measure of how
much of the correlation comes from a straight line rather than from one or two dominant sizes, is 0.98 to 0.99; the experiment's learned
model, by the reviewer's float64 computation, shows the sign asymmetry too but with slopes about twice its coupling and a chi that varies
from 0.07 to 0.997 across operators, where SYK's is uniform. Across N (the reviewer's float64, seed 1, psi_1, t = 2.8, to be replaced by the
machines' exact values): N m* = 2.31, 2.64, 2.98, 3.40, 3.36 at N = 8, 10, 12, 14, 16, crossing the coupling's 3 near N = 12.

## The machines' runs

<!-- pair-begin -->
- N = 8, three seeds, t = 2.8: started 06:59:28 AM, both machines line for line at 06:59:48 AM; protocol chain `1fb4f612eb0c3575`, winding chain `85e0da202e9a53c8` (both the PC build's).
- N = 10, three seeds, t = 2.8: started 07:00:08 AM, both machines line for line at 07:00:29 AM; protocol chain `bb5e6efa00d910e5`, winding chain `2bee1f5c161bfe55` (both the PC build's).
- N = 12, three seeds, t = 2.8: started 07:00:29 AM, both machines line for line at 07:00:59 AM; protocol chain `663739f23e1e1cc0`, winding chain `98c9dfd66af8675b` (both the PC build's).
- N = 14, three seeds, t = 2.8: started 07:01:00 AM, both machines line for line at 07:02:38 AM; protocol chain `10a186a9f6e5fa44`, winding chain `ab106e45624bd39a` (both the PC build's).
- the 2022 experiment's model, all seven operators, t = 2.8: started 06:59:48 AM, both machines line for line at 07:00:08 AM; protocol chain `347d411caa01caff`, winding chain `60c2a928ab4105d0` (both the PC build's).
- N = 8, three seeds, t = 1 to 6: started 06:56:59 AM, both machines line for line at 06:57:26 AM; protocol chain `6ff9e3ebd0fa0c7f`, winding chain `8834d99dd0bb4f02` (both the PC build's).
- N = 10, three seeds, t = 1 to 6: started 06:57:26 AM, both machines line for line at 06:58:00 AM; protocol chain `b101fd5a26c155fe`, winding chain `8137da08c0f42a68` (both the PC build's).
- N = 12, three seeds, t = 1 to 6: started 06:58:00 AM, both machines line for line at 06:59:28 AM; protocol chain `1b590e73c01e1f3c`, winding chain `b089fc2bab21c535` (both the PC build's).
- N = 14, three seeds, t = 1 to 6: started 07:02:38 AM, both machines line for line at 07:10:52 AM; protocol chain `9d2559bc3d7810dd`, winding chain `40835081d0b3c8ac` (both the PC build's).

**KSY's fit numbers from the machines' exact q(l), beta = 4, t = 2.8** (COMPUTED by `code/wind_summary.py` from the lines above;
every seed and both operators; the machine's own 2|C(mu)| beside the fit):

| N | seed | t | op | W | r | rbar | C* | m* | N m* | chi | 2abs(C(-12))/C* | 2abs(C(+12))/C* | machine 2abs(C(-12)), 2abs(C(+12)) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 8 | 1 | 2.8 | 1 | 0.852058831 | 0.952440 | 0.6785 | 0.949124 | +0.2892 | +2.314 | 0.9670 | 0.9908 | 0.5258 | 0.940387, 0.499004 |
| 8 | 1 | 2.8 | 2 | 0.833036532 | 0.949196 | 0.6957 | 0.941122 | +0.3296 | +2.637 | 0.9305 | 0.9978 | 0.5358 | 0.939026, 0.504228 |
| 8 | 2 | 2.8 | 1 | 0.903193015 | 0.949709 | 0.4805 | 0.941159 | +0.1902 | +1.522 | 0.8162 | 0.9631 | 0.6769 | 0.906396, 0.637086 |
| 8 | 2 | 2.8 | 2 | 0.791685752 | 0.925321 | 0.6415 | 0.922653 | +0.3704 | +2.963 | 0.9800 | 1.0000 | 0.5037 | 0.922632, 0.464732 |
| 8 | 3 | 2.8 | 1 | 0.753199442 | 0.892553 | 0.5646 | 0.890564 | +0.3972 | +3.178 | 0.9857 | 0.9995 | 0.4992 | 0.890114, 0.444562 |
| 8 | 3 | 2.8 | 2 | 0.740707061 | 0.913208 | 0.6653 | 0.909531 | +0.4306 | +3.445 | 0.9787 | 0.9967 | 0.4431 | 0.906567, 0.403006 |
| 10 | 1 | 2.8 | 1 | 0.831261941 | 0.926210 | 0.5627 | 0.922471 | +0.2641 | +2.641 | 0.9606 | 0.9981 | 0.5994 | 0.920736, 0.552893 |
| 10 | 1 | 2.8 | 2 | 0.827557190 | 0.935412 | 0.6255 | 0.928159 | +0.2848 | +2.848 | 0.9327 | 0.9997 | 0.5973 | 0.927863, 0.554360 |
| 10 | 2 | 2.8 | 1 | 0.796957650 | 0.926786 | 0.6394 | 0.918454 | +0.3382 | +3.382 | 0.9358 | 0.9983 | 0.5834 | 0.916870, 0.535846 |
| 10 | 2 | 2.8 | 2 | 0.874571311 | 0.945731 | 0.5673 | 0.935939 | +0.2170 | +2.170 | 0.8624 | 0.9904 | 0.6639 | 0.926924, 0.621396 |
| 10 | 3 | 2.8 | 1 | 0.810334944 | 0.911733 | 0.5346 | 0.910304 | +0.2795 | +2.795 | 0.9859 | 0.9994 | 0.5822 | 0.909746, 0.529963 |
| 10 | 3 | 2.8 | 2 | 0.756551817 | 0.904425 | 0.6074 | 0.898089 | +0.3519 | +3.519 | 0.9572 | 0.9964 | 0.5391 | 0.894832, 0.484149 |
| 12 | 1 | 2.8 | 1 | 0.813727163 | 0.917956 | 0.5596 | 0.914356 | +0.2484 | +2.980 | 0.9655 | 1.0000 | 0.6105 | 0.914352, 0.558231 |
| 12 | 1 | 2.8 | 2 | 0.809618427 | 0.917305 | 0.5656 | 0.915163 | +0.2501 | +3.001 | 0.9801 | 1.0000 | 0.5959 | 0.915163, 0.545331 |
| 12 | 2 | 2.8 | 1 | 0.822302898 | 0.918031 | 0.5387 | 0.915058 | +0.2408 | +2.890 | 0.9689 | 0.9998 | 0.6270 | 0.914919, 0.573781 |
| 12 | 2 | 2.8 | 2 | 0.817575280 | 0.921375 | 0.5690 | 0.916898 | +0.2506 | +3.007 | 0.9569 | 1.0000 | 0.6191 | 0.916898, 0.567662 |
| 12 | 3 | 2.8 | 1 | 0.772167841 | 0.900466 | 0.5631 | 0.892875 | +0.2850 | +3.420 | 0.9408 | 0.9978 | 0.5881 | 0.890950, 0.525122 |
| 12 | 3 | 2.8 | 2 | 0.802756713 | 0.909518 | 0.5413 | 0.904465 | +0.2501 | +3.001 | 0.9527 | 1.0000 | 0.6074 | 0.904465, 0.549337 |
| 14 | 1 | 2.8 | 1 | 0.786902343 | 0.901104 | 0.5359 | 0.898487 | +0.2432 | +3.404 | 0.9771 | 0.9982 | 0.6157 | 0.896825, 0.553194 |
| 14 | 1 | 2.8 | 2 | 0.797399710 | 0.909096 | 0.5513 | 0.906496 | +0.2419 | +3.386 | 0.9767 | 0.9984 | 0.6234 | 0.905006, 0.565122 |
| 14 | 2 | 2.8 | 1 | 0.825067007 | 0.919250 | 0.5384 | 0.915853 | +0.2139 | +2.994 | 0.9639 | 1.0000 | 0.6492 | 0.915853, 0.594544 |
| 14 | 2 | 2.8 | 2 | 0.791902792 | 0.905853 | 0.5476 | 0.902376 | +0.2404 | +3.365 | 0.9695 | 0.9985 | 0.6159 | 0.901009, 0.555761 |
| 14 | 3 | 2.8 | 1 | 0.791587220 | 0.908538 | 0.5612 | 0.905850 | +0.2433 | +3.407 | 0.9770 | 0.9981 | 0.6093 | 0.904134, 0.551951 |
| 14 | 3 | 2.8 | 2 | 0.791731595 | 0.909453 | 0.5652 | 0.905049 | +0.2496 | +3.494 | 0.9626 | 0.9974 | 0.6199 | 0.902665, 0.561013 |
| 8 | 1 | 2.8 | 1 | 0.671952823 | 0.967658 | 0.9014 | 0.966769 | +0.6588 | +5.270 | 0.9970 | 0.8995 | 0.4558 | 0.869606, 0.440638 |
| 8 | 1 | 2.8 | 2 | 0.751398424 | 0.955065 | 0.8192 | 0.953664 | +0.5364 | +4.292 | 0.9931 | 0.9553 | 0.5535 | 0.911000, 0.527833 |
| 8 | 1 | 2.8 | 3 | 0.622718737 | 0.906206 | 0.7514 | 0.901507 | +0.5863 | +4.690 | 0.9801 | 0.9204 | 0.3946 | 0.829727, 0.355767 |
| 8 | 1 | 2.8 | 4 | 0.828425060 | 0.978633 | 0.8755 | 0.916832 | +0.5199 | +4.159 | 0.0725 | 0.9822 | 0.7602 | 0.900524, 0.697012 |
| 8 | 1 | 2.8 | 5 | 0.609474913 | 0.901046 | 0.7466 | 0.899575 | +0.6396 | +5.117 | 0.9938 | 0.8982 | 0.4205 | 0.808043, 0.378309 |
| 8 | 1 | 2.8 | 6 | 0.732588475 | 0.966148 | 0.8734 | 0.963182 | +0.7877 | +6.302 | 0.9546 | 0.8984 | 0.6057 | 0.865319, 0.583398 |
| 8 | 1 | 2.8 | 7 | 0.842254091 | 0.973397 | 0.8314 | 0.900320 | +0.3727 | +2.982 | 0.1778 | 0.9978 | 0.7759 | 0.898332, 0.698601 |

The slope against N (seed 1, the first operator, t = 2.8): N m* = +2.314 at N = 8, +2.641 at N = 10, +2.980 at N = 12, +3.404 at N = 14, against the coupling's exact 3.

Per size, N = 14, three seeds, t = 1 to 6, the first operator, from the machines' lines (the exact values rounded to five places):

l = 1 p 0.57560 |q| 0.57422 arg-arg1 +0.000 r_l 0.9976
- l = 3 p 0.16726 |q| 0.16636 arg-arg1 +1.552 r_l 0.9946
- l = 5 p 0.23016 |q| 0.18110 arg-arg1 +0.175 r_l 0.7868
- l = 7 p 0.02020 |q| 0.01699 arg-arg1 +1.789 r_l 0.8411
- l = 9 p 0.00669 |q| 0.00484 arg-arg1 +0.223 r_l 0.7234
- l = 11 p 0.00007 |q| 0.00005 arg-arg1 +1.843 r_l 0.7372
- l = 13 p 0.00001 |q| 0.00000 arg-arg1 -0.027 r_l 0.9460

The raw WINDING and C(mu) lines of one machine: `results/2_syk_winding_board1.txt`.
<!-- pair-end -->
Follow below as they land (this night's queue). Each run both machines line for line; the protocol chain equal to version 5's PC build
for the same arguments; the winding chain printed beside it.

## Files

`code/wh_syk6.c` (the program, with the derivation in its header), `code/wind_summary.py` (KSY's fit numbers from the printed lines),
`results/2_syk_winding_board1.txt` (the raw output of one machine, when it lands).
