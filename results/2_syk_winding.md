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
  zeros, and the program skips them (or, with `windcheck=1`, computes them and reports the largest, which must be 0).
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

Printed with every operator and time: SUM_l P(l) - 1 in units of 2^-60 (must be a few units), the largest even-size sum (must be 0), W =
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

Follow below as they land (this night's queue). Each run both machines line for line; the protocol chain equal to version 5's PC build
for the same arguments; the winding chain printed beside it.

## Files

`code/wh_syk6.c` (the program, with the derivation in its header), `code/wind_summary.py` (KSY's fit numbers from the printed lines),
`results/2_syk_winding_board1.txt` (the raw output of one machine, when it lands).
