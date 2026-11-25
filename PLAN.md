# The plan: three computations, exact

Preston Ericson, 2026-11-11. The first line, again: the programs here compute metrics (descriptions of how space and time are shaped)
and quantum protocols; computing a metric does not create it. Each of the three computations below is a published equation or claim
that its authors left uncomputed, left as "future work", or that is in dispute, and each result is a number others can check with the
code in `code/`.

**In plain words.** This is the plan of the work: what each of three published pieces of physics left open, what each computation must
deliver, in what order, and, at the end, where each stands. Computation 1 is about the "humanly traversable" wormhole of Maldacena and
Milekhin, whose mouths are wide enough that a person would survive its tides (gravity's stretching of a body): first the shape of its
extra dimension (a fifth dimension that their model adds to the four of space and time; 1a), then the full geometry the paper did not
write down (1b, the research part, still under way). Computation 2 is the wormhole-teleportation protocol (a recipe that inserts a
quantum bit on one side of two linked systems and reads it out on the other) behind the 2022 quantum-computer experiment, computed
exactly at sizes past the experiment's, to test the experiment's disputed signature (the pattern it reported as a sign of wormhole
behaviour). Computation 3 checks two warp-drive papers exactly: how much negative energy the classic (Alcubierre) design needs, and
whether the 2024 shell of positive-energy matter passes the energy conditions (standard tests that matter behaves like ordinary
matter) and up to what speed.

## What "exact" means here

Fixed-point integer arithmetic in a declared format, deterministic rounding, every constant written as an integer, so that two
machines or two runs give the same bits and anyone with the code gets the same bits. Error bounds derived on paper beside the code.
Floating point is used to cross-check on a PC and in the few places named under [Why exact](#why-exact) below; no C program's printed result depends on it.
Every number is labelled: MEASURED (printed by the two machines, the run named), COMPUTED (derived, or printed by a PC build as a
labelled development run, the inputs named) or ESTIMATE (and what it would take to compute).

## Where it runs

The programs are plain C computing in 64-bit integers with no math library (the few floating-point lines are named under [Why exact](#why-exact)
below). Each result is run on two separate ARM Cortex-A9 boards at the same time, and the two outputs must agree line for line before
it is called measured (PC builds run first, as labelled development runs); the same source compiled on any other machine prints the
same digits. That is the whole point of the method: the result is a property of the
equations and the declared arithmetic, not of one computer.

## 1. The Maldacena-Milekhin solution left uncomputed (arXiv:2008.06618)

Their words (section 3.2): "We have argued that the solution exists, but one could imagine trying to find it numerically in a more
explicit form." Two parts, honestly sized:

- **1a (days): the flow.** Their eq. (3.19), two ODEs in rho for lambda(rho) and sigma(rho) with the scale l_B:

      0 = -6 + 3 l_B^-4 e^{-4 sigma} + lambda'^2 + 4 lambda' sigma' + sigma'^2
      0 =  4 - 4 l_B^-4 e^{-4 sigma} - 2 lambda' sigma' - 2 sigma'^2 - sigma''

  UV: lambda = sigma = rho (AdS5); IR: sigma -> const (e^{-sigma} = l_B), lambda' = sqrt 3 (AdS3 x R2). Their Appendix A solves it
  with l_B = 1: sigma = 0, lambda = sqrt 3 rho is a solution; linearize around it, take the growing exponential as the IR boundary
  condition, integrate out to the UV; their Figure 3 shows sigma and lambda. Then eq. (A.30), the same flow on a two-sphere of radius
  1 with q flux units (the extra -e^{-2 sigma} and +e^{-2 sigma} terms), whose IR radius is e^{-sigma} = l_B + l_B^3/16 + ... (A.31).
  DELIVERABLE: both flows integrated bit-exact, the IR constants to the printed digits, (A.31) checked against the numerics, Figure 3
  reproduced as a table. This is the warm-up and the check that the arithmetic reaches their equations.
  PROGRAM: `code/wh_flow.c`, now version 3 (one rounding per update with the remainder carried, each flow on its own unstable manifold,
  the gauge invariants; `results/1a_flow.md`), and `code/wh_flow60.c`, version 4 in Q.56. Version 2, as first run: Q.40 fixed point in
  int64 with 128-bit intermediates, products and divisions rounded to nearest; exp by
  range reduction on ln 2 and 24 Taylor terms; sqrt exact to the floor; RK4 with step h = 2^-12; both flows started on the IR growing
  mode at rho = -8 with amplitude 2^-10 (a gauge: any other amplitude gives the same curve shifted in rho); tables every 0.25 in rho;
  the UV asymptotes; a printed self-test; a chain hash of every step. The first version (Q.28, floor rounding, amplitude 2^-20) was
  withdrawn before any run: its first increments in the deep IR were a sixth of one unit, which the floor made zero, so the flow would
  have stood still. Sizing the first step on paper caught it; the fixed-point checks alone had not. Checked on paper before the first run: the UV fixed point has sigma' = 1 and lambda' =
  -2 + sqrt 9 = 1 (their lambda = sigma = rho + const); the sphere flow's IR fixed point solves 4 + u - 4u^2 = 0 with u = e^{-2 sigma},
  u = (1 + sqrt 65)/8 = 1.1328, so e^{-sigma} = 1.0643 against (A.31)'s 1.0625 at l_B = 1 (their expansion is for small l_B).

- **1b (weeks, research): the patched five-dimensional geometry.** The full solution they did not write down: the AdS2 x S2 throat
  (their 2.7) glued to two extremal magnetic black holes (2.4) at rho = +/- rho_c, the fifth dimension flowing from AdS5 near the
  Planck brane to AdS3 x S2 at distance l_B (3.25), the exterior region of length d/(pi ell) where the field weakens; their section
  3.2 sketches it on R x S3 with the mouths antipodal and says "It is possible to describe the metric more explicitly and find the
  solution numerically, but we will not do it here." DELIVERABLE: that explicit metric, numerically, with the Einstein equations'
  residuals printed; the parameters of their (3.27) (re = 1.5 x 10^7 m, R5 = 50 micrometres, c2 ~ 7 x 10^72, ell ~ 3,000 light
  years) reproduced from the solution rather than from the estimates. Needs: the 5D Einstein-Maxwell equations with the S2 symmetry
  reduced to a 2D elliptic problem, a relaxation solver in fixed point, and a careful write-up. The result would be the first explicit
  humanly-traversable-wormhole geometry; it still "has no plausible mechanism for formation" (their conclusion), and we say so.

## 2. The 2022 "wormhole on a quantum processor", done honestly

The experiment (Jafferis et al., Nature 612, 51, 2022): two copies of a machine-learned SPARSE SYK model (N = 7 Majoranas a side),
the thermofield-double state, a qubit inserted on the left at time -t, the coupling e^{i g V} applied, the qubit read on the right at
+t; "wormhole teleportation" when the coupling has the negative-energy sign, not the positive; 9 qubits, 164 two-qubit gates on
Sycamore. The comment (Kobrin, Schuster and Yao, arXiv:2302.07897) argues the learned sparse Hamiltonian is too simple to carry
gravitational physics. The protocol itself: Gao and Jafferis, "A traversable wormhole teleportation protocol in the SYK model"
(arXiv:1911.07416).
DELIVERABLE: the FULL (not sparsified) q = 4 SYK protocol computed exactly for N = 7, 8, 10, 12, 14 Majoranas a side (dense exact
state vectors up to about 2^16 amplitudes, sparse Hamiltonian application), the teleportation fidelity and mutual information against
g, t and beta for both signs of the coupling, averaged over disorder samples with a declared integer seed; the "size winding" signature
checked as the 2022 paper and Kobrin, Schuster and Yao state it (the protocol paper defines none); then the SAME computation on the
learned sparse model of the 2022 paper (its couplings are in
the supplement) to see whether the signature stands or falls, which is the live dispute. Exactness: fixed-point complex amplitudes with
a declared format and a norm check every step; the random couplings from an xorshift generator with a declared seed so the ensemble is
reproducible bit for bit. Floating point on a PC only to cross-check the smallest N.

## 3. The energy conditions of the warp metrics

Bobrick and Martire (2021, arXiv:2102.06824, read whole): a general warp-drive spacetime class; "the first general model for
subluminal positive-energy, spherically symmetric warp drives"; Alcubierre optimizations "that decrease the negative energy
requirements by two orders of magnitude". Fuchs et al. (2024, arXiv:2405.02709, Class. Quantum Grav. 41, 095013, read whole): "a
constant-velocity subluminal warp drive that satisfies all of the energy conditions"; the metric is generated numerically with their
public toolkit Warp Factory (github.com/NerdsWithAttitudes/WarpFactory), which we read as the definition of their construction.

### 3a. The Alcubierre integral, 1D exact
Bobrick-Martire eq. (3): T^00 = -(1/8 pi) (rho^2 v^2 / 4 r^2) f'^2, negative everywhere; their (A.7): E = -(v^2/12) INT r^2 f'^2 dr.
For their optimal wall f = min(r0/r, 1) this is E = -v^2 r0/12 exactly (f' = -r0/r^2 and INT r0^2/r^2 dr = r0). For Alcubierre's tanh
wall (their eq. 2) with a thin wall it is -v^2 R^2 sigma/36 (sech^4 integrates to 4/3), so the ratio between the two is R sigma / 3,
not a constant 3: their "about a factor of three" is the value at R sigma = 9. DELIVERABLE: the integral in Q.40 for R sigma = 4, 6,
8, 9, 10, 12, 16 by Simpson's rule in fixed point on a declared grid (2^14 intervals on [0, 4]), the ratio tabulated, the asymptote
checked; the flattening theorem E -> E/alpha_X (their 4.2) checked on a 2D (x, rho) grid for alpha_X = 2 and 8 as a test of the
arithmetic. Their "two orders of
magnitude" is flattening by 10 times the wall factor: the table says what it is at each wall thickness.

### 3b. The 2024 warp shell: their construction reproduced in integers, then judged exactly
Their code (metricGet_WarpShellComoving.m, TOVconstDensity.m, alphaNumericSolver.m, compactSigmoid.m; the example W1_Warp_Shell.mlx)
fixes every number the paper left out:
- R1 = 10 m, R2 = 20 m, M = (R2 c^2 / 2G) / 3 = 4.489 x 10^27 kg (the paper's 4.49 x 10^27), shift 0.02 c, Rbuff = 0, sigma = 0;
- the radial sample: 10^5 points on [0, 1.2 x 42.29 m] (step 5.07 x 10^-4 m); density = M / (4/3 pi (R2^3 - R1^3)) between R1 and
  R2; pressure = the closed-form constant-density TOV solution with P(R2) = 0;
- MATLAB smooth() = a moving average, FOUR passes, span 4000 samples (2.03 m) for the pressure and 1.79 x 4000 = 7160 samples
  (3.63 m) for the density (the paper's text says 1.72; the code says 1.79), with MATLAB's shrinking windows at the ends;
- the mass M(r) by a cumulative trapezoid of 4 pi r^2 times the smoothed density; B = (1 - 2GM(r)/r c^2)^-1; a(r) by a cumulative
  trapezoid of (GM/c^2 + 4 pi G r^3 P/c^4)/(r^2 - 2GMr/c^2), pinned to 1/2 ln(1 - 2GM/r c^2) at the last sample; g_tt = -e^{2a};
- the shift S(r) = the compact sigmoid of their eq. (28) smoothed twice with span 4000; g_tx = -0.02 S(r).
Note that four passes of a 3.63 m average smear the density from about 2.7 m to about 27.3 m (each pass widens the support by half a
window), so "R1 = 10, R2 = 20" names the unsmoothed shell only; the exact profile says where the matter is.
Their Einstein tensor: the profiles interpolated (3rd-order Lagrange) onto a 300 x 300 x 5 Cartesian grid of 0.2 m, derivatives by
4th-order central differences, only the middle z-plane fully stencilled, time derivatives zero (one time slice: stationary). Their
observers: 100 directions on a Fibonacci sphere (golden-ratio spiral) times 10 speeds linspace(0,1) (so one "timelike" shell is
null), transformed to the Eulerian frame by a Cholesky factor of the metric; NEC and WEC are the minimum of T(k,k) over them; SEC with
the trace; DEC as -(sign sqrt |eta(Y,Y)|) of Y = -T V, which checks non-spacelike only (future-pointing is not checked). Their floor:
values ~10^39 J/m^3, noise ~10^34 (their Appendix A), i.e. 10^-5 to 10^-6 relative.

OUR VERSION: the same profiles computed in int64 fixed point (the trapezoids and the moving averages are finite integer sums, so they
reproduce bit for bit); then the Einstein tensor NOT by Cartesian differences but from the profiles themselves on an (r, theta) grid
with the shift direction as the polar axis. The metric is axisymmetric: g_tr = -beta cos theta, g_t theta = beta r sin theta, g_rr = B,
g_theta theta = r^2, g_phi phi = r^2 sin^2 theta, g_tt = -e^{2a}; the angular derivatives are analytic and the radial derivatives come
from the profiles' own definitions (a' is the integrand, M' = 4 pi r^2 rho), so no truncation error enters the curvature. The Eulerian
frame by an exact Cholesky in fixed point; the four conditions over THEIR exact 1000 observers (the Fibonacci angles in fixed-point
trigonometry) AND over a denser set; the minimum of each condition at every r and theta printed with its sign. Then three numbers:
1. The sign below their floor: are the values they call noise (their Figure A1, 10^30 to 10^35 outside the shell) zeros, or small
   violations? Vacuum must give exactly zero, so any residue is the interpolation and the stencil, and we can say so.
2. THE SHIFT LIMIT: the shift scanned (0.02, 0.05, 0.1, 0.2, 0.3, ...) at fixed shell until the first condition fails, and where
   (r, theta, which observer). The paper: "there is an upper limit ... This upper limit is a future direction of work."
3. Their light-ray test: along the axis a null ray has dt/dx = (-beta +- sqrt(beta^2 + e^{2a} B)) / e^{2a}, so the round-trip
   difference is 2 INT beta e^{-2a} dx / c; their 7.6 ns recomputed exactly from the same profiles, and the Shapiro delay against flat
   space, which they say stays a delay, not an advance.
SIZE: the profiles are 10^5-point 1D sums (seconds); the (r, theta) grid 2000 x 181 with 1000 observers times 16 components is about
6 x 10^9 multiply-adds per condition set, minutes per shift value in int64 on one Cortex-A9 core per machine (the measured runs took
25 to 137 minutes a grid), the two machines agreeing line for line as in 1a. Floating-point cross-check of the smallest grid on a PC only.
HONEST LINE: a positive-energy shell of 2.4 Jupiter masses in a 20 m radius, moving at a constant 0.02 c with no way to accelerate it,
is not a drive; the paper says so itself ("the question of how to make physical and efficient acceleration is one of the foremost
problems"). What we add is exactness, numbers they left open, and a public check.

## Order

1a first (it tests the arithmetic on real equations in days), then 3 (two papers, mechanical, high value: a published claim checked or
not), then 2 (the biggest, a live dispute), then 1b (the research project). Each result is a number others can check, written up with
its format, its bounds and its seed. None of them builds a door, and each write-up says so in its first line.

## Status (2026-11-12, 9:11 AM)

- 1a: DONE 2026-11-11 4:26 PM, `results/1a_flow.md`: both machines identical, 124 lines; the flow from the throat to AdS5, printed to 12
  places and good to about 7 (the review note on the page; version 3, run on the pair (the two boards) at 10:17 PM, brings the
  invariants to 2.5e-11).
- 3a: DONE 2026-11-11 4:28 PM, `results/3a_alcubierre.md`: both machines identical, 18 lines; the law I = R sigma/3 + (pi^2-6)/(36 R sigma). Error bound on paper: Simpson truncation below 8e-8 absolute at R sigma = 16 (3e-10 at 4), rounding
  below 8e-8 (each term carries about 20 roundings of 2^-40 on values up to 1024, and the weights integrate to the range 4): I to 7 digits.
- 3b part 1: DONE 2026-11-11 7:36 PM, `results/3b_shell.md`: both machines identical, 148 lines; the shell as published passes all four
  energy conditions (margins 0.014-0.023 inside the matter; zeros to the arithmetic's floor outside); its smoothing raises the mass by
  1.9 percent and spreads the matter to 2.7-27 m. (The plan's earlier line calling the pressure-kept cavity "the code's variant" was a
  misreading of Warp Factory's `TOVconstDensity.m`, which multiplies by the step density: the code zeroes the cavity pressure as the paper
  does. The pressure-kept case is kept as HYPOTHETICAL.) Part 2 (`wh_warp.c`, the shift, axisymmetric; the shift limit) is written and
  did not pass its first review (a division guard in the shared arithmetic, fixed in `whq.h` v2; a rounding floor on the sign tests; the
  cells next to the axis; a self-test with a constant shift), was rewritten with all four, and is DONE 2026-11-11 9:21 PM,
  `results/3b_warp.md`: both machines identical, 44 lines; the four conditions hold to beta0 = 0.03033 and fail together just above it,
  at 12.3 m on the equator (0.030330352784 to 0.030330429078 on the 0.2 m radial grid; 0.030325622559 to 0.030325698853 on the 0.05 m grid,
  both on the pair, 11:58 PM, the fine grid's chain equal to the PC build's). CORRECTION 12:05 AM from review: a truncated dr^2 in
  the kernel's second derivative lowered every crossing by 3.25e-5; the corrected limits are 0.030362854005 to 0.030362930299 (0.2 m) and
  0.030358123781 to 0.030358200075 (0.05 m), COMPUTED by the reviewer and printed digit for digit by version 3 on the pair (the 0.2 m
  grid from 1:52 AM, the 0.05 m grid from 2:27 AM, 11-12; chains 09f6b75f837adef2 and 4f27e30f55bd65b8, the PC build's;
  `results/3b_warp_v3_board1.txt`). The review confirmed the general relativity (the Ricci convention, the
  Eulerian frame, the four margins) against an independent derivation.
- 4a: DONE 2026-11-11 7:39 PM (v2) and REDONE 9:36 PM (v4), `results/4a_casimir.md`: both machines identical, 21 lines; the lattice law
  1 + (5 pi^2/24)/N^2 + O(1/N^4) toward the continuum coefficient pi^2/1440 for one scalar, exact to N = 32 (v2's Q.40 frequencies had
  lost digits past N = 16 and its table said 2.1/N^2; the review caught it the same evening. The first draft had used the electromagnetic
  pi^2/720, corrected in review before the first run.)
- 2: FIRST RUN DONE 2026-11-11 7:46 PM, `results/2_syk.md`: the full q = 4 SYK protocol exact at N = 8 over three seeds, both signs of
  the coupling and mu = 0 (the exact zero), and the 2022 experiment's five-term model at its own parameters; both machines identical,
  64 + 8 lines; a double-precision cross-check agrees to 1e-11. 9:26 PM, from review: the density-matrix formulas (2.63)-(2.66) are the
  paper's ensemble-averaged forms; the page is labelled "evaluated per seed", and v3 of the program (`code/wh_syk3.c`) computes the exact
  six-element matrix (GJ 2.15-2.20) beside them, verified against an independent float implementation to 1e-12. 9:36-9:43 PM: v2
  (`code/wh_syk2.c`, the factorized layout: left Majoranas on one half of the qubits, right on the other, D = 2^(N/2) a side,
  e^{-iHt} Psi = U_L Psi U_R^T) run on the pair at N = 8, 10 and 12 in 28 s, 63 s and 5.5 min, both machines identical and the N = 8
  chain equal to the PC build's. v3 ran on the pair at N = 8-14 (`results/2_syk_exact.md`, pushed 9:57 PM); v4 (Q.60 end to end, twelve
  digits verified against a float implementation, the time steps built once per seed) ran on the pair at N = 8-16 and on the 2022 model
  from 11:59 PM to 12:25 AM (`results/2_syk_exact_v4_board1.txt`); v4.1 derives the right-hand operators exactly and v5
  (`code/wh_syk5.c`) fits N = 18 in the machines' memory (a PC development run
  at N = 18 is on the page, labelled). The size-winding signature (derived in review; version 6) ran on the two machines at N = 8 to 16
  and for the 2022 model on 11-12, 6:56 to 7:32 AM, with the PC build's chains (`results/2_syk_winding.md`). Next: the pair's N = 18
  (N = 16 ran on the pair from 12:08 to 12:24 AM, 16.4 min for one seed, the PC's digits), the N-dependence of the sign
  asymmetry written up.
- 1a: a review note added to `results/1a_flow.md` (a rounding bias of 4-7e-8 on the UV constants from ties-away rounding, the sphere
  flow's starting mode, the gauge invariants to publish). v3 of `wh_flow.c` BUILT 10:05 PM from the reviewer's specification and
  reference program: one rounding per update with the remainder carried, each flow on its own unstable manifold from hard-coded start
  integers, stop at sigma > 16, the gauge invariants printed; the PC build's 160 lines are byte-identical to the reference's and its two
  chains equal; K_sigma, K_lambda agree with an independent 40-digit integration to 2.5e-11 (the arithmetic's whole error now). The pair
  ran it at 10:17:38 to 10:18:09 PM, 160 lines identical, chains equal to the PC build's (`results/1a_flow_v3_board1.txt`); its section
  is on the 1a page.
- 4a again: the reviewer certified the exact values (three routes to 25 digits) and found v4 still +1.7e-9 / +2.5e-9 at N = 32 / 24 from
  three roundings; v5 (10:10 PM) removes them and lands within 8.4e-12 of the certified ratio on the PC build; the law's second
  coefficient is 9.3762 (closed form on the page), not the fitted 9.6. v5 ran on the pair at 10:18:15 to 10:20:25 PM, chain
  708f274a7b9e6043, and the page's table is v5's.
- 1b: STARTED 2026-11-11 10:50 PM. The reviewer derives the reduced system in two verified stages (the ansatz and gauge, the 5D
  Einstein-Maxwell equations as a 2D elliptic system with the constraints separate, the residuals, the boundary conditions edge by edge;
  then the (3.27) parameters as formulas in the solved functions, the sizing, the 1a checks). The fixed-point relaxation kernel is built
  and verified ahead of the equations: `code/wh_relax.c`, a Q.60 red-black SOR with 128-bit residuals, tested on Poisson with an exact
  quartic solution (for which the five-point Laplacian is exact): the residual reaches its rounding floor, a few G^2 units of 2^-60, in
  about 6 G sweeps as the spectral radius predicts, and the solution is within 4e-17, 6e-16 and 8e-15 of exact at G = 32, 128 and 512
  (PC build); the pair ran version 3 at 12:28 to 12:32 AM (below). Stage 1 (11-12, 1:18 AM, checked independently in nearly every step): only
  the rotationally symmetric reductions are well-posed two-dimensional problems (the paper's flat-space pair has axial symmetry only and needs
  three coordinates and a strut); the paper's fine-tuned circle of mouths has no static solution without both extra positive energy and a positive cosmological
  constant on the brane, or matter with rho + 3p below zero (`results/1b_nogo.md`);
  the order chosen: the
  solver test, then the near-brane throat (an ODE matched to the flow of 1a; the flow has no inner end, so this is an approximation, not
  a closed boundary-value problem), then the paper's own rotationally symmetric example with antipodal mouths, which it called findable
  numerically; the no-go for the fine-tuned circle holds where the mouths' spacing is large against re/(2 ln(L_IR/R5)). The kernel is at version 4: a Robin boundary option eliminated into the first interior row
  and grids stretched toward a boundary with exact integer stencil weights, both landing on the known solution to tens of units of 2^-60;
  an adaptive relaxation parameter was tried and dropped (it misreads SOR's early transient). The pair ran version 3 at G = 128 and 256
  and the nonlinear G = 128 (12:28-12:32 AM, chains equal to the PC's); version 4 has not yet run on the machines. The system kernel
  (`code/wh_relaxs.c`, 2:00 AM) relaxes three coupled fields with nonlinear sources, the Newton step over the local 3 x 3 Jacobian: on the exact test it lands
  within 3 to 172 units of 2^-60 at G = 32 to 128 in about 8 G sweeps; where the source's own derivative is a third of the Laplacian's
  diagonal, the Newton step (block or diagonal) still lands exactly while the lagged plain update overshoots: its source passes the Q.60
  ceiling within two sweeps (with the source formed in 128 bits the program stops; the earlier 2e-2 "floor" was int64 wrap-around).
  Neither has run on the machines yet. The kernel in the geometry's own rows (`code/wh_relaxt.c`, 4:12 AM): rows uniform in t = ln(-y) from the cut
  to the deep end, the operator F_tt - F_t + e^{2t} F_xx normalized by its row diagonal into per-row weights, Dirichlet at the cut, the centred
  ghost form of Neumann on the other sides, a fixed relaxation parameter 2 - 90/N with the change guard (the reconciled sizing: never a field
  alone). Verified on an exact polynomial test (within 4e-14 to 5e-13 of exact at N = 64 to 256) and on the AdS3 x S2 reference carrying the
  reviewer's three P0 sources (P0: the reduced problem's simplest case, with no dependence along the circle) in deviation variables:
  with the reference's discretization defect subtracted, zero deviation is the exact discrete
  solution (zero corrections at sweep 1) and bumps of 0.01 and 0.1 relax back to it in 3,921 sweeps at N = 128 and 5,284 at N = 256; with the
  defect left in, the end deviation is the reference's discretization error, 4.8e-3 and 2.1e-2 (U and B) at N = 64, 1.2e-3 and 5.3e-3 at 128,
  3.0e-4 and 1.3e-3 at 256: 4.0x per doubling on every row including the Neumann one, while V's defect is exactly zero (its bracket vanishes
  for the reference by construction). MEASURED on the pair at 6:14 to 6:18 AM: the four configurations (the exact test at N = 128, the
  perturbed reference at N = 128, the defect left in at N = 128 and 64) each identical on the two machines with the PC build's chains. The real
  P0 takes this form with the brane's cut data (the reviewer's series in the conformal coordinate, with the v^4 log term) at Q2 = 6, where
  the reference's w solves 12 w^2 - 3 w - 12 = 0, and rows in Q.56 for the deep cut (t = -11 passes the Q.60 ceiling); the reviewer's spec
  (4:00 AM) is being checked independently. DONE on the PC at 4:41 AM (`results/1b_p0.md`): the data cannot be imposed at once (every direct start
  leaves the basin inside 63 to 128 sweeps while a quarter of the data converges), so the brane radius is ramped over ten stages (nine
  doublings, a = 3.9e-5 to 0.02); every stage
  converges in 2,559 to 3,274 sweeps; the centre value U_c = -0.594115 (N = 128), -0.594579 (N = 256), extrapolating to -0.594734 =
  the reviewer's independent plain-row extrapolation; the error falls 4.00x per doubling (N = 512: U_c = -0.594695). MEASURED on the pair
  at N = 128 at 6:12 AM (191 identical lines, the PC's chain 46e8b5ae9ae59aa2) and at N = 256 at 6:58 AM (30,608 sweeps, the PC's chain
  5a89a7f4c5f17234); CORRECTED 7:45 AM by the review: a 64-bit sum of two neighbours in the point update wrapped on the cut rows
  (U high by about 3.5e-7); fixed, confirmed to every decimal (chains 09d5ac7eddb94cd4, 7b497cb99b5103b9), Richardson -0.594734016 against the
  continuum -0.594734119; the pair reruns the fixed kernel;
  the reviewer's own defect-corrected run and a direct check at alpha = 1/8 against its specification's continuum values follow (the check,
  6:26 AM: V and B agree to 1e-8 and 2e-7; U agrees as well, to 6e-8 (COMPUTED by the reviewer), once the kernel's int64 wrap
  (found in review) was fixed). P3, the
  near-brane throat, DONE on the PC at 6:25 AM (`results/1b_p3.md`, `code/wh_p3.c`): the reviewer's verified specification coded as an
  integer shot in Q.56, bit-identical to its mimic in all three cases (the raw ells, the halvings, the event depths to the step); the
  machines have not yet run the three cases. P1's kernel (P1: the paper's own rotationally symmetric example with antipodal
  mouths, which depends on the circle direction; `code/wh_dtk.c`, 6:58 AM): the Einstein-DeTurck form's five fields on a periodic
  circle coordinate by a depth coordinate, the principal part h^ab d_a d_b with the reference metric's inverse per point (the weights and
  the row scale formed once), the cross-derivative term over the four diagonal neighbours, four colours with all fields from one state; on a
  manufactured problem with a non-diagonal reference metric the largest error of each field falls 3.77x to 4.15x per doubling from
  Ge = N = 16 to 32 and 3.94x to 4.01x from 32 to 64 (COMPUTED, PC development runs; from 16 to 32, 3.83x to 4.08x without the cross term
  and 3.77x to 4.14x with a lagged coupling lam = 0.1); its production edges pass too (7:29 AM): deviation variables
  from a reference table that is not a solution, with the reference's own discrete defect optionally subtracted (the zero start then stops
  at sweep 1), two reflection lines with a parity per field (Neumann for the even fields, zero for the odd), and a Neumann deep end by the
  mirror ghost, each 4.0x per doubling alone and together. Two more Q.60 ceiling crossings were found in the tests on the way (the shared
  cosine, which squared arguments past 2.83 after reducing them only to [-pi, pi] and now reduces to [-pi/2, pi/2], and a test curvature
  of exactly 8), neither in the kernel. The equations and the reference metric are to come from the reviewer's P1 specification; the
  boundary data are built (`code/wh_p1data.c`, 7:44 AM): the paper's
  SO(3) example as one conformal class in the cylinder frame, the throat, mouth and S3 pieces blended in ln rho_o by the smooth step, the
  local invariant w, the Fefferman-Graham coefficients and the flux log term at the cutoff z = lhat/30, all in fixed point with the circle
  coordinate in Q.56; at the first grid (gamma = 1e3, a = ell, lhat = 1/4) it reproduces the reviewer's independent table of rho_o and w
  at the edges and centres of the blend intervals to the table's precision (PC development run). The yardstick (`results/1b_cost.md`,
  7:36 AM; COMPUTED on the PC by `code/wh_cost.py`): the paper's own chain evaluated at every mouth size with our constant, the
  human-tides condition dropped; it says what each size costs (charge, field, negative energy, crossing time) and that no size is cheap:
  the negative energy grows as the mouth shrinks, from 5 mm to 1 AU the dark magnetic charge is 5e31 to 2e45 Dirac units (for g4 = e) and
  the field at the mouth at least 2e7 T in energy density, and the paper's own wormhole takes about ten thousand years to cross from
  outside (10,230 in the chain). Any proposal to scale the thing down is measured against that page.

## Details from the README

These sections stood on the README until 2026-11-12, when the README was shortened to a front page. They are moved here unchanged,
except for one reworded sentence each in "How it was built" and in the SYK list.

### Why exact

Every C program here computes in 64-bit integers in a declared fixed-point format (40, 56 or 60 fractional bits, with 128-bit
intermediates built from 32-bit halves; [`code/whq.h`](code/whq.h), [`code/whq60.h`](code/whq60.h), [`code/whq56.h`](code/whq56.h)) and uses no math library.
Floating point appears in three places, and no C program's printed result depends on it:
the 2022 model's five decimal couplings pass once through a double on their way to integers (the same integers as exact decimal rounding);
two relaxation kernels (`wh_relax.c`, `wh_relaxs.c`) show their integer residual to three digits on their progress lines; and three
helpers are not integer programs ([`code/wind_summary.py`](code/wind_summary.py), float64 fits of the printed exact sums; [`code/wh_cost.py`](code/wh_cost.py), the yardstick page's
formula chain in 30-digit mpmath; [`code/p2_nogo/`](code/p2_nogo/), the reviewer's sympy and mpmath checks).

Exponentials come from a range reduction on ln 2 and a Taylor series with a declared number of terms; square roots from integer Newton
iteration or a digit-by-digit 128-bit root; differential equations from a fixed-step Runge-Kutta scheme whose updates are rounded once
with the remainder carried; matrix exponentials from powers accumulated in 128 bits and rounded once. The same source compiled anywhere
prints the same digits, so a result is a fact about the equations, not about one machine's rounding.

A result is run on two separate ARM boards, and the two outputs must agree line for line before it is called measured; until then a PC
build's numbers appear only as labelled development runs. Error bounds are derived on paper beside each program, and every number on a
page is labelled: MEASURED (printed by the two machines, with the run named, or a published value, named), COMPUTED (derived, or printed
by a PC build as a development run, with the inputs named) or ESTIMATE.

### How it was built, and how it was checked

- Each program was sized on paper before its first run (the format, the step, where overflow could sit), and its first run on the
  machines goes to both at once (a PC development build runs first).
- An independent reviewer (another AI session working from copies, with its own derivations, its own float and 128-bit float builds and
  its own symbolic checks) reads every program and every page. In the first evening it found a division the shared header never guarded,
  a coefficient for the wrong field, a published Casimir table with digits lost past N = 16, a second one off by one unit in the ninth
  decimal, a flow started off its manifold, and the paper's density-matrix formulas being the ensemble-averaged forms rather than the
  per-realization ones. Every one of those is on its page, with the reviewer credited, the corrected version run again on the two
  machines, and the lineage of versions kept. After 10 PM on 11-11 it found a truncated dr^2 in the warp kernel worth 3.25e-5 in the
  published shift limit (the correction went up within a minute of its reading, before the fix; the fixed program then ran on the machines),
  derived the size-winding observable for our layout, and found a step of the 1b plan that was a matched approximation rather than a
  closed problem; and an early claim about how the first night's runs were split did not match the program's own output, which is now
  written where that claim stood. At 3:45 AM on 11-12 it found the shared division routine's guard a factor of two too weak in all seven copies (six
  programs and [`code/whq60.h`](code/whq60.h); it stopped at sixteen times the divisor instead of eight; in between, a quotient with the wrong sign came
  back silently): every copy now carries two checks, and every published chain of the programs that carry it reruns identical from the
  fixed sources, so no published run
  had passed through the range. The same morning the P1 kernel's own tests found that the shared header's cosine ([`code/whq60.h`](code/whq60.h)) squared its argument after reducing it only to
  [-pi, pi], so past 2.83 the square passed the format's ceiling; it now reduces to [-pi/2, pi/2] with a sign flip, and no program had
  called it past pi/2 before. And at 6:45 AM, solving the P0 kernel's own form exactly in its own code, the reviewer found a 64-bit sum of two
  neighbours in the kernel's point update that wrapped on the rows nearest the cut and left U high by about 3.5e-7 (1.4e-5 at a larger brane);
  fixed in the three kernels, confirmed to every printed decimal, the machines rerun it. The record of what was wrong is as public as the
  record of what is right.
- A result moves from "computed on a PC" to "measured" only when the two machines have run it and agreed; PC runs on this page are labelled
  as development runs.
- The code is written to be read: every file opens with what it computes, from which equation of which paper, in what arithmetic, with
  what self-tests, and with its own history of versions.

### The programs by computation

  - 1a, the Maldacena-Milekhin flow: `wh_flow.c` (Q.40) and `wh_flow60.c` (version 4, the state in 56 fraction bits on `whq60.h`: the
    flat flow's invariants to fourteen digits, the sphere flow's to twelve; a PC development run so far).
  - 1b, the patched five-dimensional geometry: the solver kernels built ahead of the equations, each verified on a problem with a known
    answer: `wh_relax.c` (one field: red-black SOR, a Robin boundary option, grids stretched toward a boundary with exact integer weights),
    `wh_relaxs.c` (several coupled nonlinear fields, a Newton step over the local Jacobian block, including the regime where a lagged
    source overshoots), `wh_relaxt.c` (the kernel in the geometry's own rows, uniform in t = ln(-y) from the cut to the deep end, carrying the
    reviewer's three sources in deviation variables around the AdS3 x S2 reference, whose discretization defect it measures and subtracts),
    `wh_p3.c` with `whq56.h` (the near-brane throat as an integer shot: three rescaled variables inward from the brane in a 56-bit format,
    the shot bisected on exact event tests), `wh_dtk.c` (the five-field kernel for the two-dimensional problems in the Einstein-DeTurck
    form: the variable-coefficient principal part with its cross term, per-point weights formed once, on a manufactured problem),
    `wh_p1data.c` (P1's boundary data in the cylinder frame: the paper's SO(3) example as a conformal class, blended and expanded to the
    Fefferman-Graham cutoff, checked against the reviewer's table), and
    `wh_ode.c` (the flow's exact integrator for a vector state); `wh_cost.py` (the paper's construction at every mouth size, the yardstick page);
    and `p2_nogo/`, the reviewer's scripts behind the first 1b result, the no-go for the paper's circle of mouths.
  - 2, the SYK teleportation protocol, seven versions kept as a lineage: `wh_syk.c` (the paper's layout), `wh_syk2.c` (the factorized
    layout that made N = 14 minutes instead of weeks), `wh_syk3.c` (the exact per-realization density matrix), `wh_syk4.c` (the same in
    Q.60 end to end, twelve digits against a float implementation), `wh_syk5.c` (the memory layout for N = 18, written to split every
    product in two through `par2.h`; so far the two halves run one after the other, with the same bits), `wh_syk6.c` (the size-winding observable the 2022 experiment reported, with the summary
    numbers of Kobrin, Schuster and Yao's comment, from the states version 5 already holds; `wind_summary.py` turns its printed per-size
    sums into their fit numbers), and `wh_syk7.c` (every matrix product through the reviewer's three-limb kernel `limb_q60.c` and
    `limb_tile.c` (build: `cc -O2 -std=gnu11 -o whsyk7 wh_syk7.c limb_q60.c limb_tile.c`): 24-bit limbs, int64 accumulators, the exact
    128-bit sums rebuilt, one rounding, bit-equal to the 128-bit kernel for every input inside its proved range (no INT64_MIN entry, at
    most 16,384 terms a sum, which every product here satisfies); 2.3 to 3.1 times fewer A9 instructions per complex multiply-accumulate
    in the reviewer's qemu counts on the PC, and an estimated four to five times fewer at N = 18).
  - 3a, the Alcubierre negative-energy integral and the wall-thickness law: `wh_alcub.c`.
  - 3b, the 2024 warp shell: `wh_shell.c` (the construction and its energy conditions) and `wh_warp.c` with `wh_prof.h` (the shift and the
    shift limit the paper left open; version 3 carries the reviewer's full fix).
  - 4a, the Casimir lattice sums: `wh_casimir.c`.

### Status (2026-11-12, 9:11 AM)

- **1a, the Maldacena-Milekhin flow** ([`results/1a_flow.md`](results/1a_flow.md)): their eq. (3.19) and (A.30) integrated exactly from the AdS3 x R2 throat
  to AdS5, their Figure 3 as a table; then, after review, a third version with one rounding per step, each flow on its own unstable
  manifold, and the gauge invariants, which agree with an independent 40-digit integration to 2.5e-11; then a fourth version in a
  56-bit format (a PC development run; the machines have not yet run it) that carries the flat flow's invariants to fourteen digits and
  prints a new gauge-invariant constant of the flat flow, Xi =
  0.121789029475802, which fixes the constant inside the paper's eq. (3.24) that it had estimated to one figure (kappa_flat = 1/4 - Xi =
  0.128210970524, verified by four independent routes and printed by both boards to twelve digits on 11-12); with the solved flow's logarithm (the
  constant plus a term from eq. 3.25 that the paper's estimate left out) the paper's derived numbers move by 4 to 9 percent, the constant alone by
  under 1 percent.
- **3a, the Alcubierre negative-energy integral** ([`results/3a_alcubierre.md`](results/3a_alcubierre.md)): the exact integral for seven wall thicknesses and the law
  behind "a factor of about three": I = R sigma / 3 + (pi^2 - 6)/(36 R sigma) plus exponentially small terms.
- **3b, the 2024 warp shell** ([`results/3b_shell.md`](results/3b_shell.md), [`results/3b_warp.md`](results/3b_warp.md)): the construction reproduced as integer sums, its four energy
  conditions judged exactly with a written floor under every verdict; then its shift added and the limit the paper left open computed:
  the conditions hold to a shift of about 0.03035 c and fail together just above it, 12.3 m from the centre on the equator (on the two boards,
  the two radial grids first printed 0.03033; a review found a truncated dr^2 in the kernel worth 3.25e-5; version 3 with the reviewer's full fix
  ran on the two boards the same night and printed 0.03036 on both grids, the PC's chains; the correction and the rerun are on the page).
- **4a, the Casimir sum** ([`results/4a_casimir.md`](results/4a_casimir.md)): one scalar field between plates as a finite lattice sum with no divergent subtraction,
  and the exact law by which the lattice approaches the continuum, 1 + (5 pi^2/24)/N^2 + 9.376/N^4 + ..., the ratio to ten decimals at
  N = 32 after two corrections from review.
- **2, the SYK teleportation protocol** ([`results/2_syk.md`](results/2_syk.md), [`results/2_syk_exact.md`](results/2_syk_exact.md)): the paper's averaged formulas in exact
  arithmetic at N = 8 for three seeds and the 2022 experiment's model; then the exact per-realization density matrix at N = 8 to 16 on
  the machines (18 on the PC, labelled), where the paper's averaged formulas are off by up to 14 percent at N = 8, 4 percent at N = 12
  and 1.6 percent at N = 14 (one seed), the wormhole sign wins over the wrong sign at every size
  by 1.4 to 1.6 at the peak, and the paper's own layout and ours agree to 7e-11 across 270 numbers at N = 10; then the size-winding signature (part 3,
  [`results/2_syk_winding.md`](results/2_syk_winding.md)), computed exactly from the same states and verified on the PC against the reviewer's independent float64 at
  N = 14 (the teleporting sign reaches 99.8 percent of the best possible line, the wrong sign 62 percent); the two machines ran it at
  N = 8 to 16 (N = 16 for one seed) and for the 2022 model between 6:56 and 7:32 AM on 11-12, with the PC build's chains
  ([`results/2_syk_winding_board1.txt`](results/2_syk_winding_board1.txt)).
- **1b, the patched five-dimensional geometry** ([`results/1b_nogo.md`](results/1b_nogo.md), [`results/1b_p0.md`](results/1b_p0.md), [`results/1b_p3.md`](results/1b_p3.md)): the reviewer derived the
  reduced Einstein-Maxwell system in verified stages, and it is being solved piece by piece. Done so far: the no-go for the paper's circle of
  mouths on its fine-tuned brane (a derivation with its checking scripts); the relaxation kernels, each verified on a problem with a known
  answer (PC builds; the machines ran `wh_relax.c` version 3 and `wh_relaxt.c`'s four test configurations with the PC's chains, raw
  output not yet in `results/`); the first solution of the reduced problem's simplest case, the one with no dependence along the
  circle (P0): the centre value U_c = -0.594116, -0.594579 and -0.594695 (the kernel as corrected at 7:45 AM) at N = 128, 256 and 512 rows, falling 4.00x per doubling toward
  -0.594734, the reviewer's independent extrapolation to six digits (after the review's correction of a wrapped sum, the
  two-grid Richardson is -0.594734016 against the continuum's -0.594734119), and a central charge 5.3e-5 above the fixed point's where the
  reviewer's law gives 5.26e-5; MEASURED on the two machines at N = 128 and 256 (26,011 and 30,608 sweeps each, identical lines, the
  PC's chains), N = 512 the PC's. And the near-brane throat (P3): the paper's log ln(L_IR/R5) as an integer shot at
  three charges, bit-identical to the reviewer's independent integer mimic in all three; at the paper's charge L_IR = 1.13679 l_B (its
  "L_IR ~ l_B" made exact) and the AdS2 radius exceeds the S2 radius at the brane by 1.77 percent; PC development runs; the machines
  have not yet run it (the three chains to match are on the page). The circle-direction problem (P1, the mouths) is next: its kernel
  in the Einstein-DeTurck form, five fields with a
  variable-coefficient principal part and a cross-derivative term, is built ([`code/wh_dtk.c`](code/wh_dtk.c)) and second order on all five fields of a
  manufactured problem (PC development runs), with its reflection lines, its Neumann deep end and its deviation form, ahead of the equations. And a reading of
  the paper's construction at every mouth size ([`results/1b_cost.md`](results/1b_cost.md), [`code/wh_cost.py`](code/wh_cost.py), computed on the PC from the paper's chain with our
  constant): the one condition that fixed its size was a human surviving the tides; without it the mouth can be much smaller (the table
  starts at a hundred times the fifth dimension, where the construction's approximations turn marginal; ESTIMATE), but every size from
  5 mm to 1 AU needs a dark magnetic charge of 5e31 to 2e45 Dirac units (for g4 = e), a dark field at the mouth with the energy density
  of 2e7 T or more, and a binding energy that grows as the mouth shrinks (minus 6.5e19 kg at 5 mm against minus 4.8e9 kg at the paper's
  size with our constant, where the paper estimated minus 5e9 kg; the negative Casimir energy that holds the throat open is twice the
  binding energy, the paper's 2.11 and 2.13); even the paper's wormhole takes about ten thousand years to cross as seen from outside
  (10,230 in the paper's leading-order chain with our constant). The table is the yardstick any proposal to shrink it has to beat.
