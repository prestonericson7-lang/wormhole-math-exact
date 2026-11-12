# The plan: three computations, exact

Preston Ericson, 2026-11-11. The first line, again: the programs here compute metrics and quantum protocols; computing a metric does
not create it. Each of the three computations below is a published equation or claim that its authors left uncomputed, left as
"future work", or that is in dispute, and each result is a number others can check with the code in `code/`.

## What "exact" means here

Fixed-point integer arithmetic in a declared format, deterministic rounding, every constant written as an integer, so that two
machines or two runs give the same bits and anyone with the code gets the same bits. Error bounds derived on paper beside the code.
Floating point is used only to cross-check small cases on a PC, never as the result. Every number is labelled: COMPUTED (which
program, which md5, which format and step) or ESTIMATE (and what it would take to compute).

## Where it runs

The programs are plain C using only 64-bit integers: no floating point, no math library. We run each one on two separate ARM
Cortex-A9 boards at the same time and require the two outputs to agree line for line before a result is written down; the same
source compiled on any other machine prints the same digits. That is the whole point of the method: the result is a property of the
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
  PROGRAM: `code/wh_flow.c` (v2). Q.40 fixed point in int64 with 128-bit intermediates, products and divisions rounded to nearest; exp by
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
checked as the protocol paper defines it; then the SAME computation on the learned sparse model of the 2022 paper (its couplings are in
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
8, 9, 10, 12, 16 by a fixed-point trapezoid on a declared grid, the ratio tabulated, the asymptote checked; the flattening theorem
E -> E/alpha_X (their 4.2) checked on a 2D (x, rho) grid for alpha_X = 2 and 10 as a test of the arithmetic. Their "two orders of
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
6 x 10^9 multiply-adds per condition set, minutes per shift value in int64 on two Cortex-A9 cores, which split theta and must agree
line for line as in 1a. Floating-point cross-check of the smallest grid on a PC only.
HONEST LINE: a positive-energy shell of 2.4 Jupiter masses in a 20 m radius, moving at a constant 0.02 c with no way to accelerate it,
is not a drive; the paper says so itself ("the question of how to make physical and efficient acceleration is one of the foremost
problems"). What we add is exactness, numbers they left open, and a public check.

## Order

1a first (it tests the arithmetic on real equations in days), then 3 (two papers, mechanical, high value: a published claim checked or
not), then 2 (the biggest, a live dispute), then 1b (the research project). Each result is a number others can check, written up with
its format, its bounds and its seed. None of them builds a door, and each write-up says so in its first line.

## Status (2026-11-11, updated in the evening)

- 1a: DONE 2026-11-11 4:26 PM, `results/1a_flow.md`: both machines identical, 124 lines; the flow from the throat to AdS5 to 12 places.
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
  0.030358123781 to 0.030358200075 (0.05 m), COMPUTED by the reviewer; version 3 of the program reruns on the pair tonight. The review confirmed the general relativity (the Ricci convention, the
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
  digits verified against a float implementation, the time steps built once per seed) runs on the pair at N = 8-16 and the 2022 model
  tonight; v4.1 derives the right-hand operators exactly and v5 (`code/wh_syk5.c`) fits N = 18 in the machines' memory (a PC development run
  at N = 18 is on the page, labelled). Next: the pair's N = 18 (N = 16 ran on the pair at 12:08 AM, 16.4 min for one seed, the PC's digits), the size-winding signature (derived in review; version 6 computes it, verified on the PC at N = 14 against the reviewer's float64, on the machines tonight; `results/2_syk_winding.md`), the N-dependence of the sign
  asymmetry written up.
- 1a: a review note added to `results/1a_flow.md` (a rounding bias of 4-7e-8 on the UV constants from ties-away rounding, the sphere
  flow's starting mode, the gauge invariants to publish). v3 of `wh_flow.c` BUILT 10:05 PM from the reviewer's specification and
  reference program: one rounding per update with the remainder carried, each flow on its own unstable manifold from hard-coded start
  integers, stop at sigma > 16, the gauge invariants printed; the PC build's 160 lines are byte-identical to the reference's and its two
  chains equal; K_sigma, K_lambda agree with an independent 40-digit integration to 2.5e-11 (the arithmetic's whole error now). The pair
  runs it in tonight's queue; its section is on the 1a page.
- 4a again: the reviewer certified the exact values (three routes to 25 digits) and found v4 still +1.7e-9 / +2.5e-9 at N = 32 / 24 from
  three roundings; v5 (10:10 PM) removes them and lands within 8.6e-12 of the certified ratio on the PC build; the law's second
  coefficient is 9.3762 (closed form on the page), not the fitted 9.6. The pair runs v5 tonight; the page's table is replaced then.
- 1b: STARTED 2026-11-11 10:50 PM. The reviewer derives the reduced system in two verified stages (the ansatz and gauge, the 5D
  Einstein-Maxwell equations as a 2D elliptic system with the constraints separate, the residuals, the boundary conditions edge by edge;
  then the (3.27) parameters as formulas in the solved functions, the sizing, the 1a checks). The fixed-point relaxation kernel is built
  and verified ahead of the equations: `code/wh_relax.c`, a Q.60 red-black SOR with 128-bit residuals, tested on Poisson with an exact
  quartic solution (for which the five-point Laplacian is exact): the residual reaches its rounding floor, a few G^2 units of 2^-60, in
  about 6 G sweeps as the spectral radius predicts, and the solution is within 4e-17, 6e-16 and 8e-15 of exact at G = 32, 128 and 512
  (PC build); the pair runs it tonight for bit-equality. Stage 1 (11-12, 1:18 AM, verified in nearly all of it by the reviewer's second builder): only
  the rotationally symmetric reductions are well-posed two-dimensional problems (the paper's flat-space pair has axial symmetry only and needs
  three coordinates and a strut); the paper's fine-tuned circle of mouths has no static solution without extra positive energy on the brane;
  the order chosen: the
  solver test, then the near-brane throat (an ODE matched to the flow of 1a; the flow has no inner end, so this is an approximation, not
  a closed boundary-value problem), then the paper's own rotationally symmetric example with antipodal mouths, which it called findable
  numerically; the no-go for the fine-tuned circle holds where the mouths' spacing is large against re/(2 ln(L_IR/R5)). The kernel is at version 4: a Robin boundary option eliminated into the first interior row
  and grids stretched toward a boundary with exact integer stencil weights, both landing on the known solution to tens of units of 2^-60;
  an adaptive relaxation parameter was tried and dropped (it misreads SOR's early transient). The pair ran version 3 at G = 128 and 256
  and the nonlinear G = 128 (12:28-12:32 AM, chains equal to the PC's) and runs version 4 tonight. The system kernel (`code/wh_relaxs.c`,
  2:00 AM) relaxes three coupled fields with nonlinear sources, the Newton step over the local 3 x 3 Jacobian: on the exact test it lands
  within 3 to 172 units of 2^-60 at G = 32 to 128 in about 8 G sweeps; where the source's own derivative is a third of the Laplacian's
  diagonal, the Newton step (block or diagonal) still lands exactly while the lagged plain update stalls at a false floor with a 2e-2 error.
  The pair runs both tonight.
