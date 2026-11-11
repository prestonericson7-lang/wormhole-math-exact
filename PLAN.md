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
  PROGRAM: `code/wh_flow.c`. Q.28 fixed point in int64; exp by range reduction on ln 2 and 24 Taylor terms; sqrt by integer Newton;
  RK4 with step h = 2^-12; both flows started on the IR growing mode at rho = -8 with amplitude 2^-20; tables every 0.25 in rho; the
  UV asymptotes; a chain hash of every step. Checked on paper before the first run: the UV fixed point has sigma' = 1 and lambda' =
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
not a constant 3: their "about a factor of three" is the value at R sigma = 9. DELIVERABLE: the integral in Q.28 for R sigma = 4, 6,
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

## Status (2026-11-11)

- 1a: `code/wh_flow.c` built; first run pending; its table and the two machines' agreement line will appear in `results/`.
- 3: recipe written (above); program not yet started.
- 2, 1b: planned.
