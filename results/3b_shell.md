# Result 3b, part 1: the 2024 "warp shell" without its shift, reproduced in integers and judged exactly, 2026-11-11

**First line:** a shell of positive-energy matter around a flat cavity is not a warp drive and moves nothing; this result reproduces the
2024 construction as published, declares precisely what object it is, and reports its energy conditions with the arithmetic and its floor
stated. Part 2 (the shift that makes it "move") is the next program, and is under review before its first run.

## What was computed

Fuchs, Helmerich, Bobrick, Sellers, Melcher and Martire, "Constant velocity physical warp drive solution", Class. Quantum Grav. 41,
095013 (2024), arXiv:2405.02709, section 3, build a spherical shell of constant density between R1 = 10 m and R2 = 20 m with the total
mass M = R2 c^2 / (6 G) (2.4 Jupiter masses), a TOV pressure for constant density set to zero inside R1, MATLAB's `smooth` applied four
times to the density (span 1.79 x 4000) and to the pressure (span 4000) on a grid of 100,000 points out to rmax = 1.2 x sqrt(29.9^2 +
29.9^2 + 0.4^2) m, the mass function by cumulative trapezoid, and the lapse by cumulative trapezoid of (m + 4 pi r^3 P)/(r^2 - 2 m r)
pinned to (1/2) ln(1 - 2m/r) at the grid's end. Their public code (WarpFactory: `metricGet_WarpShellComoving.m`, `TOVconstDensity.m`,
`alphaNumericSolver.m`) fixes the numbers the paper leaves out, and does the same as the paper at the one place we first thought it
differed: `TOVconstDensity.m` multiplies by the step density array, which is zero inside R1, so the cavity pressure is zero in the code
too. (An earlier version of this program's comments and of the plan called the pressure-kept case "the code's variant"; that was a
misreading, caught in review on 2026-11-11 before this page was written. The pressure-kept case is computed anyway and labelled
HYPOTHETICAL: it is neither the paper nor its code.)

What "exactly" means for a construction built from smoothing and trapezoids: every one of their steps is a finite sum, so it is done as
one in Q.40 (2^-40 = 9.1e-13) with 128-bit intermediates and one rounding per stored value; then the metric is DECLARED to be the
piecewise interpolant of the resulting tables (m and a piecewise quadratic and C^1 between grid points, since their integrands are taken
linear there), and the Einstein tensor of THAT object is evaluated in closed form at every interval's midpoint (the static spherically
symmetric formulas), where its first and second derivatives are both defined. Nothing is differentiated numerically.

- Units: R2 = 1 (20 m), G = c = 1; m_total = 1/6; density 1/(7 pi) = 0.045472840883; energies and pressures in c^4 / (G R2^2) =
  3.0256e41 J/m^3 (the paper's ~1e39 J/m^3 is ~3.3e-3 in these units).
- Program: `code/wh_shell.c` with `code/wh_prof.h` on `code/whq.h`; binary md5 8a874b79 (ARM build).
- Arithmetic: Q.40, 128-bit intermediates, round to nearest; the smoothing means and the trapezoids in 128-bit accumulators, one rounding
  per stored value; the TOV pressure from the closed form with sqrt exact to the floor and the exact mass 1/6 (their code uses the
  trapezoid of the unsmoothed step density for that mass, about 4e-5 relative away from 1/6: a known deviation, to be matched in v2).
- Run: 7:36 PM, both machines at once; 148 output lines each; **identical line for line**, chain 375780b722176169. The raw output:
  `3b_shell_board1.txt` (board 2's lines are byte-identical).

## The numbers

**The construction, exactly.** Four smoothing passes of span 7159 (3.63 m) spread the shell's matter over r = 0.137442 .. 1.362548, that
is from 2.75 m to 27.25 m: the "shell between 10 and 20 m" is, after the smoothing the paper applies, a shell between 2.7 and 27 m. The
smoothed total mass is 0.169802392148 against the nominal 1/6 = 0.166666666666: the smoothing raises the mass by 1.88 percent, because a
symmetric average on a 4 pi r^2 measure moves matter outward on net. The smoothed pressure is nonzero over r = 0.297796 .. 1.199429 and
exactly zero at the centre; a(0) = -0.272157027890, a(rmax) = -0.071848331002 = (1/2) ln(1 - 2m/r) there.

**The energy conditions, as published** (the minimum of each margin over each region; in units of 3.0256e41 J/m^3):

| region | rho + P_r (NEC) | rho + P_t (NEC) | rho (WEC) | rho + P_r + 2 P_t (SEC) | rho - abs(P_r) (DEC) | rho - abs(P_t) (DEC) |
|---|---|---|---|---|---|---|
| cavity r < 0.5 | 0 | 0 | 0 | -1.3e-11 at r = 0.139 | 0 | -1.6e-8 at r = 0.148 |
| shell 0.5..1 | 0.022996 at r = 1.0 | 0.022626 at r = 1.0 | 0.022736 at r = 1.0 | 0.022774 at r = 1.0 | 0.020736 at r = 0.5 | 0.014277 at r = 0.5 |
| exterior r > 1 | 0 | -1.2e-9 at r = 2.39 | 0 | -2.4e-9 at r = 2.39 | -1e-12 at r = 1.36 | -1.2e-9 at r = 2.43 |

Inside the matter every condition holds with a margin of 0.014 to 0.023, which is 4 to 7 x 10^39 J/m^3. In the cavity and the exterior
the margins are zero to within the arithmetic's floor: the exterior residual (the maximum of |P_r| and |P_t| for r > 1.52, where the
metric must be exactly Schwarzschild) is 1.22e-9, and the cavity's -1.6e-8 sits at the inner edge of the smoothed density, where the
second derivative of the lapse is a difference of table values over dr = 2.5e-5 and amplifies the stored rounding by 4e4 (bound below).
In SI those floors are 3.7e32 and 5e33 J/m^3, inside the 1e30-1e35 J/m^3 band the paper itself calls numerical noise in its Figure A1. So
the verdict the paper claims for its unshifted shell stands exactly: all four energy conditions hold, and where the program prints
VIOLATED outside the matter it is printing its own floor, not the shell's physics. The raw lines carry those VIOLATED tags; this page is
where they are read.

**The hypothetical variant** (the TOV pressure kept inside R1, which the paper rules out and its code never does): P(0) = 0.005757, and the
cavity then has pressure without density, which fails the dominant energy condition genuinely: rho - |P_r| and rho - |P_t| are about -0.0057
across the cavity (the table rows at r = 0.05 and 0.10 read rho = 0, P_r = 0.00573, P_t = 0.00571). The printed cavity minima of that
variant (-0.0369, -0.0073 at r = 0.0002) are the 1/r terms of the formulas at the first midpoints and overstate it; the honest number is
-P(0). That is presumably why the paper zeroes the cavity pressure, and the exact computation agrees with the choice.

## Error bound, on paper

Each stored table value carries one rounding (+-2^-41); each smoothing pass is a mean of up to 7159 such values and adds at most one
more; four passes: below 1e-11 absolute on rho~ and P~. The trapezoids sum up to 10^5 terms in 128-bit: below 1e-7 absolute worst case on
m and a, 1e-9 as a random walk. The Einstein-tensor formulas use a'' = (da[i+1] - da[i])/dr with dr = 2.54e-5, which turns a 2^-41
rounding of da into 3.6e-8 in a'', and divide by r^2 and (r^2 - 2 m r): the floor on the printed rho_T, P_r, P_t is therefore about
4e-8 near r = 0.14 and 1e-9 or better at r > 1. Every "VIOLATED" minimum outside the shell is inside that floor; every margin inside the
shell is a million times above it. The two machines' bit-equality checks the arithmetic, not this bound.

## What it says and does not say

A positive-energy shell at rest is ordinary general relativity, and this one obeys all four energy conditions once its construction is
taken as the declared object. Two exact facts the paper does not state come with it: the smoothing raises the mass by 1.9 percent and
moves the matter to 2.7-27 m. The paper's claim is about the shifted version; part 2 computes it with the same declared-object method,
after its review (a division guard in the shared arithmetic, a rounding floor on the sign tests, and the cells next to the axis) is done.
