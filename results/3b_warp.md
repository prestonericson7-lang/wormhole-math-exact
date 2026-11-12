# Result 3b, part 2: the 2024 warp shell with its shift, judged exactly, and the shift limit its authors left open, 2026-11-11

**First line:** a shell of positive-energy matter with a shift of three percent of the speed of light inside it is a solution of
Einstein's equations on paper and a 2.4-Jupiter-mass object nobody can build; nothing here changes that. This result computes, exactly,
how far the shift of the 2024 construction can be raised before the first energy condition breaks, a number the paper calls "a future
direction of work", and finds it just above 0.03.

## What was computed

Fuchs, Helmerich, Bobrick, Sellers, Melcher and Martire ("Constant velocity physical warp drive solution", Class. Quantum Grav. 41,
095013, 2024, arXiv:2405.02709, section 4) add to the shell of part 1 a shift vector along one axis, g_tx = -beta0 S(r), with S the
compact sigmoid of their eq. (27)-(28) smoothed twice, and state that all four energy conditions hold at beta0 = 0.02 (0.02 c), "likely
not an upper limit". The upper limit is computed here.

- The metric, with the shift direction as the polar axis: g_tt = -e^{2a}, g_tr = -beta cos theta, g_t theta = beta r sin theta,
  g_rr = 1/(1 - 2m/r), g_theta theta = r^2, g_phi phi = r^2 sin^2 theta, beta = beta0 S(r); the profiles m, a, S are the declared
  interpolants of the construction's tables (part 1), S from the exact triangle-kernel sums of its two smoothing passes.
- The Einstein tensor in the warped-product form g = h + f^2 dphi^2, f = r sin theta: the three-dimensional Ricci tensor of the
  (t, r, theta) block by the generic Christoffel code, the Hessian of f on it, and R_ab = R[h]_ab - (nabla_a nabla_b f)/f,
  R_phiphi = -f Laplacian f. Nothing divides by sin^4 theta, which an earlier draft did and which overflowed the fixed point next to
  the axis (caught in review before any run).
- The four conditions in the Eulerian frame over 100 Fibonacci directions x 10 speeds in the coordinate triad, refined by a local
  descent; the dominant energy condition with a 128-bit norm.
- A written floor under every verdict: F = 3 x 2^-40 [0.1/r^4 + 10 + 2/(r sin theta)] + 3 x 2^-41 / (8 pi dr), the curvature code's
  rounding envelope (measured in review against 50-digit evaluations) and the tables' rounding seen through a''. A margin is VIOLATED
  only below -F. Cells where the declared matter is exactly zero are reported apart as residuals of the declared object.
- Self-tests, all printed: beta0 = 0 against the shell's closed form (agreement to 1-2 units in the twelfth place); flat space with a
  constant shift (every frame component zero to 2e-12); the warped-product code against the earlier four-index code off the axis.
- Program: `code/wh_warp.c` with `code/wh_prof.h` on `code/whq.h` (v2); binary md5 48e33083 (ARM build).
- Runs: the pair, both machines at once, at r every 400 table midpoints (0.2 m = 0.01 R2) and 31 theta cells, 4278 cells, 8:56 to
  9:21 PM on 11-11 (25 minutes), the two outputs identical line for line, chain cd537a1c30e7ddc9 (`results/3b_warp_board1.txt`); and
  the pair again on the fine grid, r every 100 midpoints (0.05 m) and 31 cells, 17,112 cells, 10:22 PM to 11:58 PM (96 minutes), identical
  line for line and equal to the PC build's chain, ce224896a73d0e49 (`results/3b_warp_fine_board1.txt`); a coarse PC run (0.2 m, 15 cells)
  before both.

## The numbers

On the pair (MEASURED; r every 0.2 m, 31 theta cells): all four conditions hold at beta0 = 0, 0.02, 0.025 and 0.03 within the floor, with
no matter cell below it (the worst margin at 0.03 is -2.7e-10 against a floor of 3.6e-9, at r = 0.14 R2 next to the axis); at 0.035 all
four fail together in 61 cells, worst -0.006872 at r = 0.6176 R2 = 12.35 m on the equator. The bisection over the 154 candidate cells (16
halvings) puts the first failure between

    beta0 = 0.030330352784 and 0.030330429078   at r = 0.6176 R2 = 12.35 m, theta = pi/2,

the dominant energy condition first by its |Y| term, the other three within 1e-6 of it. The 527 cells where the declared matter is
exactly zero show residuals of up to 1.05e-8 against their own floor of 5e-9 to 6e-9: the declared object's own rounding, reported apart,
not a violation.

On the fine grid (MEASURED on the pair; r every 0.05 m, 31 cells): the limiting cell moves to r = 0.61506 R2 = 12.30 m and the bracket to

    beta0 = 0.030325622559 to 0.030325698853   at r = 0.61506 R2 = 12.30 m, theta = pi/2,

570 candidate cells of 17,112; chain ce224896a73d0e49, the PC build's exactly. So the radial grid moves the limit by 4.8e-6, downward, as
a finer grid can only do, and the fine grid's value is the one to quote: 0.030326. An independent review with its own derivation of the curvature and its own reduced grids
found 0.03035 at 12.32 m for its own copy of the object (its total mass taken as exactly 1/6 where the tables' trapezoid sum gives
0.16666023, and its own smoothing); the review is now computing the acceptance for this object, and that difference of 2.4e-5 is the
object, not the grid.

## The light rays along the axis

The round-trip difference 4 INT beta e^{-2a} dr between a ray with the shift and one against it: 6.839932 ns in coordinate time at
beta0 = 0.02 (R2 = 20 m) on the pair's grid (6.839926 on the coarse grid: the integral's own grid dependence, in the sixth digit),
6.365734 ns in the proper time of an emitter at the grid's edge; both are delays, and both scale linearly in beta0 to 34.2 and 31.8 ns
at beta0 = 0.1. The paper's 7.6 ns at 0.02 c is
reproduced by neither definition (its Table 1 is captioned 0.04 c; the same definition does reproduce the paper's Alcubierre row of
8.0 ns): an open discrepancy, stated as such.

## Error bound, on paper

The floor F is the bound: at the limiting cell (r = 0.62 R2, equator) it is 2.2e-9 in the margin, and the margin's slope in beta0 there
is about 1.5 (from -0.00687 at 0.035 to zero at 0.03033), so the floor contributes about 1.5e-9 to beta0; the bisection's bracket is
7.6e-8; the radial grid is the largest term: from 0.2 m to 0.05 m spacing the limiting cell moved by 0.0025 R2 and the limit by 4.8e-6.
The two grids' limits, 0.030330 and 0.030326, differ by the grid's own term, 4.8e-6; the fine grid's 0.0303256 is good to its fifth
decimal, and a still finer grid would lower it by less than that. The two machines' bit-equality checks the arithmetic, not this bound.

## What it says and does not say

The construction of the 2024 paper holds its four energy conditions up to a shift of about 0.0303 c and fails them, all together, just
above it, at 12.3 m from the centre on the equator, for a null observer moving along the shift. That is the paper's open question
answered with the paper's own object, exactly. It does not make the object buildable, and it says nothing about acceleration, which the
paper does not have either.
