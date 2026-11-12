# Result 1a: the Maldacena-Milekhin flow, solved exactly on two machines, 2026-11-11

**First line:** this computes a geometry on paper; it does not create one. The flow below is the shape of the fifth dimension in the
"humanly traversable wormhole" of Maldacena and Milekhin (arXiv:2008.06618), a configuration they themselves call allowed by the
equations with no plausible mechanism for its formation.

## What was computed

Their eq. (3.19) with l_B = 1 (Appendix A, eq. A.28-A.29), the flat flow:

    0 = -6 + 3 e^{-4 sigma} + lambda'^2 + 4 lambda' sigma' + sigma'^2
    0 =  4 - 4 e^{-4 sigma} - 2 lambda' sigma' - 2 sigma'^2 - sigma''

and their eq. (A.30), the same flow on a two-sphere of unit radius (the extra -e^{-2 sigma} and +e^{-2 sigma} terms). Both integrated
from deep in the IR (AdS3 x R2, where sigma is at its fixed point and lambda' = sqrt 3) to the UV (AdS5, where sigma and lambda grow
like rho), starting on the growing linearized mode, as the paper prescribes.

- Program: `code/wh_flow.c` on `code/whq.h`; binary md5 e72cb881 (ARM build, -O2).
- Arithmetic: Q.40 fixed point in 64-bit integers with 128-bit intermediates, products and divisions rounded to nearest, exp by
  range reduction on ln 2 with 24 Taylor terms, square roots exact to the floor. RK4 with step h = 2^-12 in rho.
- Gauge: the linear mode has amplitude 2^-10 at rho = -8 (sigma - sigma* = 2^-10 e^{k (rho + 8)}, k = sqrt 19 - sqrt 3 =
  2.626848135972). Any other amplitude gives the same curve shifted in rho.
- Run: 4:25-4:26 PM, both machines at once; 124 output lines each; **identical line for line**, both chain hashes equal
  (flat 5642434347e87c07, sphere 96c768410f49dff7). The raw output: `1a_flow_board1.txt` (board 2's lines are byte-identical).

## The numbers

| quantity | value (12 places as printed) | note |
|---|---|---|
| flat flow, IR fixed point | sigma* = 0, lambda' = 1.732050807568 | sqrt 3, as the paper says |
| flat flow, UV | sigma - rho -> 5.769688646998, lambda - rho -> -4.196210084188 | constant to all printed places from rho = 1 on |
| flat flow, UV | lambda - sigma -> -9.965898731186 | the UV constant of the flat flow, this gauge |
| flat flow, UV slopes at the end | sigma' = 0.999999999612, lambda' = 1.000000000386 | both -> 1; the residual 4e-10 is the rounding floor (bound below) |
| sphere flow, IR fixed point | e^{-sigma*} = 1.064322422264, sigma* = -0.062338373459 | the exact root of 4 + u - 4 u^2 = 0, u = e^{-2 sigma*} = (1 + sqrt 65)/8 |
| sphere flow, IR fixed point | lambda' = 1.811959037795, sigma'' = 0.000000000021 | the residual 2e-11 is the Newton solve for sigma* and rounding |
| their (A.31) at l_B = 1 | e^{-sigma} = l + l^3/16 + ... = 1.0625 | the two-term expansion is low by 0.17% at l_B = 1; the expansion is for small l_B |
| sphere flow, UV | sigma - rho -> 5.920889751429, lambda - rho -> -4.797406486407 | constant to all printed places from rho = 4.5 on |
| sphere flow, UV | lambda - sigma -> -10.718296237836 | the UV constant of the sphere flow, this gauge |
| steps | flat 58,288 (rho from -8 to 6.23), sphere 57,669 (to 6.08) | stopped where sigma passed 12 (e^{-4 sigma} < 2^-69) |

The flows, every 1.0 in rho (the full tables at every 0.25 are in the raw file):

    flat:   rho | sigma | sigma' | lambda' | lambda - rho | sigma - rho
    -8 | 0.000976562500 | 0.002565281382 | 1.730298961471 | -5.856406460545 | 8.000976562500
    -7 | 0.013308192824 | 0.034396029157 | 1.708585620423 | -5.132763687648 | 7.013308192824
    -6 | 0.153955295356 | 0.336522256317 | 1.499316034586 | -4.497334135803 | 6.153955295356
    -5 | 0.803697426936 | 0.897518719022 | 1.085262583789 | -4.225389707028 | 5.803697426936
    -4 | 1.770940793482 | 0.995624003966 | 1.003962441724 | -4.197357997238 | 5.770940793482
    -3 | 2.769723160749 | 0.999873517183 | 1.000118770466 | -4.196242671191 | 5.769723160749
    -2 | 3.769689490030 | 0.999996835225 | 1.000003023409 | -4.196210893291 | 5.769689490030
    -1 | 4.769688665538 | 0.999999926576 | 1.000000070833 | -4.196210102871 | 5.769688665538
     0 | 5.769688646998 | 0.999999998436 | 1.000000001515 | -4.196210084188 | 5.769688646998
     1 | 6.769688646998 | 0.999999999611 | 1.000000000386 | -4.196210084188 | 5.769688646998

    sphere: rho | sigma | sigma' | lambda' | lambda - rho | sigma - rho
    -7.5 | -0.058418749705 | 0.011030379028 | 1.804028003005 | -6.091774045492 | 7.441581250294
    -6.5 | -0.000670600086 | 0.161554925635 | 1.694706809966 | -5.321556955528 | 6.499329399913
    -5.5 |  0.466761145177 | 0.808054708558 | 1.192499518193 | -4.858916440632 | 5.966761145177
    -4.5 |  1.415846203668 | 1.003229395149 | 1.004838019178 | -4.796818952046 | 5.915846203668
    -3.5 |  2.419705376134 | 1.002151880223 | 0.999135784817 | -4.796872978696 | 5.919705376134
    -2.5 |  3.420714949903 | 1.000343966245 | 0.999833572381 | -4.797320582382 | 5.920714949903
    -1.5 |  4.420865730791 | 1.000047905774 | 0.999976179426 | -4.797394514176 | 5.920865730791
    -0.5 |  5.420886491619 | 1.000006515214 | 0.999996745416 | -4.797404859891 | 5.920886491619
     0.5 |  6.420889309737 | 1.000000882384 | 0.999999558913 | -4.797406267675 | 5.920889309737
     1.5 |  7.420889691511 | 1.000000119367 | 0.999999940355 | -4.797406458039 | 5.920889691511

What the table shows, in words: the transition from the AdS3 x R2 throat (lambda' = sqrt 3, sigma flat) to AdS5 (both slopes 1)
takes about two units of rho and is centred near rho = -5 in this gauge. In the flat flow sigma' rises to 1 from below and stays
below; in the sphere flow sigma' overshoots 1 by 0.3% near rho = -4.5 and lambda' undershoots by 0.09% near -3.5 before both settle.
Each UV approach is exponential: sigma' - 1 shrinks by a factor of about 45 per unit of rho (e^{-4 sigma}, with sigma = rho + const),
until it reaches the arithmetic's floor of a few 10^-10. This is their Figure 3 as digits.

## Self-test lines (printed by the program, in the run log)

    sqrt3 = 1.732050807568 (1.732050807569), sqrt19 = 4.358898943540 (4.358898943541), k = 2.626848135972 (2.626848135972),
    e^ln2 = 2.000000000000, e^-ln2 = 0.500000000000, 1/3 = 0.333333333333, 3*5 = 15.000000000000, sqrt4 = 2.000000000000,
    e^0 = 1.000000000000, sqrt3^2 = 2.999999999997 (3 - at most 2^-38), e^-48 = 0.000000000000

The square roots print the floor at 2^-40, so their last place can read one below the rounded value; that is the declared behaviour.

## Error bound (on paper)

- RK4 truncation on the growing mode sigma ~ eps e^{k rho}: about (k h)^5/120 = (2.627 x 2^-12)^5/120 = 9 x 10^-19 relative per step,
  4 x 10^-15 per unit of rho, 5 x 10^-14 over the whole flow. Negligible.
- Rounding: each RK4 update rounds once to nearest, +-0.5 unit = +-4.5 x 10^-13 on sigma, sigma' and lambda per step; over about
  58,000 steps the worst case (every rounding the same way) is 1.4 x 10^-8 absolute, the random-walk expectation about 10^-10. The
  observed UV residuals (4 x 10^-10 in the slopes) sit at that expectation.
- So the printed values are solutions of the ODE to about 8 decimal places (10^-9 relative in the UV where sigma ~ 12; in the deep IR
  where sigma ~ 10^-3 the worst case is a gauge shift of rho by at most 5 x 10^-6). The digits past the 8th belong to the declared
  arithmetic. Two machines agreeing bit for bit checks the arithmetic, not this bound; the bound is the paper above.

## Against the paper

Their Appendix A solves the flat flow numerically and plots sigma and lambda (Figure 3); the digits above are that figure as a table,
with the IR and UV asymptotes they state (lambda' = sqrt 3, then sigma' = lambda' = 1) reproduced. Their (A.31) gives the sphere flow's
IR radius as an expansion in l_B; at l_B = 1 the exact fixed point is 1.0643 against the expansion's 1.0625. Nothing here conflicts
with the paper; it fills in what they left as "one could imagine trying to find it numerically in a more explicit form" for the flow
itself. The patched five-dimensional geometry (their section 3.2) remains the open computation 1b.

## Review note, 2026-11-11 evening (an independent review of this program; the numbers above are left as printed)

A second, independent reviewer (another AI session working from copies, with its own derivations and runs) raised three points on this
result, recorded here before the program is rerun:

1. **Rounding bias.** The RK4 update rounds its increments with `rdiv`, which rounds ties away from zero; over the flow's ~60,000 steps
   that is a one-sided bias of about 1/12 unit per step, amplified in the deep IR, and it puts a 4-7e-8 bias on the UV constants. The
   bound on this page (1.4e-8 worst case, 1e-10 random walk) assumed symmetric rounding and is therefore too tight: the UV constants are
   good to about 7 decimals, not 12, and the error grows as the step shrinks. A remainder-carrying update (the reviewer's estimate: 2.5e-11)
   is the fix for v3.
2. **The sphere flow's starting mode.** The sphere flow was started on the flat flow's growing mode; its own linearized exponent is 2.8301,
   not 2.627. The sphere flow's UV constants shift by 0.006-0.014 when it is started on its own mode. The flat flow is unaffected.
3. **What to publish.** The gauge-dependent combination sigma - rho should be replaced by the gauge invariants (the reviewer's values for
   the flat flow: K_sigma 0.40780011419034, K_lambda -0.27169993536961, from its own run; ours follow with v3).

None of this changes the two machines' bit-equality, which checks the arithmetic, not the bound; it changes how many of the printed digits
are the equation's. v3 of `wh_flow.c` carries the remainder, starts the sphere flow on its own mode and prints the invariants.

## Lineage

A first build (Q.28, floor rounding, amplitude 2^-20) was withdrawn before any run: its first increments in the deep IR were a sixth of
a unit and floored to zero, so it would have printed a flat table. The first attempt to run this build failed because the target's
file store was full; the run script now clears it first and refuses to compare empty outputs. Both are in the commit history.
