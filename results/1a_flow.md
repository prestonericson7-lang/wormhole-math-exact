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

- Program: `code/wh_flow.c` on `code/whq.h`; binary md5 e72cb881 (ARM build, -O2) for this first run; the file in `code/` is now version 3
  (its section below), whose first-run-equivalent lines differ from these by the gauge and the sphere's mode, as that section explains.
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

## Version 3, 2026-11-11 10:05 PM: the review's three points answered

The reviewer wrote the specification and a reference program; version 3 of `code/wh_flow.c` is that specification, four changes and none to
the arithmetic header:

1. **One rounding per update.** Each variable's RK4 increment h X / 6 (X = k1 + 2 k2 + 2 k3 + k4, an exact int64) is one exact int64
   floor division by 6 x 2^12, and the remainder is carried to the next step in an int64 per variable. After n steps the state is exactly
   the starting value plus the running sum of the increments rounded once (proved by induction, checked by the reviewer with a 128-bit
   accumulator at every one of the 74,674 + 74,023 steps): the accumulated rounding error is in [-1/2, 1/2) unit at every step and never
   grows. The UV slopes now settle at exactly 1.000000000000: the dead zone of the first run is gone.
2. **Each flow starts on its own unstable manifold** at rho = -8 with the linear coordinate x = 2^-10, from the manifold's series to
   order 4 (the first order whose Q.40 integers are final), lambda in the canonical gauge; the start integers are hard-coded (COMPUTED by
   the reviewer in 60-digit arithmetic and re-derived by its verifier). The sphere flow's own growing exponent is 2.830105157, not the flat
   flow's 2.627.
3. **Stop at sigma > 16** (was 12), so the sphere flow's sigma - rho has converged too.
4. **The gauge invariants are printed:** K_sigma = (sigma - rho)_end + GAMMA_sigma and K_lambda = (lambda - rho)_end + GAMMA_lambda, the UV
   constants carried to the gauge in which the IR mode is e^{k rho} with unit amplitude; GAMMA are the reviewer's exact gauge constants
   (COMPUTED; the same in eight gauges to 2e-15). sigma - rho alone depends on where rho = 0 is put; K does not.

**Run (MEASURED):** the pair, both machines at once, 10:17:38 to 10:18:09 PM on 11-11 (31 seconds), 160 output lines each, identical line
for line; chains 24d2f870a2abea59 (flat) and 85a3df31c75d8935 (sphere), equal to the PC build's and to the reviewer's reference program's,
so three builds of the arithmetic on two architectures agree bit for bit. The raw output of one machine: `results/1a_flow_v3_board1.txt`
(`1a_flow_board1.txt` stays the first run's). Binary md5 5c7a26b9 (ARM build).

| | flat flow | sphere flow |
|---|---|---|
| start at rho = -8 (sigma, sigma', lambda) | 0.000975399492, 0.002559176054, -13.857071005204 | -0.061363013201, 0.002756975087, -14.496374133760 |
| steps to sigma > 16 | 74,674 (rho = 10.2310) | 74,023 (rho = 10.0720) |
| UV slopes at the end, sigma' and lambda' | 1.000000000000, 1.000000000000 | 1.000000000000, 1.000000000000 |
| sigma - rho, lambda - rho, lambda - sigma at the end | 5.769097209892, -4.196441803900, -9.965539013793 | 5.928170409940, -4.804059361524, -10.732229771464 |
| K_sigma, K_lambda (this program) | 0.407800114165, -0.271699935356 | 0.377362726326, -0.297030895771 |
| K_sigma, K_lambda (the reviewer, 40 digits) | 0.40780011419036, -0.27169993536963 | 0.37736272634195, -0.29703089577475 |
| this program minus the 40-digit value | -2.5e-11, +1.4e-11 | -1.6e-11, +0.4e-11 |

For a direct comparison with the paper's Figure 3: the reviewer read the figure's gauge off the plot and finds this repository's rho equal
to the figure's rho minus 5.3619 (COMPUTED), so the first run's table can be laid over the figure with that shift.

The differences from the first run's table are the gauge and the mode, not errors: the flat flow's sigma - rho moves by -5.9e-4 because
the first run's start was off the manifold at second order by that much in rho, and the sphere's constants move by +0.0073 (sigma - rho),
-0.0067 (lambda - rho) and -0.0139 (lambda - sigma) because it now starts on its own mode. The invariants are the numbers to quote: they
agree with an independent 40-digit integration to 2.5e-11, which is this arithmetic's whole error (the reviewer's error table: the v3 rule
leaves at most 3.7e-11 in the UV constants and the error no longer depends on the step).

## Lineage

A first build (Q.28, floor rounding, amplitude 2^-20) was withdrawn before any run: its first increments in the deep IR were a sixth of
a unit and floored to zero, so it would have printed a flat table. The first attempt to run this build failed because the target's
file store was full; the run script now clears it first and refuses to compare empty outputs. Both are in the commit history.

## Version 4 (11-12, 3:10 AM): the flow in Q.56, the invariants to twelve digits

The same two flows and the same RK4 (`code/wh_flow60.c` on `code/whq60.h`), with the state in a fixed-point format of 56 fraction bits instead
of 40, the exponentials and roots from the Q.60 header rounded once, and every update rounded once with its remainder carried. From version
3's start integers (exact to 2^-40; the start's error is a translation of the flow, which the invariants ignore) the PC build printed
the invariants to eleven and twelve digits (3:10 AM); with the reviewer's exact manifold point at rho = -8 as 56-bit integers and its gauge
constants in closed form (Gamma_sigma = -8 + 10 ln 2/k, Gamma_lambda = (1 - sqrt3) Gamma_sigma) it prints (MEASURED, 3:16 AM) K_sigma =
0.407800114190352 and K_lambda = -0.271699935369623 for the flat flow, against the reviewer's 30-digit integration 0.407800114190356327 and
-0.271699935369626228: fourteen digits, where version 3 had ten; the sphere flow, still from its version-3 start, prints 0.377362726342 and
-0.297030895775, its twelve. The same program integrates a fourth variable J along the flat flow and prints the gauge-invariant constant
Xi = lim [sigma - e^{2(sigma - lambda)} J] = 0.121789029475802, the reviewer's 0.121789029475803078 to fifteen digits (`code/wh_ode.c` holds the
Q.40 version, which reaches 0.121789029377, the format's floor). What Xi means (COMPUTED by the reviewer; verified 3:37 AM by an adversarial rebuild
with its own code, and by four routes in all: the linearized brane matching, a shooting of the closed mirror-brane throat at fourteen
flux values extrapolating to within 4e-9, an independent agent's integration, and this program): on the paper's fine-tuned brane the
throat obeys re^2/Q2 = ln(re/S0) + 1/4 - Xi with S0 = (Q2 L^2/6)^(1/4), so the constant inside their eq. (3.24), which the paper estimated
to one figure, is exact: kappa_flat = 1/4 - Xi = 0.128210970524, and L_IR = e^(1/4 - Xi) l_B = 1.1368 l_B. With it, the paper's eq. (3.27)
numbers move by 4 to 9 percent, each staying within its one-figure estimate (the paper's value in brackets; COMPUTED by the reviewer, with
the paper's four-dimensional energetics used as printed and re read as the throat's S2 radius on the brane; if the paper's re is meant as
the AdS2 radius instead, c2 is 1.7 percent lower): ln(L_IR/R5) = 14.457 (13.21), c2 = 6.709e72 (7.02e72), ell = 3,257 light years (3,113),
gamma = 2.054e12 (1.964e12), E_bin = -4.79e9 kg (-5.24e9), l_B = 83.58 m (81.72), and the throat's AdS2 radius exceeds its S2 radius by
1.775 percent (the paper's own 1.947). Only c2 and the logarithm are the solution's; ell, gamma and E_bin rest on the paper's balance
equations as before. The machines run version 4 tonight.
