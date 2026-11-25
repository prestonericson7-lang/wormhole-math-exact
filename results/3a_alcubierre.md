# Result 3a: the Alcubierre negative-energy integral in integer fixed-point arithmetic, bit-identical on two machines, and the law behind "a factor of about three", 2026-11-11

**In plain words.** The Alcubierre "warp drive" is a proposed shape of space and time (a metric: the description of how space and time
are shaped) in which a bubble of space of radius R moves at a speed v, and the bubble's wall needs negative energy (energy below that of
empty space). A 2021 paper by Bobrick and Martire showed that a differently shaped wall needs less of it, "by about a factor of 3". This
page computes the amount for Alcubierre's own wall at several thicknesses in fixed-point integer arithmetic (every calculation done with
whole numbers in a set format, so any computer prints the same digits; the two machines, two small ARM computers, printed identical
output) and finds the simple law behind that factor: it is sR/3 (s the wall's steepness, R the bubble's radius) plus a small correction,
so, with the optimal wall's size set equal to the bubble's radius, the "factor of about 3" is its value at sR = 9 (3.012), Alcubierre's
own plotted wall (sR = 8) gives 2.68 (both MEASURED), and the factor is not a constant: thinner walls give proportionally more. Nothing
here makes a warp drive possible: the negative energy is still required, and for a 100 m bubble at the speed of light it is about 16
Jupiter masses of negative mass with Alcubierre's plotted wall and about 6 with the optimal one (COMPUTED).

**First line:** the Alcubierre metric needs negative energy density everywhere its wall is; nothing here changes that. This result
computes how much for the published wall shapes, in integer arithmetic printed to ten decimals (seven guaranteed by the error bound
below), and finds the closed-form law that the optimal-wall comparison follows.

## What was computed

Bobrick and Martire ("Introducing physical warp drives", Class. Quantum Grav. 38, 105009, 2021, arXiv:2102.06824) give, in G = c = 1
units, the Eulerian energy density of the Alcubierre drive, their eq. (3),

    T^00 = -(1/8 pi) (rho^2 v^2 / 4 r^2) f'(r)^2,

and its total, their eq. (A.7), E = -(v^2/12) INT_{r0}^inf r^2 f'(r)^2 dr (the optimal wall's f' vanishes below r0; for the tanh
wall this page integrates from 0). They show by variation that the wall f = min(r0/r, 1)
minimizes |E| (giving INT = r0 exactly) and write that it "decreases the needed (absolute value of) negative energy of the Alcubierre
solution by about a factor of 3" compared with Alcubierre's tanh wall, f = [tanh(s(r+R)) - tanh(s(r-R))] / (2 tanh(sR)). They also
show that flattening the bubble along its motion by a factor a divides E by a.

- Program: `code/wh_alcub.c` on `code/whq.h`; binary md5 11944bed (ARM build, -O2).
- Arithmetic: Q.40 fixed point, 128-bit intermediates, round to nearest; Simpson's rule with 2^14 intervals on [0, 4] (h = 2^-12) in
  one dimension and 2^10 x 2^10 in two, every step an exact power of two; sums in 128-bit accumulators. Units R = 1, v = 1.
- Run: 4:26-4:28 PM, both machines at once; 18 output lines each; **identical line for line**, chain 69b34aea7baff8af. The raw
  output: `3a_alcubierre_board1.txt` (board 2's lines are byte-identical).

## The numbers

I(sR) = INT_0^4 r^2 f'^2 dr for the tanh wall (I_opt = r0 = 1 for the optimal wall, so the second column is also the ratio):

| sR | I(sR) = ratio to the optimal wall | sR/3 | sR (I - sR/3) |
|---|---|---|---|
| 4 | 1.362028350131 | 1.333333333333 | 0.114780067193 |
| 6 | 2.017964426844 | 2.000000000000 | 0.107786561067 |
| 8 | 2.680103999487 | 2.666666666666 | 0.107498662568 |
| 9 | 3.011943406954 | 3.000000000000 | 0.107490662592 |
| 10 | 3.344082262017 | 3.333333333333 | 0.107489286847 |
| 12 | 4.008957418211 | 4.000000000000 | 0.107489018533 |
| 16 | 5.340051396534 | 5.333333333333 | 0.107489011221 |

**The law.** The last column converges to a constant, and the constant is identifiable: write the sech^2 bump at r = R in the variable
u = s (r - R); then INT r^2 f'^2 dr = (s/4) INT (1 + u/s)^2 sech^4(u) du + (exponentially small in s). The three moments are exact:
INT sech^4 = 4/3, INT u sech^4 = 0, INT u^2 sech^4 = (pi^2 - 6)/9. Hence

    I(sR) = sR/3 + (pi^2 - 6)/(36 sR) + O(sR e^{-2 sR}),      (pi^2 - 6)/36 = 0.107489011141...
          = coth^2(sR) [sR/3 + (pi^2 - 6)/(36 sR)] + O((sR)^2 e^{-4 sR})   (the e^{-2 sR} part exactly)

The boards' column gives 0.107489011221 at sR = 16 (the difference, 8 x 10^-11, is almost all arithmetic: 7.5e-11, against a
remainder of 4.3e-12) and already 0.1074890185 at sR = 12. With the optimal wall's r0 set equal to R (the paper leaves r0 free;
half-height matching, r0 = R/2, would move the factor 3 to sR = 4.42), the "factor of about 3" is the value of this law at sR = 9
(3.012); at Alcubierre's own plotted wall, sR = 8, it is 2.68; thin walls give proportionally more. The ratio is not a constant of the problem but R sigma / 3
plus a small universal correction.

**Checks printed by the program.**
- The optimal wall's own Simpson integral over [1, 4]: 0.749999999999 against the exact 1 - 1/4 = 0.75 (the tail to infinity is 1/4
  exactly, so I_opt = 1).
- Flattening at sR = 8: J(a) = INT_0^{4/a} dx INT_0^4 rho (df/drho)^2 drho for a = 1, 2, 8 gives a J(a) = 1.786735999639,
  1.786735999639, 1.786735999638: on the power-of-two grid the three 2D sums are one and the same integer (a x runs over the same
  points for every a), so this checks the code, and the last digit at a = 8 is the rounding of the final scaling; the theorem itself is a
  change of variable. J(1) against the one-dimensional value (2/3) I(8) = 1.786735999658: the 1.9 x 10^-11 gap is mostly the 2D rule's
  end term at rho = 0 (COMPUTED by the reviewer: -2.1 x 10^-11).
- Self-test: tanh 1 = 0.761594155956 (exact 0.761594155956), sech^2 1 = 0.419974341613 (0.419974341614), sqrt 2 = 1.414213562372
  (1.414213562373; square roots print the floor).

**Floating-point cross-check (PC, mpmath, not the result):** I(4) = 1.3620283501297, I(8) = 2.6801039994811, I(16) = 5.3400513965299.
The exact integer results agree with these to 1.3 x 10^-12, 6 x 10^-12 and 4 x 10^-12.

## Error bound (on paper)

Simpson truncation is at most (4/180) h^4 max|f''''| with h = 2^-12; the integrand's fourth derivative scales as s^4 x 1024 x 16, so the
bound is 8 x 10^-8 absolute at sR = 16 and 3 x 10^-10 at sR = 4. Rounding: each term carries about 20 roundings of 2^-40 on values up
to 1024, so at most 2 x 10^-8 per term, and the Simpson weights integrate to the range 4, so at most 8 x 10^-8 on the integral. Total
below 1.6 x 10^-7 absolute, 3 x 10^-8 relative at sR = 16: seven digits of I are the ODE-free integral's, and the cross-check above
shows the actual error is four orders smaller than the bound.

## In SI units (COMPUTED from the table; c = 299792458 m/s exact, G = 6.67430 x 10^-11 m^3 kg^-1 s^-2, CODATA 2018)

E = -(v^2/12) I R c^4/G with c^4/G = 1.2103 x 10^44 J per metre (the run's SI line printed 1.2104e44, a rounding slip, since corrected in the program). For a bubble of R = 100 m at v = c (the papers' usual scale):
- tanh wall with sR = 8: I R = 268.01 m, E = -2.703 x 10^45 J, equivalent to -3.01 x 10^28 kg, about 15.8 Jupiter masses of
  negative mass;
- the optimal wall with r0 = 100 m: I R = 100 m, E = -1.009 x 10^45 J, -1.12 x 10^28 kg, about 5.9 Jupiter masses.
Flattening by a divides both by a. None of these is a quantity anyone can make; they are the exact prices of the published metrics.

## Against the paper

Their statements hold: the optimal wall reduces |E| by a factor near three for walls of the thickness they had in mind, and flattening
divides E by the flattening factor exactly. What this adds is the exact law for the factor, R sigma / 3 + (pi^2 - 6)/(36 R sigma), and
the numbers to ten places (the largest error, at sR = 12, is 1.1e-11) that anyone can regenerate from `code/wh_alcub.c` on any machine.
