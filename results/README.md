# Results

**In plain words.** Each page below is one result: what was computed, the numbers, the error bound (how far the printed digits can be
from the exact answer) worked out on paper and the comparison with the paper it checks. A run on the two machines (two small ARM
computer boards) comes with the raw output of one of them (`*_board1.txt`); the other machine's output was identical line for line, as
each page states. Every number carries a label: MEASURED (printed identically by the two machines, or a published value, named),
COMPUTED (derived, or printed by a PC development run) or ESTIMATE. PC runs are development runs and are labelled so on their pages
(the one raw PC file is named `_pc_`); [`1b_nogo.md`](1b_nogo.md) is a derivation and ran on no machine. Nothing here builds
anything: these are checkable numbers about published equations.

- **The shape of the extra dimension in Maldacena and Milekhin's "humanly traversable" wormhole, solved step by step in integer
  arithmetic.**
  [`1a_flow.md`](1a_flow.md),
  [`1a_flow_board1.txt`](1a_flow_board1.txt) (the first run, the program's version 2), [`1a_flow_v3_board1.txt`](1a_flow_v3_board1.txt)
  (version 3): the Maldacena-Milekhin flow (arXiv:2008.06618 eq. 3.19 and A.30), 2026-11-11; version 3 prints the gauge invariants
  K_sigma and K_lambda, which agree with an independent 40-digit integration to 2.5e-11; version 4 (Q.56, a PC development run so far)
  prints the constant Xi.
- **What the same wormhole would cost at other sizes, from 5 mm to 1 AU, the distance from the Earth to the Sun (a formula
  evaluation, not a new solution).**
  [`1b_cost.md`](1b_cost.md): the paper's wormhole construction (arXiv:2008.06618) evaluated at mouth sizes from 5 mm to 1 AU with our
  exact constant in place of its estimate, and what each size costs (charge, field, negative energy, crossing time); COMPUTED on a PC
  from the paper's formula chain by [`code/wh_cost.py`](../code/wh_cost.py), 2026-11-12.
- **A derivation showing that the paper's circle of mouths (the loop its magnetic field lines make, through the wormhole and back
  outside it) has no static (unchanging) solution.** [`1b_nogo.md`](1b_nogo.md): the paper's
  circle of mouths has no static solution on its fine-tuned brane; the reviewer's derivation, with its scripts in
  [`code/p2_nogo/`](../code/p2_nogo/), 2026-11-12.
- **The first solution of the explicit wormhole geometry (its full shape, which the paper did not write down), in its simplest
  case.** [`1b_p0.md`](1b_p0.md),
  [`1b_p0_board1.txt`](1b_p0_board1.txt) (N = 128), [`1b_p0_n256_board1.txt`](1b_p0_n256_board1.txt) (N = 256): the first solution of
  the reduced 1b problem, its case with no dependence along the circle (P0), on the two machines at N = 128 and 256; N = 512 and the
  alpha = 1/8 check are PC development runs, 2026-11-12.
- **The wormhole's throat (the tube between its two mouths) next to the brane (the edge of the extra dimension, where ordinary matter
  lives), computed in integers.**
  [`1b_p3.md`](1b_p3.md), [`1b_p3_pc_6_2pow28.txt`](1b_p3_pc_6_2pow28.txt): the near-brane throat (P3) as an integer shot; PC
  development runs so far, 2026-11-12.
- **The wormhole-teleportation protocol (a recipe that inserts a quantum bit on one side of two linked systems of particles and reads
  it out on the other): its formulas, which average over the random strengths of the particles' interactions (the couplings),
  evaluated in exact arithmetic.** [`2_syk.md`](2_syk.md),
  [`2_syk_board1.txt`](2_syk_board1.txt): the SYK traversable-wormhole teleportation protocol's averaged formulas in exact arithmetic at
  N = 8 (three seeds, both signs of the coupling, mu = 0 as the exact zero) and the 2022 experiment's five-term model at its own
  parameters (arXiv:1911.07416, 2302.07897), 2026-11-11; the formulas are the paper's ensemble-averaged forms evaluated per seed.
- **The same protocol computed exactly for one draw of its random couplings; the averaged formulas are off by up to 14 percent
  (COMPUTED from the machines' MEASURED values).**
  [`2_syk_exact.md`](2_syk_exact.md), [`2_syk_exact_board1.txt`](2_syk_exact_board1.txt) (version 3),
  [`2_syk_exact_v4_board1.txt`](2_syk_exact_v4_board1.txt) (version 4): the exact per-realization density matrix, N = 8 to 16 on the
  machines (N = 18 a PC development run), 2026-11-11 and 11-12.
- **The "size winding" pattern the 2022 experiment reported (a pattern in how an inserted signal spreads over the particles, expected
  when the system behaves like a wormhole), computed exactly from the same states.**
  [`2_syk_winding.md`](2_syk_winding.md), [`2_syk_winding_board1.txt`](2_syk_winding_board1.txt): the size-winding signature, on the two
  machines at N = 8 to 16 (N = 16 for one seed) and for the 2022 model, 2026-11-12.
- **The negative energy the classic (Alcubierre) warp drive needs, in integer arithmetic, and the law behind a published "factor of
  about three" (how much a better wall shape lowers that energy).**
  [`3a_alcubierre.md`](3a_alcubierre.md), [`3a_alcubierre_board1.txt`](3a_alcubierre_board1.txt): the Alcubierre negative-energy
  integral and the law behind the optimal-wall factor (arXiv:2102.06824 eq. 3, A.7), 2026-11-11.
- **The 2024 warp shell without its shift (the part that makes it "move"), rebuilt in integers, with its four energy conditions
  judged and the rounding floor (the size of rounding error each verdict must clear) stated.** [`3b_shell.md`](3b_shell.md),
  [`3b_shell_board1.txt`](3b_shell_board1.txt): the 2024 warp shell without its shift, the construction reproduced as integer sums, the
  energy conditions judged with their floor (arXiv:2405.02709 section 3), 2026-11-11.
- **The 2024 warp shell's speed limit, which its paper left open: about 0.03035 c (COMPUTED; the machines' finest grid puts it no
  higher than 0.030358 c, MEASURED).** [`3b_warp.md`](3b_warp.md),
  [`3b_warp_board1.txt`](3b_warp_board1.txt) and [`3b_warp_fine_board1.txt`](3b_warp_fine_board1.txt) (version 2),
  [`3b_warp_v3_board1.txt`](3b_warp_v3_board1.txt) (version 3): the 2024 warp shell with its shift and the limit the paper left open
  (arXiv:2405.02709 section 4): about 0.03035 c, 12.32 m from the centre on the equator; version 3's 0.05 m grid bounds it from above at
  0.030358; version 2's value was 3.25e-5 low from a truncated dr^2 found in review, 2026-11-11 and 11-12.
- **The Casimir energy (the energy of empty space between two close surfaces, which can be below zero) as a finite sum, with nothing
  infinite subtracted.**
  [`4a_casimir.md`](4a_casimir.md), [`4a_casimir_board1.txt`](4a_casimir_board1.txt): the Casimir energy of one scalar field as a finite
  lattice sum, no divergent subtraction, and the lattice law 1 + (5 pi^2/24)/N^2 + 9.376/N^4 + ..., the ratio to ten decimals at N = 20
  to 32 (version 5 of the program; version 2 had lost digits past N = 16 and version 4 one unit in the ninth decimal, both caught in
  review), 2026-11-11.
