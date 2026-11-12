# Wormhole math, exact

**Three computations the wormhole and warp-drive papers left undone or in dispute, done in exact integer arithmetic so that anyone
with a C compiler gets the same bits.** By Preston Ericson, working with Claude (Anthropic's AI) as research partner.

The first line, before anything else: **none of this builds a wormhole, a warp drive or a time machine.** Computing a metric does not
create it. The papers we check say the same about their own results. What this repository does is narrower and checkable: it takes
equations that were published with "one could imagine finding the solution numerically" or "this upper limit is a future direction of
work" or a disputed experimental signature, and computes them, exactly, in the open, with every constant written as an integer and
every result reproducible bit for bit on any machine.

## Why exact

Every program here uses only 64-bit integers in a declared fixed-point format (Q.40: 40 fractional bits, with 128-bit intermediates
built from 32-bit halves; `code/whq.h`). No floating point, no math library. Exponentials come from a range reduction on ln 2 and a Taylor series with a declared number of terms; square
roots from integer Newton iteration; differential equations from a fixed-step Runge-Kutta scheme. The same source compiled anywhere
prints the same digits, so a result is a fact about the equations, not about one machine's rounding. We run every program on two
separate ARM boards and require the two outputs to agree line for line before a result is written down. Error bounds are derived on
paper beside each program.

## What is here

- `PLAN.md`: the three computations, what each paper left open, what we deliver, in what order, and what "exact" means here.
- `THE_GLUE.md`: the reading that led here: what holds a traversable wormhole open in the post-2016 physics, why those wormholes
  cannot be time machines by construction, and where the energy scale meets the collider.
- `code/`: the programs. `wh_flow.c` solves the Maldacena-Milekhin flow equations (computation 1a); `wh_alcub.c` computes the
  Alcubierre negative-energy integral, the wall-thickness law behind the "factor of three", and the flattening check (computation 3a);
  `wh_shell.c` reproduces the 2024 warp shell's matter shell and judges its energy conditions exactly (computation 3b, part 1);
  `wh_warp.c` with `wh_prof.h` adds the shift and finds the shift limit (3b, part 2); `wh_casimir.c` computes the Casimir lattice sums
  (4a); `wh_syk.c`, `wh_syk2.c` and `wh_syk3.c` run the SYK teleportation protocol (2) in the paper's layout, in the factorized layout
  that makes N = 14-16 a matter of minutes, and with the exact per-seed density matrix. All of them sit on `whq.h`, the fixed-point header.
- `results/`: each result as a plain-text table with the program's md5, the format, the step, and the two boards' agreement line.

## Status

Six results are in `results/`, each run on two separate ARM machines at once with byte-identical output (2026-11-11, 4:26 PM to 9:43 PM):

- **1a, the Maldacena-Milekhin flow** (`results/1a_flow.md`): their eq. (3.19) and (A.30) integrated exactly from the AdS3 x R2
  throat to AdS5; their Figure 3 as a table to twelve places; the sphere flow's exact IR radius 1.0643 against their expansion's 1.0625.
- **3a, the Alcubierre negative-energy integral** (`results/3a_alcubierre.md`): the exact integral for seven wall thicknesses, the
  flattening theorem checked in two dimensions, and the law behind "a factor of about three": I = R sigma / 3 + (pi^2 - 6)/(36 R sigma)
  plus exponentially small terms, with the constant identified to ten digits.
- **3b, the 2024 warp shell** (`results/3b_shell.md`, `results/3b_warp.md`): the construction reproduced as integer sums and its four
  energy conditions judged exactly with a written floor under every verdict; then its shift added and the limit the paper left open
  computed: the conditions hold to a shift of 0.03033 c and fail together just above it, at 12.35 m on the equator.
- **4a, the Casimir sum** (`results/4a_casimir.md`): one scalar field between plates as a finite lattice sum, no divergent subtraction,
  and the exact law by which the lattice approaches the continuum, 1 + (5 pi^2/24)/N^2 + O(1/N^4), checked to N = 32.
- **2, the SYK teleportation protocol** (`results/2_syk.md`): exact at N = 8 for three seeds, both signs of the coupling and mu = 0 as the
  exact zero, and the 2022 experiment's five-term model; the density-matrix formulas used are the paper's ensemble-averaged forms
  evaluated per seed (labelled so after review); the exact per-seed matrix is version 3 of the program, on the machines tonight.

Next: the SYK protocol's exact per-seed density matrix at N = 8 to 16 (version 3, beside the averaged forms), the warp shell's limit on
the fine radial grid on the pair, then 1b (the patched five-dimensional geometry). See `PLAN.md`.

## Licenses

Code: MIT (`code/LICENSE`). Text and tables: CC BY 4.0 (`LICENSE-TEXT`). Please cite as: Preston Ericson, "Wormhole math, exact",
2026, this repository.
