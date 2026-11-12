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
  `wh_shell.c` reproduces the 2024 warp shell's matter shell and judges its energy conditions exactly (computation 3b, part 1).
- `results/`: each result as a plain-text table with the program's md5, the format, the step, and the two boards' agreement line.

## Status

Two results are in `results/`, each run on two separate ARM machines at once with byte-identical output (2026-11-11, 4:26-4:28 PM):

- **1a, the Maldacena-Milekhin flow** (`results/1a_flow.md`): their eq. (3.19) and (A.30) integrated exactly from the AdS3 x R2
  throat to AdS5; their Figure 3 as a table to twelve places; the sphere flow's exact IR radius 1.0643 against their expansion's 1.0625.
- **3a, the Alcubierre negative-energy integral** (`results/3a_alcubierre.md`): the exact integral for seven wall thicknesses, the
  flattening theorem checked in two dimensions, and the law behind "a factor of about three": I = R sigma / 3 + (pi^2 - 6)/(36 R sigma)
  plus exponentially small terms, with the constant identified to ten digits.

Next: 3b (the 2024 warp shell's energy conditions, exactly, with the shift limit its authors left open), then 2 (the SYK protocol), then
4 (the Casimir sums), then 1b (the patched five-dimensional geometry). See `PLAN.md`.

## Licenses

Code: MIT (`code/LICENSE`). Text and tables: CC BY 4.0 (`LICENSE-TEXT`). Please cite as: Preston Ericson, "Wormhole math, exact",
2026, this repository.
