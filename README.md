# Wormhole math, exact

**Three computations the wormhole and warp-drive papers left undone or in dispute, done in exact integer arithmetic so that anyone
with a C compiler gets the same bits.** By Preston Ericson, working with Claude (Anthropic's AI) as research partner.

The first line, before anything else: **none of this builds a wormhole, a warp drive or a time machine.** Computing a metric does not
create it. The papers we check say the same about their own results. What this repository does is narrower and checkable: it takes
equations that were published with "one could imagine finding the solution numerically" or "this upper limit is a future direction of
work" or a disputed experimental signature, and computes them, exactly, in the open, with every constant written as an integer and
every result reproducible bit for bit on any machine.

## Why exact

Every program here uses only 64-bit integers in a declared fixed-point format (for example Q.28: 28 fractional bits). No floating
point, no math library. Exponentials come from a range reduction on ln 2 and a Taylor series with a declared number of terms; square
roots from integer Newton iteration; differential equations from a fixed-step Runge-Kutta scheme. The same source compiled anywhere
prints the same digits, so a result is a fact about the equations, not about one machine's rounding. We run every program on two
separate ARM boards and require the two outputs to agree line for line before a result is written down. Error bounds are derived on
paper beside each program.

## What is here

- `PLAN.md`: the three computations, what each paper left open, what we deliver, in what order, and what "exact" means here.
- `THE_GLUE.md`: the reading that led here: what holds a traversable wormhole open in the post-2016 physics, why those wormholes
  cannot be time machines by construction, and where the energy scale meets the collider.
- `code/`: the programs. `wh_flow.c` solves the Maldacena-Milekhin flow equations (computation 1a).
- `results/`: each result as a plain-text table with the program's md5, the format, the step, and the two boards' agreement line.

## Status

Computation 1a (the Maldacena-Milekhin flow, arXiv:2008.06618 eq. 3.19 and A.30) is built and waiting for its first run. Nothing in
`results/` yet. Computation 3 (the energy conditions of the published warp metrics) has its exact recipe written in `PLAN.md`.

## Licenses

Code: MIT (`code/LICENSE`). Text and tables: CC BY 4.0 (`LICENSE-TEXT`). Please cite as: Preston Ericson, "Wormhole math, exact",
2026, this repository.
