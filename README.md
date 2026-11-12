# Wormhole math, exact

**The computations that the traversable-wormhole and warp-drive papers left undone or in dispute, done in exact integer arithmetic on two
small ARM boards, so that anyone with a C compiler gets the same bits.** By Preston Ericson, working with Claude (Anthropic's AI) as
research partner. Preston directs the work, owns and runs the machines, and decides what is published; Claude writes the programs and
the pages with him; a second, independent AI session reviews every program and every page, before and after it goes up, and its
corrections are recorded here under its own name.

The first line, before anything else: **none of this builds a wormhole, a warp drive or a time machine.** Computing a metric does not
create it. The papers we check say the same about their own results. What this repository does is narrower and checkable: it takes
equations that were published with "one could imagine finding the solution numerically" or "this upper limit is a future direction of
work" or a disputed experimental signature, and computes them, exactly, in the open, with every constant written as an integer and
every result reproducible bit for bit on any machine.

## Check it yourself in two minutes

    cd code && cc -O2 -std=gnu11 -o whflow wh_flow.c && ./whflow | tail -n 4

The last lines print two chain hashes, `24d2f870a2abea59` for the flat flow and `85a3df31c75d8935` for the sphere flow. They are the hashes
of every step of the integration; the two ARM machines printed the same two, and so did the reviewer's own build of the same source on a
third machine. Every program here ends the same way: a chain hash and a `done` line, and every results page quotes the chain its numbers
came from. If your build prints a different hash, something is different, and we want to know what.

## Why exact

Every program here uses only 64-bit integers in a declared fixed-point format (40 or 60 fractional bits, with 128-bit intermediates built
from 32-bit halves; `code/whq.h`). No floating point, no math library. Exponentials come from a range reduction on ln 2 and a Taylor series
with a declared number of terms; square roots from integer Newton iteration or a digit-by-digit 128-bit root; differential equations from
a fixed-step Runge-Kutta scheme whose updates are rounded once with the remainder carried; matrix exponentials from powers accumulated in
128 bits and rounded once. The same source compiled anywhere prints the same digits, so a result is a fact about the equations, not about
one machine's rounding. We run every program on two separate ARM boards and require the two outputs to agree line for line before a result
is written down. Error bounds are derived on paper beside each program, and every number on a page is labelled: MEASURED (printed by a
program, with the run named), COMPUTED (derived, with the inputs named) or ESTIMATE.

## How it was built, and how it was checked

- Each program was sized on paper before its first run (the format, the step, where overflow could sit), and the first run of anything
  new goes to the two machines at once.
- An independent reviewer (another AI session working from copies, with its own derivations, its own float and 128-bit float builds and
  its own symbolic checks) reads every program and every page. In the first evening it found a division the shared header never guarded,
  a coefficient for the wrong field, a published Casimir table with digits lost past N = 16, a second one off by one unit in the ninth
  decimal, a flow started off its manifold, and the paper's density-matrix formulas being the ensemble-averaged forms rather than the
  per-realization ones. Every one of those is on its page, with the reviewer credited, the corrected version run again on the two
  machines, and the lineage of versions kept. In the second night it found a truncated dr^2 in the warp kernel worth 3.25e-5 in the
  published shift limit (the correction went up within a minute of its reading, before the fix; the fixed program then ran on the machines),
  the paper's density-matrix formulas' sign-winding counterpart derived for our layout, and a step of the 1b plan that was a matched
  approximation rather than a closed problem; the machines' firmware giving the process one core was found by the program's own line and is
  written where the two-core claim stood. The record of what was wrong is as public as the record of what is right.
- A result moves from "computed on a PC" to "measured" only when the two machines have run it and agreed; PC runs on this page are labelled
  as development runs.
- The code is written to be read: every file opens with what it computes, from which equation of which paper, in what arithmetic, with
  what self-tests, and with its own history of versions.

## What is here

- `PLAN.md`: the computations, what each paper left open, what we deliver, in what order, what "exact" means here, and the status of each.
- `THE_GLUE.md`: the reading that led here: what holds a traversable wormhole open in the post-2016 physics, why those wormholes
  cannot be time machines by construction, and where the energy scale meets the collider.
- `results/`: each result as a page with the numbers, the error bound on paper, the comparison with the paper it checks, the program's md5,
  the chain, and the two machines' agreement line, plus the raw output of one machine.
- `code/`: the programs, all on `whq.h`, the fixed-point header (Q.40 and Q.60, 128-bit intermediates). By computation:
  - 1a, the Maldacena-Milekhin flow: `wh_flow.c`.
  - 1b, the patched five-dimensional geometry: the solver kernels built ahead of the equations, each verified on a problem with a known
    answer: `wh_relax.c` (one field: red-black SOR, a Robin boundary option, grids stretched toward a boundary with exact integer weights),
    `wh_relaxs.c` (several coupled nonlinear fields, a Newton step over the local Jacobian block, including the regime where a lagged
    source stalls), and `wh_ode.c` (the flow's exact integrator for a vector state, so the near-brane throat is a right-hand side there);
    and `p2_nogo/`, the reviewer's scripts behind the first 1b result, the no-go for the paper's circle of mouths.
  - 2, the SYK teleportation protocol, seven versions kept as a lineage: `wh_syk.c` (the paper's layout), `wh_syk2.c` (the factorized
    layout that made N = 14 minutes instead of weeks), `wh_syk3.c` (the exact per-realization density matrix), `wh_syk4.c` (the same in
    Q.60 end to end, twelve digits against a float implementation), `wh_syk5.c` (the memory layout for N = 18, written to split every
    product across the machines' two cores through `par2.h`; on the first night the firmware gave the process one core, so those runs
    fell back to one core with the same bits, and the cause is being traced), `wh_syk6.c` (the size-winding observable of Kobrin, Schuster
    and Yao from the states version 5 already holds; `wind_summary.py` turns its printed per-size sums into their fit numbers), and
    `wh_syk7.c` (every matrix product through the reviewer's three-limb kernel `limb_q60.c` and `limb_tile.c`: 24-bit limbs, int64
    accumulators, the exact 128-bit sums rebuilt, one rounding, bit-equal to the 128-bit kernel for every input, four to five times fewer
    instructions on the machines).
  - 3a, the Alcubierre negative-energy integral and the wall-thickness law: `wh_alcub.c`.
  - 3b, the 2024 warp shell: `wh_shell.c` (the construction and its energy conditions) and `wh_warp.c` with `wh_prof.h` (the shift and the
    shift limit the paper left open; version 3 carries the reviewer's full fix).
  - 4a, the Casimir lattice sums: `wh_casimir.c`.

## Status (2026-11-11, one day)

- **1a, the Maldacena-Milekhin flow** (`results/1a_flow.md`): their eq. (3.19) and (A.30) integrated exactly from the AdS3 x R2 throat
  to AdS5, their Figure 3 as a table; then, after review, a third version with one rounding per step, each flow on its own unstable
  manifold, and the gauge invariants, which agree with an independent 40-digit integration to 2.5e-11.
- **3a, the Alcubierre negative-energy integral** (`results/3a_alcubierre.md`): the exact integral for seven wall thicknesses and the law
  behind "a factor of about three": I = R sigma / 3 + (pi^2 - 6)/(36 R sigma) plus exponentially small terms.
- **3b, the 2024 warp shell** (`results/3b_shell.md`, `results/3b_warp.md`): the construction reproduced as integer sums, its four energy
  conditions judged exactly with a written floor under every verdict; then its shift added and the limit the paper left open computed:
  the conditions hold to a shift of about 0.0304 c and fail together just above it, 12.3 m from the centre on the equator (the pair's
  two radial grids printed 0.03033; a review found a truncated dr^2 in the kernel worth 3.25e-5, so the corrected value is 0.03036, and
  version 3 reruns on the pair tonight; the correction is on the page).
- **4a, the Casimir sum** (`results/4a_casimir.md`): one scalar field between plates as a finite lattice sum with no divergent subtraction,
  and the exact law by which the lattice approaches the continuum, 1 + (5 pi^2/24)/N^2 + 9.376/N^4 + ..., the ratio to ten decimals at
  N = 32 after two corrections from review.
- **2, the SYK teleportation protocol** (`results/2_syk.md`, `results/2_syk_exact.md`): exact at N = 8 for three seeds and the 2022
  experiment's model; then the exact per-realization density matrix at N = 8 to 16 on the machines (18 on the PC, labelled), where
  the paper's averaged formulas are off by 14 percent at N = 8 and 1.6 at N = 14, the wormhole sign wins over the wrong sign at every size
  by 1.4 to 1.6 at the peak, and the paper's own layout and ours agree to 7e-11 across 270 numbers at N = 10; then the size-winding signature (part 3,
  `results/2_syk_winding.md`), computed exactly from the same states and verified on the PC against the reviewer's independent float64 at
  N = 14 (the teleporting sign reaches 99.8 percent of the best possible line, the wrong sign 62 percent), on the machines tonight.
- **1b, the patched five-dimensional geometry**: started. The reviewer derives the reduced Einstein-Maxwell system in verified stages (its
  first stage finds that only the rotationally symmetric reductions are well-posed two-dimensional problems, and that the paper's fine-tuned
  circle of mouths has no static solution without both extra positive energy and a positive cosmological constant on the brane, or
  matter with rho + 3p below zero: `results/1b_nogo.md`, the first result of 1b, a derivation with its checking scripts); the relaxation kernel that will solve it is built and
  verified on a problem with a known answer, now with a Robin boundary option and grids stretched toward a boundary, both exact.

Next: the pair's N = 18 and its size-winding runs, 1b's equations into the kernel. See `PLAN.md`.

## Licenses

Code: MIT (`code/LICENSE`). Text and tables: CC BY 4.0 (`LICENSE-TEXT`). Please cite as: Preston Ericson, "Wormhole math, exact",
2026, this repository.
