# Wormhole math, exact

[![Watch the video](https://img.shields.io/badge/Watch-the_video-FF0000?logo=youtube&logoColor=white)](https://www.youtube.com/watch?v=4lfyN1AOj1o) [![Watch the companion video](https://img.shields.io/badge/Watch-the_companion_video-FF0000?logo=youtube&logoColor=white)](https://www.youtube.com/watch?v=wmKaPpUYgok) [![YouTube channel](https://img.shields.io/badge/YouTube-%40LittlebrainsBigmess-FF0000?logo=youtube&logoColor=white)](https://www.youtube.com/@LittlebrainsBigmess) [![ORCID](https://img.shields.io/badge/ORCID-0009--0009--9123--2731-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0009-9123-2731)

**Related:** [the-winter-that-wasnt](https://github.com/prestonericson7-lang/the-winter-that-wasnt) · [monroe-stargate-record](https://github.com/prestonericson7-lang/monroe-stargate-record) · [all our work](https://github.com/prestonericson7-lang)

---

**Computations that the traversable-wormhole and warp-drive papers left undone or in dispute, done in exact integer arithmetic and checked
on two small ARM boards, so that anyone with a C compiler gets the same bits.** The largest, the explicit wormhole geometry (1b), is under
way; results not yet run on the boards are labelled as PC development runs.

By Preston Ericson, working with Mira (Claude, Anthropic's AI) as research partner. Preston directs the work, owns and runs the machines, and
decides what is published; Claude writes the programs and the pages with him; a second, independent AI session reviews every program and
every page, before and after it goes up, and its corrections are recorded here under its own name.

The first line, before anything else: **none of this builds a wormhole, a warp drive or a time machine.** Computing a metric (the
description of how space and time are shaped) does not create it. The papers we check say the same about their own results. What this
repository does is narrower and checkable: it takes equations that were published with "one could imagine trying to find it
numerically" or "this upper limit is a future direction of work" or a disputed experimental signature, and computes them, exactly, in
the open, with every constant written as an integer and every result reproducible bit for bit on any machine.

**In plain words.** Physics papers since 2016 describe traversable wormholes (tunnels between two places that something could pass
through) and warp drives that the equations allow on paper, and some leave numbers open or in dispute: a speed limit left as future
work, a solution its authors left with "one could imagine trying to find it numerically", a quantum-computer experiment whose wormhole
signature other physicists dispute. The C programs here compute those numbers with fixed-point integer arithmetic (every calculation
done with whole numbers in a set format, so any computer prints the same digits), and a result counts as measured only when two
separate small ARM computer boards print it identically. Each C program ends with a chain hash (a fingerprint of a run's results, so
two runs can be compared exactly), so anyone can rerun it and compare the hash. The three results below show what came out.
Nothing here builds anything: these are checkable numbers about published equations.

## The results to read first

- **The 2024 warp shell's speed limit, which its paper left open** ([`results/3b_warp.md`](results/3b_warp.md)). Fuchs and co-authors
  described on paper a warp drive built from a shell of positive-energy matter, moving at a constant 2 percent of the speed of light
  (0.02 c; the paper gives no way to accelerate it), with all four energy conditions holding (standard tests that the matter behaves like
  ordinary matter, for example that no observer sees negative energy), and called the upper limit on that speed "a future direction of
  work". The limit is about **0.03035 c** (COMPUTED by the reviewer for the shell as defined); the boards' finest grid (their most closely
  spaced set of points) puts it no higher than **0.030358 c** (MEASURED). Just above it all four conditions fail together, 12.3 m from
  the shell's centre.
- **The teleportation protocol behind the 2022 quantum-computer wormhole experiment, computed exactly**
  ([`results/2_syk_exact.md`](results/2_syk_exact.md)). The protocol (a recipe that inserts a quantum bit on one side of two linked
  systems and reads it out on the other) runs on the SYK model (a model of many interacting particles that physicists use as a stand-in
  for a wormhole), whose couplings (the strengths of the interactions) are random. The protocol paper's formulas are averages over those
  couplings; for one draw of the couplings they are off from the exact result by up to **14 percent at N = 8** particles a side and
  4 percent at N = 12. And teleportation works better when the coupling that links the two sides has the wormhole sign (the one that
  gives negative energy) than with the opposite sign, at every size from N = 8 to 16, by **1.4 to 1.6 times** at the peak (both
  COMPUTED from the boards' MEASURED values).
- **The shape of the extra dimension of the "humanly traversable" wormhole, solved step by step in integer arithmetic**
  ([`results/1a_flow.md`](results/1a_flow.md)). Maldacena and Milekhin give equations (their eq. 3.19 and A.30) for the shape of the
  fifth dimension (one that their model adds to the four of space and time) in their wormhole, whose mouths are wide enough that a
  person would survive its tides (gravity's stretching of a body). The two boards print the same chain hashes (MEASURED), and so do a PC
  build and the reviewer's separately written program (COMPUTED): builds on two kinds of chip (ARM and x86) agreeing bit for bit. The
  numbers that do not depend on where the coordinate's zero is placed (the gauge invariants) agree with an independent 40-digit
  calculation to **2.5e-11** (COMPUTED; 2.5e-11 is 0.000000000025). A later version computes a constant the paper had estimated to one
  digit, **kappa_flat = 0.128210970524** (MEASURED: both boards printed it to twelve digits, line for line, on 11-12; confirmed four independent
  ways). With the logarithm the solved flow gives (this constant plus a term from the paper's eq. 3.25 that its estimate left out), the paper's
  derived numbers move by 4 to 9 percent, each within the paper's one-digit estimate; the constant alone moves them by under 1 percent
  (COMPUTED by the reviewer).

The rest, including the explicit wormhole geometry still under way, is in [`PLAN.md`](PLAN.md#details-from-the-readme), under "Details from the README".

**Why exactness matters.** Every C program computes with integers in a declared format and no math library (square roots and
exponentials come from the repository's own integer code), so the same source prints the same digits on any machine: a result is a
fact about the equations, not about one computer's rounding, and two runs either agree to the last bit or they do not. That is what
lets two small boards, a PC and the reviewer's own programs check one another line for line ([Why exact](PLAN.md#why-exact), in `PLAN.md`, has the details).

**How it was checked.** Every number is labelled: MEASURED (printed identically by the two boards, or a published value, named),
COMPUTED (derived, or printed by a PC development run) or ESTIMATE; error bounds (how far each printed number can be from the exact
answer) are worked out on paper beside each program. The independent reviewer, the second AI named above, checks every program and
every page from copies, with its own derivations and programs, and real errors it found have been fixed: a truncated term in
the warp program that put the speed limit 3.25e-5 too low (COMPUTED); a 64-bit sum in a geometry solver that overflowed next to the
boundary and moved a result by about 3.5e-7 (COMPUTED); a division check a factor of two too weak in seven copies (no published run had
reached it); a table of the Casimir energy (the energy of empty space between two close surfaces, which can be below zero) that had
lost digits; and averaged formulas used where the exact single-draw ones were needed. Each correction is written on its page, and the
earlier numbers are kept ([How it was built, and how it was checked](PLAN.md#how-it-was-built-and-how-it-was-checked), in `PLAN.md`).

## Check it yourself in two minutes

These commands are for a Unix shell: Linux, macOS, WSL (Ubuntu) on Windows, or Git Bash. Windows PowerShell does not accept `&&`, and Windows has no `cc` until a C compiler is installed. From PowerShell with WSL installed, the one-liner is `wsl bash -c "cd code && cc -O2 -std=gnu11 -o whflow wh_flow.c && ./whflow | grep -oE 'chain [0-9a-f]+|WHFLOW: done'"`.

```sh
git clone https://github.com/prestonericson7-lang/wormhole-math-exact.git
cd wormhole-math-exact/code
cc -O2 -std=gnu11 -o whflow wh_flow.c && ./whflow | grep -oE 'chain [0-9a-f]+|WHFLOW: done'
```

On a little-endian machine (x86, ARM) the run prints three lines: the two chain hashes, `24d2f870a2abea59` for the flat flow and
`85a3df31c75d8935` for the sphere flow, and `WHFLOW: done`. They are the hashes of every step of the integration; the two ARM machines
printed the same two, and so did a PC build of this source and the reviewer's own reference program, written separately from the same
specification and built for x86 and for 32-bit ARM under an emulator. Every C program here ends the same way: a chain hash and a `done`
line, and every results page for a run of one of them quotes the chain its numbers came from. If your build prints a different hash,
something is different, and we want to know what.

## What is here

- [`PLAN.md`](PLAN.md): the computations, what each paper left open, what we deliver, in what order, what "exact" means here, and the
  status of each.
- [`THE_GLUE.md`](THE_GLUE.md): the reading that led here: what holds a traversable wormhole open in the post-2016 physics, why those
  wormholes cannot be time machines by construction, and where the energy scale meets the collider.
- [`results/`](results/README.md): each result as a page with the numbers, the error bound on paper and the comparison with the paper it
  checks; for a run on the machines, the program, its chain, the two machines' agreement and the raw output of one machine
  (`*_board1.txt`). PC development runs and the one derivation (`1b_nogo.md`) are labelled as such; `1b_cost.md` is a formula
  evaluation on the PC.
- [`code/`](code/): the programs, on the fixed-point headers `whq.h` (Q.40), `whq60.h` (Q.60, built on `whq.h`) and `whq56.h` (Q.56), with
  128-bit intermediates built from 32-bit halves; SYK versions 4 to 7, the kernels `wh_relax.c` and `wh_relaxs.c`, and `wh_casimir.c`
  carry their Q.60 routines in their own source, and the limb kernel has `limb_common.h`; the Python helpers (`wh_cost.py`,
  `wind_summary.py` and the scripts in `p2_nogo/`) are not on the headers. The programs by computation are listed in [`PLAN.md`](PLAN.md#the-programs-by-computation).

## Licenses

Code: MIT ([`code/LICENSE`](code/LICENSE)). Text and tables: CC BY 4.0 ([`LICENSE-TEXT`](LICENSE-TEXT)). Please cite as: Preston Ericson, "Wormhole math, exact",
2026, this repository.
