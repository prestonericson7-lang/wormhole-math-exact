/* par2.h (the public copy, 11-12-2026) -- the two-core split of a matrix product or a relaxation sweep: par2(fn, a0, a1) runs fn(a0)
 * and fn(a1) and returns when both have ended. On the two ARM machines a version of this header hands fn(a1) to the second core
 * through the machines' own system software and falls back to this one-core form wherever that is unavailable; that version is
 * specific to those machines and is not published. The two halves write disjoint rows and share nothing else, so the result is the
 * same bits either way: every chain on the results pages was produced with the split on, and this file reproduces it on one core. */
#ifndef PAR2_H
#define PAR2_H
static inline void par2(void (*fn)(void *), void *a0, void *a1) { fn(a0); fn(a1); }
#endif
