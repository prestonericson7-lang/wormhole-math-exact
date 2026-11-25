/* par2.h (the public copy, 11-12-2026) -- the two-core split of a matrix product or a relaxation sweep: par2(fn, a0, a1) runs fn(a0)
 * and fn(a1) and returns when both have ended. On the two ARM machines a machine-specific version of this header runs fn(a1) on the
 * second core where one is available, and both halves on one core otherwise, as here. The two halves write disjoint rows and share
 * nothing else, so the result is the same bits either way: this one-core copy reproduces the published chains of the programs that
 * include it. */
#ifndef PAR2_H
#define PAR2_H
static inline void par2(void (*fn)(void *), void *a0, void *a1) { fn(a0); fn(a1); }
#endif
