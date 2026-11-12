/* limb_common.h (11-11-2026) -- shared by limb_q40 and limb_q60.
 *
 * lq_s128   the layout of the main program's s128 (wh_syk5.c:31): the value hi 2^64 + lo, two's complement modulo 2^128.
 * lq_c64    the layout of the main program's c60 (wh_syk5.c:30) and cq (wh_syk.c:43): one complex fixed-point entry {re, im}.
 * lq_round_sh(a, sh)   the main program's s128_round_sh (wh_syk5.c:70-75; wh_syk.c:65-73 is its sh = 40 case), bit for bit, for 1 <= sh <= 63:
 *     neg = (a < 0); U = |a| (mod 2^128); r = floor(U / 2^sh) mod 2^64; r += (U mod 2^sh >= 2^(sh-1)); result = neg ? -r : r (mod 2^64).
 *   That is round to nearest, ties away from zero, whenever the rounded value fits an int64. Written here in unsigned arithmetic, so it
 *   has no signed overflow (the main program's copy negates hi in int64, which overflows only for a = -2^127; the bits are the same wherever
 *   theirs is defined). sh = 64 is NOT supported here (wh_syk5.c:73 shifts a uint64 by 64 there, which C leaves undefined).
 * lq_bits(v)   the smallest B >= 0 with |v| < 2^B (64 for INT64_MIN): the operand bound the contracts are stated in.
 */
#ifndef LIMB_COMMON_H
#define LIMB_COMMON_H
#include <stdint.h>

typedef struct { int64_t hi; uint64_t lo; } lq_s128;
typedef struct { int64_t re, im; } lq_c64;

static inline int64_t lq_round_sh(lq_s128 a, int sh)
{
    const int neg = a.hi < 0;
    uint64_t lo = a.lo, hi = (uint64_t)a.hi;
    if (neg) { lo = ~lo + 1; hi = ~hi + (lo == 0); }
    uint64_t r = (lo >> sh) | (hi << (64 - sh));
    if ((lo & (((uint64_t)1 << sh) - 1)) >= ((uint64_t)1 << (sh - 1))) r++;
    return (int64_t)(neg ? (uint64_t)0 - r : r);
}

static inline int lq_bits(int64_t v)
{
    const uint64_t u = v < 0 ? (uint64_t)0 - (uint64_t)v : (uint64_t)v;
#if defined(__GNUC__)
    return u ? 64 - __builtin_clzll(u) : 0;
#else
    uint64_t w = u; int b = 0; while (w) { b++; w >>= 1; } return b;
#endif
}

/* the limb product of the kernels: int32 x int32 -> int64, added to a SIGNED int64 accumulator: one SMLAL on the A9.
   The accumulators are signed on purpose. With unsigned (wrapping) sums gcc 11 -O2 reassociates each k term's products into a separate
   SMULL/SMLAL tree added with ADDS/ADC, and the A9 code grew from 73 to 147 instructions a complex MAC (Q.60, measured, count_a9.out
   of the first build); signed sums stay as in-place SMLAL chains. The contracts keep every signed accumulator inside int64, so no
   overflow happens under them. (For the same reason do not build the kernel files with -fwrapv: correct, but slower.) */
#define LQ_P(x, y) ((int64_t)(x) * (int64_t)(y))

#endif
