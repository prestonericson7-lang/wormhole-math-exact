/* whq60.h (11-12-2026, 3:05 AM) -- the Q.60 fixed-point header for the 1b programs and the twelve-digit flow: int64 with 60 fraction bits
 * (values below 8 in size), every product and quotient through 128 bits and rounded once to nearest, ties away from zero. The pieces are the
 * ones proved in the SYK programs (wh_syk4.c to wh_syk7.c): the 128-bit accumulator, mul60, div60, the digit-by-digit 128-bit square root,
 * ln by range reduction and artanh, cos and sin by Taylor; new here: exp60 by ln 2 reduction and a Taylor series with a declared number of
 * terms (as whq.h's expq, in Q.60), a range-reduced cossin60 for any angle, parse60 for decimal arguments, and a self-test. No floating
 * point, no libm: the same bits on every machine. Requires whq.h beside it for mul128 and div128 (its 64 x 64 -> 128 product and 128/64
 * division) and fnv. */
#ifndef WHQ60_H
#define WHQ60_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "whq.h"

#define F60 60
static const int64_t ONE60 = (int64_t)1 << F60;
typedef struct { int64_t hi; uint64_t lo; } s128;
static inline void s128_add(s128 *a, int64_t hi, uint64_t lo) { uint64_t l = a->lo + lo; a->hi += hi + (l < lo); a->lo = l; }
static inline void s128_addv(s128 *a, int64_t v) { s128_add(a, v < 0 ? -1 : 0, (uint64_t)v); }
static inline void mulprod(int64_t a, int64_t b, int64_t *hi, uint64_t *lo)
{
    int neg = (a < 0) ^ (b < 0); uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b, h, l;
    mul128(ua, ub, &h, &l); if (neg) { l = ~l + 1; h = ~h + (l == 0); } *hi = (int64_t)h; *lo = l;
}
static inline void s128_addp(s128 *a, int64_t x, int64_t y) { int64_t h; uint64_t l; mulprod(x, y, &h, &l); s128_add(a, h, l); }
static inline int64_t s128_round_sh(s128 a, int sh)   /* round(a / 2^sh) to nearest, ties away from zero; 1 <= sh <= 63 */
{
    int neg = a.hi < 0; if (neg) { a.lo = ~a.lo + 1; a.hi = ~a.hi + (a.lo == 0); }
    uint64_t r = (a.lo >> sh) | ((uint64_t)a.hi << (64 - sh)); if ((a.lo & ((((uint64_t)1) << sh) - 1)) >= (((uint64_t)1) << (sh - 1))) r++;
    return neg ? -(int64_t)r : (int64_t)r;
}
static inline int64_t mul60(int64_t a, int64_t b) { s128 s = {0, 0}; s128_addp(&s, a, b); return s128_round_sh(s, F60); }
static inline int64_t rdiv64(int64_t x, int64_t n) { return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n); }
static int64_t div60(int64_t a, int64_t b)             /* round(a 2^60 / b), b > 0, |a| < 8 b; both conditions checked (the quotient must fit an int64) */
{
    if (b <= 0) { printf("WHQ60: div60 by a non-positive value -- stop\n"); exit(1); }
    int neg = a < 0; uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a; const uint64_t hi = ua >> 4, lo = ua << 60;
    if (hi >= (uint64_t)b) { printf("WHQ60: div60 overflow -- stop\n"); exit(1); }
    uint64_t rem; uint64_t qq = div128(hi, lo, (uint64_t)b, &rem); if (2 * rem >= (uint64_t)b) qq++;
    if (qq > (uint64_t)INT64_MAX) { printf("WHQ60: div60 overflow (|a| at or above 8 b) -- stop\n"); exit(1); }   /* the check above alone stops only at 16 b; between 8 b and 16 b the result came back with the wrong sign (found by the reviewer, 4:00 AM, 11-12) */
    return neg ? -(int64_t)qq : (int64_t)qq;
}
static int64_t div128s(s128 num, int64_t den)         /* round(num / den), num signed 128-bit, 0 < den < 2^63, the quotient an int64 */
{
    int neg = num.hi < 0; if (neg) { num.lo = ~num.lo + 1; num.hi = ~num.hi + (num.lo == 0); }
    if ((uint64_t)num.hi >= (uint64_t)den) { printf("WHQ60: a quotient does not fit -- stop\n"); exit(1); }
    uint64_t rem; uint64_t q = div128((uint64_t)num.hi, num.lo, (uint64_t)den, &rem); if (2 * rem >= (uint64_t)den) q++;
    return neg ? -(int64_t)q : (int64_t)q;
}
static uint64_t isqrt128(uint64_t hi, uint64_t lo)     /* floor(sqrt(hi:lo)), digit by digit */
{
    uint64_t rhi = 0, rlo = 0, root = 0;
    for (int i = 0; i < 64; i++) {
        const int top = i < 32 ? (int)((hi >> (62 - 2 * i)) & 3) : (int)((lo >> (126 - 2 * i)) & 3);
        rhi = (rhi << 2) | (rlo >> 62); rlo = (rlo << 2) | (uint64_t)top;
        const uint64_t thi = root >> 62, tlo = (root << 2) | 1;
        root <<= 1;
        if (rhi > thi || (rhi == thi && rlo >= tlo)) { const uint64_t nlo = rlo - tlo; rhi = rhi - thi - (rlo < tlo); rlo = nlo; root |= 1; }
    }
    return root;
}
static int64_t sqrt60(int64_t x)                       /* floor(sqrt(x)) in Q.60 for x >= 0: sqrt(x 2^60) with x 2^60 held in 128 bits */
{
    if (x <= 0) return 0;
    return (int64_t)isqrt128((uint64_t)x >> 4, (uint64_t)x << 60);
}
static int64_t artanh2_60(int64_t z)                  /* 2 artanh z for |z| < 1/3 or so: the series with 40 terms */
{
    const int64_t z2 = mul60(z, z); int64_t term = z, sum = 0;
    for (int k = 0; k < 40; k++) { const int64_t add = rdiv64(term, 2 * k + 1); if (!add && k > 2) break; sum += add; term = mul60(term, z2); if (!term) break; }
    return 2 * sum;
}
static int64_t LN2_60_cache = 0;
static int64_t ln2_60(void) { if (!LN2_60_cache) LN2_60_cache = artanh2_60(div60(ONE60, 3 * ONE60)); return LN2_60_cache; }   /* ln 2 = 2 artanh(1/3) */
static int64_t ln60(int64_t x)                         /* ln x for x > 0: x = 2^k m, m in [1/2, 1); ln m = 2 artanh((m-1)/(m+1)) */
{
    if (x <= 0) { printf("WHQ60: ln of a non-positive value -- stop\n"); exit(1); }
    int k = 0; int64_t m = x;
    while (m < (ONE60 >> 1)) { m <<= 1; k--; }
    while (m >= ONE60) { m = rdiv64(m, 2); k++; }
    const int64_t z = div60(m - ONE60, m + ONE60);
    return artanh2_60(z) + (int64_t)k * ln2_60();
}
static int EXP60_TERMS = 0;                            /* the largest number of Taylor terms any exp60 took (printed by the programs) */
static int64_t exp60n(int64_t n, int64_t f)            /* e^{n + f} in Q.60, n an integer in [-62, 2], f a Q.60 value: x = n + f = k ln2 + r, |r| <= ln2/2,
                                                          e^r by Taylor to a unit, then the shift by k. The sum n 2^60 + f is formed in 128 bits, since a
                                                          Q.60 number itself holds no value of size 8 or more (the 1b exponents reach -60). */
{
    const int64_t L = ln2_60();
    s128 x = { n >> 4, (uint64_t)n << 60 };            /* n 2^60 in 128 bits (the arithmetic shift keeps the sign) */
    s128_addv(&x, f);
    int64_t k = div128s(x, L);                         /* the nearest multiple of ln 2 */
    s128 rr = x; s128_addp(&rr, -k, L); int64_t r = s128_round_sh(rr, 0);
    if (r > L / 2) { r -= L; k++; }
    if (r < -L / 2) { r += L; k--; }
    int64_t term = ONE60, sum = ONE60; int it;
    for (it = 1; it <= 40; it++) { term = rdiv64(mul60(term, r), it); if (!term) break; sum += term; }
    if (it > EXP60_TERMS) EXP60_TERMS = it;
    if (k > 3 || (k == 3 && sum >= ONE60)) { printf("WHQ60: exp60 above the Q.60 ceiling -- stop (argument %lld + %lld units of 2^-60, k = %lld)\n", (long long)n, (long long)f, (long long)k); exit(1); }   /* e^x must stay below 8: k = 3 is allowed while e^r < 1 (x below 3 ln 2 = 2.079) */
    if (k <= -62) return 0;
    if (k >= 0) return sum << k;
    return (sum + ((int64_t)1 << (-k - 1))) >> (-k);
}
static int64_t exp60(int64_t x) { return exp60n(0, x); }   /* e^x for a Q.60 x (|x| < 7) */
static const int64_t PI60 = 3622009729038561421ll;    /* pi 2^60 = 3622009729038561421.3 */
static void cossin60(int64_t x, int64_t *c, int64_t *s)   /* any x: reduce to [-pi, pi] by 2 pi (one rounding of 2 pi), then to [-pi/2, pi/2] by a
                                                             sign flip (cos(x - pi) = -cos x: x^2 would pass the Q.60 ceiling at |x| > 2.83; found 6:50 AM 11-12
                                                             by the first caller past pi/2), then Taylor */
{
    const int64_t twopi = 2 * PI60;
    int64_t k = x / twopi; x -= k * twopi;
    if (x > PI60) x -= twopi; if (x < -PI60) x += twopi;
    int sign = 1; if (x > PI60 / 2) { x -= PI60; sign = -1; } else if (x < -PI60 / 2) { x += PI60; sign = -1; }
    const int64_t x2 = mul60(x, x); int64_t tc = ONE60, ts = x, sc = ONE60, ss = x;
    for (int kk = 1; kk < 40; kk++) {
        tc = -rdiv64(mul60(tc, x2), (int64_t)(2 * kk - 1) * (2 * kk)); ts = -rdiv64(mul60(ts, x2), (int64_t)(2 * kk) * (2 * kk + 1));
        sc += tc; ss += ts; if (!tc && !ts) break;
    }
    *c = sign * sc; *s = sign * ss;
}
static int64_t parse60(const char *s)                  /* a decimal into Q.60, round(v 2^60); |v| < 7 */
{
    int neg = 0; if (*s == '-') { neg = 1; s++; }
    uint64_t ip = 0; while (*s >= '0' && *s <= '9') ip = ip * 10 + (uint64_t)(*s++ - '0');
    uint64_t fp = 0, scale = 1; if (*s == '.') { s++; while (*s >= '0' && *s <= '9' && scale < 1000000000000ull) { fp = fp * 10 + (uint64_t)(*s++ - '0'); scale *= 10; } }
    if (ip >= 7) { printf("WHQ60: a value at or above 7 cannot be held in Q.60\n"); exit(1); }
    int64_t v = (int64_t)(ip << F60);
    if (fp) { uint64_t h, l, rem; mul128(fp, (uint64_t)1 << F60, &h, &l); uint64_t qq = div128(h, l, scale, &rem); if (2 * rem >= scale) qq++; v += (int64_t)qq; }
    return neg ? -v : v;
}
static void fmt60(int64_t v, char *out)                /* twelve decimals */
{
    int neg = v < 0; uint64_t u = neg ? -(uint64_t)v : (uint64_t)v;
    uint64_t ip = u >> F60, frac = u & ((((uint64_t)1) << F60) - 1), h, l;
    mul128(frac, 1000000000000ull, &h, &l);
    uint64_t dg = (l >> F60) | (h << (64 - F60)); if ((l >> (F60 - 1)) & 1) dg++;
    if (dg >= 1000000000000ull) { dg -= 1000000000000ull; ip++; }
    sprintf(out, "%s%llu.%012llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dg);
}
static void fmt60_18(int64_t v, char *out)             /* eighteen decimals: the full Q.60 content (2^-60 = 8.7e-19) */
{
    int neg = v < 0; uint64_t u = neg ? -(uint64_t)v : (uint64_t)v;
    uint64_t ip = u >> F60, frac = u & ((((uint64_t)1) << F60) - 1), h, l;
    mul128(frac, 1000000000000000000ull, &h, &l);
    uint64_t dg = (l >> F60) | (h << (64 - F60)); if ((l >> (F60 - 1)) & 1) dg++;
    if (dg >= 1000000000000000000ull) { dg -= 1000000000000000000ull; ip++; }
    sprintf(out, "%s%llu.%018llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dg);
}
/* the self-test: prints the known values beside the computed ones; returns the largest deviation in units */
static int64_t whq60_selftest(const char *tag)
{
    char a[48], b[48]; int64_t worst = 0;
    struct { const char *name; int64_t got; const char *want; } t[8];
    int n = 0;
    t[n].name = "ln 2"; t[n].got = ln2_60(); t[n++].want = "0.693147180559945309";
    t[n].name = "e"; t[n].got = exp60(ONE60); t[n++].want = "2.718281828459045235";
    t[n].name = "e^-1"; t[n].got = exp60(-ONE60); t[n++].want = "0.367879441171442322";
    t[n].name = "sqrt 2"; t[n].got = sqrt60(2 * ONE60); t[n++].want = "1.414213562373095049";
    t[n].name = "ln 3"; t[n].got = ln60(3 * ONE60); t[n++].want = "1.098612288668109691";
    { int64_t c, s; cossin60(PI60 / 6, &c, &s); t[n].name = "sin(pi/6)"; t[n].got = s; t[n++].want = "0.500000000000000000";
      cossin60(PI60 / 3, &c, &s); t[n].name = "cos(pi/3)"; t[n].got = c; t[n++].want = "0.500000000000000000"; }
    t[n].name = "e^-20"; t[n].got = exp60n(-20, 0); t[n++].want = "0.000000002061153622";
    for (int i = 0; i < n; i++) {
        fmt60_18(t[i].got, a);
        /* the deviation in units from the 18-decimal reference */
        const char *w = t[i].want; uint64_t ip = 0, fp = 0; while (*w != '.') ip = ip * 10 + (uint64_t)(*w++ - '0'); w++; while (*w) fp = fp * 10 + (uint64_t)(*w++ - '0');
        uint64_t h, l, rem; mul128(fp, (uint64_t)1 << F60, &h, &l); uint64_t qq = div128(h, l, 1000000000000000000ull, &rem); if (2 * rem >= 1000000000000000000ull) qq++;
        const int64_t want = (int64_t)(ip << F60) + (int64_t)qq, dev = t[i].got - want; const int64_t ad = dev < 0 ? -dev : dev; if (ad > worst) worst = ad;
        sprintf(b, "%lld", (long long)dev);
        printf("%s selftest %s = %s (%s; %s units)\n", tag, t[i].name, a, t[i].want, b);
    }
    return worst;
}
#endif
