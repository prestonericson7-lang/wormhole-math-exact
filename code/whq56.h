/* whq.h (11-11-2026) -- exact fixed-point arithmetic shared by the wormhole-math programs (wh_flow.c, wh_alcub.c, ...).
 *
 * Every quantity is an int64 in Q.40 (40 fractional bits, 2^-40 = 9.1e-13; magnitudes below 2^23). No floating point, no libm.
 *   mulq(a, b)  round(a b / 2^40) to nearest, ties away from zero, from the full 128-bit product (four 32x32 halves); needs |a b| < 2^103.
 *   divq(a, b)  round(a 2^40 / b) to nearest for b > 0, by a 128-by-64 long division; needs |a / b| < 2^23.
 *   sqrtq(x)    floor(sqrt(x)) exactly, by integer Newton with a 128-bit fix-up.
 *   rdiv(x, n)  round(x / n) to nearest for a small positive integer n.
 *   expq(x)     e^x by range reduction on ln 2 and a 24-term Taylor series (|r| <= ln2/2, so the 25th term is below 2^-40 / 10^5).
 *   fmt(v, s)   the value as a decimal with 12 places, exact (two 6-digit stages).
 *   fnv(h,p,n)  the FNV-1a chain hash used to fingerprint every step of a run.
 * The same bits on every machine and every run: the result IS the printed digits.
 */
/* whq56.h (11-12-2026, 6:30 AM) -- whq.h's fixed-point primitives at 56 fraction bits (values below 128), generated from whq.h by
 * mk_whq56.py with only the format and its two constants changed: FB 56, ln 2 = 49946518145322874, pi = 226375608064910089 (both
 * round(x 2^56)). mulq and divq round to nearest, sqrtq is the floor root, expq saturates at 64 (k >= 6), lnq by range reduction and
 * artanh. This is the arithmetic the reviewer's P3 mimic used ("whq.h's semantics at 56 bits"), so wh_p3.c can match its integers. */
#ifndef WHQ56_H
#define WHQ56_H
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define FB 56
#define ONE ((int64_t)1 << FB)
#define MASK (ONE - 1)
typedef int64_t q;

static void mul128(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
    const uint64_t a1 = a >> 32, a0 = a & 0xffffffffull, b1 = b >> 32, b0 = b & 0xffffffffull;
    const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const uint64_t mid = p01 + p10;
    const uint64_t c1 = mid < p01;
    const uint64_t l = p00 + (mid << 32);
    const uint64_t c2 = l < p00;
    *lo = l;
    *hi = p11 + (mid >> 32) + (c1 << 32) + c2;
}
static q mulq(q a, q b)
{
    const int neg = (a < 0) != (b < 0);
    const uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b;
    uint64_t hi, lo; mul128(ua, ub, &hi, &lo);
    uint64_t r = (lo >> FB) | (hi << (64 - FB));
    if ((lo & MASK) >= ((uint64_t)1 << (FB - 1))) r++;
    return neg ? -(q)r : (q)r;
}
static uint64_t div128(uint64_t hi, uint64_t lo, uint64_t b, uint64_t *rem)   /* (hi:lo) / b, needs hi < b */
{
    uint64_t qt = 0;
    for (int i = 63; i >= 0; i--) {
        const uint64_t top = hi >> 63;
        hi = (hi << 1) | (lo >> 63); lo <<= 1;
        if (top || hi >= b) { hi -= b; qt |= (uint64_t)1 << i; }
    }
    *rem = hi; return qt;
}
static q divq(q a, q b)                                 /* a / b rounded to nearest; b of either sign (v2 11-11 7:35 PM: v1 assumed b > 0 and
                                                           returned a/2^24-sized garbage for b < 0; 1a, 3a and the shell never divided by a negative) */
{
    if (b < 0) { a = -a; b = -b; }
    const int neg = a < 0;
    const uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a;
    uint64_t rem;
    uint64_t r = div128(ua >> (64 - FB), ua << FB, (uint64_t)b, &rem);
    if (2 * rem >= (uint64_t)b) r++;
    return neg ? -(q)r : (q)r;
}
static q sqrtq(q x)
{
    if (x <= 0) return 0;
    const uint64_t X = (uint64_t)x;
    uint64_t s = 0; for (int bit = 31; bit >= 0; bit--) { const uint64_t t = s | ((uint64_t)1 << bit); if (t * t <= X) s = t; }
    uint64_t R = (s + 1) << (FB / 2);
    for (int i = 0; i < 8; i++) {
        uint64_t rem; const uint64_t qq = div128(X >> (64 - FB), X << FB, R, &rem);
        const uint64_t nR = (R + qq) >> 1;
        if (nR >= R) break;
        R = nR;
    }
    for (;;) {
        uint64_t hi, lo; mul128(R, R, &hi, &lo);
        const uint64_t Xhi = X >> (64 - FB), Xlo = X << FB;
        if (hi > Xhi || (hi == Xhi && lo > Xlo)) { R--; continue; }
        mul128(R + 1, R + 1, &hi, &lo);
        if (hi < Xhi || (hi == Xhi && lo <= Xlo)) { R++; continue; }
        break;
    }
    return (q)R;
}
static q rdiv(q x, int n) { return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n); }
static const q LN2 = 49946518145322874LL;                    /* ln 2 * 2^40 = 762123384786.05... */
static q expq(q x)
{
    int64_t k = x / LN2;
    q r = x - k * LN2;
    if (r > LN2 / 2) { r -= LN2; k++; }
    if (r < -LN2 / 2) { r += LN2; k--; }
    q term = ONE, sum = ONE;
    for (int n = 1; n <= 24; n++) { term = rdiv(mulq(term, r), n); sum += term; }
    if (k >= 63 - FB - 1) return (q)1 << 62;
    if (k <= -62) return 0;
    if (k >= 0) return sum << k;
    return (sum + ((q)1 << (-k - 1))) >> (-k);
}
static const q PI = 226375608064910089LL;                    /* pi * 2^40 = 3454217652357.64 */
static q lnq(q x)                                     /* ln x for x > 0: x = 2^k y with y in [1, 2); ln y = 2 artanh((y-1)/(y+1)), 30 terms */
{
    if (x <= 0) return -((q)1 << 62);
    int k = 0; q y = x;
    while (y >= 2 * ONE) { y >>= 1; k++; }
    while (y < ONE) { y <<= 1; k--; }
    const q t = divq(y - ONE, y + ONE), t2 = mulq(t, t);
    q term = t, sum = 0;
    for (int n = 1; n <= 59; n += 2) { sum += rdiv(term, n); term = mulq(term, t2); }
    return (q)k * LN2 + 2 * sum;
}
static q cosq(q x);
static q sinq(q x)                                    /* sin x: reduce to [-pi, pi], Taylor to the 21st power (|x|^23/23! < 2^-40 for |x| <= pi) */
{
    const q twopi = 2 * PI;
    q k = x / twopi; x -= k * twopi;
    if (x > PI) x -= twopi; if (x < -PI) x += twopi;
    const q x2 = mulq(x, x);
    q term = x, sum = x;
    for (int n = 1; n <= 11; n++) { term = -rdiv(mulq(term, x2), (2 * n) * (2 * n + 1)); sum += term; }
    return sum;
}
static q cosq(q x)
{
    const q twopi = 2 * PI;
    q k = x / twopi; x -= k * twopi;
    if (x > PI) x -= twopi; if (x < -PI) x += twopi;
    const q x2 = mulq(x, x);
    q term = ONE, sum = ONE;
    for (int n = 1; n <= 11; n++) { term = -rdiv(mulq(term, x2), (2 * n - 1) * (2 * n)); sum += term; }
    return sum;
}
static void fmt40_only(q v, char *out)   /* whq.h's twelve-decimal formatter; its 10^6 stage overflows at 56 bits: unused here */
{
    const int neg = v < 0;
    const uint64_t a = (uint64_t)(neg ? -v : v);
    const uint64_t ip = a >> FB, fp = a & MASK;
    const uint64_t t1 = fp * 1000000ull;
    const uint64_t d1 = t1 >> FB, t2 = (t1 & MASK) * 1000000ull, d2 = t2 >> FB;
    sprintf(out, "%s%llu.%06llu%06llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)d1, (unsigned long long)d2);
}
static uint64_t fnv(uint64_t h, const void *p, size_t n) { const uint8_t *b = p; for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 0x100000001b3ull; } return h; }

/* a 128-bit unsigned accumulator for long Simpson sums */
typedef struct { uint64_t hi, lo; } acc128;
static inline void acc_add(acc128 *a, uint64_t v) { a->lo += v; if (a->lo < v) a->hi++; }
static inline uint64_t acc_shr(const acc128 *a, int sh) { return (a->lo >> sh) | (a->hi << (64 - sh)); }   /* needs hi < 2^sh */

#endif
