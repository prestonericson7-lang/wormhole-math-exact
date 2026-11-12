/* wh_flow.c (11-11-2026; v2 4:10 PM) -- WORMHOLE MATH 1a: the Maldacena-Milekhin flow, exact.
 *
 * Maldacena & Milekhin, "Humanly traversable wormholes", arXiv:2008.06618, eq. (3.19) with l_B = 1 (their Appendix A, eq. A.28-A.29):
 *   0 = -6 + 3 e^{-4s} + L'^2 + 4 L' s' + s'^2            (so  L' = -2 s' + sqrt(3 s'^2 + 6 - 3 e^{-4s}),  the root with L' = sqrt3 at s = s' = 0)
 *   0 =  4 - 4 e^{-4s} - 2 L' s' - 2 s'^2 - s''           (so  s'' = 4 - 4 e^{-4s} - 2 L' s' - 2 s'^2)
 * where s = sigma(rho) and L = lambda(rho). IR (rho -> -inf): s -> 0, L' -> sqrt3 (AdS3 x R2); UV (rho -> +inf): s ~ L ~ rho + const (AdS5).
 * Their boundary condition: the growing solution of the linearized equation around s = 0, s'' + 2 sqrt3 s' - 16 s = 0, i.e. s = eps e^{k rho},
 * k = sqrt(19) - sqrt(3) = 2.6268..., integrated toward the UV. Their Figure 3 shows sigma and lambda; this program prints them as a table.
 * And the flow on a two-sphere, their eq. (A.30) with l_B-hat = 1:
 *   0 = -6 - e^{-2s} + 3 e^{-4s} + L'^2 + 4 L' s' + s'^2 ;  0 = 4 + e^{-2s} - 4 e^{-4s} - 2 L' s' - 2 s'^2 - s''
 * whose IR fixed point solves 4 + u - 4 u^2 = 0, u = e^{-2s}; their (A.31) says e^{-s} = l + l^3/16 + ... (l = 1: 1.0625).
 *
 * THE GAUGE: rho is defined so that the linear mode has amplitude eps = 2^-10 at rho = -8 (s - s* = eps e^{k (rho + 8)}). Any other amplitude
 * gives the same curve shifted in rho; the decaying mode excited by starting with s' = k (s - s*) dies as e^{-(sqrt19 + sqrt3) rho}.
 *
 * EXACT: every quantity is an int64 in Q.40 (40 fractional bits, 2^-40 = 9.1e-13); products by mulq() (the full 128-bit product from four
 * 32x32 halves, shifted right by 40 with round-to-nearest); divisions by divq() (a 128-bit numerator over a 64-bit divisor, long division,
 * round-to-nearest); exp by range reduction on ln2 and a 24-term Taylor series; sqrt exact to the floor (integer Newton with a 128-bit check);
 * RK4 with a fixed step h = 2^-12, its sixths rounded to nearest. No floating point anywhere, no libm: the same bits on every machine and
 * every run. The result IS the printed digits. v1 (Q.28, floor) is withdrawn: its first increments in the deep IR were a sixth of a unit,
 * floored to zero, so the flow would have stood still; found 4:08 PM before any run.
 * Build:  cc -O2 -std=gnu11 -o whflow wh_flow.c
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define FB 40
#define ONE ((int64_t)1 << FB)
#define MASK (ONE - 1)
typedef int64_t q;

/* the full 128-bit product of two unsigned 64-bit numbers, from four 32x32 halves */
static void mul128(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo)
{
    const uint64_t a1 = a >> 32, a0 = a & 0xffffffffull, b1 = b >> 32, b0 = b & 0xffffffffull;
    const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    uint64_t mid = p01 + p10;
    const uint64_t c1 = mid < p01;                    /* carry out of the middle sum */
    const uint64_t l = p00 + (mid << 32);
    const uint64_t c2 = l < p00;
    *lo = l;
    *hi = p11 + (mid >> 32) + (c1 << 32) + c2;
}
static q mulq(q a, q b)                               /* round(a*b / 2^FB) to nearest, ties away from zero; |a*b| < 2^103 */
{
    const int neg = (a < 0) != (b < 0);
    const uint64_t ua = a < 0 ? -(uint64_t)a : (uint64_t)a, ub = b < 0 ? -(uint64_t)b : (uint64_t)b;
    uint64_t hi, lo; mul128(ua, ub, &hi, &lo);
    uint64_t r = (lo >> FB) | (hi << (64 - FB));
    if ((lo & MASK) >= ((uint64_t)1 << (FB - 1))) r++;
    return neg ? -(q)r : (q)r;
}
static uint64_t div128(uint64_t hi, uint64_t lo, uint64_t b, uint64_t *rem)   /* (hi:lo) / b, needs hi < b; restoring long division */
{
    uint64_t qt = 0;
    for (int i = 63; i >= 0; i--) {
        const uint64_t top = hi >> 63;
        hi = (hi << 1) | (lo >> 63); lo <<= 1;
        if (top || hi >= b) { hi -= b; qt |= (uint64_t)1 << i; }
    }
    *rem = hi; return qt;
}
static q divq(q a, q b)                               /* round(a * 2^FB / b) to nearest for b > 0; |a/b| < 2^23 */
{
    const int neg = a < 0;
    const uint64_t ua = neg ? -(uint64_t)a : (uint64_t)a;
    uint64_t rem;
    uint64_t r = div128(ua >> (64 - FB), ua << FB, (uint64_t)b, &rem);
    if (2 * rem >= (uint64_t)b) r++;
    return neg ? -(q)r : (q)r;
}
static q sqrtq(q x)                                   /* floor(sqrt(x)) in Q.FB exactly: R = floor(sqrt(X 2^FB)) with a 128-bit check */
{
    if (x <= 0) return 0;
    const uint64_t X = (uint64_t)x;
    /* integer sqrt of X, then Newton on the 128-bit radicand X << FB */
    uint64_t s = 0; for (int bit = 31; bit >= 0; bit--) { const uint64_t t = s | ((uint64_t)1 << bit); if (t * t <= X) s = t; }
    uint64_t R = (s + 1) << (FB / 2);                 /* an upper start: (s+1)^2 > X */
    for (int i = 0; i < 8; i++) {
        uint64_t rem; const uint64_t qq = div128(X >> (64 - FB), X << FB, R, &rem);
        const uint64_t nR = (R + qq) >> 1;
        if (nR >= R) break;
        R = nR;
    }
    for (;;) {                                        /* fix-up: R^2 <= X 2^FB < (R+1)^2, compared in 128 bits */
        uint64_t hi, lo; mul128(R, R, &hi, &lo);
        const uint64_t Xhi = X >> (64 - FB), Xlo = X << FB;
        if (hi > Xhi || (hi == Xhi && lo > Xlo)) { R--; continue; }
        mul128(R + 1, R + 1, &hi, &lo);
        if (hi < Xhi || (hi == Xhi && lo <= Xlo)) { R++; continue; }
        break;
    }
    return (q)R;
}
static q rdiv(q x, int n)                             /* round(x / n) to nearest, ties away from zero, for a small positive n */
{
    return x >= 0 ? (x + n / 2) / n : -((-x + n / 2) / n);
}
static const q LN2 = 762123384786;                    /* ln 2 * 2^40 = 762123384786.05... */
static q expq(q x)                                    /* e^x in Q.FB: 2^k * e^r, |r| <= ln2/2, Taylor 24 terms; our arguments stay in [-48, 1] */
{
    int64_t k = x / LN2;
    q r = x - k * LN2;
    if (r > LN2 / 2) { r -= LN2; k++; }
    if (r < -LN2 / 2) { r += LN2; k--; }
    q term = ONE, sum = ONE;
    for (int n = 1; n <= 24; n++) { term = rdiv(mulq(term, r), n); sum += term; }
    if (k >= 63 - FB - 1) return (q)1 << 62;          /* out of our range: never reached */
    if (k <= -62) return 0;
    if (k >= 0) return sum << k;
    return (sum + ((q)1 << (-k - 1))) >> (-k);        /* round to nearest */
}
static void fmt(q v, char *out)                       /* a Q.FB value as a decimal with 12 places, exact (no float; two 6-digit stages) */
{
    const int neg = v < 0;
    const uint64_t a = (uint64_t)(neg ? -v : v);
    const uint64_t ip = a >> FB, fp = a & MASK;
    const uint64_t t1 = fp * 1000000ull;              /* fp < 2^40, times 1e6 < 2^60 */
    const uint64_t d1 = t1 >> FB, t2 = (t1 & MASK) * 1000000ull, d2 = t2 >> FB;
    sprintf(out, "%s%llu.%06llu%06llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)d1, (unsigned long long)d2);
}

/* the two flows: sphere = 0 is (3.19)/(A.28-29), sphere = 1 is (A.30) */
static q Lprime(q s, q sp, int sphere)
{
    const q e4 = expq(-4 * s), e2 = sphere ? expq(-2 * s) : 0;
    const q under = 3 * mulq(sp, sp) + 6 * ONE - 3 * e4 + e2;      /* (2s')^2 - (s'^2 + 3e^-4s + e^-2s - 6) */
    return -2 * sp + sqrtq(under);
}
static q spp(q s, q sp, int sphere)
{
    const q e4 = expq(-4 * s), e2 = sphere ? expq(-2 * s) : 0;
    const q lp = Lprime(s, sp, sphere);
    return 4 * ONE + e2 - 4 * e4 - 2 * mulq(lp, sp) - 2 * mulq(sp, sp);
}
static void rk4(q *s, q *sp, q *L, q h, int sphere)   /* one step of (s, s', L) with s'' = spp and L' = Lprime */
{
    const q k1s = *sp, k1p = spp(*s, *sp, sphere), k1L = Lprime(*s, *sp, sphere);
    const q h2 = h / 2;
    const q s2 = *s + mulq(h2, k1s), p2 = *sp + mulq(h2, k1p);
    const q k2s = p2, k2p = spp(s2, p2, sphere), k2L = Lprime(s2, p2, sphere);
    const q s3 = *s + mulq(h2, k2s), p3 = *sp + mulq(h2, k2p);
    const q k3s = p3, k3p = spp(s3, p3, sphere), k3L = Lprime(s3, p3, sphere);
    const q s4 = *s + mulq(h, k3s), p4 = *sp + mulq(h, k3p);
    const q k4s = p4, k4p = spp(s4, p4, sphere), k4L = Lprime(s4, p4, sphere);
    *s += rdiv(mulq(h, k1s + 2 * k2s + 2 * k3s + k4s), 6);
    *sp += rdiv(mulq(h, k1p + 2 * k2p + 2 * k3p + k4p), 6);
    *L += rdiv(mulq(h, k1L + 2 * k2L + 2 * k3L + k4L), 6);
}

static uint64_t fnv(uint64_t h, const void *p, size_t n) { const uint8_t *b = p; for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 0x100000001b3ull; } return h; }

int main(void)
{
    char a[48], b[48], c[48], d[48], e[48], f[48];
    const q h = ONE >> 12;                            /* rho step 2^-12 */
    const q sqrt3 = sqrtq(3 * ONE), sqrt19 = sqrtq(19 * ONE), k = sqrt19 - sqrt3;
    printf("WHFLOW: Q.%d fixed point, int64 with 128-bit intermediates, round to nearest; exp by ln2 reduction + 24 Taylor terms; sqrt exact; RK4 h = 2^-12; gauge eps = 2^-10 at rho = -8\n", FB);
    /* self-test, printed so the run log carries it: the arithmetic on known values */
    fmt(sqrt3, a); fmt(sqrt19, b); fmt(k, c); fmt(expq(LN2), d); fmt(expq(-LN2), e); fmt(divq(ONE, 3 * ONE), f);
    printf("WHFLOW selftest: sqrt3 = %s (1.732050807569), sqrt19 = %s (4.358898943541), k = %s (2.626848135972), e^ln2 = %s, e^-ln2 = %s, 1/3 = %s\n", a, b, c, d, e, f);
    fmt(mulq(3 * ONE, 5 * ONE), a); fmt(sqrtq(4 * ONE), b); fmt(expq(0), c); fmt(mulq(sqrt3, sqrt3), d); fmt(expq(-48 * ONE), e);
    printf("WHFLOW selftest: 3*5 = %s, sqrt4 = %s, e^0 = %s, sqrt3^2 = %s (3 - at most 2^-38), e^-48 = %s\n", a, b, c, d, e);
    for (int sphere = 0; sphere < 2; sphere++) {
        uint64_t hash = 0xcbf29ce484222325ull;
        /* the IR fixed point: sphere 0: s* = 0 (e^-4s = 1); sphere 1: u = e^-2s solves 4 + u - 4u^2 = 0, u = (1 + sqrt65)/8 */
        q sstar = 0;
        if (sphere) {
            const q u = divq(ONE + sqrtq(65 * ONE), 8 * ONE);     /* e^{-2 s*} */
            const q ems = sqrtq(u);                                /* e^{-s*} */
            q s0 = 0; for (int i = 0; i < 60; i++) { const q ex = expq(-s0); s0 += divq(ex - ems, ex); }   /* s* = -ln(ems) by Newton */
            sstar = s0;
            fmt(ems, a); fmt(sstar, b); fmt(Lprime(sstar, 0, 1), c); fmt(spp(sstar, 0, 1), d);
            printf("WHFLOW sphere: IR fixed point e^-s* = %s (their A.31 expansion 1 + 1/16 = 1.0625), s* = %s, L' there = %s, s'' there = %s (0)\n", a, b, c, d);
        } else {
            fmt(Lprime(0, 0, 0), a); fmt(spp(0, 0, 0), b);
            printf("WHFLOW flat: IR fixed point s* = 0, L' there = %s (sqrt3), s'' there = %s (0)\n", a, b);
        }
        /* start deep in the IR on the growing linearized mode: s - s* = eps e^{k (rho + 8)}, s' = k (s - s*), eps = 2^-10 */
        q rho = -8 * ONE, s = sstar + (ONE >> 10), sp = mulq(k, ONE >> 10), L = mulq(Lprime(sstar, 0, sphere), rho);
        int steps = 0;
        const char *name = sphere ? "sphere" : "flat";
        printf("WHFLOW %s: rho | sigma | sigma' | lambda' | lambda - rho | sigma - rho\n", name);
        while (rho < 12 * ONE) {
            if ((steps & 1023) == 0) {                /* every 0.25 in rho */
                const q lp = Lprime(s, sp, sphere);
                fmt(rho, a); fmt(s, b); fmt(sp, c); fmt(lp, d); fmt(L - rho, e); fmt(s - rho, f);
                printf("WHFLOW %s: %s | %s | %s | %s | %s | %s\n", name, a, b, c, d, e, f);
            }
            rk4(&s, &sp, &L, h, sphere);
            hash = fnv(hash, &s, 8); hash = fnv(hash, &sp, 8); hash = fnv(hash, &L, 8);
            rho += h; steps++;
            if (s > 12 * ONE) break;                  /* deep in the UV: e^-4s has vanished (2^-69), nothing more to learn */
        }
        const q lp = Lprime(s, sp, sphere);
        fmt(rho, a); fmt(sp, b); fmt(lp, c); fmt(L - s, d);
        printf("WHFLOW %s: END at rho = %s after %d steps: sigma' = %s (UV: 1), lambda' = %s (UV: 1), lambda - sigma = %s (UV: a constant); chain %016llx\n",
               name, a, steps, b, c, d, (unsigned long long)hash);
    }
    printf("WHFLOW: done\n");
    return 0;
}
