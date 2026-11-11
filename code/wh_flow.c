/* wh_flow.c (11-11-2026) -- WORMHOLE MATH 1a: the Maldacena-Milekhin flow, exact.
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
 * EXACT: every quantity is an int64 in Q.28 (28 fractional bits); products by mulq() (a 64x64 -> 64 product in two halves, floor);
 * exp by range reduction on ln2 and a 24-term Taylor series in Q.28; sqrt by integer Newton (floor); RK4 with a fixed step h = 2^-12.
 * No floating point anywhere, no libm: the same bits on every machine and every run. The result IS the printed digits.
 * Build:  cc -O2 -std=gnu11 -o whflow wh_flow.c   (any C compiler, any machine: only 64-bit integers are used, so the output is the same everywhere)
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define FB 28
#define ONE ((int64_t)1 << FB)
#define MASK (ONE - 1)
typedef int64_t q;

static q mulq(q a, q b)                               /* floor(a*b / 2^FB), exact in two halves (no 128-bit type on the A9) */
{
    const q bh = b >> FB, bl = b & MASK;              /* b = bh*2^FB + bl, bl in [0, 2^FB) */
    const q p1 = a * bh;                              /* |a| < 2^40 and |bh| < 2^12 here: fits */
    const q p2 = a * bl;                              /* |a| < 2^35 needed: all our magnitudes are below 2^12 (4096.0) */
    return p1 + (p2 >> FB);                           /* >> on a negative value floors: deterministic */
}
static q divq(q a, q b)                               /* floor(a*2^FB / b) with a, b in Q.FB; |a| < 2^35 */
{
    return (a << FB) / b;                             /* C division truncates toward zero; fine for our positive divisors */
}
static const q LN2 = 186065279;                       /* ln 2 * 2^28 = 186065279.36... */
static q expq(q x)                                    /* e^x in Q.FB for x in about [-40, 20]: 2^k * e^r, |r| <= ln2/2, Taylor 24 terms */
{
    int64_t k = x / LN2;                              /* truncates toward zero */
    q r = x - k * LN2;
    if (r > LN2 / 2) { r -= LN2; k++; }
    if (r < -LN2 / 2) { r += LN2; k--; }
    q term = ONE, sum = ONE;
    for (int n = 1; n <= 24; n++) { term = mulq(term, r) / n; sum += term; }
    if (k >= 34) return (q)1 << 62;                   /* out of our range: never reached (our arguments stay in [-48, 1]) */
    if (k <= -62) return 0;                           /* e^x below 2^-62 of one: zero in Q.28 (a shift by 63 or more is undefined in C) */
    if (k >= 0) return sum << k;
    return sum >> (-k);
}
static q sqrtq(q x)                                   /* floor(sqrt(x)) in Q.FB: Newton on the integer sqrt of x*2^FB */
{
    if (x <= 0) return 0;
    uint64_t n = (uint64_t)x << FB;                   /* x < 2^35 so n < 2^63 */
    uint64_t r = (uint64_t)1 << 31;
    for (int i = 0; i < 64; i++) { uint64_t nr = (r + n / r) >> 1; if (nr >= r) { if (nr - r <= 1) break; } r = nr; }
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return (q)r;
}
static void fmt(q v, char *out)                       /* a Q.FB value as a decimal with 9 places, exact (no float) */
{
    const int neg = v < 0;
    uint64_t a = (uint64_t)(neg ? -v : v);
    uint64_t ip = a >> FB, fp = a & MASK;
    uint64_t dec = (fp * 1000000000ull) >> FB;        /* fp < 2^28, times 1e9 < 2^58: fits */
    sprintf(out, "%s%llu.%09llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dec);
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
    *s += mulq(h, k1s + 2 * k2s + 2 * k3s + k4s) / 6;
    *sp += mulq(h, k1p + 2 * k2p + 2 * k3p + k4p) / 6;
    *L += mulq(h, k1L + 2 * k2L + 2 * k3L + k4L) / 6;
}

static uint64_t fnv(uint64_t h, const void *p, size_t n) { const uint8_t *b = p; for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 0x100000001b3ull; } return h; }

int main(void)
{
    char a[40], b[40], c[40], d[40];
    const q h = ONE >> 12;                            /* rho step 2^-12 */
    const q sqrt3 = sqrtq(3 * ONE), sqrt19 = sqrtq(19 * ONE), k = sqrt19 - sqrt3;
    fmt(sqrt3, a); fmt(k, b);
    printf("WHFLOW: Q.%d fixed point, int64, exp by ln2 reduction + 24 Taylor terms, sqrt by Newton, RK4 h = 2^-12; sqrt3 = %s, k = sqrt19 - sqrt3 = %s\n", FB, a, b);
    for (int sphere = 0; sphere < 2; sphere++) {
        uint64_t hash = 0xcbf29ce484222325ull;
        /* the IR fixed point: sphere 0: s* = 0 (e^-4s = 1); sphere 1: u = e^-2s solves 4 + u - 4u^2 = 0, u = (1 + sqrt65)/8 */
        q sstar = 0;
        if (sphere) {
            const q u = divq(ONE + sqrtq(65 * ONE), 8 * ONE);     /* e^{-2 s*} */
            const q ems = sqrtq(u);                                /* e^{-s*} */
            /* s* = -ln(ems): Newton on exp */
            q s0 = 0; for (int i = 0; i < 60; i++) { const q e = expq(-s0); s0 += divq(e - ems, e); }
            sstar = s0;
            fmt(ems, a); fmt(sstar, b); fmt(Lprime(sstar, 0, 1), c);
            printf("WHFLOW sphere: IR fixed point e^-s* = %s (their A.31 expansion 1 + 1/16 = 1.0625), s* = %s, L' there = %s\n", a, b, c);
        } else {
            fmt(Lprime(0, 0, 0), a);
            printf("WHFLOW flat: IR fixed point s* = 0, L' there = %s (sqrt3)\n", a);
        }
        /* start deep in the IR on the growing linearized mode: s - s* = eps e^{k rho}, s' = k (s - s*), with eps = 2^-20 at rho0 = -8 */
        q rho = -8 * ONE, s = sstar + (ONE >> 20), sp = mulq(k, ONE >> 20), L = mulq(Lprime(sstar, 0, sphere), rho);
        int steps = 0, printed = 0;
        const char *name = sphere ? "sphere" : "flat";
        printf("WHFLOW %s: rho | sigma | sigma' | lambda' | lambda - rho | sigma - rho\n", name);
        while (rho < 10 * ONE) {
            if ((steps & 1023) == 0) {                /* every 0.25 in rho */
                const q lp = Lprime(s, sp, sphere);
                fmt(rho, a); fmt(s, b); fmt(sp, c); fmt(lp, d);
                char e[40], f[40]; fmt(L - rho, e); fmt(s - rho, f);
                printf("WHFLOW %s: %s | %s | %s | %s | %s | %s\n", name, a, b, c, d, e, f);
                printed++;
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
