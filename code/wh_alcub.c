/* wh_alcub.c (11-11-2026) -- WORMHOLE MATH 3a: the Alcubierre negative-energy integral, exact.
 *
 * Bobrick & Martire, "Introducing physical warp drives", arXiv:2102.06824 (Class. Quantum Grav. 38, 105009, 2021), in G = c = 1 units:
 *   eq. (3)    T^00 = -(1/8pi) (rho^2 v^2 / 4 r^2) f'(r)^2          (negative wherever f' != 0)
 *   eq. (A.7)  E    = -(v^2/12) INT_0^inf r^2 f'(r)^2 dr            (the total energy seen by Eulerian observers)
 *   their optimal wall f = min(r0/r, 1): f' = -r0/r^2 for r > r0, so INT = r0 exactly and E = -v^2 r0/12.
 *   Alcubierre's wall (their eq. 2): f = [tanh(s(r+R)) - tanh(s(r-R))] / (2 tanh(sR)),  f' = s [sech^2(s(r+R)) - sech^2(s(r-R))] / (2 tanh(sR)).
 *   Their text: the optimal wall "decreases the needed (absolute value of) negative energy of the Alcubierre solution by about a factor of 3".
 *   Thin-wall estimate on paper: INT r^2 f'^2 dr ~ R^2 s/3 (sech^4 integrates to 4/3), so the factor is R s / 3, not a constant: 3 at R s = 9.
 *   Their section 4.2: flattening the bubble along x by a factor a (f(r) -> f(sqrt((a x)^2 + rho^2))) divides E by a exactly (a change of variable).
 *
 * This program prints, in units R = 1 and v = 1 (every quantity dimensionless; multiply INT by R to get meters):
 *   I(s) = INT_0^4 r^2 f'^2 dr by Simpson's rule for sR = 4, 6, 8, 9, 10, 12, 16; the ratio I/I_opt with I_opt = r0 = 1; the thin-wall
 *   asymptote s/3 and the scaled remainder s (I - s/3); the optimal wall's own Simpson integral over [1, 4] against its exact 1 - 1/4;
 *   the 2D check J(a) = INT_0^{4/a} dx INT_0^4 rho (df/drho)^2 drho for a = 1, 2, 8 (df/drho = f'(r) rho / r), which must give J(a) a = J(1)
 *   and J(1) = (2/3) I (their E = -(v^2/16) INT dx INT rho drho (df/drho)^2 over all x equals -(v^2/12) I, so the half space gives 2/3 I).
 *
 * EXACT: every quantity an int64 in Q.28 (28 fractional bits); products by mulq() (a 64x64 -> 64 product in two halves, round to nearest);
 * exp by range reduction on ln 2 and a 24-term Taylor series; sqrt by integer Newton (floor); Simpson's rule with N = 2^14 intervals on
 * [0, 4] (h = 2^-12 exactly) in 1D and 2^10 x 2^10 in 2D (steps 2^-8, 2^-9, 2^-11 by 2^-8, all exact powers of two, so every product by a
 * step is a shift). No floating point anywhere, no libm: the same bits on every machine and every run. The result IS the printed digits.
 * Build:  cc -O2 -std=gnu11 -o whalcub wh_alcub.c
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define FB 28
#define ONE ((int64_t)1 << FB)
#define HALF ((int64_t)1 << (FB - 1))
#define MASK (ONE - 1)
typedef int64_t q;

static q mulq(q a, q b)                               /* round(a*b / 2^FB) to nearest (half up), exact in two halves; needs |a| < 2^35 */
{
    const q bh = b >> FB, bl = b & MASK;
    return a * bh + ((a * bl + HALF) >> FB);
}
static q divq(q a, q b)                               /* round(a*2^FB / b) to nearest for b > 0; needs |a| < 2^35 */
{
    const q n = a << FB;
    return n >= 0 ? (n + b / 2) / b : -((-n + b / 2) / b);
}
static const q LN2 = 186065279;                       /* ln 2 * 2^28 = 186065279.36... */
static q expq(q x)                                    /* e^x in Q.FB for x in about [-40, 20]: 2^k * e^r, |r| <= ln2/2, Taylor 24 terms */
{
    int64_t k = x / LN2;
    q r = x - k * LN2;
    if (r > LN2 / 2) { r -= LN2; k++; }
    if (r < -LN2 / 2) { r += LN2; k--; }
    q term = ONE, sum = ONE;
    for (int n = 1; n <= 24; n++) { term = mulq(term, r) / n; sum += term; }
    if (k >= 34) return (q)1 << 62;
    if (k <= -62) return 0;
    if (k >= 0) return sum << k;
    return sum >> (-k);
}
static q sqrtq(q x)                                   /* floor(sqrt(x)) in Q.FB: Newton on the integer sqrt of x*2^FB, from a bit-length guess */
{
    if (x <= 0) return 0;
    uint64_t n = (uint64_t)x << FB;
    int bits = 0; for (uint64_t t = n; t; t >>= 1) bits++;
    uint64_t r = (uint64_t)1 << ((bits + 1) / 2);
    for (int i = 0; i < 64; i++) { uint64_t nr = (r + n / r) >> 1; if (nr >= r) break; r = nr; }
    while (r * r > n) r--;
    while ((r + 1) * (r + 1) <= n) r++;
    return (q)r;
}
static void fmt(q v, char *out)                       /* a Q.FB value as a decimal with 9 places, exact (no float) */
{
    const int neg = v < 0;
    uint64_t a = (uint64_t)(neg ? -v : v);
    uint64_t ip = a >> FB, fp = a & MASK;
    uint64_t dec = (fp * 1000000000ull) >> FB;
    sprintf(out, "%s%llu.%09llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dec);
}
static uint64_t fnv(uint64_t h, const void *p, size_t n) { const uint8_t *b = p; for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 0x100000001b3ull; } return h; }

/* sech^2(x) = 4 e^{-2|x|} / (1 + e^{-2|x|})^2, written so that the exponent is never positive */
static q sech2(q x)
{
    if (x < 0) x = -x;
    const q e = expq(-2 * x);
    const q d = ONE + e;
    return divq(4 * e, mulq(d, d));
}
static q tanhq(q x)                                   /* x >= 0: (1 - e^{-2x}) / (1 + e^{-2x}) */
{
    const q e = expq(-2 * x);
    return divq(ONE - e, ONE + e);
}
/* f'(r) for Alcubierre's wall with R = 1: s [sech^2(s(r+1)) - sech^2(s(r-1))] / (2 tanh s) */
static q fprime(q r, q s, q two_tanh_s)
{
    const q a = sech2(mulq(s, r + ONE)), b = sech2(mulq(s, r - ONE));
    return divq(mulq(s, a - b), two_tanh_s);
}

int main(void)
{
    char A[40], B[40], C[40], D[40], E[40];
    uint64_t hash = 0xcbf29ce484222325ull;
    printf("WHALCUB: Q.%d fixed point, int64, products rounded to nearest, exp by ln2 reduction + 24 Taylor terms, Simpson N = 2^14 on [0, 4] (1D), 2^10 x 2^10 (2D); R = 1, v = 1\n", FB);

    /* the optimal wall over [1, 4]: INT 1/r^2 dr = 1 - 1/4 exactly; the tail to infinity is 1/4 */
    {
        const int N = 1 << 14; const q h = (3 * ONE) / N;            /* 3/16384 is not a power of two: the step is this integer, declared */
        q sum = 0;
        for (int i = 0; i <= N; i++) {
            const q r = ONE + (q)i * h;
            const q t = divq(ONE, mulq(r, r));                         /* r^2 (r0/r^2)^2 = 1/r^2 */
            const int w = (i == 0 || i == N) ? 1 : (i & 1) ? 4 : 2;
            sum += w * t;
        }
        const q I = sum >> 14;                                          /* sum * h / 3 with h = 3 * 2^14 in Q.28: (sum * 3 * 2^14 / 3) >> 28 = sum >> 14, exact */
        fmt(I, A); fmt(ONE - ONE / 4, B);
        printf("WHALCUB optimal wall f = min(1/r, 1): Simpson INT_1^4 r^2 f'^2 dr = %s (exact 1 - 1/4 = %s; the tail INT_4^inf = 1/4; so I_opt = 1 exactly)\n", A, B);
    }

    /* Alcubierre's wall: I(s) for s R = 4, 6, 8, 9, 10, 12, 16 */
    const int SV[7] = {4, 6, 8, 9, 10, 12, 16};
    q Ival[7];
    printf("WHALCUB tanh wall: sR | I = INT_0^4 r^2 f'^2 dr | I / I_opt | s/3 (thin wall) | s (I - s/3)\n");
    for (int si = 0; si < 7; si++) {
        const q s = (q)SV[si] * ONE, tt = 2 * tanhq(s);
        const int N = 1 << 14; const q h = ONE >> 12;                 /* 4 / 2^14 = 2^-12 exactly */
        q sum = 0;
        for (int i = 0; i <= N; i++) {
            const q r = (q)i * h;
            const q g = mulq(r, fprime(r, s, tt));                     /* r f' , |g| <= 32 */
            const q t = mulq(g, g);                                    /* r^2 f'^2 <= 1024 */
            const int w = (i == 0 || i == N) ? 1 : (i & 1) ? 4 : 2;
            sum += w * t;
            hash = fnv(hash, &t, 8);
        }
        const q I = (sum >> 12) / 3;                                   /* sum * h / 3 with h = 2^-12: a shift, then the division */
        Ival[si] = I;
        const q third = s / 3, rem = mulq(s, I - third);
        fmt(I, A); fmt(I, B); fmt(third, C); fmt(rem, D);
        printf("WHALCUB tanh wall: %2d | %s | %s | %s | %s\n", SV[si], A, B, C, D);
    }

    /* the 2D flattening check at sR = 8: J(a) = INT_0^{4/a} dx INT_0^4 rho (f'(r) rho / r)^2 drho, a = 1, 2, 8 */
    {
        const q s = 8 * ONE, tt = 2 * tanhq(s);
        const int AV[3] = {1, 2, 8}, SH[3] = {8, 9, 11};              /* hx = 2^-8, 2^-9, 2^-11 (4/a over 2^10 points); hrho = 2^-8 */
        const int N = 1 << 10;
        q J1 = 0;
        printf("WHALCUB flattening (sR = 8): a | J(a) | a J(a) | J(1) | (2/3) I(8)\n");
        for (int ai = 0; ai < 3; ai++) {
            const q hx = ONE >> SH[ai], hr = ONE >> 8;
            int64_t S = 0;                                              /* sum of wi wj rho (df/drho)^2, in Q.28 units: < 2^60 */
            for (int i = 0; i <= N; i++) {
                const q x = (q)i * hx, ax = (q)AV[ai] * x;
                const int wi = (i == 0 || i == N) ? 1 : (i & 1) ? 4 : 2;
                int64_t Sr = 0;
                for (int j = 0; j <= N; j++) {
                    const q rho = (q)j * hr;
                    q t = 0;
                    if (j > 0) {
                        const q r = sqrtq(mulq(ax, ax) + mulq(rho, rho));
                        const q fp = fprime(r, s, tt);
                        const q u = divq(rho, r);                         /* rho / r <= 1 */
                        const q dfr = mulq(fp, u);                        /* df/drho */
                        t = mulq(rho, mulq(dfr, dfr));                    /* rho (df/drho)^2 <= 4 * 16 */
                    }
                    const int wj = (j == 0 || j == N) ? 1 : (j & 1) ? 4 : 2;
                    Sr += wj * t;
                }
                S += wi * (Sr >> 8);                                    /* times hrho = 2^-8 */
                hash = fnv(hash, &Sr, 8);
            }
            const q J = (S >> SH[ai]) / 9;                              /* times hx, over 9 */
            if (ai == 0) J1 = J;
            fmt(J, A); fmt(J * AV[ai], B); fmt(J1, C); fmt((2 * Ival[2]) / 3, D);
            printf("WHALCUB flattening (sR = 8): %d | %s | %s | %s | %s\n", AV[ai], A, B, C, D);
        }
    }
    (void)E;
    printf("WHALCUB: chain %016llx\n", (unsigned long long)hash);
    printf("WHALCUB: SI: E = -(v^2/12) I R c^4/G with c = 299792458 m/s exactly and G = 6.67430e-11 m^3 kg^-1 s^-2 (CODATA 2018): c^4/G = 1.2104e44 J per meter of I R\n");
    printf("WHALCUB: done\n");
    return 0;
}
