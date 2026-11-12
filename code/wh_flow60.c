/* wh_flow60.c (11-12-2026, 3:10 AM) -- WORMHOLE MATH 1a, THE FLOW IN Q.56: wh_flow.c version 3's two Maldacena-Milekhin flows (flat (3.19)/(A.28-29),
 * sphere (A.30)) and the constant Xi of the flat flow (defined in wh_ode.c; its meaning for the paper is under review), integrated in a fixed-point format with 56 fraction bits instead of 40,
 * so that the gauge invariants K_sigma, K_lambda and Xi come out to twelve digits instead of ten. The state (sigma, sigma', lambda, J) is held
 * in Q.56 (values below 128: sigma reaches 16 at the stop, past Q.60's ceiling of 8); every product is a 128-bit product rounded once to
 * nearest; the exponentials e^{-4 sigma}, e^{-2 sigma} and e^{2(lambda - sigma)} come from whq60.h's exp60n (ln 2 reduction, Taylor to a
 * unit of 2^-60, an integer-plus-fraction entry since the exponents reach -64) and are rounded once to Q.56; the square root in Lprime is the
 * 128-bit digit-by-digit root of the Q.112 radicand; RK4 with the step 2^-12, each component's update h (k1 + 2 k2 + 2 k3 + k4)/6 formed in
 * 128 bits and rounded once with its remainder carried (wh_flow.c v3's rule). The start: wh_flow.c v3's manifold integers at rho = -8
 * (Q.40, the review's) scaled to Q.56, and J's exact tail integer (the review's, 3:08 AM) likewise, until the review supplies Q.56-exact ones;
 * with the Q.40-exact start the K's inherit its 2^-40 error, Xi (translation- and gauge-invariant) does not.
 * ACCEPTANCE: the review's 30-digit values K_sigma = 0.4078001141903564, K_lambda = -0.2716999353696263 (flat), 0.377362726342 /
 * -0.297030895775 (sphere), Xi = 0.12178902947580307791; the Q.40 flow's values are reproduced to its ten digits.
 * Usage: whflow60 [h=12]      Build: cc -O2 -std=gnu11 -o whflow60 wh_flow60.c   (whq.h and whq60.h beside it) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "whq60.h"

#define FS 56                                           /* the state's fraction bits */
static const int64_t ONE56 = (int64_t)1 << FS;
static inline int64_t mul56(int64_t a, int64_t b) { s128 s = {0, 0}; s128_addp(&s, a, b); return s128_round_sh(s, FS); }
static inline int64_t q60to56(int64_t v) { s128 s = { v < 0 ? -1 : 0, (uint64_t)v }; return s128_round_sh(s, 4); }
static int64_t exp56_of(int64_t a, int64_t v)        /* e^{a v} for an integer a and a Q.56 v, as a Q.56 value: the exponent a v in 128 bits, split */
{
    s128 X = {0, 0}; s128_addp(&X, a, v);               /* a v in Q.56, 128 bits */
    /* to Q.60 integer + fraction: n = floor(X / 2^56), f60 = (X - n 2^56) 2^4 */
    const int64_t n = (int64_t)(((uint64_t)X.hi << 8) | (X.lo >> 56));
    const int64_t f60 = (int64_t)((X.lo & (((uint64_t)1 << 56) - 1)) << 4);
    if (n < -62) return 0;
    if (n > 2) { printf("WHFLOW60: an exponential above the ceiling -- stop\n"); exit(1); }
    return q60to56(exp60n(n, f60));
}
static int64_t sqrt56_of128(s128 X)                   /* sqrt of a Q.112 value (two Q.56 factors' product), as Q.56 */
{
    if (X.hi < 0) { printf("WHFLOW60: a negative radicand -- stop\n"); exit(1); }
    return (int64_t)isqrt128((uint64_t)X.hi, X.lo);
}
/* the two flows: sphere = 0 is (3.19)/(A.28-29), sphere = 1 is (A.30) with l_B-hat = 1 */
static int64_t Lprime(int64_t s, int64_t sp, int sphere)
{
    const int64_t e4 = exp56_of(-4, s), e2 = sphere ? exp56_of(-2, s) : 0;
    s128 under = {0, 0};                                /* (2s')^2 - (s'^2 + 3 e^-4s + e^-2s - 6) = 3 s'^2 + 6 - 3 e4 + e2, in Q.112 */
    s128_addp(&under, 3 * sp, sp); s128_addp(&under, 6 * ONE56 - 3 * e4 + e2, ONE56);
    return -2 * sp + sqrt56_of128(under);
}
static int64_t spp(int64_t s, int64_t sp, int sphere)
{
    const int64_t e4 = exp56_of(-4, s), e2 = sphere ? exp56_of(-2, s) : 0;
    const int64_t lp = Lprime(s, sp, sphere);
    s128 r = {0, 0};                                    /* 4 + e2 - 4 e4 - 2 lp sp - 2 sp^2, in Q.112 */
    s128_addp(&r, 4 * ONE56 + e2 - 4 * e4, ONE56); s128_addp(&r, -2 * lp, sp); s128_addp(&r, -2 * sp, sp);
    return s128_round_sh(r, FS);
}
#define XI_SHIFT (10 * ONE56)
static int HS = 12;
static int64_t REM[4];
static int64_t upd(s128 X, int64_t *R)                  /* the increment h X / 6 for h = 2^-HS, rounded once (half up), the remainder kept in *R */
{
    const int64_t D = (int64_t)6 << HS;
    s128 T = X; s128_addv(&T, *R + D / 2);
    int neg = T.hi < 0; if (neg) { T.lo = ~T.lo + 1; T.hi = ~T.hi + (T.lo == 0); }
    uint64_t rem; uint64_t qq = div128((uint64_t)T.hi, T.lo, (uint64_t)D, &rem);
    int64_t inc = (int64_t)qq;
    if (neg) { inc = -inc; if (rem) { inc -= 1; rem = (uint64_t)D - rem; } }
    *R = (int64_t)rem - D / 2;
    return inc;
}
static void rhs(const int64_t *y, int64_t *dy, int sphere)
{
    dy[0] = y[1]; dy[1] = spp(y[0], y[1], sphere); dy[2] = Lprime(y[0], y[1], sphere);
    s128 e = {0, 0}; s128_addp(&e, 2, y[2] + XI_SHIFT - y[0]);   /* 2 (lambda + 10 - sigma) */
    const int64_t n = (int64_t)(((uint64_t)e.hi << 8) | (e.lo >> 56)), f60 = (int64_t)((e.lo & (((uint64_t)1 << 56) - 1)) << 4);
    dy[3] = n < -62 ? 0 : q60to56(exp60n(n, f60));
}
static void rk4(int64_t *y, int64_t h, int sphere)
{
    int64_t k1[4], k2[4], k3[4], k4[4], y2[4], y3[4], y4[4]; const int64_t h2 = h / 2;
    rhs(y, k1, sphere);
    for (int i = 0; i < 4; i++) y2[i] = y[i] + mul56(h2, k1[i]);
    rhs(y2, k2, sphere);
    for (int i = 0; i < 4; i++) y3[i] = y[i] + mul56(h2, k2[i]);
    rhs(y3, k3, sphere);
    for (int i = 0; i < 4; i++) y4[i] = y[i] + mul56(h, k3[i]);
    rhs(y4, k4, sphere);
    for (int i = 0; i < 4; i++) { s128 X = { k1[i] < 0 ? -1 : 0, (uint64_t)k1[i] }; s128_addv(&X, 2 * k2[i]); s128_addv(&X, 2 * k3[i]); s128_addv(&X, k4[i]); y[i] += upd(X, &REM[i]); }
}
/* the v3 manifold integers at rho = -8 (Q.40, the review's) and J's exact tail (Q.40, the review's, gauge lambda + 10), scaled to Q.56 */
static const int64_t START40[2][3] = { { 1072463084LL, 2813843829LL, -15236010697140LL }, { -67469346530LL, 3031326166LL, -15938931920661LL } };
static const int64_t J0_40 = 141633791LL;
static const int64_t GAMMA40[2][2] = { { -5894808496714LL, 4315299320483LL }, { -6103177591682LL, 4955530204813LL } };
static void fmt56(int64_t v, char *out)                 /* twelve decimals of a Q.56 value */
{
    int neg = v < 0; uint64_t u = neg ? -(uint64_t)v : (uint64_t)v;
    uint64_t ip = u >> FS, frac = u & ((((uint64_t)1) << FS) - 1), h, l;
    mul128(frac, 1000000000000ull, &h, &l);
    uint64_t dg = (l >> FS) | (h << (64 - FS)); if ((l >> (FS - 1)) & 1) dg++;
    if (dg >= 1000000000000ull) { dg -= 1000000000000ull; ip++; }
    sprintf(out, "%s%llu.%012llu", neg ? "-" : "", (unsigned long long)ip, (unsigned long long)dg);
}
int main(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) if (!strncmp(argv[i], "h=", 2)) HS = atoi(argv[i] + 2);
    if (HS < 4 || HS > 30) { printf("WHFLOW60: h= must be 4..30\n"); return 1; }
    const int64_t h = ONE56 >> HS;
    char a[48], b[48], c[48], d[48], e[48], f[48];
    printf("WHFLOW60: Q.%d state (int64, 128-bit products rounded once), exponentials and roots from whq60.h (Q.60, rounded once to Q.56); RK4 with the step 2^-%d, each update rounded once with its remainder carried; start = v3's Q.40 manifold integers at rho = -8 scaled, J0 = the review's exact tail\n", FS, HS);
    const int64_t worst = whq60_selftest("WHFLOW60"); printf("WHFLOW60 selftest: largest deviation %lld units of 2^-60\n", (long long)worst);
    for (int sphere = 0; sphere < 2; sphere++) {
        uint64_t hash = 0xcbf29ce484222325ull;
        int64_t y[4] = { START40[sphere][0] << 16, START40[sphere][1] << 16, START40[sphere][2] << 16, J0_40 << 16 }, rho = -8 * ONE56;
        memset(REM, 0, sizeof REM);
        const char *name = sphere ? "sphere" : "flat";
        int steps = 0;
        while (rho < 12 * ONE56) {
            rk4(y, h, sphere);
            for (int i = 0; i < 4; i++) hash = fnv(hash, &y[i], 8);
            rho += h; steps++;
            if (y[0] > 16 * ONE56) break;
        }
        fmt56(rho, a); fmt56(y[1], b); fmt56(Lprime(y[0], y[1], sphere), c); fmt56(y[0] - rho, d); fmt56(y[2] - rho, e); fmt56(y[2] - y[0], f);
        printf("WHFLOW60 %s: END at rho = %s after %d steps: sigma' = %s (UV: 1), lambda' = %s (UV: 1), sigma - rho = %s, lambda - rho = %s, lambda - sigma = %s; chain %016llx\n",
               name, a, steps, b, c, d, e, f, (unsigned long long)hash);
        fmt56(y[0] - rho + (GAMMA40[sphere][0] << 16), a); fmt56(y[2] - rho + (GAMMA40[sphere][1] << 16), b);
        printf("WHFLOW60 %s: invariants K_sigma = %s (%s), K_lambda = %s (%s)   [the GAMMA constants and the start are Q.40-exact until the review's Q.56 ones arrive]\n",
               name, a, sphere ? "0.377362726342" : "0.407800114190", b, sphere ? "-0.297030895775" : "-0.271699935370");
        if (!sphere) {
            const int64_t xi = y[0] - mul56(exp56_of(2, y[0] - y[2] - XI_SHIFT), y[3]);
            fmt56(xi, a); fmt56(y[3], b); fmt56((ONE56 / 4) - xi, c);
            printf("WHFLOW60 flat: Xi = %s (the review's 30-digit value 0.121789029475803), J (gauge lambda + 10) = %s, kappa_flat = 1/4 - Xi = %s (0.128210970524)\n", a, b, c);
        }
    }
    printf("WHFLOW60: exp60 took at most %d Taylor terms\n", EXP60_TERMS);
    printf("WHFLOW60: done\n");
    return 0;
}
