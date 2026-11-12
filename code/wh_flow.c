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
 * EXACT: the Q.40 arithmetic of whq.h (int64 with 128-bit intermediates, products and divisions rounded to nearest, exp by ln2 reduction
 * and 24 Taylor terms, sqrt exact to the floor); RK4 with a fixed step h = 2^-12, its sixths rounded to nearest. No floating point anywhere, no libm: the same bits on every machine and
 * every run. The result IS the printed digits. v1 (Q.28, floor) is withdrawn: its first increments in the deep IR were a sixth of a unit,
 * floored to zero, so the flow would have stood still; found 4:08 PM before any run.
 * Build:  cc -O2 -std=gnu11 -o whflow wh_flow.c      (whq.h beside it)
 */
#include "whq.h"

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
