/* wh_flow.c (11-11-2026; v2 4:10 PM; v3 10:05 PM) -- WORMHOLE MATH 1a: the Maldacena-Milekhin flow, exact.
 *
 * Maldacena & Milekhin, "Humanly traversable wormholes", arXiv:2008.06618, eq. (3.19) with l_B = 1 (their Appendix A, eq. A.28-A.29):
 *   0 = -6 + 3 e^{-4s} + L'^2 + 4 L' s' + s'^2            (so  L' = -2 s' + sqrt(3 s'^2 + 6 - 3 e^{-4s}),  the root with L' = sqrt3 at s = s' = 0)
 *   0 =  4 - 4 e^{-4s} - 2 L' s' - 2 s'^2 - s''           (so  s'' = 4 - 4 e^{-4s} - 2 L' s' - 2 s'^2)
 * where s = sigma(rho) and L = lambda(rho). IR (rho -> -inf): s -> 0, L' -> sqrt3 (AdS3 x R2); UV (rho -> +inf): s ~ L ~ rho + const (AdS5).
 * Their boundary condition: the growing solution of the linearized equation around s = 0, s'' + 2 sqrt3 s' - 16 s = 0, i.e. s = eps e^{k rho},
 * k = sqrt(19) - sqrt(3) = 2.6268..., integrated toward the UV. Their Figure 3 shows sigma and lambda; this program prints them as a table.
 * And the flow on a two-sphere, their eq. (A.30) with l_B-hat = 1:
 *   0 = -6 - e^{-2s} + 3 e^{-4s} + L'^2 + 4 L' s' + s'^2 ;  0 = 4 + e^{-2s} - 4 e^{-4s} - 2 L' s' - 2 s'^2 - s''
 * whose IR fixed point solves 4 + u - 4 u^2 = 0, u = e^{-2s}; their (A.31) says e^{-s} = l + l^3/16 + ... (l = 1: 1.0625). Its own growing
 * exponent is sqrt(19 + 9u/4) - sqrt(3 + u/4) = 2.8301 (v2 had started it on the flat flow's 2.6268: 97.8% growing, 2.2% decaying).
 *
 * v3 (the review's spec of 10:00 PM, four changes, none to whq.h; its reference program reproduced here line for line):
 *  (1) THE UPDATE ROUNDS ONCE. Each variable's RK4 increment h X / 6 (X = k1 + 2 k2 + 2 k3 + k4, an exact int64) is one exact int64 floor
 *      division by D = 6 x 2^12, and the remainder is carried to the next step in an int64 per variable. After n steps the state is exactly
 *      v0 + floor((X_1 + ... + X_n + D/2) / D): the running sum of the increments rounded once, an error in [-1/2, 1/2) unit at every step
 *      that never grows. v2 rounded twice per update (mulq by h, then rdiv by 6) and rdiv's ties-away rule gave +1/12 unit per step on
 *      average, a bias of 4-7e-8 on the UV constants that grew as h shrank (the review, MEASURED against binary128); v3's error is at most
 *      3.7e-11 and does not depend on h. The dead zone is gone too: small increments accumulate in the remainder instead of rounding to 0,
 *      so the UV slopes settle at exactly 1.
 *  (2) EACH FLOW STARTS ON ITS OWN UNSTABLE MANIFOLD: the point (sigma, sigma', lambda) at rho = -8 with the linear coordinate x = 2^-10,
 *      from the manifold's series to order 4 (the first order whose Q.40 integers are final; a second-order start still misses by 1.1e-6
 *      in rho), lambda in the canonical gauge (lambda - Lambda* rho -> 0 as rho -> -inf). The integers are hard-coded (START): COMPUTED
 *      by the review in 60-digit arithmetic and re-derived by its verifier. v2's start (s* + 2^-10, k 2^-10) was off the manifold at second
 *      order, a 5.9e-4 shift in rho.
 *  (3) stop at sigma > 16 (v2: 12), so the sphere's e^{-2 sigma}/6 < 2.1e-15 and its sigma - rho has converged too.
 *  (4) THE GAUGE INVARIANTS: K_sigma = (sigma - rho)_end + GAMMA_sigma and K_lambda = (lambda - rho)_end + GAMMA_lambda, the UV constants
 *      carried to the gauge in which the IR mode is e^{k rho} with unit amplitude (GAMMA: the review's exact gauge constants, COMPUTED, the
 *      same in 8 gauges to 2e-15). sigma - rho alone depends on where rho = 0 is put; K does not.
 *
 * EXACT: the Q.40 arithmetic of whq.h (int64 with 128-bit intermediates, products and divisions rounded to nearest, exp by ln2 reduction
 * and 24 Taylor terms, sqrt exact to the floor); RK4 with a fixed step h = 2^-12, the stages as in v2, the update as above. No floating
 * point anywhere, no libm: the same bits on every machine and every run. v1 (Q.28, floor) is withdrawn: its first increments in the deep IR
 * were a sixth of a unit, floored to zero, so the flow would have stood still; found 4:08 PM before any run.
 * Build:  cc -O2 -std=gnu11 -o whflow wh_flow.c      (whq.h beside it)
 */
#include "whq.h"

/* the two flows: sphere = 0 is (3.19)/(A.28-29), sphere = 1 is (A.30) with l_B-hat = 1 -- unchanged from v2 */
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
static q Rs, Rp, RL;                                 /* v3: the carried remainders of sigma, sigma', lambda; 0 at the start of each flow */
static q upd(q X, q *R)                              /* v3: the increment h X / 6 for h = 2^-12, rounded once (half up), remainder kept in *R */
{
    const q D = (q)6 << 12;                          /* 24576 */
    const q T = *R + X + D / 2;                      /* exact: |X| < 2^47 */
    q inc = T / D;
    if (T % D < 0) inc--;                            /* floor, also for T < 0 */
    *R = T - D / 2 - inc * D;                        /* = R + X - inc D, in [-D/2, D/2) */
    return inc;
}
static void rk4(q *s, q *sp, q *L, q h, int sphere)   /* one step of (s, s', L) with s'' = spp and L' = Lprime; stages as v2 */
{
    const q k1s = *sp, k1p = spp(*s, *sp, sphere), k1L = Lprime(*s, *sp, sphere);
    const q h2 = h / 2;
    const q s2 = *s + mulq(h2, k1s), p2 = *sp + mulq(h2, k1p);
    const q k2s = p2, k2p = spp(s2, p2, sphere), k2L = Lprime(s2, p2, sphere);
    const q s3 = *s + mulq(h2, k2s), p3 = *sp + mulq(h2, k2p);
    const q k3s = p3, k3p = spp(s3, p3, sphere), k3L = Lprime(s3, p3, sphere);
    const q s4 = *s + mulq(h, k3s), p4 = *sp + mulq(h, k3p);
    const q k4s = p4, k4p = spp(s4, p4, sphere), k4L = Lprime(s4, p4, sphere);
    *s += upd(k1s + 2 * k2s + 2 * k3s + k4s, &Rs);   /* v3 (v2: *s += rdiv(mulq(h, ...), 6)) */
    *sp += upd(k1p + 2 * k2p + 2 * k3p + k4p, &Rp);
    *L += upd(k1L + 2 * k2L + 2 * k3L + k4L, &RL);
}
/* v3: (sigma, sigma', lambda) at rho = -8 on each flow's unstable manifold, x = 2^-10, units of 2^-40 (the review's integers, 10:00 PM) */
static const q START[2][3] = { { 1072463084LL, 2813843829LL, -15236010697140LL },        /* flat */
                               { -67469346530LL, 3031326166LL, -15938931920661LL } };    /* sphere */
/* v3: K_sigma = (sigma - rho)_end + GAMMA[.][0], K_lambda = (lambda - rho)_end + GAMMA[.][1], units of 2^-40 (the review's constants) */
static const q GAMMA[2][2] = { { -5894808496714LL, 4315299320483LL }, { -6103177591682LL, 4955530204813LL } };

int main(void)
{
    char a[48], b[48], c[48], d[48], e[48], f[48];
    const q h = ONE >> 12;                            /* rho step 2^-12 */
    const q sqrt3 = sqrtq(3 * ONE), sqrt19 = sqrtq(19 * ONE), k = sqrt19 - sqrt3;
    printf("WHFLOW: Q.%d fixed point, int64 with 128-bit intermediates, round to nearest; exp by ln2 reduction + 24 Taylor terms; sqrt exact; RK4 h = 2^-12, "
           "each update rounded once, remainder carried (v3); start on the unstable manifold, x = 2^-10 at rho = -8\n", FB);
    fmt(sqrt3, a); fmt(sqrt19, b); fmt(k, c); fmt(expq(LN2), d); fmt(expq(-LN2), e); fmt(divq(ONE, 3 * ONE), f);
    printf("WHFLOW selftest: sqrt3 = %s (1.732050807569), sqrt19 = %s (4.358898943541), k = %s (2.626848135972), e^ln2 = %s, e^-ln2 = %s, 1/3 = %s\n", a, b, c, d, e, f);
    fmt(mulq(3 * ONE, 5 * ONE), a); fmt(sqrtq(4 * ONE), b); fmt(expq(0), c); fmt(mulq(sqrt3, sqrt3), d); fmt(expq(-48 * ONE), e);
    printf("WHFLOW selftest: 3*5 = %s, sqrt4 = %s, e^0 = %s, sqrt3^2 = %s (3 - at most 2^-38), e^-48 = %s\n", a, b, c, d, e);
    for (int sphere = 0; sphere < 2; sphere++) {
        uint64_t hash = 0xcbf29ce484222325ull;
        q kown = k;                                   /* v3: the flow's own growing exponent (printed only) */
        if (sphere) {
            const q u = divq(ONE + sqrtq(65 * ONE), 8 * ONE);     /* e^{-2 s*} */
            const q ems = sqrtq(u);                                /* e^{-s*} */
            q s0 = 0; for (int i = 0; i < 60; i++) { const q ex = expq(-s0); s0 += divq(ex - ems, ex); }
            const q sstar = s0;
            fmt(ems, a); fmt(sstar, b); fmt(Lprime(sstar, 0, 1), c); fmt(spp(sstar, 0, 1), d);
            printf("WHFLOW sphere: IR fixed point e^-s* = %s (their A.31 expansion 1 + 1/16 = 1.0625), s* = %s, L' there = %s, s'' there = %s (0)\n", a, b, c, d);
            kown = sqrtq(19 * ONE + divq(9 * u, 4 * ONE)) - sqrtq(3 * ONE + divq(u, 4 * ONE));   /* sqrt(19 + 9u/4) - sqrt(3 + u/4) */
        } else {
            fmt(Lprime(0, 0, 0), a); fmt(spp(0, 0, 0), b);
            printf("WHFLOW flat: IR fixed point s* = 0, L' there = %s (sqrt3), s'' there = %s (0)\n", a, b);
        }
        q rho = -8 * ONE, s = START[sphere][0], sp = START[sphere][1], L = START[sphere][2];
        Rs = Rp = RL = 0;
        const char *name = sphere ? "sphere" : "flat";
        fmt(kown, a); fmt(s, b); fmt(sp, c); fmt(L, d);
        printf("WHFLOW %s: own mode k = %s (%s), start on its unstable manifold at rho = -8: sigma = %s, sigma' = %s, lambda = %s (raw %lld %lld %lld)\n",
               name, a, sphere ? "2.830105157339" : "2.626848135972", b, c, d, (long long)s, (long long)sp, (long long)L);
        printf("WHFLOW %s: rho | sigma | sigma' | lambda' | lambda - rho | sigma - rho\n", name);
        int steps = 0;
        while (rho < 12 * ONE) {
            if ((steps & 1023) == 0) {                /* every 0.25 in rho */
                const q lp = Lprime(s, sp, sphere);
                fmt(rho, a); fmt(s, b); fmt(sp, c); fmt(lp, d); fmt(L - rho, e); fmt(s - rho, f);
                printf("WHFLOW %s: %s | %s | %s | %s | %s | %s\n", name, a, b, c, d, e, f);
            }
            rk4(&s, &sp, &L, h, sphere);
            hash = fnv(hash, &s, 8); hash = fnv(hash, &sp, 8); hash = fnv(hash, &L, 8);
            rho += h; steps++;
            if (s > 16 * ONE) break;                  /* v3: e^{-2 sigma}/6 < 2.1e-15, so sigma - rho has converged for the sphere too */
        }
        const q lp = Lprime(s, sp, sphere);
        fmt(rho, a); fmt(sp, b); fmt(lp, c); fmt(s - rho, d); fmt(L - rho, e); fmt(L - s, f);
        printf("WHFLOW %s: END at rho = %s after %d steps: sigma' = %s (UV: 1), lambda' = %s (UV: 1), sigma - rho = %s, lambda - rho = %s, lambda - sigma = %s; chain %016llx\n",
               name, a, steps, b, c, d, e, f, (unsigned long long)hash);
        fmt(s - rho + GAMMA[sphere][0], a); fmt(L - rho + GAMMA[sphere][1], b); fmt(GAMMA[sphere][0], c); fmt(GAMMA[sphere][1], d);
        printf("WHFLOW %s: invariants K_sigma = %s (%s), K_lambda = %s (%s); K = (end value) + gauge constant %s, %s\n",
               name, a, sphere ? "0.377362726342" : "0.407800114190", b, sphere ? "-0.297030895775" : "-0.271699935370", c, d);
    }
    printf("WHFLOW: done\n");
    return 0;
}
